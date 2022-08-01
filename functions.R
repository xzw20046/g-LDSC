#rcpp calculator#####
library(Rcpp)
library(RcppEigen)
library(inline)
library(emulator)
library(rhdf5)
library(doParallel)

#
crossprodCpp<-'using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
const int m(A.rows()), n(A.cols());
MatrixXd AtA(MatrixXd(n, n).setZero().
selfadjointView<Lower>().rankUpdate(A.adjoint()));
//MatrixXd AAt(MatrixXd(m, m).setZero().
//selfadjointView<Lower>().rankUpdate(A));
return wrap(AtA);'

fcprd <- cxxfunction(signature(AA = "matrix"), crossprodCpp, "RcppEigen")

#
LDSCmatrix<-function(LDmatrix,Amatrix){
  cat<-ncol(Amatrix)
  res<-diag(LDmatrix)
  for (an in 1:cat){
    Rnew<-(LDmatrix**2)%*%as.vector(Amatrix[,an])
    res<-cbind(res,Rnew)
  }
  return(as.data.frame(res))
}

##LDSM function V2####
LDSMsum2<-function(LDmatrix,Amatrix){
  middle<-apply(Amatrix, 1, sum)
  Rnew<-LDmatrix*as.vector(sqrt(middle))
  res<-fcprd(Rnew)
  return(res)
}

#LDSM manual tau version
LDSMman<-function(LDmatrix,betavar){
  Rnew<-LDmatrix*(sqrt(betavar))
  res<-fcprd(Rnew)
  return(res)
}


##eigen cut function
eigen.cut<-function(vector,percentage){
  total<-sum(vector)
  k=1
  while (sum(vector[1:k])<percentage*total) {
    k=k+1
  }
  return(k)
}

gls.left.right<-function(i,x,y,block.names.all,LDmatrix,raw.Ntau,
                         LDSM,snp.list,cutpoint=0.99,intercept=T){
  block.name<-block.names.all[i]
  snp.list.temp<-snp.list[[block.name]]

  used.snp<-intersect(snp.list.temp,y$SNP)

  #message(paste0(i,' OK'))
  if(length(used.snp)<2){return(list(0,0))} else{
  used.loc<-match(used.snp, snp.list.temp)
  y.temp<-as.vector(y[match(used.snp,y$SNP),'Z']**2-1)
  #message(paste0(i,' Y OK'))
  x.temp<-as.matrix(x[[block.name]])[used.loc,]
  #message(paste0(i,' X OK'))
  #raw.Ntau<-mean(y.temp)/mean(apply(x.temp, 1, sum))

  if(intercept==T){x.temp<-cbind(rep(1,nrow(x.temp)),x.temp)}


  ldsm.temp<-LDSM[[block.name]][used.loc,used.loc]
  #message(paste0(i,'ldsm OK'))
  ldsm.temp2<-(raw.Ntau*ldsm.temp+LDmatrix[[block.name]][used.loc,used.loc])**2
  #message(paste0(i,'ldsm2 OK'))
  eigdi<-eigen(ldsm.temp2,symmetric = T)
  eigval<-eigdi$values
  #
  cut.point<-eigen.cut(eigval,cutpoint)
  eigval.inv<-(1/eigval[1:cut.point])

  omiga.inv<-crossprod(t(eigdi$vectors[,1:cut.point])*sqrt(eigval.inv))

  omiga.inv2<-t(x.temp)%*%omiga.inv

  gls.left<-omiga.inv2%*%x.temp

  gls.right<-as.vector(omiga.inv2%*%y.temp)

  return(list(gls.left,gls.right))
 }
}


##input left right output set of taus
gls.estimator<-function(left,right,A,M.anno,anno.total,N,intercept=T){
  #point estimator
  if(intercept==T){
    test.tau<-solve(left)%*%right

    tau.all<-as.vector(test.tau)/N
    tau<-tau.all[-1]
    constant=N*tau.all[1]+1

    est.h<-crossprod(A)%*%tau
    h.total<-est.h[1]
    e<-est.h/M.anno/h.total

    #hypothesis test
    M<-anno.total[1]
    h0.test<-(est.h/anno.total)-(h.total-est.h)/(M-anno.total)

    #variance in theory
    cov.tau.all<-solve(left)/N^2
    cov.tau<-cov.tau.all[-1,-1]
    constant.var=N^2*cov.tau.all[1,1]

    hv.var<-2*crossprod(A)%*%cov.tau%*%crossprod(A)
    h0.theory.var<-c()
    for (i in 2:length(anno.total)) {
      C<-anno.total[i]
      h0.theory.var[i]<-1/(C^2*(M-C)^2)*(C^2*hv.var[1,1]+M^2*hv.var[i,i]-2*M*C*hv.var[1,i])
    }
    P.theory<-pnorm(h0.test/sqrt(h0.theory.var), lower.tail = F)

    result<-list(tau, est.h, e, constant, sqrt(constant.var), h0.test, P.theory)
    names(result)<-c('Taus', 'Partition H2', 'Enrichment', 'intercept', 'intercept SD', 'e.stat','P')
    return(result)
  }else{
    test.tau<-solve(left)%*%right

    tau<-as.vector(test.tau)/N

    est.h<-crossprod(A)%*%tau
    h.total<-est.h[1]
    e<-est.h/M.anno/h.total

    #hypothesis test
    M<-anno.total[1]
    h0.test<-(est.h/anno.total)-(h.total-est.h)/(M-anno.total)

    #variance in theory
    cov.tau<-solve(left)/N^2

    hv.var<-2*crossprod(A)%*%cov.tau%*%crossprod(A)
    h0.theory.var<-c()
    for (i in 2:length(anno.total)) {
      C<-anno.total[i]
      h0.theory.var[i]<-1/(C^2*(M-C)^2)*(C^2*hv.var[1,1]+M^2*hv.var[i,i]-2*M*C*hv.var[1,i])
    }
    P.theory<-pnorm(h0.test/sqrt(h0.theory.var), lower.tail = F)

    result<-list(tau, est.h, e, 1, 0, h0.test, P.theory)
    names(result)<-c('Taus', 'Partition H2', 'Enrichment', 'intercept', 'intercept SD', 'e.stat','P')
    return(result)
  }
}

oneout<-function(i,raw.result,left.total,right.total,A,M.anno,anno.total,N,intercept){
  left=left.total-raw.result[[i]][[1]]
  right=right.total-raw.result[[i]][[2]]
  #the gls.generator
  if(intercept==T){
    test.tau<-solve(left)%*%right

    tau.all<-as.vector(test.tau)/N
    tau<-tau.all[-1]
    constant=tau.all[1]+1

    est.h<-crossprod(A)%*%tau
    h.total<-est.h[1]
    e<-est.h/M.anno/h.total

    #hypothesis test
    #M<-anno.total[1]
    #h0.test<-(est.h/anno.total)-(h.total-est.h)/(M-anno.total)

    #variance in theory
    #cov.tau.all<-solve(left)/N^2
    #cov.tau<-cov.tau.all[-1,-1]
    #constant.var=N^2*cov.tau.all[1,1]

    #hv.var<-2*crossprod(A)%*%cov.tau%*%crossprod(A)
    #h0.theory.var<-c()
    #for (i in 2:length(anno.total)) {
    #  C<-anno.total[i]
    #  h0.theory.var[i]<-1/(C^2*(M-C)^2)*(C^2*hv.var[1,1]+M^2*hv.var[i,i]-2*M*C*hv.var[1,i])
    #}
    #P.theory<-1-pnorm(h0.test/sqrt(h0.theory.var))

    result<-list(tau, est.h, e, constant)
    names(result)<-c('Taus', 'Partition H2', 'Enrichment')
    return(result)
  }else{
    test.tau<-solve(left)%*%right

    tau<-as.vector(test.tau)/N

    est.h<-crossprod(A)%*%tau
    h.total<-est.h[1]
    e<-est.h/M.anno/h.total

    #hypothesis test
    #M<-anno.total[1]
    #h0.test<-(est.h/anno.total)-(h.total-est.h)/(M-anno.total)

    #variance in theory
    #cov.tau<-solve(left)/N^2

    #hv.var<-2*crossprod(A)%*%cov.tau%*%crossprod(A)
    #h0.theory.var<-c()
    #for (i in 2:length(anno.total)) {
    #  C<-anno.total[i]
    #  h0.theory.var[i]<-1/(C^2*(M-C)^2)*(C^2*hv.var[1,1]+M^2*hv.var[i,i]-2*M*C*hv.var[1,i])
    #}
    #P.theory<-1-pnorm(h0.test/sqrt(h0.theory.var))

    result<-list(tau, est.h, e, 1)
    names(result)<-c('Taus', 'Partition H2', 'Enrichment', 'intercept')
    return(result)
    }
}



refpannel.par<-
function(batch,LDinfo=LDinfo,LD.path=LD.path,chrs=chrs[chr],
         anno=anno,maf.snp=maf.snp,gwas.snp=gwas.snp,
         used.a=used.a){
  batch.name<-LDinfo[batch,'name']

  #merge all snps
  LDinfo.temp<-h5read(file = paste(LD.path,chrs,sep = '/'),name = LDinfo[batch,'name'])
  snp1<-LDinfo.temp$snplist

  #snp1.intersect<-Reduce(intersect, list(snp1,anno[,'SNP'],maf.snp,gwas.snp,ldsc.snp))
  snp1.intersect<-Reduce(intersect, list(snp1,anno[,'SNP'],maf.snp,gwas.snp))

  batch.used=0
  if(length(snp1.intersect)>0){
    Anno.temp<-anno[match(snp1.intersect,anno[,'SNP']),used.a]

    Anno.temp<-as.matrix(Anno.temp)

    LD.temp<-LDinfo.temp$ldblk
    LD.temp<-LD.temp[match(snp1.intersect,snp1),match(snp1.intersect,snp1)]

    LDSM.temp<-LDSMsum2(LD.temp,Anno.temp)
    x.ldsc.temp<-LDSCmatrix(LD.temp,Anno.temp)[,-1]
    batch.used=1
  }else{snp1.intersect<-LD.temp<-LDSM.temp<-x.ldsc.temp<-batch.used<-0}
  return(list(batch.used,snp1.intersect,LD.temp,LDSM.temp,
              x.ldsc.temp))
}

refpannel.par.man<-
  function(batch,LDinfo,LD.path,chrs,
           anno,maf.snp,gwas.snp,
           used.a,betavar){
    batch.name<-LDinfo[batch,'name']

    #merge all snps
    LDinfo.temp<-h5read(file = paste(LD.path,chrs,sep = '/'),name = LDinfo[batch,'name'])

    snp1<-LDinfo.temp$snplist

    #snp1.intersect<-Reduce(intersect, list(snp1,anno[,'SNP'],maf.snp,gwas.snp,ldsc.snp))
    snp1.intersect<-Reduce(intersect, list(snp1,anno[,'SNP'],maf.snp,gwas.snp))

    batch.used=0
    if(length(snp1.intersect)>0){
      Anno.temp<-anno[match(snp1.intersect,anno[,'SNP']),used.a]

      Anno.temp<-as.matrix(Anno.temp)

      LD.temp<-LDinfo.temp$ldblk
      LD.temp<-LD.temp[match(snp1.intersect,snp1),match(snp1.intersect,snp1)]

      betavar<-as.vector(betavar)
      betavar2=betavar[match(snp1.intersect,snp1)]
      LDSM.temp<-LDSMman(LD.temp,betavar=betavar2)
      x.ldsc.temp<-LDSCmatrix(LD.temp,Anno.temp)[,-1]
      batch.used=1
    }else{snp1.intersect<-LD.temp<-LDSM.temp<-x.ldsc.temp<-batch.used<-0}
    return(list(batch.used,snp1.intersect,LD.temp,LDSM.temp,
                x.ldsc.temp))
  }

mkLDSM<-function(LD.path,anno.path,maf.path,annotation='all',snp=NULL,out,MAF,numCores=4){
  #register cores
  registerDoParallel(numCores)
  #Check reference files#####
  LD.files <- list.files(LD.path)
  if(any(grepl(x = LD.files, pattern = "ldblk_*"))){
    chrs <- LD.files[grep(x = LD.files, pattern = "ldblk_*")]
  }else{
    error.message <- "Ref pannel (LDfiles) input error"
    stop(error.message)
  }

  anno.files<- list.files(anno.path)
  if(any(grepl(x = anno.files, pattern = "*.annot.gz"))){
    annos <- anno.files[grep(x = anno.files, pattern = "*.annot.gz")]
    ldscs <- anno.files[grep(x = anno.files, pattern = "*.ldscore.gz")]
  }else{
    error.message <- "Ref pannel (ANNOfile) input error"
    stop(error.message)
  }

  maf.files<- list.files(maf.path)
  if(any(grepl(x = maf.files, pattern = "1000G.mac5eur.*"))){
    mafs <- maf.files[grep(x = maf.files, pattern = "1000G.mac5eur.*")]
  }else{
    error.message <- "Ref pannel (MAFfiles) input error"
    stop(error.message)
  }

  if(length(chrs)==length(annos)&length(chrs)==length(mafs)){
    chr.length<-length(chrs)
  }else{
    error.message <- "Ref pannel (number of chr) input error"
    stop(error.message)
  }


  #Preparation#####
  if(is.null(snp)==F){
    gwas.df<-read.table(snp,header = T)
    gwas.df<-na.omit(gwas.df)
    gwas.snp<-as.vector(gwas.df[,'SNP'])
    print(paste('load',length(gwas.snp),'SNPs from custom snplist',sep = ' '))
  }

  #Starts#####
  snp.list<-LDmatrix<-LDSM<-x.gls<-list()
  Atotal<-NA
  start_time <- Sys.time()
  for (chr in 1:chr.length) {

    chr.exact<-substr(chrs[chr],nchar('ldblk_1kg_')+1,nchar(chrs[chr])-nchar('.hdf5'))
    #read ref by chr
    LDinfo<-h5ls(file = paste(LD.path,chrs[chr], sep = '/'),recursive=F)

    file1<-gzfile(paste(anno.path,annos[chr], sep = '/'),'rt')
    file3<-gzfile(paste(maf.path,mafs[chr], sep = '/'),'rt')

    anno2<-read.table(file1, header = T,sep = '')
    maf<-read.table(file3,header = T,sep = '')

    close(file1)
    #close(file2)
    close(file3)

    #select annotation
    if(annotation[1]=='all'){used.a<-names(anno2)[-4:-1]
    }else{used.a<-annotation}

    #select chosen annotation
    anno<-anno2[,c(names(anno2)[1:4],used.a)]
    maf.snp<-as.vector(maf[maf$FRQ<(1-MAF)&maf$FRQ>MAF,'SNP'])
    #ldsc.snp<-ldsc$SNP
    #

    #generate by block
    if(is.null(snp)){gwas.snp=maf.snp}
    raw.result<-foreach (batch=1:nrow(LDinfo)) %dopar% refpannel.par(batch,LDinfo=LDinfo,LD.path=LD.path,chrs=chrs[chr],
                                                                     anno=anno,maf.snp=maf.snp,gwas.snp=gwas.snp,
                                                                     used.a=used.a)
    #record result
    for (batch in 1:length(raw.result)) {
      if(raw.result[[batch]][[1]]==1){
        batch.name<-LDinfo[batch,'name']
        snp.list[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[2]]
        LDmatrix[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[3]]
        LDSM[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[4]]
        x.gls[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[5]]

        Anno.temp<-anno[match(raw.result[[batch]][[2]],anno[,'SNP']),used.a]
        Anno.temp<-as.matrix(Anno.temp)
        Atotal<-rbind(Atotal,Anno.temp)
      }
    }
    print(paste0('LD Score Matrix of CHR ',chr,' has been complited'))
  }

  Atotal<-Atotal[-1,]
  snplist<-unlist(snp.list)
  M<-length(snplist)
  anno.total<-apply(Atotal, 2, sum)
  M.anno<-anno.total/M
  print(paste(M,'SNPs left for analysis',sep = ' '))
  print(Sys.time()-start_time)

  #output
  result<-list()
  result[['LDinfo']]<-LDinfo
  result[['LDmatrix']]<-LDmatrix
  result[['LDSM']]<-LDSM
  result[['snp.list']]<-snp.list
  result[['x.gls']]<-x.gls
  result[['anno.total']]<-anno.total
  result[['M.anno']]<-M.anno
  result[['Atotal']]<-Atotal
  result[['chrs']]<-chrs
  return(result)
  saveRDS(result,paste0(out,'.Rdata'))
}

gldsc<-function(panel,gwas,out,jackknife,intercept2,numCores=4){
  start_time=Sys.time()
  message('Analysis begin')
  ref.pannel<-readRDS(panel)
  message('load reference panel...')
  snplist.pan<-unlist(ref.pannel[['snp.list']])
  message(paste0('read ',length(snplist.pan),' SNPs from reference panel'))
  registerDoParallel(numCores)
  x.sum<-Reduce(sum,lapply(ref.pannel[['x.gls']], sum))


  #gls
  gwas.df<-read.table(gwas,header = T,fill = TRUE)
  message('read GWAS...')
  gwas.df<-na.omit(gwas.df)
  message(paste0('read ',nrow(gwas.df),' SNPs from GWAS'))
  snp.used=intersect(snplist.pan,gwas.df$SNP)
  message(paste0('After merging with reference panel ',length(snp.used),' SNPs remain'))
  #
  merged.loc<-match(snp.used,snplist.pan)

  gwas.df<-gwas.df[match(snp.used, gwas.df$SNP),]
  y.sum<-sum(gwas.df[,'Z']**2-1,na.rm = T)
  raw.Ntau=y.sum/x.sum

  gwasN<-median(gwas.df$N,na.rm = T)
  gls.left<-gls.right<-0
  raw.result.all<-list()
  for (chr in 1:length(ref.pannel[['chrs']])) {
    chr.exact<-substr(ref.pannel[['chrs']][chr],nchar('ldblk_1kg_')+1,nchar(ref.pannel[['chrs']][chr])-nchar('.hdf5'))
    block.names.all<-names(ref.pannel[['LDSM']])[grep(x = names(ref.pannel[['LDSM']]), pattern = paste0(chr.exact,'_'))]
    total.block<-length(block.names.all)

    if(total.block>0){
      raw.result<-foreach (i=1:total.block) %dopar% gls.left.right(i=i,x=ref.pannel[['x.gls']],y=gwas.df,
                                                                   block.names.all = block.names.all,
                                                                   LDmatrix = ref.pannel[['LDmatrix']],LDSM = ref.pannel[['LDSM']],
                                                                   snp.list = ref.pannel[['snp.list']],
                                                                   raw.Ntau = raw.Ntau,intercept = intercept2
      )
      for (batch in 1:length(raw.result)) {
        gls.left<-gls.left+raw.result[[batch]][[1]]
        gls.right<-gls.right+raw.result[[batch]][[2]]
      }
      raw.result.all[[chr.exact]]<-raw.result
    }
    message(paste0(chr.exact,' calculation complited'))
  }


  #first time result
  result.one<-gls.estimator(left = gls.left,right = gls.right,
                            A=ref.pannel[['Atotal']],M.anno = ref.pannel[['M.anno']],anno.total = ref.pannel[['anno.total']],N = gwasN,intercept = intercept2)
  message('point estimate completed')
  #result
  if(jackknife==F){
    message(Sys.time()-start_time)
    return(result.one)
  }else if(jackknife==T){
    message('start Jackknife estimation')
    #jackknife starts here
    jack.tau<-jack.h2<-jack.e<-jack.intersecp<-jack.stat<-NA
    jack.time=0
    for (chr in 1:length(ref.pannel[['chrs']])) {
      chr.exact<-substr(ref.pannel[['chrs']][chr],nchar('ldblk_1kg_')+1,nchar(ref.pannel[['chrs']][chr])-nchar('.hdf5'))
      raw.result.temp=raw.result.all[[chr.exact]]
      total.block2<-length(raw.result.temp)
      jack.result<-foreach (i=1:total.block2) %dopar% oneout(i,raw.result=raw.result.temp,left.total = gls.left,right.total = gls.right,
                                                             A=ref.pannel[['Atotal']],M.anno = ref.pannel[['M.anno']],
                                                             anno.total = ref.pannel[['anno.total']],N = gwasN,intercept = intercept2)
      for (i in 1:total.block2) {
        jack.tau<-rbind(jack.tau,as.vector(jack.result[[i]][[1]]))
        jack.h2<-rbind(jack.h2,as.vector(jack.result[[i]][[2]]))
        jack.e<-rbind(jack.e,as.vector(jack.result[[i]][[3]]))
        #jack.intersecp<-c(jack.intersecp,jack.result[[i]][[4]])
        #jack.stat<-rbind(jack.stat,as.vector(jack.result[[i]][[5]]))
        jack.time=jack.time+1
      }
    }
    jack.tau.var<-apply(jack.tau[-1,], 2, var)*(jack.time-1)^2/jack.time
    jack.h2.var<-apply(jack.h2[-1,], 2, var)*(jack.time-1)^2/jack.time
    jack.e.var<-apply(jack.e[-1,], 2, var)*(jack.time-1)^2/jack.time
    #jack.intersecp.var<-var(jack.intersecp[-1])*(jack.time-1)^2/jack.time
    #jack.stat.var<-apply(jack.stat[-1,], 2, var)*(jack.time-1)^2/jack.time


    result.one[['tau.var.jack']]=jack.tau.var
    result.one[['h2.var.jack']]=jack.h2.var
    result.one[['e.var.jack']]=jack.e.var
    #result.one[['inter.var.jack']]=jack.intersecp.var
    #result.one[['P']]=pnorm(result.one[['e.stat']]/sqrt(jack.stat.var),lower.tail=F)
    message(Sys.time()-start_time)
    return(result.one)
  }else{print('Jackknife should be TRUE of FALSE')}
}



