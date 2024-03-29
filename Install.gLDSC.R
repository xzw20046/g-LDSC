##### Install gLDSC required packages #####
#Rcpp
if(!require("Rcpp",character.only = TRUE)){
  try_error <- try(install.packages("Rcpp", 
                                    dependencies=TRUE), silent = TRUE)
  if(class(try_error) == "try-error"){
    try_error <- try(install.packages("Rcpp", 
                                      dependencies=TRUE, 
                                      lib = Sys.getenv("R_LIBS_USER")), silent = TRUE)
  }
}

#inline
if(!require("inline",character.only = TRUE)){
  try_error <- try(install.packages("inline", 
                                    dependencies=TRUE), silent = TRUE)
  if(class(try_error) == "try-error"){
    try_error <- try(install.packages("inline", 
                                      dependencies=TRUE, 
                                      lib = Sys.getenv("R_LIBS_USER")), silent = TRUE)
  }
}

#emulator
if(!require("emulator",character.only = TRUE)){
  try_error <- try(install.packages("emulator", 
                                    dependencies=TRUE), silent = TRUE)
  if(class(try_error) == "try-error"){
    try_error <- try(install.packages("emulator", 
                                      dependencies=TRUE, 
                                      lib = Sys.getenv("R_LIBS_USER")), silent = TRUE)
  }
}

#doParallel
if(!require("doParallel",character.only = TRUE)){
  try_error <- try(install.packages("doParallel", 
                                    dependencies=TRUE), silent = TRUE)
  if(class(try_error) == "try-error"){
    try_error <- try(install.packages("doParallel", 
                                      dependencies=TRUE, 
                                      lib = Sys.getenv("R_LIBS_USER")), silent = TRUE)
  }
}

#rhdf5
if(!require("rhdf5",character.only = TRUE)){
  try_error <- try(install.packages("rhdf5", 
                                    dependencies=TRUE), silent = TRUE)
  if(class(try_error) == "try-error"){
    try_error <- try(install.packages("rhdf5", 
                                      dependencies=TRUE, 
                                      lib = Sys.getenv("R_LIBS_USER")), silent = TRUE)
  }
}
