
# Load Required Libraries  ------------------------------------------------

library(GEOquery)
library(dplyr)
library(affy)
library(limma)
library(affyPLM)


# Create Function  --------------------------------------------------------

extract_data <- function (accession_no, workdir){

                workdir = workdir
                accession_no = accession_no
                # Fetch expression data
                # download rawfile to the directory
                getGEOSuppFiles(accession_no)
                # To open the tar file, access the folder with
                setwd(paste0(workdir,"/",accession_no))
                # Open the tar file
                untar(paste0(accession_no,"_RAW.tar"))
                # #Get a vector list of zippped GSM

                Sys.sleep(120)
                setwd(workdir)
                file_list<-list.files(accession_no, pattern  = ".CEL|.cel")
                setwd(paste0(workdir,"/",accession_no))
                # #Get a vector list of unzipped the CEL files
                sapply(file_list, gunzip)
                setwd(workdir)
                # Get a list vector of the files from the folder
                cel.list<-list.files(accession_no, pattern  = ".CEL|.cel",
                                     full.names = F)

                setwd(paste0(workdir,"/",accession_no))

                # Read the CEl files into R, normalize with RMA and get the
                #expression data set

                exp.set=justRMA()

                expr.data<-exprs(exp.set)

                #Fetch phenotype data

                datadir<-paste0(workdir,"/",accession_no)
                pdata.class <- getGEO(accession_no, GSEMatrix = T,AnnotGPL = T,
                                      destdir = datadir)
                pdata<- pdata.class[[1]]  %>%
                  pData()

                extract_list <- list(expr.data, pdata)
                return(extract_list)

    }
