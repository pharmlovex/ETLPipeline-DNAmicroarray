

# Load Required Libraries  ------------------------------------------------

library(GEOquery)
library(dplyr)
library(affy)
library(limma)
library(affyPLM)


# Extract Assay and Phenotype data -------------------------------------------

extract_data <- function (accession_no, workdir){

                workdir = workdir
                accession_no = accession_no
                setwd(workdir)
                # Fetch rawfile to the directory
                getGEOSuppFiles(accession_no)
                # To open the tar file, access the folder with
                setwd(paste0(workdir,"/",accession_no))
                # Open the tar file
                untar(paste0(accession_no,"_RAW.tar"))
                setwd(workdir)
                # Get a vector list of zippped GSM
                file_list<-list.files(accession_no, pattern  = ".CEL|.cel")
                setwd(paste0(workdir,"/",accession_no))
                # Get a vector list of unzipped the CEL files
                sapply(file_list, gunzip)

                # Read the CEl files into R, normalize with RMA and get the
                #expression data set

                exp.set=justRMA()

                expr.data<-exprs(exp.set)
                expr.df = as.data.frame(expr.data)

                #Fetch phenotype data

                datadir<-paste0(workdir,"/",accession_no)
                pdata.class <- getGEO(accession_no, GSEMatrix = T,AnnotGPL = T,
                                      destdir = datadir)
                pdata<- pdata.class[[1]]  %>%
                  pData()
                # Gather the results together
                extract_list <- list(expr.df, pdata)
                names(extract_list) = c("assaydata", "phenotype")
                setwd(workdir)
                return(extract_list)

    }


# Extract Assay  ----------------------------------------------------------


extract.data <- extract_data(accession_no = "GSE158643",
                         workdir = "C:/Users/dell/Downloads/DNAmicroarray")


# Transform  --------------------------------------------------------------

df = extract.data$assaydata
transform_df <- transformMatrixCol(expMat = df,
                                   Sep = "_",
                                   nSep = 6,
                                   pos = 1)
