########################################
## Name: 	PackageBuilder.R
## Author:	JDB
## Date:	May 2023
## 
## Purpose: Build R package in Dropbox/Github/ROC, then push to Git
## Package: ROC (or ROCtools)
########################################

#### Helpful links for building packages
## https://kbroman.org/pkg_primer/pages/docs.html
## https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/

##### libraries
library(devtools)

#### Load package from Github
devtools::install_github("statevidence/ROC")
library(ROC)
sessionInfo()


#### Assemble Package

#### For JDB; set directory
getwd()
setwd("/Users/Jeffrey/Dropbox/Github/ROC")

## libraries
library(devtools)

## Assembly Steps

build()

install()

document()

## Test
library(ROC)

?fit.roc

## Check documentation for Compliance
check()

#### Extras

## packageVersion("ROC")
## remove.packages("ROC")

###
##
#