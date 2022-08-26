## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  # Install current version
#  install.packages("chemodiv")
#  
#  # Install developmental version
#  install.packages("devtools") # Install devtools if not already installed
#  library("devtools")
#  install_github("hpetren/chemodiv")

## -----------------------------------------------------------------------------
# Load chemodiv
library(chemodiv)

## -----------------------------------------------------------------------------
data("alpinaSampData")

head(alpinaSampData)[,1:5]

## -----------------------------------------------------------------------------
data("alpinaCompData")

head(alpinaCompData)

## -----------------------------------------------------------------------------
chemoDivCheck(compoundData = alpinaCompData, sampleData = alpinaSampData)

## -----------------------------------------------------------------------------
data("alpinaPopData")

table(alpinaPopData)

## ---- eval = FALSE------------------------------------------------------------
#  alpinaNPC <- NPCTable(compoundData = alpinaCompData)
#  
#  alpinaNPC[1,] # Classification of the first compound in dataset

## ---- echo = FALSE------------------------------------------------------------
data("alpinaNPCTable")
alpinaNPC <- alpinaNPCTable
rm(alpinaNPCTable)
alpinaNPC[1,]

## ---- eval = FALSE------------------------------------------------------------
#  alpinaCompDis <- compDis(compoundData = alpinaCompData,
#                           type = "PubChemFingerprint")
#  
#  alpinaCompDis$fingerDisMat[1:4, 1:4] # Part of compound dissimilarity matrix

## ---- echo = FALSE, message = FALSE-------------------------------------------
data("alpinaCompDis")
alpinaCompDisMat <- alpinaCompDis
rm(alpinaCompDis)
alpinaCompDis <- list()
alpinaCompDis[["fingerDisMat"]] <- alpinaCompDisMat
alpinaCompDis$fingerDisMat[1:4, 1:4]

## -----------------------------------------------------------------------------
alpinaDiv <- calcDiv(sampleData = alpinaSampData, 
                     compDisMat = alpinaCompDis$fingerDisMat,
                     type = "FuncHillDiv",
                     q = 1)

head(alpinaDiv)

## -----------------------------------------------------------------------------
alpinaDivProf <- calcDivProf(sampleData = alpinaSampData,
                             compDisMat = alpinaCompDis$fingerDisMat,
                             type = "FuncHillDiv")

head(alpinaDivProf$divProf)[,1:5] # Part of the diversity profile data frame

## -----------------------------------------------------------------------------
alpinaBetaDiv <- calcBetaDiv(sampleData = alpinaSampData,
                             compDisMat = alpinaCompDis$fingerDisMat,
                             type = "FuncHillDiv")

alpinaBetaDiv

## ---- message = FALSE---------------------------------------------------------
alpinaSampDis <- sampDis(sampleData = alpinaSampData,
                         compDisMat = alpinaCompDis$fingerDisMat,
                         type = "GenUniFrac")

alpinaSampDis$GenUniFrac[1:4, 1:4] # Part of sample dissimilarity matrix

## ---- message = FALSE---------------------------------------------------------
alpinaNetwork <- molNet(compDisMat = alpinaCompDis$fingerDisMat,
                        npcTable = alpinaNPC,
                        cutOff = 0.75)

summary(alpinaNetwork)

## ---- fig.width = 12, fig.height = 8, out.width = "95%"-----------------------
molNetPlot(sampleData = alpinaSampData,
           networkObject = alpinaNetwork$networkObject,
           npcTable = alpinaNPC,
           plotNames = TRUE)

## ---- fig.width = 12, fig.height = 8, out.width = "95%"-----------------------
chemoDivPlot(compDisMat = alpinaCompDis$fingerDisMat,
             divData = alpinaDiv,
             divProfData = alpinaDivProf,
             sampDisMat = alpinaSampDis$GenUniFrac,
             groupData = alpinaPopData)

## ---- eval = FALSE------------------------------------------------------------
#  quickChemoDiv(compoundData = alpinaCompData,
#                sampleData = alpinaSampData,
#                groupData = alpinaPopData,
#                outputType = "plots") # Not run

