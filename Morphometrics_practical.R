## -------------------------------------- ##
## ----    INITIALIZE ENVIRONMENT   ----- ##
## -------------------------------------- ##

#\\ Clean environment
rm(list=ls())

#\\ Load needed packages 
library(maptools)
library(geomorph)

#\\ Set environment
#
# Option 1)
# You can either type in the following commented line
# ("wd <- "~/Workshop/") the absolute path
# to the directory you have saved the data from the GitHub repository
# (do not forget the last "/" !) and uncomment this line and 
# the following to set your environment
#
# wd <- "~/Workshop/"
# setwd(wd)
#
# Option 2)
# Or you can leave the previous lines commented and find your working
# directory with the next commands. If you are not going to use this
# option and prefer option 1), do not run the following commands or 
# just uncomment them
filename <- "Morphometrics_practical.R"
filepath <- file.choose()  # browse and select your_file.R in the window
wd = substr(filepath, 1, nchar(filepath)-nchar(filename))
# I run the next command in case it has been run in Windows
wd <- gsub(pattern="[\\]", replace="/", x=wd)
setwd(wd)

## -------------------------------------- ##
## ************************************** ##


## -------------------------------- ##
## ----    CREATE VARIABLES    ---- ##
## -------------------------------- ##

#\\ Load data set - 9 species x (3 info columns + 144 coordinates)
df <- read.table(file=paste(wd,"data/Triturus_and_Calotriton_lmk_reduced.csv", sep=""), header=T,sep=",", dec=".", stringsAsFactors=F)

# The first 3 columns do not contain the coordinates 
# Get from the 4th to the nth column, where the data is 
mm <- df[,4:(dim(df)[2])] # 9sp x 144coords (144/3=48 lmk)
# Get name of specimens as row names
rownames(mm) <- df[,1]

#\\ Create an empty array to store the coordinates in the format 
#\\ p x k x n (num.coords x coor3D x ns)
num.coords <- dim(mm)[2]
coor3D <- 3
ns <- 9
ma <- array(dim=c(num.coords/coor3D, coor3D, ns)) # 48 lmks, 3D, 9 specimens

#\\ Select the positions from 1 to 144 that correspond to 
#\\ x coords, y coords, and z coords
xi <- seq(from=1, to=num.coords, by=coor3D); yi <- xi+1; zi <- yi+1

## -------------------------------------- ##
## ************************************** ##


## -------------------------------- ##
## ----    PREPARE DATA SET    ---- ##
## -------------------------------- ##

#\\ Fill in array p x k x n
#\\ p landmarks, k dimensions, n specimens
#
# Note that we need to first unlist the objects
# "mm[i,xi]", "mm[i,yi]", and "mm[i,zi]"
# because each coordinate belongs to one list of vectors  
# from the data.frame "mm"
# E.g. * mm$astl is a column vector
#      * mm[1,1] is coordinate x from column one 
#         and is equal to mm$astl[1]
# Uncomment to the following command to prove it
# all.equal(mm$astl[1], mm[1,1])
#
# When we get only coordinates "x", for instance,
# we need to unlist the 102 column vectors (mm$astl, 
# mm$astl.1, ...) so we can extract one coordinate 
# "x" per species to put it in the array "ma"
for (i in 1:ns) {
  ma[,1,i] <- unlist(mm[i,xi])
  ma[,2,i] <- unlist(mm[i,yi])
  ma[,3,i] <- unlist(mm[i,zi])
}

## -------------------------------------- ##
## ************************************** ##


## ------------------------------- ##
## ----  PROCRUSTES ANALYSIS  ---- ##
## ------------------------------- ##

#\\ Get procrustes analysis done
ma.paln <- geomorph::gpagen(ma)

#\\ Convert the array form with the coordinates
#\\ previously superimposed in "ma.paln" into 
#\\ a data.frame format with dimensions 
#\\ p x n (sp x coords)
faln <- matrix(0, nrow=ns, ncol=num.coords)
for (i in 1:ns) {
  faln[i,xi] <- ma.paln$coords[,1,i]
  faln[i,yi] <- ma.paln$coords[,2,i]
  faln[i,zi] <- ma.paln$coords[,3,i]
}

#\\ Uncomment the following if you want to save the resulting
#\\ coordinates in a csv file
#
# Save only coordinates
# write.table(faln, file="PA/Triturus_and_Calotriton_after_PA_clean.csv",
#             row.names=F, col.names=F, quote=F, sep=",")
# Save coordinates + 3 columns of info about specimens
# faln.rn <- cbind(df[,1:3], faln)
# colnames(faln.rn) <- colnames(df)
# write.table(faln.rn, file="PA/Triturus_and_Calotriton_after_PA.csv",
#             row.names=F, quote=F, sep=",")

## -------------------------------------- ##
## ************************************** ##



## ----------------- ##
## ----  PLOTS  ---- ##
## ----------------- ##

#\\ Plot the data set before PA
#
# Uncomment the two commented lines if you want to save
# this plot
#
# pdf("plots/Raw_data_2D.pdf", paper="a4", width=15, height=15)
plotAllSpecimens(ma[,1:2,])
# dev.off()

#\\ Plot the data set after PA
# Uncomment the two commented lines if you want to save
# this plot
#
# pdf("plots/Data_after_PA_2D.pdf", paper="a4", width=15, height=15)
plotAllSpecimens(ma.paln$coords[,1:2,]) 
# dev.off()

#\\ Plot a 2D procrustes plot with lines x=0 and y=0 and mean shape
# 
# Here we have decided to add to the plot the mean of every group of 
# landmarks (one landmark per specimen) that fall to the same point 
# in the plot forming the skull shape --> "points()"
# 
# If you want to plot the labels of the landmarks of the 1st 
# and the 2nd specimen, just uncomment the command that calls the function
# "maptools::pointLabel". You can change the numbers if you want to see 
# other landmarks or you can also leave this uncommented
#
# NOTE: Uncomment the lines with "pdf()" and "dev.off()" if you want to save
#       this plot
#
#pdf("plots/Data_after_PA_2D_limits.pdf", paper="a4", width=15, height=15)
msh <- mshape(ma.paln$coords[,1:2,])
xl <- range(faln[,xi]); yl <- range(faln[,yi])
matplot(faln[1:ns,xi], faln[1:ns,yi], pch='+', cex=.5, asp=1, col="darkgrey", xlim=xl, ylim=yl, xlab="X", ylab="Y", las=1)
points(msh[,1], msh[,2],cex=1, pch=19,col="black")
# maptools::pointLabel(faln[1:2,xi], faln[1:2,yi], labels = paste(1:dim(msh)[1]),cex=1,col=c("black", "red"),
#                      method = c("SANN"),
#                      allowSmallOverlap = FALSE,
#                      trace = FALSE,
#                      doPlot = TRUE
# )
abline(h=0, v=0, lty=2, lwd=.25)
#dev.off()

## ----------------- ##
## ***************** ##