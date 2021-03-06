
# This example obtains the LSE of the model/method described on
# Maitra & Riddles (2010)

# We load the required packages
library(SynMRI)
library(oro.nifti)

#Let's first read the data
dat <- readNIfTI("ZHRTS1.nii")
data <- oro.nifti::img_data(dat)

# we save the dimension and total number of images
MM <- 12
N1 <- 256
N2 <- 256
N3 <- 20

# here we save the indexes of the sub-image to work with
n1_ini <- 65
n1_fin <- n1_ini + 127
n2_ini <- 65
n2_fin <- n2_ini + 127
n3_ini <- 8
n3_fin <- n3_ini + 5

# we compute the new dimensions
n1 <- n1_fin - n1_ini + 1
n2 <- n2_fin - n2_ini + 1
n3 <- n3_fin - n3_ini + 1

# let's save the data on a list of 3-dimensional arrays
index <- 1:MM
R_all <- lapply(index, function(i, data){return(data[n1_ini:n1_fin,n2_ini:n2_fin,n3_ini:n3_fin,i])}, data)

# let's save the values of the design parameters
# TE in millisecons
TE_all <- c(10, 15, 20, 10, 30, 40, 10, 40, 80, 10, 60, 100)
# TR in seconds
TR_all <- c(600, 600, 600, 1000, 1000, 1000, 2000, 2000, 2000, 3000, 3000, 3000)

# Here let's save the subset to be used
# 1, 10 and 12 are the three images we focus on
subset <- c(1, 10, 12)

R <- R_all[subset]
TE <- TE_all[subset]
TR <- TR_all[subset]

# And here we fit parameters with LSE
fit <- fit_model_SynMRI(R, TE, TR)

# Now we will use the fitted values to make predictions at the design parameter
# values we have available
R_all_pred <- predict_SynMRI(fit, TE_all, TR_all)

# Each element of the list is a syntetic MRI

