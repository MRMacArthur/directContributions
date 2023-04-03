library(readxl)
library(dplyr)
library(tidyverse)
library(pracma)

# Load raw data, in tidy format
dcData <- read_excel("Direction Contribution Data_Normalized.xlsx",
                     sheet = "Combined")

# Filter data to include only young 
# Replace "Young" with "Aged" to include only aged
# or bypass lines 12 and 13 to include all samples
dcData <- dcData %>%
  filter(Group == "Young")

# Define the squared residual function for optimization
residualFunction <- function(x){
  Y2 = as.matrix(M_ij) %*% as.matrix(x)
  return(sum((Y2-Y_ij)^2))
}

# Generate summary dataframe
dcData_summary <- dcData %>%
  group_by(Infusate, Analyte) %>%
  summarize(meanFrac = mean(Serum_Labeled_Frac),
            sdFrac = sd(Serum_Labeled_Frac)) %>%
  mutate(meanFrac = ifelse(Analyte == Infusate, 1, meanFrac))

# Split summary dataframe into separate mean and sd dataframes
meanMatrix <- dcData_summary %>%
  select(Infusate, Analyte, meanFrac) %>%
  pivot_wider(names_from = Analyte, values_from = meanFrac) %>%
  column_to_rownames('Infusate')

sdMatrix <- dcData_summary %>%
  select(Infusate, Analyte, sdFrac) %>%
  pivot_wider(names_from = Analyte, values_from = sdFrac) %>%
  column_to_rownames('Infusate')

# Definte # of monte carlo iterations
n_iter = 100

# Initiate empty lists to store final outputs
outputList_mean <- list()
outputList_sd <- list()

# Start for loop to perform optimization for each infusate
for(i in 1:length(sdMatrix)){
  
  # Define rows/columns to drop (current infusate)
  keepIndex <- setdiff(1:length(meanMatrix), i)
  
  # Mean and sd matrices with infusate removed
  M <- meanMatrix[keepIndex, keepIndex]
  dM <- sdMatrix[keepIndex, keepIndex]
  
  # Mean and sd vectors of infusate enrichments 
  # (with self-enrichment removed)
  Y <- meanMatrix[-i,i]
  dY <- sdMatrix[-i,i]
  
  # Initiate empty list for monte carlo outputs
  outputList_i <- list()
  
  # Do monte carlo
  for(j in 1:n_iter) {
    # Add normally distributed noise
    Y_ij = Y + rnorm(length(Y)) * dY
    
    M_ij = M + matrix(rnorm(
      nrow(M) * nrow(M)
    ), nrow = nrow(M)) * dM
    
    # Perform optimization with upper and lower bounds
    ub <- rep(1, length(Y))
    lb <- rep(0, length(Y))
    X0 <- solve(M_ij, Y_ij)
    
    X <- fmincon(x0 = X0,
                 fn = residualFunction,
                 lb = lb,
                 ub = ub)
    
    outputList_i[[j]] <- unlist(X[1])
    
  }
  
  # Summarize monte carlo results (mean and sd)
  outputList_mean[[i]] <- colMeans(do.call(rbind, outputList_i))
  
  outputList_sd[[i]] <- apply(
    do.call(rbind, outputList_i),
    2, sd
  )
  
}

# Add the diagonal zeros to the mean output
outputMean_frame <- do.call(rbind, outputList_mean)

outputMean_frame2 <- `diag<-`(matrix(ncol = ncol(outputMean_frame) + 1, 1,
                                     nrow = nrow(outputMean_frame)), 0)

outputMean_frame2[outputMean_frame2 == 1] <- outputMean_frame

# Add the diagonal zeros to the sd output
outputSd_frame <- do.call(rbind, outputList_sd)

outputSd_frame2 <- `diag<-`(matrix(ncol = ncol(outputSd_frame) + 1, 1,
                                     nrow = nrow(outputSd_frame)), 0)

outputSd_frame2[outputSd_frame2 == 1] <- outputSd_frame





