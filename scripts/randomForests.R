# This script takes an expression matrix and outputs a list of "important" genes
# for some p-value cutoff according to a permutation-like test.
library(randomForest)   # CRAN

### Parameters ###
# expression.matrix     - expression matrix (or data frame) in the form of samplesXgenes (rowsXcolumns)
# response              - vector containing the response (ex. Tissue Type) for the samples in expression.matrix
# p.value.cutoff        - p-value for which you want the significantly "important" genes to be selected
# num.rf.trees          - the number of decision trees to be used in the random forest (500 is typical, anything above 10,000 probably doesn't add new info)
# random.seed           - seed for reproducible results

randomForestImpGenes <- function(expression.matrix, response, p.value.cutoff, num.rf.trees, random.seed){
  set.seed(random.seed)
  
  # building the model
  rf.mod = randomForest(x = expression.matrix, y = factor(response), ntree = num.rf.trees, importance = TRUE, scale = FALSE)
  
  # obtaining importance measure (type = 1 returns "MeanDecreaseAccuracy" importance)
  impt.rf = importance(rf.mod, type = 1, scale = TRUE)
  
  # finding the number of genes that have importance measure greater than the 1-{alpha} quantile of a standard normal distribution
  num.sig.impt.genes = sum(impt.rf[, "MeanDecreaseAccuracy"] >= qnorm(1-p.value.cutoff, mean = 0, sd = 1))
  
  # obtaining the list of significantly important genes
  if(num.sig.impt.genes > 1){sig.genes = impt.rf[order(impt.rf[, "MeanDecreaseAccuracy"], decreasing = TRUE), ][1:num.sig.impt.genes]}
  else if(num.sig.impt.genes == 1){sig.genes = impt.rf[order(impt.rf[, "MeanDecreaseAccuracy"], decreasing = TRUE), ][1]}
  
  return(sig.genes)
}
