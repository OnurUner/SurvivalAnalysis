source("make.R")
make()
library(R.matlab)
scores <- readMat('C:\\r workspace\\MultiSconES\\data\\scores.mat')
dataset <- readMat('C:\\r workspace\\MultiSconES\\data\\breast_dataset.mat')
ppi_graph <- graph_from_adjacency_matrix(dataset$ppi.network)
X <- dataset$data
Y <- dataset$labels
C <- c()
task.names <- dataset$y.columns

for (name in task.names){
  task <- gsub("_", "\\.", name)
  for (i in 1:length(scores)){
    if(grepl(task, names(scores)[i])){
      print(names(scores)[i])
      C <- c(C, c(scores[[i]]))
    }
  }
} 

res <- mscones(g = ppi_graph, X = X, Y = Y, C = C)
features <- c()
for (feature_num in res[[1]]){
  features <- c(features, feature_num)
}
write.csv(features, file = "C:\\r workspace\\MultiSconES\\data\\selected_features.csv",row.names=FALSE, col.names = FALSE)
