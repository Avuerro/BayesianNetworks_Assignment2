library('dagitty')
library( 'bayesianNetworks' )
library('lavaan')
library('bnlearn')
library('graphviz')
library('NetworkDistance')
#library('igraph')
library('glue')
# Import Data
d <- read.csv("explored_forestfires.csv", colClasses=c("numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

# Show head of data
head(d)

wl <- data.frame(from=c("ISI","FFMC","DC","DMC") ,to=c("area","area","area","area"))
bl <- data.frame(from=c("area","month","FFMC","FFMC", "DMC", "DMC", "DMC","temp","DC"), to=c("month","area","temp","rain","month","temp","rain","month","month"))
structurez <- tabu(d)

structurez
graphviz.plot(structurez, shape = "circle", render = TRUE)

original_model_string <- 'dag{
  bb="0,0,1,1"
  area [outcome,pos="0.5,0.8"]
  ISI [pos="0.2,0.6"]
  FFMC [pos="0.4,0.6"]
  DMC [pos="0.6,0.6"]
  DC [pos="0.8,0.6"]
  month [pos="0.7,0.25"]
  RH [pos="0.4,0.4"]
  rain [pos="0.8,0.4"]
  temp [pos="0.6,0.4"]
  wind [pos="0.2,0.4"]
  DC -> area
  DMC -> area
  FFMC -> area
  FFMC -> ISI
  ISI -> area
  month -> rain
  month -> temp
  month -> DC
  RH -> DMC
  RH -> FFMC
  rain -> DC
  rain -> DMC
  rain -> FFMC
  temp -> DC
  temp -> DMC
  temp -> FFMC
  temp -> RH
  wind -> FFMC
  wind -> ISI
  DC -> DMC
}'

M <- lavCor(d)
best_alpha <- 0.1
load_original = function(model_string){
  g <- dagitty(model_string)
  original <- model2network(toString(g,"bnlearn"))
  original
}

original <- load_original(original_model_string)

"
  Converts the strings to dag model, bnlearn objects and adjacency matrices
"
prepare_data <- function(original_graph_string,param_list, TABU=TRUE) {
  
  n <-length(param_list) + 1 # we add original graph as well...
  print(glue("the number of parameters {n}"))
  graph_strings <- vector(mode="list",length=n) 
  graph_dag <- vector(mode="list",length=n) 
  graph_objects <- vector(mode="list",length=n) 
  adjacency_matrices <- vector(mode="list", length=n)
  algo_name <- ""
  param_name <- ""
  
  graph_dag[[1]] <-dagitty(original_graph_string)
  graph_objects[[1]] <- original
  adjacency_matrices[[1]] <- amat(original)
  for(i in 2:n){
    param_val = param_list[[i-1]]
    
    if(TABU) {
      graph_objects[[i]] = tabu(d, max.iter = param_val)
      algo_name <- "tabu"
      param_name <- "max.tabu"
    } else{
      graph_objects[[i]] = si.hiton.pc(d, alpha=param_val, undirected=FALSE)
      algo_name <-"s_hiton_pc" 
      param_name <- "alpha"
    }
    
    name = paste(algo_name,param_val,".png",sep="")
    png(file=name)
    plot(graph_objects[[i]])
    title(param_val)
    dev.off()
    rows = nrow(graph_objects[[i]]$arcs)
    graph_strings[[i]] = paste("dag{")
    for (row in 1:rows) {
      graph_strings[[i]] = paste(graph_strings[[i]], graph_objects[[i]]$arcs[row,1], " -> ", graph_objects[[i]]$arcs[row,2])
    }
    
    graph_strings[[i]] = paste(graph_strings[[i]],"}")
    graph_dag[[i]] <- dagitty(graph_strings[[i]])
    
    #print(paste("For param :",max_iter_value))
    print(glue("Algorithm {algo_name} and parameter name {param_name} with value {param_val} "))
    adjacency_matrices[[i]] <- amat(graph_objects[[i]])
    graph_data <- list("dags"=graph_dag, "objects"=graph_objects, "adjacency_matrices"=adjacency_matrices)
  }
  return(graph_data)
  
}

max_iter_values <- seq(2,20,2)
tabu_graph_data <- prepare_data(original_model_string, max_iter_values )

alphas <- list(0.005, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99)
si_hiton_graph_data <- prepare_data(original_model_string, alphas, TABU=FALSE)

#results <- matrix(nrow=n,ncol=10)
"
  Calculates the betweenness, degree for the nodes in a graph
  as well as the hamming distance to the original graph
  returns a list of metric objects
"
calculate_metrics <- function(data) {
  betweenness <- nd.centrality(data[['adjacency_matrices']],mode="Between",out.dist=FALSE,directed = TRUE)$features
  degree<- nd.centrality(data[['adjacency_matrices']],mode="Degree",out.dist=FALSE,directed = TRUE)$features
  nr_of_objects <- length(data[['objects']])
  hammingd <- vector(mode="list",length=(nr_of_objects-1) )
  
  results <- matrix(nrow=nr_of_objects-1,ncol=6)
  #eigenvector_centrality <- vector(mode="list", length=nr_of_objects)
  
  
  original_object = data[['objects']][[1]]
  #original_graphmodel = graph_from_adjacency_matrix(data[['adjacency_matrices']][[1]])
  #eigenvector_centrality[[1]] <- eigen_centrality(original_graphmodel, directed=TRUE)$vector
  for(i in 2:nr_of_objects) {
    hammingd[[i-1]] <- hamming(original_object,data[['objects']][[i]])
    dcompare <- compare(data[['objects']][[i]],original_object)
    dlocaltests <- localTests(data[['dags']][[i]], sample.cov = M, sample.nobs=nrow(d) )
    
    #graphmodel = graph_from_adjacency_matrix(data[['adjacency_matrices']][[i]])
    #eigenvector_centrality[[i]]<- eigen_centrality(graphmodel, directed=TRUE)$vector
    results[[i-1,1]] <- dcompare[[1]]
    results[[i-1,2]] <- dcompare[[2]]
    results[[i-1,3]] <- dcompare[[3]]
    results[[i-1,4]] <- sum(dlocaltests[,][1] >0.1)
    results[[i-1,5]] <- sum(dlocaltests[,][1] >0.2)
    results[[i-1,6]] <- sum(dlocaltests[,][1] >0.3)
    
  }
  
  return_data <- list("betweenness"=betweenness,"degree"=degree,"hammingd"=hammingd,"other_results"=results) #, "eigen_centrality"=eigenvector_centrality)
  
  return(return_data)
  
}


tabu_metrics <- calculate_metrics(tabu_graph_data)
hiton_metrics <- calculate_metrics(si_hiton_graph_data)


"
  Plots the evaluation metrics for the given data
"
plot_graphs <- function(data,algo_name,param_values,param_name, plot_data) {
  hamming_data <- data[['hammingd']]
  n=length(param_values)
  sum_betweenness <- vector(mode="list",length=n)
  mean_betweenness <- vector(mode="list",length=n)
  max_degree <- vector(mode="list",length=n)
  mean_degree <- vector(mode="list",length=n)
  results <- data[['other_results']]
  
  for(i in 1:n) {
    max_degree[[i]] <- max(data[['degree']][i+1,])
    mean_degree[[i]] <- mean(data[['degree']][i+1,])
    sum_betweenness[[i]] <- sum(data[['betweenness']][i+1,])
    mean_betweenness[[i]] <- mean(data[['betweenness']][i+1,])
  }
  
  png(file=glue("{algo_name}_Res_Hamming.png"), width=500, height=350)
  plot(param_values,hamming_data, main = 'Hamming distance to the original network',ylab = 'Hamming distance',xlab = glue('{param_name}'),cex = 2,type="b")
  dev.off()
  
  png(file=glue("{algo_name}_Res_Betweenness.png"), width=500, height=350)
  plot(param_values,sum_betweenness, main = 'Betweenness',ylim = c(0, max(unlist(sum_betweenness))),ylab = 'Value',xlab = glue('{param_name}'),cex = 2,type="b")
  lines(param_values,mean_betweenness,type="b",col="blue")
  abline(h = sum(data[['betweenness']][1,]),col="red")
  legend("topleft", legend=c("sum Betweenness", "Mean Betweenness","Original Betweenness"),col=c("black","blue","red"), cex=0.8, lty=1,bg = "transparent")
  dev.off()
  
  png(file=glue("{algo_name}_Res_Degree.png"), width=500, height=350)
  plot(param_values,max_degree, main = 'Graph Degree',ylim = c(min(data[['degree']]), max(data[['degree']])),ylab = 'nº',xlab = glue('{param_name}'),cex = 2,type="b")
  lines(param_values,mean_degree,type="b",col="blue")
  abline(h = max(data[['degree']][1,]),col="red")
  legend("bottomright", legend=c("sum Degree", "Mean Degree","Original Degree"),col=c("black","blue","red"), cex=0.8, lty=1,bg = "transparent")
  dev.off()
  
  png(file=glue("{algo_name}_Res_TP.png"), width=500, height=350)
  plot(param_values,results[,1], main = 'number of arcs in the original network\n also present in the new network\n(True Positives)',xlab = glue('{param_name}'),ylab = 'nº of true positives', type="b",cex = 2)
  lines(param_values,results[,1])
  dev.off()
  
  png(file=glue("{algo_name}_Res_FP.png"), width=500, height=350)
  plot(param_values,results[,2], main = 'number of arcs in the original network\n not present in the new network\n(False Positives)',xlab = glue('{param_name}'),ylab = 'nº of false positives', type="b",cex = 2)
  lines(param_values,results[,2])
  dev.off()
  
  png(file=glue("{algo_name}_Res_FN.png"), width=500, height=350)
  plot(param_values,results[,3], main = 'number of new arcs in the new network\n(False Negatives)',ylab = 'nº of false negatives',xlab = glue('{param_name}') ,type="b",cex = 2)
  lines(param_values,results[,3])
  dev.off()
  
  png(file=glue("{algo_name}_Correlations.png"), width=500, height=350)
  plot(param_values,results[,4], main = 'Missed correlations',xlab = glue('{param_name}'),ylim = c(0, max(unlist(results[,4]))),ylab = 'nº',type="h",lwd = 2)
  lines(param_values,results[,4],col="blue",lwd = 2)
  lines(param_values,results[,5],col="red",lwd = 2)
  lines(param_values,results[,6],col="green",lwd = 2)
  legend("topright", title="Higher than:",legend=c("0.1", "0.2","0.3"),col=c("blue","red","green"), cex=0.8, lty=1,bg = "transparent")
  dev.off()
  
  for(i in 1:n) {
    
    name = paste(glue("{algo_name}_{param_values[[i]]}.png"),sep="")
    png(file=name)
    plot(plot_data[[i+1]])
    
    title(param_values[[ i]])
    dev.off()
  }
  
}

plot_graphs(hiton_metrics,'hiton',param_values=alphas,'alphas', si_hiton_graph_data[['objects']])
plot_graphs(tabu_metrics,'tabu',param_values=seq(2,20,2),'max iterations', tabu_graph_data[['objects']])

" 
  Code for evaluating a graph generated by the SI-HITON-PC algorithm
  and Tabu Search algorithm
    best tabu: 10 iterations(5), hamming dist. = 14
    best hiton: alpha=0.1(3), hamming dist. = 17
"
tabu_graph_data[['objects']][[1]] = si_hiton_graph_data[['objects']][[4]]
si_hiton_graph_data[['objects']][[1]] = tabu_graph_data[['objects']][[6]]

tabu_metrics <- calculate_metrics(tabu_graph_data)
hiton_metrics <- calculate_metrics(si_hiton_graph_data)
hiton_metrics[["degree"]]
length(alphas)
png(file=glue("Hamming_tabu_hiton.png"), width=500, height=350)
plot(seq(2,20,2),tabu_metrics[['hammingd']], main = 'Hamming distance of tabu-generated graphs \n to best hiton (alpha=0.1)',ylab = 'Hamming distance',xlab = glue('max iterations'),cex = 2,type="b")
dev.off()
png(file=glue("Hamming_hiton_tabu.png"), width=500, height=350)
plot(alphas,hiton_metrics[['hammingd']], main = 'Hamming distance of hiton-generated graphs \n to best tabu (max. iterations=10)',ylab = 'Hamming distance',xlab = glue('alphas'),cex = 2,type="b")
dev.off()

