## Hardness measures adopted in Smith et al. 
## An Instance Level Analysis of Data Complexity - Machine Learning

##  kDN_measure(dataset,k) --- k-Disagreeing Neighbors (kDN)

##  TD_measure(dataset) --- Tree Depth (TD)
##   *** TD --- pruned tree (TD P) 
##   *** TU --- unpruned tree (TD U)
     
##  CL_measure(dataset)
##  *** CL --- Class Likelihood (CL)
##  *** CLD --- Class Likelihood Difference (CLD)

##  MV_measure(dataset) 
##  *** MV --- Minority Value (MV)
##  *** CB --- Class Balance (CB)

##  DS_measure(dataset) --- Disjunct Size (DS)

##  DCP_measure(dataset) --- Disjunct Class Percentage (DCP)

## HINTS
## TD: high values, high hardness
## kDN: high values, high hardness
## CL, CLD: low values, high hardness
## MV: high values, high hardness  
## CB: low values, high hardness
## DS: low values, high hardness
## DSP: low values, high hardness

########################################################
########################################################

library(rpart)
library(partykit)
library(stringr)


# Parameters for PDF generation
PDFEPS <- 1 # 0 None, 1 PDF, 2 EPS
PDFheight= 7 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
PDFwidth= 7 # 7 by default

# This function is used to generate PDFs or EPSs for the plots
openPDFEPS <- function(file, height= PDFheight, width= PDFwidth) {
  if (PDFEPS == 1) {
    pdf(paste(file, ".pdf", sep=""), width, height)
  } else if (PDFEPS == 2) {
    postscript(paste(file, ".eps", sep=""), width, height, horizontal=FALSE)
  }
}

cleanDSfunction <- function(path = "./_Toy_/_BinaryDS/ltm-3PL/"){
  load(paste(path,"irt_parameters_mc.RData",sep="")) #item_params
  datasets <- list.files(path, pattern = "*csv$")
  
  for (i in datasets){
    index <- which(datasets == i)
    dataset <- read.csv(paste(path,i,sep=""))
    df <- as.data.frame(item_param[[index]])
    whichNeg <- which(df$Dscrmn <0)
    dataset <- dataset[-whichNeg,]
    nameDS <- datasets[index]
    nameDS <- strsplit(nameDS,"[.]")[[1]][1] #keep just the name
    
    write.csv(dataset,paste(path,nameDS,"_CleanDscrmn.csv",sep=""))
  }
}


run <- function(path = "./_Toy_/_BinaryDS/ltm-3PL/", cleanDS = F){
  
  datasets <- list.files(path, pattern = "*csv$")
  #print(datasets)
  inshardList <- list()
  load(paste(path,"irt_parameters_mc.RData",sep="")) #item_params
  load(paste(path,"results_responses_mc.RData",sep=""))
  
  
  for (i in datasets){
    
    index <- which(datasets == i)
    dataset <- read.csv(paste(path,i,sep=""))
    #print(head(dataset))
   
    print(paste("Dataset: ", i)) 
    print("Calculating:  IRT Parameter")
    
    df <- as.data.frame(item_param[[index]])
   
    print("Calculating IH")
    
    df$IH <- 1- rowMeans(results[[index]])
    
    
     
    if (cleanDS){
      whichNeg <- which(df$Dscrmn <0)
      df <- df[-whichNeg,]
      dataset <- dataset[-whichNeg,]
    }

    if (length(unique(dataset$Class)) > 1){
      print("Calculating:  Instance Hardness Measures")
      df$kDN <- kDN_measure(dataset,10) ## kDN: high values, high hardness
      df$IDS <- as.vector(DS_measure(dataset)) ## DS: low values, high hardness
      df$IDCP <- as.vector(DSP_measure(dataset)) ## DSP: low values, high hardness
      
      TD <- TD_measure(dataset) ## TD: high values, high hardness
      df$TD_P <- TD$TD
      df$TD_U <- TD$TU
      
      CLAll <- CL_measure(dataset) ## CL, CLD: low values, high hardness
      df$ICL <- CLAll$CL
      df$ICLD <- CLAll$CLD
      
      MVAll <- MV_measure(dataset)
      df$MV <- as.vector(MVAll$MV) ## MV: high values, high hardness
      df$ICB <- as.vector(MVAll$CB) ## CB: low values, high hardness
      
      
      
      df$IDS <- -(df$IDS)
      df$IDCP <- -(df$IDCP)
      df$ICL <- -(df$ICL)
      df$ICLD <- -(df$ICLD)
      df$ICB <- -(df$ICB)
      
      inshardList[[index]] <- df
    }else{
      inshardList[[index]] <- df
    }
    
    
  }
 
  if(cleanDS)
  {
    name = "insHardnessCleanDscrmn.RData"
  }else{
    name = "insHardness.RData"
  }  
  save(inshardList, file = paste(path,name, sep=""))
  
}



plotCorr <- function(data, type){
  require(reshape)
  require(ggplot2)
  
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  
  upper_tri <- get_upper_tri(data)
  melted_cor <- melt(upper_tri, na.rm = TRUE)
  
  melted_cor$X1 <- factor(melted_cor$X1, levels=c("Gussng","Dffclt","Dscrmn", "IH","kDN","IDS","IDCP","TD_P","TD_U","ICL","ICLD","MV","ICB"))
  melted_cor$X2 <- factor(melted_cor$X2, levels=c("Gussng","Dffclt","Dscrmn", "IH","kDN","IDS","IDCP","TD_P","TD_U","ICL","ICLD","MV","ICB"))
  
  
  ggheatmap1 <- ggplot(data = melted_cor, aes(X2, X1, fill = value))+
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "brown1", high = "darkolivegreen3", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name=paste(type,"\nCorrelation",sep=""), na.value = "white") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed() 
  
  gg1<-ggheatmap1 + geom_text(aes(X2, X1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
  
  return(gg1)
}


calcCorr <- function(path = "./_Toy_/_BinaryDS/ltm-3PL/", cleanDS=F){
  
  require(corrplot)
  require(xtable)
  if(cleanDS){
    load(paste(path,"insHardnessCleanDscrmn.RData", sep=""))
    
  }else{
    load(paste(path,"insHardness.RData", sep=""))
    
  }
  load(paste(path,"datasets.RData", sep=""))
  
  
  CorrMethod <- rep(rep(c("spearman","pearson","kendall"), each= 3) ,length(datasets))
  Param <- rep(rep(c("Gussng","Dffclt","Dscrmn"), 3) ,length(datasets))
  data <- rep(datasets, each =9)
  results <- data.frame(data,CorrMethod,Param)
  gg1 <- list()
  gg2 <- list()
  gg3 <- list()
  
  tempCorr<-matrix(nrow= 1, ncol=9)
  for (i in 1:length(inshardList)){
    
    nameDS <- datasets[i]
    nameDS <- strsplit(nameDS,"[.]")[[1]][1] #keep just the name
    
    corrSpearm <- round(cor(inshardList[[i]], use="complete.obs", method="spearm"),2) 
    gg1[[i]] <- plotCorr(corrSpearm, paste(nameDS,"\nSpearman",sep=""))
    
  
    corrPearson <- round(cor(inshardList[[i]], use="complete.obs", method="pearson") ,2) 
    gg2[[i]] <- plotCorr(corrPearson, paste(nameDS,"\nPearson",sep=""))
    
    corrKendall <-round(cor(inshardList[[i]], use="complete.obs", method="kendall") ,2) 
    gg3[[i]] <- plotCorr(corrKendall, paste(nameDS,"\nKendall",sep=""))
    
    if (ncol(corrSpearm)>3 ) {
      tempCorr <- rbind(tempCorr, corrSpearm[1:3, 4:ncol(corrSpearm)])
      tempCorr <- rbind(tempCorr, corrPearson[1:3, 4:ncol(corrSpearm)])
      tempCorr <- rbind(tempCorr, corrKendall[1:3, 4:ncol(corrSpearm)])
      
    }else{
      tempCorr <- rbind(tempCorr, matrix(ncol=9, nrow=3))
      tempCorr <- rbind(tempCorr, matrix(ncol=9, nrow=3))
      tempCorr <- rbind(tempCorr, matrix(ncol=9, nrow=3))
    }
 
    #print(title, vp = viewport(layout.pos.row=1, layout.pos.col=1:4))
    
    
    
    #write.csv(corrSpearm,paste(path,nameDS,"corrSpearm.csv",sep=""),row.names = T)
    #write.csv(corrSpearm,paste(path,nameDS,"corrPearson.csv",sep=""),row.names = T)
    #write.csv(corrSpearm,paste(path,nameDS,"corrKendall.csv",sep=""),row.names = T)
    
    
  }
  
  PDFwidth <<- 18
  PDFheight <<- 6*length(datasets)# 7 by default
  openPDFEPS(paste(path,"VisualCorrelation", sep=""))
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(length(datasets),3)))
  
  for(j in 1:length(datasets)){
    print(gg1[[j]], vp = viewport(layout.pos.row=j, layout.pos.col=1))
    print(gg2[[j]], vp = viewport(layout.pos.row=j, layout.pos.col=2))
    print(gg3[[j]], vp = viewport(layout.pos.row=j, layout.pos.col=3))
  }
  
  
  dev.off()
  tempCorr <- tempCorr[2:nrow(tempCorr),]
  resultsCorr <- cbind(results, tempCorr)
  print(xtable(resultsCorr))

  
  save(resultsCorr, file= paste(path,"corrResults.RData",sep=""))
}

calcCorrSimple <- function(path = "./_Toy_/_BinaryDS/ltm-3PL/", cleanDS=F){
  
  require(corrplot)
  require(xtable)
  if(cleanDS){
    load(paste(path,"insHardnessCleanDscrmn.RData", sep=""))
    
  }else{
    load(paste(path,"insHardness.RData", sep=""))
    
  }
  load(paste(path,"datasets.RData", sep=""))
  
  gg1 <- list()
  gg2 <- list()
  gg3 <- list()
  
 
  for (i in 1:length(inshardList)){
    
    nameDS <- datasets[i]
    nameDS <- strsplit(nameDS,"[.]")[[1]][1] #keep just the name
    
    corrSpearm <- round(cor(inshardList[[i]], use="complete.obs", method="spearm"),2) 

    PDFwidth <<- 7
    PDFheight <<- 7# 7 by default
    openPDFEPS(paste(path,"VisualCorrelation_",nameDS, sep=""))

    g<-plotCorr(corrSpearm, paste(nameDS,"\nSpearman",sep=""))
    print(g)
    
    dev.off()

  }
  
  
}


# run(path = "./_Toy_/_BinaryDS/ltm-3PL/", cleanDS=F)
# calcCorr(path = "./_Toy_/_BinaryDS/ltm-3PL/", cleanDS = F)
# calcCorrSimple(path = "./_Toy_/_BinaryDS/ltm-3PL/", cleanDS = F)
# run(path = "./_Toy_/_MulticlassDS/ltm-3PL/", cleanDS = F)
# calcCorr(path = "./_Toy_/_MulticlassDS/ltm-3PL/", cleanDS = F)
# calcCorrSimple(path = "./_Toy_/_MulticlassDS/ltm-3PL/", cleanDS = F)



TD_measure <- function(dataset){
## Input: a data frame containing the training dataset 
## For each instance of the dataset, return the depth of the leaf node 
## that classifies the instance. Computed for both: prunned and unprunned trees.  
## Output: 
## --- TD: TD measure using a prunned tree
## --- TU: TD measure using an unprunned tree

 numInstances = nrow(dataset)
 numAtributes = ncol(dataset)

 names=names(dataset)
 formula = as.formula(paste(names[numAtributes],"~.")) 

 ## unpruned tree to use in the TU measure
 fit = rpart(formula, data = dataset,method='class',minsplit = 2,cp=0)

 ## pruned tree to use in the TD measure
 pruned_fit = rpart(formula, data = dataset,method='class')

  
 ## find the depth of terminal nodes of the UNPRUNED model
 
  terminal_ids = nodeids(as.party(fit),terminal=TRUE)
 
  list = partykit:::.list.rules.party(as.party(fit))
 
  depth = list()

  for(i in terminal_ids){

   rule = list[[toString(i)]]
   
   depth[[i]] = str_count(rule, "&") + 1

  }
  
  ## index of terminal node for each training example
  nodes = predict(as.party(fit), type = "node") 
  
  TD_U_all = unlist(depth[nodes])

  
  ## find the depth of terminal nodes of the PRUNED model
 
  terminal_ids = nodeids(as.party(pruned_fit),terminal=TRUE)
 
  list = partykit:::.list.rules.party(as.party(pruned_fit))
 
  depth = list()

  for(i in terminal_ids){

   rule = list[[toString(i)]]
   
   depth[[i]] = str_count(rule, "&") + 1

  }
  
  ## index of terminal node for each training example
  nodes = predict(as.party(pruned_fit), type = "node") 
  
  TD_P_all = unlist(depth[nodes])

  list(TD=TD_P_all,TU=TD_U_all)

}

########################################################
########################################################

library(FNN)

kDN_measure <- function(dataset,k){
## For each instance, kDN is the percentage of the k nearest neighbors 
## (using Euclidean distance) for the instance that do
## not share its target class value.

 kDN_all = list()

 numInstances = nrow(dataset)
 numAtributes = ncol(dataset)

   for (i in 1:numInstances){

     ## for each test instance i, split the training data
 
     indTrain = setdiff(1:numInstances,i)

     Train = dataset[indTrain,-numAtributes]
     Test = dataset[i,-numAtributes]

     labelsTrain = factor(dataset[indTrain,numAtributes])
 
     labelTest = factor(dataset[i,numAtributes])
     ## levels(labelTest) = levels(dataset[,numAtributes])

     ## apply the knn
     model = knn(Train,Test,labelsTrain,k=k,prob=TRUE) 

     ## retrieve the indices of the nearest neighbors  
     indices <- attr(model, "nn.index")
     
     labels_NN = dataset[indices,numAtributes]
     
     ## proportion of nn labels that are different to the class label of the test instance
     kDN = length(which(as.vector(labels_NN) != as.vector(labelTest)))/k

     kDN_all[[i]] = kDN

   }
 
   return(unlist(kDN_all))

}


########################################################
########################################################

library(e1071)

CL_measure <- function(dataset){
## CL measure: For each instance, it computes the likelihood of the
## instance belonging to its class

## CLD measure: For each instance, it computes the difference between the 
## class likelihood of an instance and the maximum likelihood for all of the other classes

numAtributes = ncol(dataset)
numInstances = nrow(dataset)

c = as.factor(dataset[,numAtributes])
labels = levels(c)

prob_prior = cbind()
for(k in 1:length(labels))
{

  pclass = length(which(dataset[,numAtributes]==labels[k]))/numInstances
  prob_prior = cbind(prob_prior, pclass)
}

model = naiveBayes(dataset[,-numAtributes],as.factor(dataset[,numAtributes])) 

CL_all = list()
CLD_all = list()


for(j in 1:numInstances)
{
  

  # choose the likelihood taking into account the target class of the instance
  target = dataset[j,numAtributes]

  prob = predict(model,dataset[j,-numAtributes],type="raw")

  aux = prob/prob_prior
  lhood = aux/(sum(aux))

  lhood_target = lhood[1,which(colnames(lhood)==target)]

  CLD = lhood_target - max(lhood[1,which(colnames(lhood)!=target)])

  CL_all[j] = lhood_target
  CLD_all[j] = CLD

} 

  list(CL=unlist(CL_all),CLD=unlist(CLD_all))


}

########################################################
########################################################

MV_measure <- function(dataset){
## For each instance, MV is the ratio of the number of instances sharing 
## its target class value to the number of instances in the majority class
## CB also measures the skewness of the class that an instance belongs
## to and offers an alternative to MV.

numAtributes = ncol(dataset)
numInstances = nrow(dataset)

c = as.factor(dataset[,numAtributes])
labels = levels(c)

num_prior = cbind()
num_classes = length(labels)

for(k in 1:num_classes)
{

  numclass = length(which(dataset[,numAtributes]==labels[k]))
  num_prior = cbind(num_prior, numclass)
}

max_prior = max(num_prior)

MV_all = cbind()
CB_all = cbind()


for(j in 1:numInstances)
{
  
  # choose the likelihood taking into account the target class of the instance
  target = dataset[j,numAtributes]
  ind_target = which(labels==target)

  MV = 1 - num_prior[ind_target]/max_prior 

  CB = num_prior[ind_target]/numInstances - 1/num_classes 

  MV_all = cbind(MV_all,MV)  
  CB_all = cbind(CB_all,CB)  
  
}

 list(MV = MV_all,CB = CB_all)

}


########################################################
########################################################

DS_measure <- function(dataset){
## DS of an instance is the number of instances in a disjunct (set of instances
## belonging to a leaf) divided by the 
## number of instances covered by the largest disjunct in a data set.

 numInstances = nrow(dataset)
 numAtributes = ncol(dataset)

 names=names(dataset)
 formula = as.formula(paste(names[numAtributes],"~.")) 

 fit = rpart(formula, data = dataset,method='class',minsplit = 1,cp=0)

 
 ## find the length of the disjunct for each terminal nodes of the UNPRUNED model
 
  ## indices of terminal nodes
  terminal_ids = nodeids(as.party(fit),terminal=TRUE)
 
  ## index of terminal node for each training example
  nodes = predict(as.party(fit), type = "node") 

  disjunct_size = cbind()

  for(i in terminal_ids){

   node_i = i
   disjunct_size = cbind(disjunct_size, length(which(nodes==node_i)))
   
  }
  
  max_disjunct_size = max(disjunct_size)
 
  DS_all = cbind()


  for(i in 1:numInstances){

     terminal_i = nodes[i]
     ind = which(terminal_ids==terminal_i)
     
     DS_i = (disjunct_size[ind]-1)/(max_disjunct_size - 1)

     DS_all = cbind(DS_all, DS_i)
  
   
  }
  
  return(DS_all)

}

########################################################
########################################################

DSP_measure <- function(dataset){
## The DCP of an instance is the number of instances in a disjunct 
## (set of instances in a leaf of a prunned tree) belonging to its class 
## divided by the total number of instances in the disjunct.

 numInstances = nrow(dataset)
 numAtributes = ncol(dataset)

 names=names(dataset)
 formula = as.formula(paste(names[numAtributes],"~.")) 

 fit = rpart(formula, data = dataset,method='class')

 
 ## find the length of the disjunct for each terminal nodes of the PRUNED model
 
  ## indices of terminal nodes
  terminal_ids = nodeids(as.party(fit),terminal=TRUE)
 
  ## index of terminal node for each training example
  nodes = predict(as.party(fit), type = "node") 

  disjunct_size = cbind()

  for(i in terminal_ids){

   node_i = i
   disjunct_size = cbind(disjunct_size, length(which(nodes==node_i)))
   
  }
 
  DSP_all = cbind()

  for(i in 1:numInstances){

     terminal_i = nodes[i]
     ind = which(terminal_ids==terminal_i)
     
     instances_node = which(nodes==terminal_i)

     same_class = length(which(dataset[instances_node,numAtributes] == dataset[i,numAtributes])) 

     DSP_i = same_class/disjunct_size[ind]
     DSP_all = cbind(DSP_all, DSP_i)

  }
  
  return(DSP_all)

}

########################################################
########################################################
