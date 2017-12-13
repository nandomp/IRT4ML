# ORIGINAL DATASET 
# runExp(CVfolds = 5, balanceTestFold = F, UnBal_exp = F, Jitter = F, FixTooMuchNegatives = F, ds = "_Toy_/", model = "LTM", NPL = 3  )
# runExp_ds(ind_dataset = -1, balanceTestFold = F, UnBal_exp = F, Jitter = F, FixTooMuchNegatives = F, ds = "_Toy_/" , model = "LTM", NPL = 2 )
# runExp_ds(ind_dataset = -1, balanceTestFold = F, UnBal_exp = F, Jitter = F, FixTooMuchNegatives = F, ds = "_Toy_/" , model = "LTM", NPL = 1 )

# runExp(CVfolds = 5, balanceTestFold = F, UnBal_exp = F, Jitter = F, FixTooMuchNegatives = F, ds = "_Toy_/", model = "LTM", NPL = 2  )
# runExp(CVfolds = 5, balanceTestFold = F, UnBal_exp = F, Jitter = F, FixTooMuchNegatives = F, ds = "_Toy_/", model = "LTM", NPL = 1  )
# runExp(CVfolds = 5, balanceTestFold = F, UnBal_exp = F, Jitter = F, FixTooMuchNegatives = F, ds = "_Toy_/", model = "MIRT", NPL = 3  )

# BALANCED TEST FOLD
# runExp(CVfolds = 5, balanceTestFold = T, UnBal_exp = T, Jitter = F, FixTooMuchNegatives = F, ds = "_Toy_/", model = "LTM", NPL = 3 )

runExp <- function(CVfolds = 5, balanceTestFold = F, UnBal_exp = F, Jitter = F, FixTooMuchNegatives = T, ds = "_Toy_/", model = "LTM", NPL = 3){
  source("methodsResponses_v1.2.R")
  #(ds= "_Toy_/", seed = 288, controlExperiment = FALSE, CVfolds = 5, balanceTestFold = FALSE)
  # CVfolds > 1 = CV experiment
  # balanceTestFold = balance test fold 
  obtain_responses(CVfolds = CVfolds, balanceTestFold = balanceTestFold, ds = ds)
  
  source("IRT_model_fit_and_plot_v1.1.R")
  #fitPlotIRT<- function(UnBal_exp=FALSE, balTestFold=FALSE, JITTER=FALSE ){ 
  #UnBal_exp = use same PCA (for unabalanced and balanced experiments)
  # balTestFold = load the balanced datasets (bigger than the original one) if done when obtaining responses from algorithms.
  fitPlotIRT(UnBal_exp = UnBal_exp, balTestFold = balanceTestFold, JITTER = Jitter, FixTooMuchNegatives = FixTooMuchNegatives,
             ds = ds, model=model, NPL = NPL)
  
}

# ind_dataset = -1 --> All datasets
# ind_dataset = 1,2,3,4... just ds id = 1,2,3,4,...

# runExp_ds(ind_dataset = -1, balanceTestFold = F, UnBal_exp = F, Jitter = F, FixTooMuchNegatives = F, ds = "_Toy_/" , model = "LTM", NPL = 2)
# runExp_ds(ind_dataset = -1, balanceTestFold = F, UnBal_exp = F, Jitter = F, FixTooMuchNegatives = F, ds = "_Toy_/" , model = "LTM", NPL = 1)
runExp_ds <- function(ind_dataset = 1, balanceTestFold = F, UnBal_exp = F, Jitter = F, FixTooMuchNegatives = F, ds = "_Toy_/", model = "LTM", NPL = 3 ){
    
    source("IRT_model_fit_and_plot_v1.1.R")
    #fitPlotIRT<- function(UnBal_exp=FALSE, balTestFold=FALSE, JITTER=FALSE ){ 
    #UnBal_exp = use same PCA (for unabalanced and balanced experiments)
    # balTestFold = load the balanced datasets (bigger than the original one) if done when obtaining responses from algorithms.
    fitPlotIRT(UnBal_exp = UnBal_exp, balTestFold = balanceTestFold, JITTER = Jitter, 
               FixTooMuchNegatives = FixTooMuchNegatives, ds = ds, ind_ds = ind_dataset, model=model, NPL = NPL)
    
}

runExpMCC <- function(CVfolds = 5, balanceTestFold = F, UnBal_exp = F, Jitter = F, FixTooMuchNegatives = T, ds = "_Toy_/", model = "LTM", NPL = 3){
  source("methodsResponses_v1.2.R")
  #(ds= "_Toy_/", seed = 288, controlExperiment = FALSE, CVfolds = 5, balanceTestFold = FALSE)
  # CVfolds > 1 = CV experiment
  # balanceTestFold = balance test fold 
  obtain_responses(CVfolds = CVfolds, balanceTestFold = balanceTestFold, ds = ds)
  
  source("IRT_model_fit_and_plot_v1.1.R")
  #fitPlotIRT<- function(UnBal_exp=FALSE, balTestFold=FALSE, JITTER=FALSE ){ 
  #UnBal_exp = use same PCA (for unabalanced and balanced experiments)
  # balTestFold = load the balanced datasets (bigger than the original one) if done when obtaining responses from algorithms.
  print("Extract Data")
  extract_data_n(nas=FALSE, all= FALSE, FixTooMuchNegatives = FALSE, ds = ds, ind_ds = 1, model = model, NPL = NPL)
  MCC(groups = 6)
  
}

runComparison_Hangout <- function(){
  source("methodsResponses.R")
  source("3PLmodelsComparison.R")
  extract_data_comp_4hangout()
}

#arff2csv(pat="*csv$")
#arff2csv(pat="*arff$")
arff2csv <- function( ds = "datasets/", pat = "*arff$"){
 
  datasets <- list.files(ds, pattern = pat)
  for(i in 1:length(datasets)){
    datos <- read.arff(paste(ds,datasets[i],sep=""))
    colnames(datos)[ncol(datos)]<-"Class"  #Changed by Adolfo: 04/04/2016
    
    if(length(levels(datos$Class))<3){
      levels(datos$Class) <- c("0","1")
    }
    
    datosComplete <- datos[complete.cases(datos),]
    
    for (j in (1:(ncol(datosComplete)-1))){
      datosComplete[,j] <- as.numeric(datosComplete[,j])  
    }
    
    datosComplete.woCons<-datosComplete[,!apply(datosComplete, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
    
    nameDS <- datasets[i]
    nameDS <- strsplit(nameDS,"[.]")[[1]][1] #keep just the name
    write.csv(datosComplete.woCons, file=paste(ds,nameDS,"-v2.csv",sep=""),row.names = FALSE)
    print(i)
  }
  
}

cleanDS_rownames <- function(){
  ds <- "_Toy_/"
  datasets <- list.files(ds, pattern = "*csv$")
  for(i in 1:length(datasets)){
    datos <- read.csv(paste(ds,datasets[i],sep=""))
    nameDS <- datasets[i]
    nameDS <- strsplit(nameDS,"[.]")[[1]][1] #keep just the name
    write.csv(datos[,2:ncol(datos)], file=paste(ds,nameDS,".2csv",sep=""), row.names = FALSE)
  }
}