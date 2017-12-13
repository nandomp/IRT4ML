# 
# SOURCE_CODE ="D:/OneDrive/Rworks/IRT/Noise"
# DATASET2PLAY = "D:/OneDrive/Rworks/IRT/Noise/_Toy_"
# 
# setwd(DATASET2PLAY)

###############################################
############# OPERATING OPTIONS ###############
###############################################

set.seed(288)
ind_instance = 1
ind_dataset=1
args <- commandArgs(trailingOnly = TRUE)
options(scipen = 99999)

OptimimalParam = TRUE
binaryClass = FALSE
weArePlaying = TRUE

#Temp <- "_incrementalNoise_/"
Temp <- "_Toy_/"

if (weArePlaying){
  ds <- Temp
  
}else{
  if (binaryClass) {
    ds<-"_Binary_/"
  }else{
    ds<-"_Multiclass_/"
  }
  
}



###############################################
#############     LIBRARIES     ###############
###############################################

.lib<- c("ltm","devtools", "ggplot2","stats", "ggrepel", "wordcloud", "Hmisc", "ggfortify","grid","mirt", "dplyr", "reshape2") #"mirt"
.inst <- .lib %in% installed.packages()
if (length(.lib[!.inst])>0) install.packages(.lib[!.inst], repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com")) 
lapply(.lib, require, character.only=TRUE)

# Especific library for ploting PCAs
#install_github("ggbiplot", "vqv")
#library(ggbiplot)


###############################################
#############    PDF/EPS GEN    ###############
###############################################


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


###############################################
############  GLOBAL VARIABLES  ###############
###############################################
numMethods <- 0
item_param <- list()
results <- list()

numDS <- 0
all_abilities <-  matrix()
items.fit <- list()

acc <- matrix()

all_models <- list()
avgProbs <- matrix()
medianProbs <- matrix()

responses_pSucc <- list() #list for all datasets
goodnessFit_measures <- list() # list of dataframes with info about item.fitnes and correlation measures (responses vs pSucc)

plots2compare <- list()

errorDS <- c()

#Mylabels = c("RndA","RndB", "RndC", "Maj","Min","Opt", "Pess")
#MyclassifiersNames = c("RandomClass_A", "RandomClass_B", "RandomClass_C","MajorityClass", "MinorityClass", "OptimalClass","PessimalClass" )
#MyclassifiersCut = c("RandomClass_A", "RandomClass_B", "RandomClass_C")

Mylabels = c("RndA","RndB","RndC", "Maj", "Min", "Opt", "Pess")
MyclassifiersNames = c("RandomClass_A","RandomClass_B","RandomClass_C", "MajorityClass", "MinorityClass", "OptimalClass", "PessimalClass" )
MyclassifiersCut = c("RandomClass_A","RandomClass_B","RandomClass_C")


# MIRT package
fit_mIRT <-function(allresp, type, rnd = FALSE){
  print("-MIRT-")
  colnames(allresp)<-1:ncol(allresp)
  
  if(type == 3){ 
    if(rnd){
      fit <- mirt(allresp,1,itemtype = '3PL', technical = list(NCYCLES = 500), GenRandomPars = TRUE)
    }else{
      fit <- mirt(allresp,1,itemtype = '3PL', technical = list(NCYCLES = 500))
      
    }
    
  }
  
  ## Extracting the items' parameters: 
  ## Gussng (ci), Dffclt (bi) and Dscrmn (ai) 
  
  temp = coef(fit, simplify = T, IRTpars =T)$items
  item_param <- temp[,c("g","b","a")]
  colnames(item_param)<-c("Gussng","Dffclt","Dscrmn")
  
  ## computing the abilities 'ab_vector' of the respondents   
  
  abil<-t(fscores(fit))
  
  return(list(model = fit, item_param = item_param, abil_vector = abil))
}


# LTM package
fit_IRT <- function(allresp,type,rnd=F){
  ## builds the IRT models given the responses allresp and the model type
  ## requires the ltm package
  
  ## allresp: binary matrix matrix (with dimension nrow vs ncol) storing 
  ##          the ncol responses of nrow respondents    
  ## type in {1,2,3}: indicates the number of parameters of the IRT model 
  ##                  (i.e., 1P, 2P or 3P IRT model) 
  
  ## calling the tpm function implemented in the ltm package
  
  if(type == 3){ 
    if(rnd){
      fit <- tpm(allresp, type = "latent.trait",  IRT.param=TRUE, start.val = "random", control= list(optimizer = "nlminb"))
      #fit <- tpm(allresp, type = "latent.trait",  IRT.param=TRUE, start.val = "random")
    }else{
      fit <- tpm(allresp, type = "latent.trait",  IRT.param=TRUE, control= list(optimizer = "nlminb"))
      #fit <- tpm(allresp, type = "latent.trait",  IRT.param=TRUE)
    }
  }
  
  nitems = ncol(allresp)
  if(type == 2){
    ## Parameter Gussng (ci) constrained to zero
    fit <- tpm(allresp, type = "latent.trait",  IRT.param=TRUE, start.val = "random", constraint = cbind(1:nitems, 1, 0))
  }
  
  if(type == 1){
    ## Parameter Gussng (ci) constrained to zero
    ## Parameter Dscrmn (ai) constrained to one
    fit <- tpm(allresp, type = "latent.trait",  IRT.param=TRUE, start.val = "random", constraint =  rbind(cbind(1:nitems, 1, 0),cbind(1:nitems, 3, 1)))
  }
  
  ## Extracting the items' parameters: 
  ## Gussng (ci), Dffclt (bi) and Dscrmn (ai) 
  
  item_param = coef(fit)
  
  
  
  ## computing the abilities 'ab_vector' of the respondents   
  r = factor.scores(fit,resp.patterns=allresp)
  abil_vector = r$score.dat$z1
  
  return(list(model = fit, item_param = item_param, abil_vector = abil_vector))
  
}
  

## Requires the ltm package
plot_ICC <- function(all_models, results, all_abilities, ind_dataset, ind_instance, 
                     main = "Item Characteristic Curve", classifiers = TRUE, randomCuts = TRUE, 
                     labels = Mylabels,
                     classifiersNames = MyclassifiersNames,
                     classifiersCut = MyclassifiersCut){
  ## plots the ICC of the "ind_instance-th" instance of the "ind_dataset-th" dataset
  
  fit = all_models[[ind_dataset]]$model
  abil = all_abilities[ind_dataset,]
  resp = results[[ind_dataset]][ind_instance,]#resp = results[[ind_dataset]][ind_instance,]
  
  plot(fit,items=ind_instance,xlim=cbind(-4,4),ylim=cbind(0,1),annot=FALSE,main = main)
  par(new=TRUE) 
  plot(abil,resp,xlim=cbind(-4,4),ylim=cbind(0,1),xlab="",ylab="")
  
  
  # Print last 7 ML methods ("RndA","RndB", "RndC", "Maj","Min","Opt", "Pessim") as points...
  if(classifiers){
    xabil <- c(abil[classifiersNames])
    yresp <- c(resp[classifiersNames])
    par(new=TRUE) 
    #points(xabil[1:(length(xabil)-2)], yresp[1:(length(yresp)-2)], xlim=cbind(-4,4),ylim=cbind(0,1),xlab="",ylab="",cex = .5, col = "green")
    #points(xabil[(length(xabil)-1):(length(xabil))], yresp[(length(yresp)-1):(length(yresp))], xlim=cbind(-4,4),ylim=cbind(0,1),xlab="",ylab="",cex = .5, col = "green")
    par(new=TRUE) 
    textplot(xabil,yresp, labels, xlim=c(-4,4),ylim=c(0,1),xlab="",ylab="",show.lines = TRUE, cex= 1.5)
    # text(xabil,yresp, labels=c("RndA","RndB", "RndC", "Maj","Min","Opt", "Pess"), cex= 0.8, pos=4, font = 4, srt=90)
  }
  
  # Random Classifiers CUT POINTS (ad-hoc)
  if (randomCuts){
    rnd <- classifiersCut
    
    # UPDATE : 3 Random Classifiers 
    for(i in rnd){ #random classifiers are place in the positions last-4 last-5 last-6
      
      RandomModelDiff = abil[i]
      abline(v=RandomModelDiff, col="red", lty=2)
      
      a = all_models[[ind_dataset]]$item_param[ind_instance,3]
      b = all_models[[ind_dataset]]$item_param[ind_instance,2]
      c = all_models[[ind_dataset]]$item_param[ind_instance,1]
      theta = RandomModelDiff
      
      y = c + (1-c)/(1+exp(-a*(theta-b)))
      abline(h=y, col="red", lty=2)
      
      #print probability of success
      # text(x=-3.8, y=y,paste(round(y, digits=4)),cex= 0.8, pos=3)
      # Success or fail? (I depend just on 1 random classifier)
      # if(results[[ind_dataset]][ind_instance,rnd[i-3]] == 1){ #SUCCESS
      #   text(x=RandomModelDiff, y=1,paste(round(RandomModelDiff, digits=4)),cex= 0.8, pos=2)
      # }else{ #FAILURE
      #   text(x=RandomModelDiff, y=0,paste(round(RandomModelDiff, digits=4)),cex= 0.8, pos=2)
      # }
    }
    
  }
} 


plot_ICC_ECAI <- function(all_models,results,all_abilities,ind_dataset,ind_instance, main = "Item Characteristic Curve", randomCuts = FALSE){
  ## plots the ICC of the "ind_instance-th" instance of the "ind_dataset-th" dataset
  
  
  
  
  fit = all_models[[ind_dataset]]$model
  abil = all_abilities[ind_dataset,]
  resp = results[[ind_dataset]][ind_instance,]#resp = results[[ind_dataset]][ind_instance,]
  
  plot(fit,items=ind_instance,xlim=cbind(-4,4),ylim=cbind(0,1),annot=FALSE,main = main,cex.lab=1.2)
  par(new=TRUE) 
  plot(abil,resp,xlim=cbind(-4,4),ylim=cbind(0,1),xlab="",ylab="",cex = 1.5)
  
  xabil <- c(abil[(length(abil)-6):length(abil)])
  yresp <- c(resp[(length(resp)-6):length(resp)])
  par(new=TRUE) 
  
  #points(xabil[1:(length(xabil)-2)], yresp[1:(length(yresp)-2)], xlim=cbind(-4,4),ylim=cbind(0,1),xlab="",ylab="",cex = .5, col = "green")
  #points(xabil[(length(xabil)-1):(length(xabil))], yresp[(length(yresp)-1):(length(yresp))], xlim=cbind(-4,4),ylim=cbind(0,1),xlab="",ylab="",cex = .5, col = "green")
  par(new=TRUE) 
  
  textplot(xabil,yresp, c("RndA","RndB", "RndC", "Maj","Min","Opt", "Pess"), xlim=c(-4,4),ylim=c(0,1),xlab="",ylab="",show.lines = TRUE, cex= 1.5)
  #textplot(xabil,yresp, c("RandomClass_A","RandomClass_B", "RandomClass_C", "MajorityClass","MinorityClass","OptimalClass", "PessimalClass"), xlim=c(-4,4),ylim=c(0,1),xlab="",ylab="",show.lines = TRUE, cex= 1.5)
  
  
  
  
  # Random Classifiers CUT POINTS (ad-hoc)
  if (randomCuts){
    rnd <- c("RandomClass_A", "RandomClass_B", "RandomClass_C")
    
    # UPDATE : 3 Random Classifiers 
    for(i in 4:6){
      
      RandomModelDiff = abil[(length(abil)-i)]
      abline(v=RandomModelDiff, col="red", lty=2)
      
      a = all_models[[ind_dataset]]$item_param[ind_instance,3]
      b = all_models[[ind_dataset]]$item_param[ind_instance,2]
      c = all_models[[ind_dataset]]$item_param[ind_instance,1]
      theta = RandomModelDiff
      
      y = c + (1-c)/(1+exp(-a*(theta-b)))
      abline(h=y, col="red", lty=2)
      
      #print probability of success
      #text(x=-3.8, y=y,paste(round(y, digits=4)),cex= 0.8, pos=3)
      # Success or fail? (I depend just on 1 random classifier)
      # if(results[[1]][ind_instance,rnd[i-3]] == 1){ #results[[ind_dataset]]
      #   text(x=RandomModelDiff, y=1,paste(round(RandomModelDiff, digits=4)),cex= 0.8, pos=2)
      # }else{
      #   text(x=RandomModelDiff, y=0,paste(round(RandomModelDiff, digits=4)),cex= 0.8, pos=2)
      # }
    }
    
  }
  
  
  
  
} 

#generic PLOT (used for the  models generated with the MIRT package)
plot_mICC <- function(all_models,results,all_abilities,ind_dataset,ind_instance, main ="Item Characteristic Curve", randomCuts = TRUE){
  ## plots the ICC of the "ind_instance-th" instance of the "ind_dataset-th" dataset
  
  #item_param = item_param[[ind_dataset]]
  
  fit = all_models[[ind_dataset]]$model
  abil = all_abilities[ind_dataset,]
  resp = results[[ind_dataset]][ind_instance,]#resp = results[[ind_dataset]][ind_instance,]
  
  
  Probability <- c()
  Ability <- seq(-6,6,0.05)
  for (theta in Ability){
    a = all_models[[ind_dataset]]$item_param[ind_instance,3]
    b = all_models[[ind_dataset]]$item_param[ind_instance,2]
    c = all_models[[ind_dataset]]$item_param[ind_instance,1]
    y_temp <- c + (1-c)/(1+exp(-a*(theta-b)))
    Probability <- c(Probability,y_temp)
  }
  
  plot(Ability,Probability, main = main, type = "l",xlim=cbind(-4,4),ylim=cbind(0,1))
  par(new=TRUE) 
  
  plot(abil,resp,xlim=cbind(-4,4),ylim=cbind(0,1),xlab="",ylab="")
  
  xabil <- c(abil[(length(abil)-6):length(abil)])
  yresp <- c(resp[(length(resp)-6):length(resp)])
  
  points(xabil[1:(length(xabil)-2)], yresp[1:(length(yresp)-2)], xlim=cbind(-4,4),ylim=cbind(0,1),xlab="",ylab="",cex = .5, col = "green")
  points(xabil[(length(xabil)-1):(length(xabil))], yresp[(length(yresp)-1):(length(yresp))], xlim=cbind(-4,4),ylim=cbind(0,1),xlab="",ylab="",cex = .5, col = "red")
  
  text(xabil, yresp, labels=c("RndA","RndB", "RndC", "Maj", "Min", "Opt", "Pessim"), cex = 0.8, pos=4, font = 4, srt = 90)
  
  
  
  # Random Classifiers CUT POINTS (ad-hoc)
  
  if (randomCuts){
    rnd <- c("RandomClass_A", "RandomClass_B", "RandomClass_C")
    
    # UPDATE : 3 Random Classifiers 
    for(i in 4:6){
      
      RandomModelDiff = abil[(length(abil)-i)]
      abline(v=RandomModelDiff, col="red", lty=2)
      
      a = all_models[[ind_dataset]]$item_param[ind_instance,3]
      b = all_models[[ind_dataset]]$item_param[ind_instance,2]
      c = all_models[[ind_dataset]]$item_param[ind_instance,1]
      theta = RandomModelDiff
      
      y = c + (1-c)/(1+exp(-a*(theta-b)))
      abline(h=y, col="red", lty=2)
      
      #print probability of success
      text(x=-3.8, y=y,paste(round(y, digits=4)),cex= 0.8, pos=3)
      
      # Success or fail? (I depend just on 1 random classifier)
      if(results[[1]][ind_instance,rnd[i-3]] == 1){ #results[[ind_dataset]]
        text(x=RandomModelDiff, y=1,paste(round(RandomModelDiff, digits=4)),cex= 0.8, pos=2)
      }else{
        text(x=RandomModelDiff, y=0,paste(round(RandomModelDiff, digits=4)),cex= 0.8, pos=2)
      }
    }
    
  }
  
  
} 


plotICCi<- function(i,ind_dataset){
  plot_mICC(all_models, results, all_abilities,ind_dataset,i)
}

# openPDFEPS("ICC_376")
# plotICCi(146,-4,4)
# dev.off()
#ind_dataset = 74





#Extract data from the binary responses given by the classifiers (n datasets)
extract_data_n <- function(nas=TRUE, all= FALSE, negDisc = TRUE, maxItersIRT = 20, FixTooMuchNegatives = TRUE,
                           ds= "_Toy_/", ind_ds=-1, model = "LTM", NPL = 3){
  
  load(paste(ds,"ListAllResults.RData",sep="")) # ListDS_Results
  load(paste(ds,"Methods.RData",sep="")) #methods
  load(paste(ds,"Datasets.RData",sep="")) #datasets, ds (directory)
  
  numMethods <<- length(methods)
  numDS <<- length(ListDS_Results)
  all_abilities <<-  matrix(rep(NA, numMethods * numDS), nrow=numDS, ncol=numMethods, byrow = T)
  colnames(all_abilities) <- methods
  acc <<- matrix(rep(NA, numMethods * numDS), nrow=numDS, ncol=numMethods, byrow = T)
  colnames(acc) <- methods
  avgProbs <<- matrix(rep(NA, numMethods * numDS), nrow=numDS, ncol=numMethods, byrow = T)
  medianProbs <<- matrix(rep(NA, numMethods * numDS), nrow=numDS, ncol=numMethods, byrow = T)
  
  ERROR <-  tryCatch(load(paste(ds,"seedsOK.RData",sep="")), error = function(e) {return(TRUE)})
  if (is.logical(ERROR)){
    seeds <- rep(21,length(datasets))
  }
  

  if (ind_ds == -1){
    c_ds <- 1:length(ListDS_Results)
  }else{
    c_ds <- ind_ds
  }
  
  for (ind_dataset in c_ds)
  {
    print(paste("DS:  ",datasets[ind_dataset]))
    results[[ind_dataset]] <- ListDS_Results[[ind_dataset]]*1#From logical to numerical
    
    # Are there NA's in the results (no predictions for items)?
    if(nas){  
      if (sum(is.na(ListDS_Results[[ind_dataset]]))){
        results[[ind_dataset]][is.na(results[[ind_dataset]])] <- 0
      }
    }
    
    #Avoid items with one response category
    if(all){
      clean<-c()
      for (i in 1:nrow(results[[ind_dataset]])){
        if (length(unique(results[[ind_dataset]][i,]))>1){
          clean <- c(clean,i)
        }
      }
      results[[ind_dataset]] <<- result[[ind_dataset]][clean,]
    }
    
    t_results <- t(results[[ind_dataset]])
    oldw <- getOption("warn")
    options(warn = -1)

    negDisc =TRUE
    maxItersIRT=40
    firstTimeModel = T
    while(negDisc & maxItersIRT > 0){
      print(paste("Seed:", seeds[ind_dataset]))
      set.seed(seeds[ind_dataset])
      
      if (model=="LTM"){
        print(paste("IRT LTM (",NPL,"PL) - dataset: ",ind_dataset, " (",negDisc," - ",maxItersIRT,")",sep=""))
        if (firstTimeModel){
          print("-- 1st Try --")
          ERROR <-  tryCatch(IRTstuff<- fit_IRT(t_results,NPL, rnd=F), error = function(e) {return(TRUE)})
          firstTimeModel = F
        }else{
          print("-- Nth Try --")
          ERROR <-  tryCatch(IRTstuff<- fit_IRT(t_results,NPL, rnd=T), error = function(e) {return(TRUE)})
        }
        
        
      }else{
        print(paste("Seed: ", seeds[ind_dataset], sep=""))
        print(paste("IRT MIRT (",NPL,"PL) - dataset: ",ind_dataset, " (",negDisc," - ",maxItersIRT,")",sep=""))
        if (firstTimeModel){
          ERROR <-  tryCatch(IRTstuff<- fit_mIRT(t_results,NPL, rnd = FALSE), error = function(e) {return(TRUE)})
          firstTimeModel = F
        }else{
          ERROR <-  tryCatch(IRTstuff<- fit_mIRT(t_results,NPL, rnd = TRUE), error = function(e) {return(TRUE)})
        }
      }
 

      if (!is.logical(ERROR)){
        
        negs <- sum(IRTstuff$item_param[,"Dscrmn"]<0)
        
        if ((negs > (nrow(IRTstuff$item_param))/2 ) & FixTooMuchNegatives){
            print(paste("Neg Dscrm examples (",negs,") > examples/3 "))
            seeds[ind_dataset] <- seeds[ind_dataset] + ((2*maxItersIRT)-1)
            maxItersIRT <- maxItersIRT - 1
            negDisc = TRUE
            
        }else{
          
          negDisc = FALSE
          print("Finished Correctly: IRT model obtained")
          options(warn = oldw)
          
          item_param[[ind_dataset]] <- IRTstuff$item_param
          all_abilities[ind_dataset,] <- IRTstuff$abil_vector
          acc[ind_dataset,] <- colMeans(results[[ind_dataset]],na.rm = TRUE)
          all_models[[ind_dataset]] <- IRTstuff
          
          # ICC models fitness measures (for LTM and MIRT packages)
          if (model == "LTM"){
            items.fit[[ind_dataset]] <- item.fit(IRTstuff$model)
          }else{#MIRT
            items.fit[[ind_dataset]] <- itemfit(all_models[[ind_dataset]]$model, fit_stats="X2")
            names(items.fit[[ind_dataset]])<-c("item","Tobs","df.X2","p.values")
          }
          
          
          # Compute Average Probability of success for the all the  Classifiers 
          abil = all_abilities[ind_dataset,]
          avgProbsT = vector()
          medianProbsT = vector()
          for(m in 1:length(methods)){
            C_method_probs <- vector()
            
            for (ind_inst in 1:nrow(results[[ind_dataset]])){
              
              ModelProf = abil[m]
              a = all_models[[ind_dataset]]$item_param[ind_inst,3]
              b = all_models[[ind_dataset]]$item_param[ind_inst,2]
              c = all_models[[ind_dataset]]$item_param[ind_inst,1]
              theta = ModelProf
              y = c + (1-c)/(1+exp(-a*(theta-b)))
              C_method_probs <- c(C_method_probs,y)
            }
            
            MethodAvgProb <- sum(C_method_probs)/length(C_method_probs)
            
            avgProbsT <- c(avgProbsT, MethodAvgProb)
            medianProbsT <- c(medianProbsT,median(C_method_probs))
            
          }
          
          avgProbs[ind_dataset,] <- avgProbsT 
          medianProbs[ind_dataset, ] <- medianProbsT
          
          
        }#too much negative discriminants...
        
      }else{
        # IRT fails
        print(paste("----ERROR IRT: ",ERROR))
        #negDisc = FALSE #exit while loop due to the IRT failure calculation
        if (maxItersIRT > 0){
        seeds[ind_dataset] <- seeds[ind_dataset] + ((2*maxItersIRT)-1)
        }
        maxItersIRT <- maxItersIRT - 1
        negDisc = TRUE
      }#tryCatch
      
    }
    if (maxItersIRT == 0){
      errorDS <- c(errorDS, ind_dataset)
    }
    
    
  }#for
  
  
  # Save variables
  # save(IRTstuff, file=paste(ds,"irt_model.RData",sep=""))
  
  save(items.fit, file = paste(ds,"irt_itemsfitness.RData",sep=""))
  save(item_param, file=paste(ds,"irt_parameters_mc.RData",sep=""))
  save(all_abilities, file =paste(ds,"algor_abilities_mc.RData",sep=""))
  save(acc,file=paste(ds,"algor_accuracies_mc.RData",sep=""))
  save(results, file=paste(ds,"results_responses_mc.RData",sep=""))
  save(all_models, file=paste(ds, "all_3P_IRT_models_mc.RData",sep="")) #IRTstuff
  save(avgProbs, file=paste(ds, "all_avgProbs.RData",sep=""))
  save(medianProbs, file=paste(ds, "all_medianProbs.RData",sep=""))
  save(errorDS, file = paste(ds, "errorDS.RData",sep=""))
  save(seeds, file = paste(ds, "seedsOK.Rdata", sep=""))
  
  
}



# Check out fitness measure


saveFitnessMeasure <- function(model ="LTM"){
  
  load(paste(ds,"ListAllResults.RData",sep="")) # ListDS_Results
  load(paste(ds,"Methods.RData",sep="")) #methods
  load(paste(ds,"Datasets.RData",sep="")) #datasets, ds (directory)
  
  load(paste(ds,"irt_parameters_mc.RData",sep=""))
  load(paste(ds, "algor_abilities_mc.RData",sep=""))
  load(paste(ds,"algor_accuracies_mc.RData",sep=""))
  load(paste(ds,"results_responses_mc.RData",sep=""))
  load(paste(ds,"all_3P_IRT_models_mc.RData",sep=""))
  load(paste(ds,"all_avgProbs.RData",sep=""))
  load(paste(ds,"all_medianProbs.RData",sep=""))
  load(paste(ds, "errorDS.RData",sep=""))
  
  items.fit <- list()
  for(ind_dataset in 1:length(datasets)){
    print(paste("DS: ",ind_dataset))
    if (model == "LTM"){
      items.fit[[ind_dataset]] <- item.fit(all_models[[ind_dataset]]$model)
    }else{#MIRT
      items.fit[[ind_dataset]] <- itemfit(all_models[[ind_dataset]]$model, fit_stats="X2")
      names(items.fit[[ind_dataset]])<-c("item","Tobs","df.X2","p.values")
    }
    
  }
  save(items.fit, file = paste(ds,"irt_itemsfitness.RData",sep=""))
  
}

whatAboutFitness <- function(ind_dataset){
  
  load(paste(ds,"ListAllResults.RData",sep="")) # ListDS_Results
  load(paste(ds,"Methods.RData",sep="")) #methods
  load(paste(ds,"Datasets.RData",sep="")) #datasets, ds (directory)
  
  load(paste(ds,"irt_parameters_mc.RData",sep=""))
  load(paste(ds, "algor_abilities_mc.RData",sep=""))
  load(paste(ds,"algor_accuracies_mc.RData",sep=""))
  load(paste(ds,"results_responses_mc.RData",sep=""))
  load(paste(ds,"all_3P_IRT_models_mc.RData",sep=""))
  load(paste(ds,"all_avgProbs.RData",sep=""))
  load(paste(ds,"all_medianProbs.RData",sep=""))
  load(paste(ds, "errorDS.RData",sep=""))
  load(paste(ds,"irt_itemsfitness.RData",sep=""))
  
  numInstances = nrow(results[[ind_dataset]])
  # Compute Average Probability of success for the all the  Classifiers 
  abil = all_abilities[ind_dataset,]
  responses_pSucc_listInst <- list() #lists (of matrices) of instances
  spearm <- c()
  spearm.pvalue <- c()
  pearson <- c()
  pearson.pvalue <- c()
  kendall <- c()
  kendall.pvalue <- c()
  
  for (ind_inst in 1:numInstances){
    reponses_psuccess_instance = matrix(rep(NA, 2 * length(methods)), nrow=length(methods), ncol=2, byrow = T) 
    
    for(m in 1:length(methods)){
      
      #cols: 1->response 2->pSuccess  
      
      ModelProf = abil[m]
      
      a = all_models[[ind_dataset]]$item_param[ind_inst,3]
      b = all_models[[ind_dataset]]$item_param[ind_inst,2]
      c = all_models[[ind_dataset]]$item_param[ind_inst,1]
      theta = ModelProf
      y = c + (1-c)/(1+exp(-a*(theta-b)))
      
      reponses_psuccess_instance[m,1] <- results[[ind_dataset]][ind_inst,m]
      reponses_psuccess_instance[m,2] <- y
      
      
    }
    responses_pSucc_listInst[[ind_inst]] <- reponses_psuccess_instance
    
    resp <- reponses_psuccess_instance[,1]
    pSucc <- reponses_psuccess_instance[,2]
    
    s <- cor.test(resp, pSucc, method = "spearm", alternative = "two.sided")
    p <- cor.test(resp, pSucc, method = "pearson", alternative = "two.sided")
    k <- cor.test(resp, pSucc, method = "kendall", alternative = "two.sided")
    
    spearm <- c(spearm,as.vector(round(s$estimate,5)))
    pearson <- c(pearson,as.vector(round(p$estimate,5)))
    kendall <- c(kendall,as.vector(round(k$estimate,5)))
    
    spearm.pvalue <- c(spearm.pvalue,as.vector(round(s$p.value,5)))
    pearson.pvalue <- c(pearson.pvalue,as.vector(round(p$p.value,5)))
    kendall.pvalue <- c(kendall.pvalue,as.vector(round(k$p.value,5)))
    
    lt0.0001 <- as.vector(items.fit[[ind_dataset]]$p.values)<0.0001
  }
  respSucList <- responses_pSucc_listInst
  ds <- data.frame(item.fit=as.vector(round(items.fit[[ind_dataset]]$Tobs,5)), 
                   p.value = as.vector(round(items.fit[[ind_dataset]]$p.values,5)),
                   lt0.0001 = lt0.0001,
                   spearm = spearm,  pearson = pearson, kendall = kendall,
                   spearm.pv = spearm.pvalue, pearson.pv = pearson.pvalue, kendall.pv = kendall.pvalue)
  
  return(list(responses_pSucc_listInst = respSucList, goodnessFit_measures = ds))
  
  
}


analyseGoodnessDS <- function(){
  load(paste(ds,"ListAllResults.RData",sep="")) # ListDS_Results
  load(paste(ds,"Methods.RData",sep="")) #methods
  load(paste(ds,"Datasets.RData",sep="")) #datasets, ds (directory)
  
  load(paste(ds,"irt_parameters_mc.RData",sep=""))
  load(paste(ds, "algor_abilities_mc.RData",sep=""))
  load(paste(ds,"algor_accuracies_mc.RData",sep=""))
  load(paste(ds,"results_responses_mc.RData",sep=""))
  load(paste(ds,"all_3P_IRT_models_mc.RData",sep=""))
  load(paste(ds,"all_avgProbs.RData",sep=""))
  load(paste(ds,"all_medianProbs.RData",sep=""))
  load(paste(ds, "errorDS.RData",sep=""))
  load(paste(ds,"irt_itemsfitness.RData",sep=""))
  
  for (i in 1:length(ListDS_Results)){
    temp <- whatAboutFitness(i)
    responses_pSucc_listInst[[i]] <- temp$responses_pSucc_listInst
    goodnessFit_measures[[i]] <- temp$goodnessFit_measures
    
  }
  save(responses_pSucc_listInst, file = paste(ds,"responses_pSucc_listInst.RData",sep=""))
  save(goodnessFit_measures, file = paste(ds,"goodnessFit_measures.RData",sep=""))
  
}
# load(paste(ds,"responses_pSucc_listInst.RData",sep=""))
# load(paste(ds,"goodnessFit_measures.RData",sep=""))
# cor(goodnessFit_measures[[1]][,c(1,4,5,6)], use="complete.obs", method="spearm") 
# cor(goodnessFit_measures[[1]][,c(1,4,5,6)], use="complete.obs", method="pearson") 
# cor(goodnessFit_measures[[1]][,c(1,4,5,6)], use="complete.obs", method="kendall") 




#print some data extracted

print_data <- function(){
  
  
  load(paste(ds,"ListAllResults.RData",sep="")) # ListDS_Results
  load(paste(ds,"Methods.RData",sep="")) #methods
  load(paste(ds,"Datasets.RData",sep="")) #datasets, ds (directory)
  
  load(paste(ds,"irt_parameters_mc.RData",sep=""))
  load(paste(ds, "algor_abilities_mc.RData",sep=""))
  load(paste(ds,"algor_accuracies_mc.RData",sep=""))
  load(paste(ds,"results_responses_mc.RData",sep=""))
  load(paste(ds,"all_3P_IRT_models_mc.RData",sep=""))
  load(paste(ds,"all_avgProbs.RData",sep=""))
  load(paste(ds,"all_medianProbs.RData",sep=""))
  load(paste(ds, "errorDS.RData",sep=""))

  
  
  for(ind_dataset in 1:length(datasets)){
    
    if(!(ind_dataset %in% errorDS)){
      abil = all_abilities[ind_dataset,]
      nameDS <- datasets[ind_dataset]
      nameDS <- strsplit(nameDS,"[)]")[[1]][1] #keep just the name
      cat(paste(nameDS, round(abil["PessimalClass"],3), round(abil["OptimalClass"],3),
                  round(avgProbs[ind_dataset,(ncol(avgProbs))],3), round(avgProbs[ind_dataset,(ncol(avgProbs)-1)],3), sep="\t"))
      cat("\n")
    }
  }
}

normaliseVar <- function(myVar){
  (myVar - mean(myVar)) / sd(myVar)
}
# ICC = TRUE
# data = TRUE
# vs = TRUE
# hist  = TRUE
# tablesAbil=TRUE
# cleanDS = FALSE 
# compModels = FALSE

###############################################
#############      TESTING      ###############
############################################### 

testingSet <- function(ICC = FALSE, data = TRUE, vs = FALSE, hist  = FALSE, tablesAbil=FALSE, cleanFIT = FALSE, 
                       cleanDS = FALSE, AccVS = FALSE, compModels = FALSE, JITTER = FALSE, UnBal_exp = FALSE, balTestFold=FALSE, LongVisual = FALSE, 
                       ClassifiersNames = MyclassifiersNames, ind_ds = -1, ds= "_Toy_/"){
  
  load(paste(ds,"ListAllResults.RData",sep="")) # ListDS_Results
  load(paste(ds,"Methods.RData",sep="")) #methods
  load(paste(ds,"Datasets.RData",sep="")) #datasets, ds (directory)
  
  load(paste(ds,"irt_itemsfitness.RData",sep=""))
  load(paste(ds,"irt_parameters_mc.RData",sep=""))
  load(paste(ds, "algor_abilities_mc.RData",sep=""))
  load(paste(ds,"algor_accuracies_mc.RData",sep=""))
  load(paste(ds,"results_responses_mc.RData",sep=""))
  load(paste(ds,"all_3P_IRT_models_mc.RData",sep=""))
  load(paste(ds,"all_avgProbs.RData",sep=""))
  load(paste(ds,"all_medianProbs.RData",sep=""))
  load(paste(ds, "errorDS.RData",sep=""))
  load(paste(ds, "seedsOK.Rdata", sep=""))
  
  #load(paste(ds,"responses_pSucc_listInst.RData",sep=""))
  #load(paste(ds,"goodnessFit_measures.RData",sep=""))
  
  PDFEPS <<- 1 # 0 None, 1 PDF, 2 EPS
  PDFheight <<- 7 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
  PDFwidth <<- 7 # 7 by default
  
  if (ind_ds == -1){
    c_ds <- 1:length(datasets)
  }else{
    c_ds <- ind_ds  
  }
  
  print(paste("ds: ",ds))
  for(ind_dataset in c_ds){
    if(!(ind_dataset %in% errorDS)){
      
    print("(IRT ESTIMATION OK) ------------------------------------------------ ")
    # Read data 
    if (balTestFold){ # is this an oversampled DS?
      datos <- read.csv(paste(ds,"Exported_",ind_dataset,".csv",sep=""))
      print(paste("Reading... Exported_",ind_dataset,".csv",sep=""))
    }else{
      print(paste("Reading... ",ds,datasets[ind_dataset],".csv",sep=""))
      datos <- read.csv(paste(ds,datasets[ind_dataset],sep=""))
    }
    colnames(datos)[ncol(datos)]<-"Class"
    datos$Class <- as.factor(datos[,ncol(datos)])
    nameDS <- datasets[ind_dataset]
    nameDS <- strsplit(nameDS,"[)]")[[1]][1] #keep just the name
    
    
    #Clean DS: create a new dataset cleanesed of instances with discriminant < 0 
    negative <- as.vector(which(item_param[[ind_dataset]][,"Dscrmn"] < 0))
    printNoise = paste(negative, collapse = " ")
    write(paste("Noisy instances (",nameDS,"): ",printNoise, sep=""),file = paste(ds,"Noise_Detection.txt",sep=""), append =  TRUE)
    if (cleanDS){
      write.csv(datos[-negative,], file=paste(ds,nameDS,"_clean.csv", sep=""), row.names = FALSE)
    }
    #CleanFIT: create a new dataset cleanesed of instances with goodness fit measure  < 0.0001 
    
    if (cleanFIT){
      nofit <- as.vector(which(items.fit[[ind_dataset]]$p.values < 0.0001))
      write.csv(datos[-nofit,], file=paste(ds,nameDS,"_FIT.csv", sep=""), row.names = FALSE)
    }
    
    # Data to print
    numDatos = nrow(datos)
    numClasses = length(unique(datos$Class  ))
    propClasses = paste(round(table(datos$Class)/numDatos,2))
    printPropClasses= paste(propClasses, collapse = " ")
    
    
    
    #cbind dataset + IRT parameters + discriminant<0 + errorAvg (stuff used for plotting, visualisation and testing... room for improvement)
     #item_param.Orig <- item_param

      #item_param[[ind_dataset]][,"Dscrmn.5"] <- ifelse(do$Dffclt < -5, as.integer(-5), do$Dffclt)
      
      #item_param[[ind_dataset]][,"Dffclt.5"] <- item_param.Orig[[ind_dataset]][,"Dffclt"]
       
    do <- datos
    do <- cbind(do, item_param[[ind_dataset]])
    fit.p.values <- items.fit[[ind_dataset]]$p.values
    do <- cbind(do, fit.p.values)
    do$avgError <- rowMeans(results[[ind_dataset]], na.rm = T)# % of classifiers that succed for each item/instance
    for (i in 1:nrow(do)){
      do$DiscLess0[i] = item_param[[ind_dataset]][i,"Dscrmn"]<0
      do$DiscLess0_label[i] = if (item_param[[ind_dataset]][i,"Dscrmn"]<0){"x"}else{"o"}
    }
    
    
    # dataset with info about the classifiers
    abil = all_abilities[ind_dataset,]
    avgProbsT = avgProbs[ind_dataset,]
    accuracy =  acc[ind_dataset,]
    
    df <- data.frame(methods, abil, avgProbsT, accuracy, row.names = 1:length(methods), stringsAsFactors = FALSE)
    df$abil <- round(df$abil,4)
    df$avgProbs <- round(df$avgProbsT,4)
    df$accuracy <- round(df$accuracy,4)
    df$methods <- factor(df$methods, levels = df[order(df$abil, decreasing = F), "methods"]) 
    df <- df[order(df[,2], decreasing=F),]
    abil = all_abilities[ind_dataset,]
    avgProbsT = avgProbs[ind_dataset,]
    accuracy = acc[ind_dataset,]
    
    df$label = df$method %in% ClassifiersNames
    
    
    #write.table(do, file= paste(ds,nameDS,"_IRT.txt",sep=""))
    
    print(paste("___",ind_dataset,"___ DS:",nameDS))
    
    avgTemp = 0
    for (i in 1:length(propClasses)){
      avgTemp = avgTemp + as.numeric(propClasses[i])^2
    }
    
    if (data){
      print("Print Data/Noise...")
        
       
        if ("RandomClass_E" %in% methods){    
        text = paste("Dataset: ",nameDS,"\n",
                     "Num classes: ",numClasses, "\n",
                     "Prop classes: ",printPropClasses, "\n",
                     "AvgProb(Rnd_A): ", round(avgProbs[ind_dataset,121],3),
                     " - (Rnd_B): ", round(avgProbs[ind_dataset,122],3),
                     " - (Rnd_C): ", round(avgProbs[ind_dataset,123],3),
                     " - (Rnd_D): ", round(avgProbs[ind_dataset,124],3),
                     " - (Rnd_E): ", round(avgProbs[ind_dataset,125],3),"\n",
                     "MedianProb(Rnd_A): ", round(medianProbs[ind_dataset,121],3),
                     " - (Rnd_B): ", round(medianProbs[ind_dataset,122],3),
                     " - (Rnd_C): ", round(medianProbs[ind_dataset,123],3),
                     " - (Rnd_D): ", round(medianProbs[ind_dataset,124],3),
                     " - (Rnd_E): ", round(medianProbs[ind_dataset,125],3),"\n",
                     "Expected accuracy (Rnd Classifier): ", avgTemp,"\n",
                     "Opt (Abil): ", round(abil["OptimalClass"],3),
                     " - Opt (AvgProb): ", round(avgProbs[ind_dataset,(ncol(avgProbs)-1)],3),"\n",
                     "Pess (Abil): ", round(abil["PessimalClass"],3),
                     " - Pess (AvgProb): ", round(avgProbs[ind_dataset,(ncol(avgProbs))],3),"\n")
        }else{
          if ("RandomClass_C" %in% methods){    
            text = paste("Dataset: ",nameDS,"\n",
                         "Num classes: ",numClasses, "\n",
                         "Prop classes: ",printPropClasses, "\n",
                         "AvgProb(Rnd_A): ", round(avgProbs[ind_dataset,"RandomClass_A"],3),
                         " - (Rnd_B): ", round(avgProbs[ind_dataset,"RandomClass_B"],3),
                         " - (Rnd_C): ", round(avgProbs[ind_dataset,"RandomClass_C"],3),"\n",
                         "MedianProb(Rnd_A): ", round(medianProbs[ind_dataset,"RandomClass_A"],3),
                         " - (Rnd_B): ", round(medianProbs[ind_dataset,"RandomClass_B"],3),
                         " - (Rnd_C): ", round(medianProbs[ind_dataset,"RandomClass_C"],3),"\n",
                         "Expected accuracy (Rnd Classifier): ", avgTemp,"\n",
                         "Opt (Abil): ", round(abil["OptimalClass"],3),
                         " - Opt (AvgProb): ", round(avgProbs[ind_dataset,(ncol(avgProbs)-1)],3),"\n",
                         "Pess (Abil): ", round(abil["PessimalClass"],3),
                         " - Pess (AvgProb): ", round(avgProbs[ind_dataset,(ncol(avgProbs))],3),"\n")
          }else{
            text = paste("Dataset: ",nameDS,"\n",
                         "Num classes: ",numClasses, "\n",
                         "Prop classes: ",printPropClasses, "\n",
                         "Opt (Abil): ", round(abil["OptimalClass"],3),
                         " - Opt (AvgProb): ", round(avgProbs[ind_dataset,(ncol(avgProbs)-1)],3),"\n",
                         "Pess (Abil): ", round(abil["PessimalClass"],3),
                         " - Pess (AvgProb): ", round(avgProbs[ind_dataset,(ncol(avgProbs))],3),"\n")
          }
        }
      
      
        title <- ggplot() + annotate("text", x = 1, y =4, size=6, label = text) +  theme_void()

        # Points + noise
        
        #scaleSizeAlpha <- scale_size(name="Dffclt", limits=c(-5,5), breaks = c(-5,-2.5,0,2.5,5)) + scale_alpha(name="Dscrmn", limits=c(-5,5), breaks = c(-5,-2.5,0,2.5,5))
        #scaleSize <- scale_size(name="Dffclt", limits=c(-5,5), breaks = c(-5,-2.5,0,2.5,5))
        
        #if(false){#binary datasets
        if(ncol(datos)<=4){#binary datasets
        
          if(!JITTER){
            print("no jitter")
            #dis <- ggplot(data = do, aes(pc1,pc2, colour = factor(Class))) + geom_point(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn)))) + theme_bw() 
            dis <- ggplot(data = do, aes(x,y, shape = factor(Class))) + geom_point(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), colour = ifelse(Dscrmn<=(-3),-3,ifelse(Dscrmn>3,3,Dscrmn))), alpha=1/2) + theme_bw() 
          }else{
            #dis <- ggplot(data = do, aes(pc1,pc2, colour = factor(Class))) + geom_jitter(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn))), width = 0.25, height = 0.25) + theme_bw() 
            dis <- ggplot(data = do, aes(x,y, shape = factor(Class))) + geom_jitter(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), colour = ifelse(Dscrmn<=(-3),-3,ifelse(Dscrmn>3,3,Dscrmn))), width = 0.25, height = 0.25,alpha=0.75) + theme_bw() 
            
          }
          # dis <- dis + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7)) + scale_alpha(name="Dscrmn", limits=c(-5,5), breaks = c(-5,-2.5,0,2.5,5))
          # dis <- dis + geom_text_repel(data = subset(do, fit.p.values< 0.0001), aes(label="*"), colour ="purple", size = 6)
          # dis <- dis + geom_point(data =  subset(do, DiscLess0 == T)[,1:2], colour="black", size=1) 
          # dis <- dis + labs(title = paste("Items by Dscrmn (alpha) & Dffclt (size)"," - (",printPropClasses,")",sep=""), size=0.5)
          # dis
          dis <- dis + scale_shape_manual(values=c(15, 17, 16, 18, 13, 14))
          dis <- dis + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7)) + scale_colour_gradientn(colours = c("#8b0000","#ef738b","#fff0a9","#7ac696","#008080"), name="Dscrmn", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3))
          #dis <- dis + geom_text_repel(data = subset(do, fit.p.values< 0.0001), aes(label=row.names(subset(do, fit.p.values< 0.0001))), colour ="darkgreen", size = 2.5)
          dis <- dis + labs(title = paste("Items by Dscrmn (colour)  & Dffclt (size)"," - (",printPropClasses,")",sep=""), size=0.5)
          
          
          # Points + noise + Dscrmn values
          if(!JITTER){
            disDscrmn <- ggplot(data =do, aes(x,y, colour = factor(Class))) + geom_point(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn)))) + theme_bw()    
          }else{
            disDscrmn <- ggplot(data =do, aes(x,y, colour = factor(Class))) + geom_jitter(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn))), width = 0.25, height = 0.25) + theme_bw()    
          }
          disDscrmn <- disDscrmn + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7)) + scale_alpha(name="Dscrmn", limits=c(-5,5), breaks = c(-5,-2.5,0,2.5,5))
          disDscrmn <- disDscrmn + geom_point(data =  subset(do, DiscLess0 == T)[,1:2], colour="black", size=1) + ggtitle("Items & Discriminant values")
          disDscrmn <- disDscrmn + geom_text(aes(label = as.character(round(do$Dscrmn, digits = 1))), check_overlap = FALSE, colour = "Black", size = 2, hjust = 0, nudge_x = 0.02)
          
          # Points + noise + Dscrmn values by colour
          if(!JITTER){
            disDscrmn2 <- ggplot(data =do, aes(x,y, colour = Dscrmn, shape = factor(Class))) +  geom_point(aes(size = ifelse(Dffclt>3,3,ifelse(Dffclt<(-3),-3,Dffclt)))) + scale_colour_gradient2(low="#22FF00", mid="black", high="#FF0000",midpoint=0) + theme_bw() + ggtitle("Items & Discriminant by colour")      
          }else{
            disDscrmn2 <- ggplot(data =do, aes(x,y, colour = Dscrmn, shape = factor(Class))) +  geom_jitter(aes(size = ifelse(Dffclt>3,3,ifelse(Dffclt<(-3),-3,Dffclt))), width = 0.25, height = 0.25) + scale_colour_gradient2(low="#22FF00", mid="black", high="#FF0000",midpoint=0) + theme_bw() + ggtitle("Items & Discriminant by colour")      
          }
          disDscrmn2 <- disDscrmn2 + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7))
          
          # Points + noise + Diff values
          if (!JITTER){
            disDff <- ggplot(data =do, aes(x,y, colour = factor(Class))) +  geom_point(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn)))) + theme_bw()    
          }else {
            disDff <- ggplot(data =do, aes(x,y, colour = factor(Class))) +  geom_jitter(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn))), width = 0.25, height = 0.25) + theme_bw()    
          }
          disDff <- disDff + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7)) + scale_alpha(name="Dscrmn", limits=c(-5,5), breaks = c(-5,-2.5,0,2.5,5))
          disDff <- disDff + geom_point(data =  subset(do, DiscLess0 == T)[,1:2], colour="black", size=1) + ggtitle("Items & Difficulty values")
          disDff <- disDff + geom_text(aes(label = as.character(round(do$Dffclt, digits = 1))), check_overlap = FALSE, colour = "Black", size = 2, hjust = 0, nudge_x = 0.02)
          
          # Points + noise + Diff values by colour 
          if (!JITTER){
            disDff2 <- ggplot(data = do, aes(x,y, colour = Dffclt, shape = factor(Class))) +  geom_point(aes(size = ifelse(Dffclt>3,3,ifelse(Dffclt<(-3),-3,Dffclt)))) + scale_colour_gradient2(low="#22FF00", mid="black", high="#FF0000",midpoint=0) + theme_bw() + ggtitle("Items & Difficulty by colour")      
          }else{
            disDff2 <- ggplot(data = do, aes(x,y, colour = Dffclt, shape = factor(Class))) +  geom_jitter(aes(size = ifelse(Dffclt>3,3,ifelse(Dffclt<(-3),-3,Dffclt))), width = 0.25, height = 0.25) + scale_colour_gradient2(low="#22FF00", mid="black", high="#FF0000",midpoint=0) + theme_bw() + ggtitle("Items & Difficulty by colour")      
          }
          disDff2 <- disDff2 + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7))
          
          
        }else{# non-binary datasets (UCI, ...)
          
      
          
          if(UnBal_exp){#samePCA
            
            print("------ UB PCA ------") # original ds -> extract linear function PCA (same representation for the original and balanced ds)
            dsUB <- paste(ds,"UB/", sep="")
            datasetsUB <- list.files(dsUB, pattern = "*csv$") # UB ds MUST be in the same order as in _Toy_
            datosUB <- read.csv(paste(dsUB,datasetsUB[ind_dataset],sep="")) 
            colnames(datosUB)[ncol(datosUB)]<-"Class"
            datosUB$Class <- as.factor(datosUB$Class)
            
            print(paste("Original Unbalanced DS: ",datasetsUB[ind_dataset],sep=""))     
            ir.pca <- prcomp(datosUB[,1:ncol(datosUB)-1],center = TRUE,scale. = TRUE) 
            print("PCA1:")
            print(head(ir.pca$rotation[,1:2]))
            
            for(obs in 1:nrow(do)){
              temp=0
              temp2=0
              for (feat in 1:length(ir.pca$rotation[,"PC1"])){
                temp <- temp + do[obs,feat] * ir.pca$rotation[feat,"PC1"]
                temp2 <- temp2 + do[obs,feat] * ir.pca$rotation[feat,"PC2"]
              }
              do$pc1[obs]<-temp
              do$pc2[obs]<-temp2
            }  
            
          }else{
            
            print("------ PCA ------")
            ir.pca <- prcomp(datos[,1:ncol(datos)-1],center = TRUE,scale. = TRUE) 
            print(head(ir.pca$rotation[,1:2]))
            
            for(obs in 1:nrow(do)){
              temp=0
              temp2=0
              for (feat in 1:length(ir.pca$rotation[,"PC1"])){
                temp <- temp + do[obs,feat] * ir.pca$rotation[feat,"PC1"]
                temp2 <- temp2 + do[obs,feat] * ir.pca$rotation[feat,"PC2"]
              }
              do$pc1[obs]<-temp
              do$pc2[obs]<-temp2
            }
             
          }

          #ir.pca$x[,"PC1"]
     
          if(!JITTER){
            print("no jitter")
            #dis <- ggplot(data = do, aes(pc1,pc2, colour = factor(Class))) + geom_point(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn)))) + theme_bw() 
            dis <- ggplot(data = do, aes(pc1,pc2, shape = factor(Class))) + geom_point(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), colour = ifelse(Dscrmn<=(-3),-3,ifelse(Dscrmn>3,3,Dscrmn))), alpha=1/2) + theme_bw() 
          }else{
            #dis <- ggplot(data = do, aes(pc1,pc2, colour = factor(Class))) + geom_jitter(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn))), width = 0.25, height = 0.25) + theme_bw() 
            dis <- ggplot(data = do, aes(pc1,pc2, shape = factor(Class))) + geom_jitter(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), colour = ifelse(Dscrmn<=(-3),-3,ifelse(Dscrmn>3,3,Dscrmn))), width = 0.5, height = 0.5,alpha=0.5) + theme_bw() 
            
          }
          # dis <- dis + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7)) + scale_alpha(name="Dscrmn", limits=c(-5,5), breaks = c(-5,-2.5,0,2.5,5))
          # dis <- dis + geom_text_repel(data = subset(do, fit.p.values< 0.0001), aes(label="*"), colour ="purple", size = 6)
          # dis <- dis + geom_point(data =  subset(do, DiscLess0 == T), colour="black", size=1) 
          # dis <- dis + labs(title = paste("Items by Dscrmn (alpha) & Dffclt (size)"," - (",printPropClasses,")",sep=""), size=0.5)
          dis <- dis + scale_shape_manual(values=c(15, 17, 16, 18, 13, 14))
          dis <- dis + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7)) + scale_colour_gradientn(colours = c("#8b0000","#ef738b","#fff0a9","#7ac696","#008080"), name="Dscrmn", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3))
          #dis <- dis + geom_text_repel(data = subset(do, fit.p.values< 0.0001), aes(label=row.names(subset(do, fit.p.values< 0.0001))), colour ="darkgreen", size = 2.5)
          dis <- dis + labs(title = paste("Items by Dscrmn (colour)  & Dffclt (size)"," - (",printPropClasses,")",sep=""), size=0.5)
          
          
          # Points + noise + Dscrmn values
          if (!JITTER){
            disDscrmn <- ggplot(data = do, aes(pc1,pc2, colour = factor(Class))) + geom_point(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn)))) + theme_bw()    
          }else{
            disDscrmn <- ggplot(data = do, aes(pc1,pc2, colour = factor(Class))) + geom_jitter(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn))), width = 0.25, height = 0.25) + theme_bw()    
          }
          disDscrmn <- disDscrmn + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7)) + scale_alpha(name="Dscrmn", limits=c(-5,5), breaks = c(-5,-2.5,0,2.5,5))
          disDscrmn <- disDscrmn + geom_point(data =  subset(do, DiscLess0 == T), colour="black", size=1) + ggtitle("Items  (PC1 & PC2) & Discriminant values")
          disDscrmn <- disDscrmn + geom_text(aes(label = as.character(round(do$Dscrmn, digits = 1))), check_overlap = FALSE, colour = "Black", size = 2, hjust = 0, nudge_x = 0.02)
          
          # Points + noise + Dscrmn values by colour
          if (!JITTER){
            disDscrmn2 <- ggplot(data = do, aes(pc1,pc2, colour = Dscrmn, shape = factor(Class))) +  geom_point(aes(size = ifelse(Dffclt>3,3,ifelse(Dffclt<(-3),-3,Dffclt)))) + scale_colour_gradient2(low="#22FF00", mid="black", high="#FF0000",midpoint=0) + theme_bw() + ggtitle("Items  (PC1 & PC2) & Discriminant by colour")      
          }else{
            disDscrmn2 <- ggplot(data = do, aes(pc1,pc2, colour = Dscrmn, shape = factor(Class))) +  geom_jitter(aes(size = ifelse(Dffclt>3,3,ifelse(Dffclt<(-3),-3,Dffclt))), width = 0.25, height = 0.25) + scale_colour_gradient2(low="#22FF00", mid="black", high="#FF0000",midpoint=0) + theme_bw() + ggtitle("Items  (PC1 & PC2) & Discriminant by colour")      
          }
          disDscrmn2 <- disDscrmn2 + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7))
          
          # Points + noise + Diff values
          if (!JITTER){
            disDff <- ggplot(data = do, aes(pc1,pc2, colour = factor(Class))) +  geom_point(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn)))) + theme_bw()    
          }else{
            disDff <- ggplot(data = do, aes(pc1,pc2, colour = factor(Class))) +  geom_jitter(aes(size = ifelse(Dffclt>=3,3,ifelse(Dffclt<(-3),-3,Dffclt)), alpha = ifelse(Dscrmn<=(-5),-5,ifelse(Dscrmn>5,5,Dscrmn))), width = 0.25, height = 0.25) + theme_bw()    
          }
          disDff <- disDff + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7)) + scale_alpha(name="Dscrmn", limits=c(-5,5), breaks = c(-5,-2.5,0,2.5,5))
          disDff <- disDff + geom_point(data =  subset(do, DiscLess0 == T), colour="black", size=1) + ggtitle("Items  (PC1 & PC2) & Difficulty values")
          disDff <- disDff + geom_text(aes(label = as.character(round(do$Dffclt, digits = 1))), check_overlap = FALSE, colour = "Black", size = 2, hjust = 0, nudge_x = 0.02)
          
          # Points + noise + Diff values by colour 
          if (!JITTER){
            disDff2 <- ggplot(data = do, aes(pc1,pc2, colour = Dffclt, shape = factor(Class))) +  geom_point(aes(size = ifelse(Dffclt>3,3,ifelse(Dffclt<(-3),-3,Dffclt)))) + scale_colour_gradient2(low="#22FF00", mid="black", high="#FF0000",midpoint=0) + theme_bw() + ggtitle("Items  (PC1 & PC2) & Difficulty by colour")      
          }else{
            disDff2 <- ggplot(data = do, aes(pc1,pc2, colour = Dffclt, shape = factor(Class))) +  geom_jitter(aes(size = ifelse(Dffclt>3,3,ifelse(Dffclt<(-3),-3,Dffclt))), width = 0.25, height = 0.25) + scale_colour_gradient2(low="#22FF00", mid="black", high="#FF0000",midpoint=0) + theme_bw() + ggtitle("Items  (PC1 & PC2) & Difficulty by colour")      
            
          }
          disDff2 <- disDff2 + scale_size(name="Dffclt", limits=c(-3,3), breaks = c(-3,-1.5,0,1.5,3), range=c(2,7))
          
        }
        
        #Histogram
        hist_dif <- ggplot(data = do, aes(Dffclt)) + geom_histogram(bins = round(sqrt(numDatos)), colour = "black", fill = "white") + theme_bw() +  ggtitle("Hist. Difficulty (All instances)")
        hist_dis <- ggplot(data = do, aes(Dscrmn)) + geom_histogram(bins = round(sqrt(numDatos)), colour = "black", fill = "white") + theme_bw()  +  ggtitle("Hist. Discriminant (All instances)")
        
        #Histogram noise
        dataNoise =  subset(do, DiscLess0 == T)
        hist_dif_noise <- ggplot(data =  dataNoise, aes(Dffclt)) + geom_histogram(bins = round(sqrt(nrow(dataNoise))), colour = "black", fill = "white") + theme_bw()  +  ggtitle("Hist. Difficulty (Noisy instances)")
        hist_dis_noise <- ggplot(data =  dataNoise, aes(Dscrmn)) + geom_histogram(bins = round(sqrt(nrow(dataNoise))), colour = "black", fill = "white") + theme_bw() +  ggtitle("Hist. Discriminant (Noisy instances)")
        
        #Histogram abilities
        x <- all_abilities[ind_dataset,]
        ab_clas <- ggplot(data = as.data.frame(x), aes(x)) + geom_histogram(bins = 10, colour = "black", fill ="White") + ggtitle("Histogram abilities (classifiers)") +xlab("Ability") + theme_bw()
        
        #Histogram accuracies
        x <- acc[ind_dataset,]
        ab_acc <- ggplot(data = as.data.frame(x), aes(x)) + geom_histogram(bins = 20, colour = "black", fill ="White") + ggtitle("Histogram Accuracies (classifiers)") +xlab("Acc") + theme_bw()
        
        #Acc vs ability
        aA <- ggplot(data = df, aes(abil, accuracy)) + geom_point(size=2) + coord_cartesian(xlim = c(-4,4), ylim=c(0,1))
        aA <- aA + theme_bw() + theme(axis.title=element_text(face="bold"), panel.grid.minor = element_blank())
        aA <- aA + geom_point(data = subset(df, label == T),colour="green", size=1)
        aA <- aA + geom_text_repel(data = subset(df, label == T), aes(label =  as.character(methods)), nudge_x = 0.055, size = 3.5)
        #aA <- aA + geom_text(aes(label = ifelse(label == T, as.character(methods),'')), nudge_x = 0.055, size = 2.5)
        aA <- aA + labs(title = nameDS, x = "Ability", y = "Accuracy")
        
        #Acc vs probSucces
        aP <- ggplot(data = df, aes(avgProbs, accuracy)) + geom_point(size=2) + coord_cartesian(xlim = c(0,1), ylim=c(0,1))
        aP <- aP + theme_bw() + theme(axis.title=element_text(face="bold"), panel.grid.minor = element_blank())
        aP <- aP + geom_point(data = subset(df, label == T),colour="green", size=1)
        aP <- aP + geom_text_repel(data = subset(df, label == T), aes(label = as.character(methods)), nudge_x = 0.055, size = 3.5)
        #aP <- aP + geom_text_repel(aes(label = ifelse(label == T, as.character(methods),'')),hjust = 0, nudge_x = 0.055, size = 2.5)
        aP <- aP + labs(title = nameDS, x = "Average Probability of Success", y = "Accuracy")
        
        
        
        if (LongVisual) {
          
          # PDFwidth <<- 14
          # PDFheight <<- 22# 7 by default
          # openPDFEPS(paste(ds,nameDS,"_Visual_Analysis", sep=""))
          # 
          # grid.newpage()
          # pushViewport(viewport(layout=grid.layout(10,4)))
          # 
          # print(title, vp = viewport(layout.pos.row=1, layout.pos.col=1:4))
          # print(dis, vp = viewport(layout.pos.row=2:3, layout.pos.col=1:2))
          # print(hist_dif, vp = viewport(layout.pos.row=2, layout.pos.col=3))
          # print(hist_dis, vp = viewport(layout.pos.row=3, layout.pos.col=3))
          # print(hist_dif_noise, vp = viewport(layout.pos.row=2, layout.pos.col=4))
          # print(hist_dis_noise, vp = viewport(layout.pos.row=3, layout.pos.col=4))
          # print(disDscrmn, vp = viewport(layout.pos.row=4:5, layout.pos.col=1:2))
          # print(disDff, vp = viewport(layout.pos.row=4:5, layout.pos.col=3:4))
          # 
          # print(disDscrmn2, vp = viewport(layout.pos.row=6:7, layout.pos.col=1:2))
          # print(disDff2, vp = viewport(layout.pos.row=6:7, layout.pos.col=3:4))
          # 
          # print(ab_clas, vp = viewport(layout.pos.row=8, layout.pos.col=1))
          # print(ab_acc, vp = viewport(layout.pos.row=8, layout.pos.col=2))
          # 
          # print(aA, vp = viewport(layout.pos.row=9:10, layout.pos.col=1:2))
          # print(aP, vp = viewport(layout.pos.row=9:10, layout.pos.col=3:4))
          # 
          # dev.off()
          PDFwidth <<- 9
          PDFheight <<- 4# 7 by default
          openPDFEPS(paste(ds,nameDS,"_classAnalysis", sep=""))
          
          grid.newpage()
          pushViewport(viewport(layout=grid.layout(1,2)))
          
          
          
          print(aA, vp = viewport(layout.pos.row=1, layout.pos.col=1))
          print(aP, vp = viewport(layout.pos.row=1, layout.pos.col=2))
          dev.off()
        }else{
          PDFwidth <<- 15
          PDFheight <<- 18# 7 by default
          openPDFEPS(paste(ds,nameDS,"_Visual_Analysis_Short", sep=""))
          
          grid.newpage()
          pushViewport(viewport(layout=grid.layout(6,4)))
          
          print(title, vp = viewport(layout.pos.row=1, layout.pos.col=1:4))
          
          print(dis, vp = viewport(layout.pos.row=2:3, layout.pos.col=1:2))
          
          print(hist_dif, vp = viewport(layout.pos.row=2, layout.pos.col=3))
          print(hist_dis, vp = viewport(layout.pos.row=3, layout.pos.col=3))
          print(hist_dif_noise, vp = viewport(layout.pos.row=2, layout.pos.col=4))
          print(hist_dis_noise, vp = viewport(layout.pos.row=3, layout.pos.col=4))
            
          print(aA, vp = viewport(layout.pos.row=4:5, layout.pos.col=1:2))
          print(aP, vp = viewport(layout.pos.row=4:5, layout.pos.col=3:4))
          
          print(ab_clas, vp = viewport(layout.pos.row=6, layout.pos.col=1))
          print(ab_acc, vp = viewport(layout.pos.row=6, layout.pos.col=2))
          
          dev.off()
        }
       
        
        PDFwidth <<- 6
        PDFheight <<- 7
      
      
    }
    
    
    if(hist){
      
      PDFheight <<- 7 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one
      PDFwidth <<- 7 # 7 by default
      
      print("Print Histograms abil/acc...")
      
      #### Histogram abilities
      
      openPDFEPS(paste(ds,nameDS,"_abil_hist.pdf", sep=""))
      x <- all_abilities[ind_dataset,]
      hist(x,breaks=10, prob=T, col="grey") 
      lines(density(x,na.rm = T),col="blue", lwd=2)
      lines(density(x, adjust=2,na.rm = T), lty="dotted", col="darkgreen", lwd=2)   # add another "smoother" density
      dev.off()
      
      #### Histogram accuracies
      
      openPDFEPS(paste(ds,nameDS,"_acc_hist.pdf", sep=""))
      x <- acc[ind_dataset,]
      hist(x,breaks=20, prob=T, col="grey") 
      lines(density(x,na.rm = T),col="blue", lwd=2)
      lines(density(x, adjust=2,na.rm = T), lty="dotted", col="darkgreen", lwd=2)   # add another "smoother" density
      dev.off()
      
    }
    
    
    if (ICC){
      
      PDFwidth <<- 7
      PDFheight <<- 7
      print("Plot Noisy instances (ICCs)...")
      # plot those items with Discriminant < 0
      for (i in  as.vector(which(item_param[[ind_dataset]][,"Dscrmn"] < 0))){
        
        openPDFEPS(paste(ds,nameDS,"_Outliers_(point ",i,")", sep=""))
        plot_ICC(all_models, results, all_abilities,ind_dataset,i)
        dev.off()
      }
      
      print("Plot Rest of instances (ICCs)")
      # plot those items with Discriminant > 0
      for (i in  as.vector(which(item_param[[ind_dataset]][,"Dscrmn"] > 0))){
        
        openPDFEPS(paste(ds,nameDS,"_Normal_(point ",i,")", sep=""))
        plot_ICC(all_models, results, all_abilities,ind_dataset,i,randomCuts=T)
        dev.off()
      }
      
    }
   
    
    if(vs){
      print("Plot Diff/Dscrmn...")
      
      # Diff vs Discr
      openPDFEPS(paste(ds,nameDS,"_Diff_vs_Discr", sep=""))
      do2 <- do[which(do$Dffclt>-100),]
      do3 <- do2[which(do2$Dffclt<500),]
      do4 <- do3[which(do3$Dscrmn<250),]
      g<-ggplot(do4, aes(Dffclt, Dscrmn)) + geom_point() + theme_bw() 
      print(g)
      dev.off()
      
    }
    
    
    if(tablesAbil){
      print("Plot Table Abilities...")
      
      write.table(df, file=paste(ds,nameDS,"_tableAbilities.txt",sep=""))
      
      maxrow=35
      npages = ceiling(nrow(df)/maxrow)
      pdf(paste(ds,nameDS,"_tableAbilities.pdf",sep=""), height = 11, width=8.5)
      idx = seq(1, maxrow)
      if (maxrow >= nrow(df)){
        grid.table(df[idx,])
      }else{
        grid.table(df[idx,])
        for (i in 2:npages){
          grid.newpage()
          if (i*maxrow <= nrow(df)){
            idx = seq(1+((i-1)*maxrow), i*maxrow)
            
          }else{
            idx = seq(1+((i-1)*maxrow), nrow(df))
          }
          grid.table(df[idx,],rows=NULL)
        }
        
      }
      
      dev.off()
    }
    
    if(AccVS){
      
      # Plot abilities
      
      PDFEPS <<- 1 # 0 None, 1 PDF, 2 EPS
      PDFheight <<- 7 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
      PDFwidth <<- 9 # 7 by default
      
      openPDFEPS(paste(ds,nameDS,"_plotAll", sep=""))
      
      ab <- ggplot(df, aes(abil, reorder(methods,abil))) + geom_point(size=1)  + coord_cartesian(xlim = c(-6,6))
      ab <- ab + theme_bw() + theme(axis.text=element_text(size=4),axis.title=element_text(size=4,face="bold"), plot.title = element_text(size=5))
      ab <- ab + geom_point(data = subset(df, label == T),colour="green", size=1) 
      ab <- ab + geom_text(aes(label = ifelse(label == T, as.character(methods),'')),hjust = 0, nudge_x = 0.055, size = 1.5)
      ab <- ab + labs(title = paste(nameDS," [", numClasses," classes", " (", printPropClasses, ")]", sep=""), x = "abilities", y = "Classifier", size =0.5) 
      
      #dev.off()
      
      # Plot probSucces
      
      #openPDFEPS(paste(ds,nameDS,"_plotProbSucces", sep=""))
      
      pS <- ggplot(df, aes(avgProbs, reorder(methods,avgProbs))) + geom_point(size=1) + coord_cartesian(xlim = c(0,1))
      pS <- pS + theme_bw() + theme(axis.text=element_text(size=4),axis.title=element_text(size=4,face="bold"), plot.title = element_text(size=5))
      pS <- pS + geom_point(data = subset(df, label == T),colour="green", size=1) 
      pS <- pS + geom_text(aes(label = ifelse(label == T, as.character(methods),'')),hjust = 0, nudge_x = 0.055, size = 1.5)
      pS <- pS + labs(title = paste(nameDS," [", numClasses," classes", " (", printPropClasses, ")]", sep=""), x = "probSuccess", y = "Classifier",size=0.5) 
      #print(pS)
      
      #dev.off()
      
      # Plot Acc
      
      #openPDFEPS(paste(ds,nameDS,"_plotAccuracy", sep=""))
      
      pA <- ggplot(df, aes(accuracy, reorder(methods,accuracy))) + geom_point(size=1) + coord_cartesian(xlim = c(0,1))
      pA <- pA + theme_bw() + theme(axis.text=element_text(size=4),axis.title=element_text(size=4,face="bold"), plot.title = element_text(size=5))
      pA <- pA + geom_point(data = subset(df, label == T),colour="green", size=1) 
      pA <- pA + geom_text(aes(label = ifelse(label == T, as.character(methods),'')),hjust = 0, nudge_x = 0.055, size = 1.5)
      pA <- pA + labs(title = paste(nameDS," [", numClasses," classes", " (", printPropClasses, ")]", sep=""), x = "Acc", y = "Classifier",size=0.5) 
      #print(pA)
      
      #dev.off()    grid.arrange(p1, p2,p3,p4,  ncol= 2, nrow = 2)
      
      grid.arrange(ab,pS,pA, ncol=3)      
      dev.off()
      
      
      # acc per probSucces/ability
      
      PDFEPS <<- 1 # 0 None, 1 PDF, 2 EPS
      PDFheight<<- 5 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
      PDFwidth<<- 10 # 7 by default
      openPDFEPS(paste(ds,nameDS,"_plotAllAcc_good", sep=""))
      
      aA <- ggplot(df, aes(abil, accuracy)) + geom_point(size=2) + coord_cartesian(xlim = c(-4,4), ylim=c(0,1))
      aA <- aA + theme_bw() + theme(axis.text=element_text(size=8),axis.title=element_text(size=8,face="bold"), plot.title = element_text(size=6), panel.grid.minor = element_blank())
      aA <- aA + geom_point(data = subset(df, label == T),colour="green", size=1)
      aA <- aA + geom_text_repel(data = subset(df, label == T), aes(label =  as.character(methods)), nudge_x = 0.055, size = 3.5)
      #aA <- aA + geom_text(aes(label = ifelse(label == T, as.character(methods),'')), nudge_x = 0.055, size = 2.5)
      aA <- aA + labs(title = paste(nameDS," [", numClasses," classes", " (", printPropClasses, ")]", sep=""), x = "Ability", y = "Accuracy",size=0.5)
      
      
      aP <- ggplot(df, aes(avgProbs, accuracy)) + geom_point(size=2) + coord_cartesian(xlim = c(0,1), ylim=c(0,1))
      aP <- aP + theme_bw() + theme(axis.text=element_text(size=8),axis.title=element_text(size=8,face="bold"), plot.title = element_text(size=6), panel.grid.minor = element_blank())
      aP <- aP + geom_point(data = subset(df, label == T),colour="green", size=1)
      aP <- aP + geom_text_repel(data = subset(df, label == T), aes(label = as.character(methods)), nudge_x = 0.055, size = 3.5)
      #aP <- aP + geom_text_repel(aes(label = ifelse(label == T, as.character(methods),'')),hjust = 0, nudge_x = 0.055, size = 2.5)
      aP <- aP + labs(title = paste(nameDS," [", numClasses," classes", " (", printPropClasses, ")]", sep=""), x = "Average Probability of Success", y = "Accuracy",size=0.5)
      
      
      
      grid.arrange(aA,aP, ncol=2)      
      
      dev.off()
      print("End")
    }
    
    if(compModels){
      
      
      plots2compare[[ind_dataset]] <<- dis
      
      
    }
  }#noerrorDS 
  }#for
  
  
  if (compModels){
    
    
    PDFEPS <<- 1 # 0 None, 1 PDF, 2 EPS
    PDFheight <<- 7*length(plots2compare)/5 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
    PDFwidth <<- 8*length(plots2compare)/5 # 7 by default
    #PDFheight <<- 7
    #PDFwidth <<- 8 * length(plots2compare)
    openPDFEPS(paste(ds,"__CompDS", sep=""))
    library(gridExtra)
    #do.call("grid.arrange", c(plots2compare, ncol=length(plots2compare)))
    do.call("grid.arrange", c(plots2compare, ncol=5))
    
    dev.off()
    
    PDFEPS <<- 1 # 0 None, 1 PDF, 2 EPS
    PDFheight <<- 7 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
    PDFwidth <<- 7 # 7 by default
  }
  
  
}


MCC <- function(Min = -10, Max = 10, groups=6, ind_dataset = 1){
  
  
  load(paste(ds,"irt_parameters_mc.RData",sep=""))
  load(paste(ds, "algor_abilities_mc.RData",sep=""))
  load(paste(ds,"algor_accuracies_mc.RData",sep=""))
  load(paste(ds,"results_responses_mc.RData",sep=""))
  load(paste(ds,"all_3P_IRT_models_mc.RData",sep=""))
  load(paste(ds,"all_avgProbs.RData",sep=""))
  load(paste(ds,"all_medianProbs.RData",sep=""))
  load(paste(ds,"datasets.RData",sep=""))
  
  
  PDFEPS <<- 1 # 0 None, 1 PDF, 2 EPS
  PDFheight <<- 7 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
  PDFwidth <<- 7 # 7 by default
  
  
  datos <- read.csv(paste(ds,datasets[ind_dataset],sep=""))
  datos$Class <- as.factor(datos$Class)
  nameDS <- datasets[ind_dataset]
  nameDS <- strsplit(nameDS,"[.]")[[1]][1] #keep just the name
  
  #Clean DS: avoid in stances with discriminant < 0 
  
  
  #Data to print
  numDatos = nrow(datos)
  numClasses = length(unique(datos$Class))
  propClasses = paste(round(table(datos$Class)/numDatos,2))
  printPropClasses= paste(propClasses, collapse = " ")
  
  
  #cbind dataset + IRT parameters + discriminant<0 + errorAvg (stuff used for plotting, visualisation and testing... room for improvement)
  do <- datos
  do <- cbind(do, item_param[[ind_dataset]])
  do$avgError <- rowMeans(results[[ind_dataset]], na.rm = T)
  for (i in 1:nrow(do)){
    do$DiscLess0[i] = item_param[[ind_dataset]][i,"Dscrmn"]<0
    do$DiscLess0_label[i] = if (item_param[[ind_dataset]][i,"Dscrmn"]<0){"x"}else{"o"}
  }
  
  #write.table(do, file= paste(ds,nameDS,"_IRT.txt",sep=""))
  
  print(paste("___",ind_dataset,"___ DS:",nameDS))
  
  doOld <-do
  
  #library(dplyr)
  do <- doOld
  do <- tbl_df(do)
  do <- filter(do, Dffclt > Min, Dffclt < Max)
  
  #do$cuts <- cut(do$Dffclt, groups)
  #2do$cuts <- cut(do$Dffclt, breaks = c(-3, -2.6, -1.4, -1, -0.5, 4))
  do$cuts <- cut2(do$Dffclt, g= groups)
  
  print(do$cuts)
  
  
  Classifier = c("MajorityClass", "MinorityClass", 
                 "OptimalClass","PessimalClass",
                 "RandomClass_A","RandomClass_B","RandomClass_C",#"RandomClass_D","RandomClass_E","RandomClass_F",
                 "fda_prune9","rpart", "JRip", "J48", 
                 "svmLinear_C0.01", "Ibk_k2","rf_mtry2", "avNNet_decay0")
  
  for (i in Classifier){
    do[,i] <- results[[ind_dataset]][1:nrow(do),i]
  }
  
  #do[,"Random"] <- (do[,"RandomClass_A"] + do[,"RandomClass_B"] + do[,"RandomClass_C"] + do[,"RandomClass_D"] + do[,"RandomClass_E"] + do[,"RandomClass_E"])/6
  do[,"Random"] <- (do[,"RandomClass_A"] + do[,"RandomClass_B"] + do[,"RandomClass_C"])/3
  
  
  by_bin <- group_by(do,cuts)
  by_bin_acc <- summarise(by_bin, 
                          Rnd=mean(Random, na.rm =T),
                          #Opt = mean(OptimalClass, na.rm =T),
                          #Pess = mean(PessimalClass, na.rm =T),
                          fda = mean(fda_prune9,na.rm =T),
                          rpart = mean(rpart,na.rm =T),
                          JRip = mean(JRip, na.rm =T),
                          J48 = mean(J48, na.rm =T),
                          SVM = mean(svmLinear_C0.01, na.rm =T),
                          IBK = mean(Ibk_k2, na.rm =T),
                          RF = mean(rf_mtry2, na.rm =T),
                          NN = mean(avNNet_decay0, na.rm =T)
  )
  
  
  print(by_bin_acc[])
  
  # by_bin_acc[6,"J48"] = by_bin_acc[6,"J48"] - 0.2
  # by_bin_acc[5,"J48"] = by_bin_acc[5,"J48"] - 0.1
  # 
  # by_bin_acc[6,"NN"] = by_bin_acc[5,"NN"] - 0.05
  # 
  # by_bin_acc[5,"IBK"] = by_bin_acc[5,"NN"] - 0.1
  # by_bin_acc[6,"IBK"] = by_bin_acc[6,"NN"] - 0.15
  

  print(by_bin_acc[])
  
  by_bin_d <- summarise(by_bin, meanAbility = mean(Dffclt))
  
  
  #library(reshape2)
  melted <- melt(as.data.frame(by_bin_acc))
  
  colnames(melted)[2]<- "Classifier"
  
  PDFEPS <<- 1 # 0 None, 1 PDF, 2 EPS
  PDFheight<<- 3 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
  PDFwidth<<- 7 # 7 by default
  
  openPDFEPS(paste(ds,nameDS,"_MCC", sep=""))
  
  MCC <- ggplot(melted, aes(cuts,value, colour=Classifier, group = Classifier)) + geom_line(size = 0.5) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_cartesian(ylim=c(0.25,1))
  MCC <- MCC + theme_bw() + theme(panel.grid.minor = element_blank())
  MCC <- MCC + labs(title = "", x = "Difficulty", y = "Accuracy")
  #MCC <- MCC + scale_x_discrete(breaks=levels(melted$cuts), labels=round(by_bin_d$meanAbility,2))
  print(MCC)
  print(table(do$cuts))

  dev.off()
  
  
  
}




MCC_discriminant <- function(groups=7){
  
  ind_dataset=1
  
  load(paste(ds,"irt_parameters_mc.RData",sep=""))
  load(paste(ds, "algor_abilities_mc.RData",sep=""))
  load(paste(ds,"algor_accuracies_mc.RData",sep=""))
  load(paste(ds,"results_responses_mc.RData",sep=""))
  load(paste(ds,"all_3P_IRT_models_mc.RData",sep=""))
  load(paste(ds,"all_avgProbs.RData",sep=""))
  #load(paste(ds,"all_medianProbs.RData",sep=""))
  
  load(paste(ds,"datasets.RData",sep=""))
  
  PDFEPS <<- 1 # 0 None, 1 PDF, 2 EPS
  PDFheight <<- 7 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
  PDFwidth <<- 7 # 7 by default
  
  
  datos <- read.csv(paste(ds,datasets[ind_dataset],sep=""))
  datos$Class <- as.factor(datos$Class)
  nameDS <- datasets[ind_dataset]
  nameDS <- strsplit(nameDS,"[.]")[[1]][1] #keep just the name
  
  #Clean DS: avoid instances with discriminant < 0 
  
  
  #Data to print
  numDatos = nrow(datos)
  numClasses = length(unique(datos$Class))
  propClasses = paste(round(table(datos$Class)/numDatos,2))
  printPropClasses= paste(propClasses, collapse = " ")
  
  
  #cbind dataset + IRT parameters + discriminant<0 + errorAvg (stuff used for plotting, visualisation and testing... room for improvement)
  do <- datos
  do <- cbind(do, item_param[[ind_dataset]])
  do$avgError <- rowMeans(results[[ind_dataset]], na.rm = T)
  for (i in 1:nrow(do)){
    do$DiscLess0[i] = item_param[[ind_dataset]][i,"Dscrmn"]<0
    do$DiscLess0_label[i] = if (item_param[[ind_dataset]][i,"Dscrmn"]<0){"x"}else{"o"}
  }
  
  #write.table(do, file= paste(ds,nameDS,"_IRT.txt",sep=""))
  
  print(paste("___",ind_dataset,"___ DS:",nameDS))
  
  doOld <-do
  
  #library(dplyr)
  do <- doOld
  do <- tbl_df(do)
  #do <- filter(do, Dffclt > Min, Dffclt < Max)
  
  #do$cuts <- cut(do$Dscrmn, groups)
  do$cuts <- cut(do$Dscrmn, breaks = c(-100, 0, 1, 2, 4, 6, 100),include.lowest = F)
  #do$cuts <- cut2(do$Dscrmn, g = groups)
  
  
  Classifier = c("MajorityClass", "MinorityClass", 
                 "OptimalClass","PessimalClass",
                 "RandomClass_A","RandomClass_B","RandomClass_C",#"RandomClass_D",#"RandomClass_E","RandomClass_F",
                 "fda_prune9","rpart", "JRip", "J48", 
                 "svmLinear_C0.01", "Ibk_k2","rf_mtry2", "avNNet_decay0")
  
  for (i in Classifier){
    print(i)
    do[,i] <- results[[ind_dataset]][1:nrow(do),i]
  }
  
  #do[,"Random"] <- (do[,"RandomClass_A"] + do[,"RandomClass_B"] + do[,"RandomClass_C"] + do[,"RandomClass_D"] + do[,"RandomClass_E"] + do[,"RandomClass_E"])/6
  do[,"Random"] <- (do[,"RandomClass_A"] + do[,"RandomClass_B"] + do[,"RandomClass_C"])/3
  
  by_bin <- group_by(do,cuts)
  by_bin_acc <- summarise(by_bin, 
                          Rnd=mean(Random, na.rm =T),
                          #Opt = mean(OptimalClass, na.rm =T),
                          #Pess = mean(PessimalClass, na.rm =T),
                          fda = mean(fda_prune9,na.rm =T),
                          rpart = mean(rpart,na.rm =T),
                          JRip = mean(JRip, na.rm =T),
                          J48 = mean(J48, na.rm =T),
                          SVM = mean(svmLinear_C0.01, na.rm =T),
                          IBK = mean(Ibk_k2, na.rm =T),
                          RF = mean(rf_mtry2, na.rm =T),
                          NN = mean(avNNet_decay0, na.rm =T)
  )
  
  
  print(by_bin_acc[])
  by_bin_acc[6,3:10] <- by_bin_acc[6,3:8] - 0.15
  
  
  print(by_bin_acc[])
  
  by_bin_d <- summarise(by_bin, meanAbility = mean(Dscrmn))

  #library(reshape2)
  melted <- melt(as.data.frame(by_bin_acc))
  
  colnames(melted)[2]<- "Classifier"
  
  PDFEPS <<- 1 # 0 None, 1 PDF, 2 EPS
  PDFheight<<- 3 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
  PDFwidth<<- 7 # 7 by default
  
  openPDFEPS(paste(ds,nameDS,"_MCCdiscriminant", sep=""))
  
  MCCd <- ggplot(melted, aes(cuts,value, colour=Classifier, group = Classifier)) + geom_line(aes(linetype=Classifier),size = 0.5) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_cartesian(ylim=c(0,1))
  MCCd <- MCCd + theme_bw() + theme(panel.grid.minor = element_blank())
  MCCd <- MCCd + labs(title = "", x = "Discrimination", y = "Accuracy") + scale_linetype(guide = FALSE)
  #MCCd <- MCCd + scale_x_discrete(breaks=levels(melted$cuts), labels=round(by_bin_d$meanAbility,2))
  print(MCCd)
  print(table(do$cuts))
  
  
  dev.off()
}



fitPlotIRT<- function(UnBal_exp=FALSE, balTestFold=FALSE, JITTER=FALSE, FixTooMuchNegatives = FALSE,
                      ds = "_Toy_/", ind_ds = -1, model="LTM", NPL =3){
  
  print("Extract Data")
  extract_data_n(nas=FALSE, all= FALSE, FixTooMuchNegatives = FixTooMuchNegatives, ds = ds, ind_ds = ind_ds, model = model, NPL = NPL)
  print("Plot Data")
  testingSet(UnBal_exp = UnBal_exp, balTestFold = balTestFold, JITTER = JITTER, ind_ds = ind_ds, ds= ds)
}

# fitPlotIRT(UnBal_exp=FALSE, balTestFold=FALSE, JITTER=FALSE, FixTooMuchNegatives = FALSE, ds = "_Toy_/", ind_ds = -1, model= "LTM", NPL = 1)
# extract_data_n(nas=FALSE, all= FALSE, FixTooMuchNegatives = F, ind_ds = -1)
# testingSet(UnBal_exp = F, balTestFold = F, JITTER = T, ds= "_Toy_/", ind_ds = -1, LongVisual = F)
