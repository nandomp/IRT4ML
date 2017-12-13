
#abils <- data.frame(Classifier = c("Rnd","fda","rpart","JRip","J48","SVM","IBK","RF","NN"),
#                    Ability = c(-2.54881738, -0.64994470, -0.57730877, -0.27572872, -0.04333066, 1.09320626, 0.38542027, 0.98248523, 1.70960798))


abils <- data.frame(Classifier = c("Rnd","fda","rpart","JRip","J48","SVM","IBK","RF","NN"),
                    Ability = c(-3.54881738, 
                                -0.64994470, 
                                -0.47730877, 
                                -0.27572872, 
                                0.04333066, 
                                1.09320626, 
                                -0.18542027, 
                                0.98248523, 
                                1.70960798))


  #1PL
  myf <- function(x) {1 / (1 + exp(x-abils["Rnd"])) }
  plot(myf,xlim=c(-3.32,0.54), ylim=c(0.25,1.0), col="red",xlab= "Difficulty", ylab= "Probability",lwd=3)
  
  myf <- function(x) {1 / (1 + exp(x-abils["fda_prune9"])) }
  plot(myf,xlim=c(-3.32,0.54), ylim=c(0.25,1.0), col="orange",lwd=3, add = T)
  
  myf <- function(x) {1 / (1 + exp(x-abils["rpart"])) }
  plot(myf,xlim=c(-3.32,0.54), ylim=c(0.25,1.0), col="darkolivegreen4",lwd=3, add = T)
  
  myf <- function(x) {1 / (1 + exp(x-abils["JRip"])) }
  plot(myf,xlim=c(-3.32,0.54), ylim=c(0.25,1.0), col="limegreen",lwd=3, add = T)
  
  myf <- function(x) {1 / (1 + exp(x-abils["J48"])) }
  plot(myf,xlim=c(-3.32,0.54), ylim=c(0.25,1.0), col="mediumturquoise",lwd=3, add = T)
  
  myf <- function(x) {1 / (1 + exp(x-abils["svmLinear_C0.01"])) }
  plot(myf,xlim=c(-3.32,0.54), ylim=c(0.25,1.0), col="lightskyblue",lwd=3, add = T)
  
  #IBK
  myf <- function(x) {1 / (1 + exp(x-abils["Ibk_k2"])) }
  plot(myf,xlim=c(-3.32,0.54), ylim=c(0.25,1.0), col="royalblue",lwd=3, add = T)
  
  myf <- function(x) {1 / (1 + exp(x-abils["rf_mtry2"])) }
  plot(myf,xlim=c(-3.32,0.54), ylim=c(0.25,1.0), col="maroon1",lwd=3, add = T)
  
  myf <- function(x) {1 / (1 + exp(x-abils["avNNet_decay0"])) }
  plot(myf,xlim=c(-3.32,0.54), ylim=c(0.25,1.0), col="purple",lwd=3, add = T)
  
  legend('bottomleft', c("Rnd","fda","rpart","JRip","J48","SVM","IBK","RF","NN"), 
         lty=1, lwd=3,col=c('red', 'orange', 'darkolivegreen4', 'limegreen',' mediumturquoise',' lightskyblue',' royalblue',' purple', 'maroon1'), bty='n', cex=0.5)
  
  
  plotsCCC_diffProb()  
  
plotsCCC_diffProb <- function(model = "3"){
  model = model
  ds = "_Toy_/"
  load(paste(ds,"irt_parameters_mc.RData",sep=""))
  load(paste(ds, "algor_abilities_mc.RData",sep=""))
  load(paste(ds,"algor_accuracies_mc.RData",sep=""))
  load(paste(ds,"results_responses_mc.RData",sep=""))
  load(paste(ds,"all_3P_IRT_models_mc.RData",sep=""))
  load(paste(ds,"all_avgProbs.RData",sep=""))
  load(paste(ds,"all_medianProbs.RData",sep=""))  
  
  #Si puedes calcular para cada dificultad cuál es la discriminación media y el guess medio por 
  # cada valor de dificultad para los datos y con eso sacas cada punto, igual se ve mejor, 
  #porque si no las curvas no se cruzan, que es lo que tiene más gracia.
  
  Classifier = c("MajorityClass", "MinorityClass", 
                 "OptimalClass","PessimalClass",
                 "RandomClass_A","RandomClass_B","RandomClass_C","RandomClass_D","RandomClass_E","RandomClass_F",
                 "fda_prune9","rpart", "JRip", "J48", 
                 "svmLinear_C0.01", "Ibk_k2","rf_mtry2", "avNNet_decay0")
  
  abils<-all_abilities[1,Classifier]
  abils["Rnd"] <- sum(abils[5:9])/6 
  abils.fin <- abils[c("Rnd",
                       "fda_prune9","rpart", "JRip", "J48", 
                       "svmLinear_C0.01", "Ibk_k2","rf_mtry2", "avNNet_decay0")]  
  

  params <- as.data.frame(all_models[[1]]$item_param)
  getA <- function(b){
    t <- filter(params, Dffclt>=(round(b,2)-0.5) & Dffclt<=(round(b,2)+0.5))
    return(mean(t$Dscrmn))
  }
  getC <- function(b){
    t <- filter(params, Dffclt>=(round(b,2)-0.5) & Dffclt<=(round(b,2)+0.5))
    return(mean(t$Gussng))
    
  }
  
  Probability <- function(b,theta, mod = model) {
    if(mod == "3"){
      a =getA(b)
      c = getC(b)
    }else{
      if(mod == "2"){
        a =getA(b)
        c = 0
      }else{
        a=1
        c=0
      }
    }

    #theta = abils["Rnd"]
    return(c + ((1-c)/(1+exp(-a*(theta-b)))))
    #1 / (1 + exp(x-abils["Rnd"])) 
  }
  
  #plot(Probability(b = c(-3.32:0.54), theta = abils["Rnd"]),xlim=c(-3.32,0.54), ylim=c(0.25,1.0), col="red",xlab= "Difficulty", ylab= "Probability",lwd=3)
  
  #plot(Probability(b = seq(-7,7,0.1), theta = abils["Rnd"]))
  c("Rnd",
    "fda_prune9","rpart", "JRip", "J48", 
    "svmLinear_C0.01", "Ibk_k2","rf_mtry2", "avNNet_decay0")
  
  
  df <- data.frame(difficulty=seq(-3.32,0.54,0.001))
  for(i in 1:nrow(df)){
    df[i,"Rnd"] = Probability(df[i,1], abils["Rnd"])
    df[i,"fda"] = Probability(df[i,1], abils["fda_prune9"])
    df[i,"rpart"] = Probability(df[i,1], abils["rpart"])
    df[i,"JRip"] = Probability(df[i,1], abils["JRip"])
    df[i,"J48"] = Probability(df[i,1], abils["J48"])
    df[i,"SVM"] = Probability(df[i,1], abils["svmLinear_C0.01"])
    df[i,"IBK"] = Probability(df[i,1], abils["Ibk_k2"])
    df[i,"RF"] = Probability(df[i,1], abils["rf_mtry2"])
    df[i,"NN"] = Probability(df[i,1], abils["avNNet_decay0"])
    #print(paste("diff: ",i,"y:",y))
    
  }
  
  #df <- df[,!names(df)%in%c("Prob")]
  df.m <- melt(df, id.vars = c("difficulty"))
  
  
  PDFEPS <<- 1 # 0 None, 1 PDF, 2 EPS
  PDFheight<<- 3 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
  PDFwidth<<- 7 # 7 by default
  
  openPDFEPS(paste("DiffProbClassifiers", sep=""))
  
  ggplot(df.m, aes(difficulty, value, colour = variable)) + geom_line(aes(linetype=variable),size = 0.5) + theme_bw() + xlim(c(-3.32,0.6)) +
    xlab("Difficulty") + ylab("Prob. Correct Response") + scale_color_discrete("Classifier") + scale_linetype(guide = FALSE)
  
  dev.off()
  
}
  
  