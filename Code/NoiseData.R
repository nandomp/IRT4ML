options( java.parameters = "-Xmx6g" )

.lib<- c("mlbench","ggplot2","dplyr")
.inst <- .lib %in% installed.packages()
if (length(.lib[!.inst])>0) install.packages(.lib[!.inst], repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com")) 
lapply(.lib, require, character.only=TRUE)

set.seed(288)

PDFEPS <- 1 # 0 None, 1 PDF, 2 EPS
PDFheight= 7 # 7 by default, so 14 makes it double higher than wide, 5 makes letters bigger (in proportion) for just one 
PDFwidth= 9 # 7 by default

# This function is used to generate PDFs or EPSs for the plots
openPDFEPS <- function(file, height= PDFheight, width= PDFwidth) {
  if (PDFEPS == 1) {
    pdf(paste(file, ".pdf", sep=""), width, height)
  } else if (PDFEPS == 2) {
    postscript(paste(file, ".eps", sep=""), width, height, horizontal=FALSE)
  }
}

# type == 1 => cassini 
# type == 2 => 2dnormals

genDS <- function(type=1, examples=1000, noise=100, numC=3){
  
  set.seed(288)

  if (type == 1){
    x <- mlbench.cassini(examples)
    t <- "Cassini"
  }else if(type == 2){
    x <- mlbench.2dnormals(examples,numC)
    t <- "2dNormals"
  }
  
  dataNoise <- data.frame(x=x$x[,1], y = x$x[,2], Class=x$classes)
  dataNoise <- dataNoise[sample(nrow(dataNoise)),]

  for (i in sample(nrow(dataNoise),noise)){
    newClass = sample(unique(dataNoise$Class),1)
    
    while (newClass == dataNoise$Class[i]){
      newClass = sample(unique(dataNoise$Class),1)
    }
    dataNoise$Class[i] <- newClass
  }
  
  openPDFEPS(paste(t,"_",c,"c_", examples,"e_",noise,"n", sep=""))
  #plot(dataNoise$x,dataNoise$y,col=dataNoise$Class, pch=16)
  p<-ggplot(dataNoise, aes(x, y, colour = Class)) + geom_point(size=4) + theme_light()
  print(p)
  dev.off()
  
  write.csv(dataNoise,file=paste(t,"_",c,"c_", examples,"e_",noise,"n.csv",sep=""), row.names = FALSE)
}




genDS_incNoise <- function(type=1, examples=200, numC=3, extra = 1){
  
  if (type == 1){
    x<-mlbench.cassini(examples)
  }else if(type == 2){
    x <- mlbench.2dnormals(examples,numC)
  }
  
  dataNoise <- data.frame(x=x$x[,1], y = x$x[,2], Class=x$classes)
  dataNoise <- dataNoise[sample(nrow(dataNoise)),]
  print(paste("DS: ", nrow(dataNoise)))
  
  noiseSample = sample(nrow(dataNoise),examples/5)
  print(paste("Noise: ", length(noiseSample)))
  
  for (i in seq(0,100,10)){
    up2 = (examples/100 * i)/5
    sample <- noiseSample[0:up2]  
    print(paste("Samples: ",length(sample)))
    
    printsample = paste(sample, collapse = " ")
    write(paste("Noisy instances (",i/5,"%): ",printsample, sep=""),file = "Noise_information.txt", append =  TRUE)
    
    
    for (j in sample){
      newClass = sample(unique(dataNoise$Class),1)
      while (newClass == dataNoise$Class[j]){
        newClass = sample(unique(dataNoise$Class),1)
      }
      dataNoise$Class[j] <- newClass
    }
    
    if (type == 1){t <- "Cassini"}
    if (type == 2){t <- "2dNormals"}
    openPDFEPS(paste(extra,"__",t,"_",c,"c_", examples,"e_",i/5,"n", sep=""))
    
    p<- ggplot(dataNoise,aes(x,y, colour= factor(Class), label= row.names(dataNoise))) + geom_point(size = 5.5) + geom_text(check_overlap = F ,size=4, hjust = 0, nudge_x = 0.055)
    print(p)
    dev.off()
    
    write.csv(dataNoise,file=paste(extra,"__",t,"_",c,"c_", examples,"e_",i/5,"n.csv",sep=""), row.names = FALSE)
  }
    
}
  
run4avg <- function(){
  for(i in 1:10){
    genDS_incNoise(1,200,3,letters[i]) 
  }
  
}
  
  
  
  


runGenDS <- function(type = 1, instances = 200){
  
  for(i in 1:10){
    noise <- i * 5
    genDS(type,instances,noise)
  }
  
}




Noisify <- function(data) {
  
  if (is.vector(data)) {
    noise <- runif(length(data), -0.00001, 0.00001)
    noisified <- data + noise
  } else {
    length <- dim(data)[1] * dim(data)[2]
    noise <- matrix(runif(length, -0.0001, 0.00001), dim(data)[1])
    noisified <- data + noise
  }
  return(noisified)
}


imbalanceDF_donut_reduceProp <- function(nEx = 500){

  set.seed(288)
  
  c=c(rep(1,90),rep(2,10))
  ds=data.frame(x,y,c)
  ggplot(ds,aes(x,y,colour=c))+geom_point()
    x = rnorm(nEx,0,0.5) 
    y = rnorm(nEx,0,0.5)
    df = data.frame(x, y, Class =rep(1,nEx))
    
    for (range in seq(0.25,0,-0.05)){
      
      dfnew <- df
      dfnew[(x < (0 + range) & x > (0 - range) & y < (0 + range) & y > (0 - range)),"Class"] <- "2"
      dfnew$Class <- as.factor(dfnew$Class)
      
      openPDFEPS(paste("ImbalanceDS_range_0",range*100,sep=""))
      p <- ggplot(dfnew,aes(x,y,colour=Class)) + geom_point()
      print(p)
      dev.off()
      write.csv(dfnew,file=paste("ImbalanceDS_range_0",range*100,".csv",sep=""), row.names = FALSE)
    
    }
 
}

imbalanceDF_donut_moveMean <- function(nEx = 500){
  
  set.seed(288)
  
  c=c(rep(1,90),rep(2,10))
  ds=data.frame(x,y,c)
  ggplot(ds,aes(x,y,colour=c))+geom_point()
  x = rnorm(nEx,0,0.5) 
  y = rnorm(nEx,0,0.5)
  df = data.frame(x, y, Class =rep(1,nEx))
  
  for (meanShift in seq(0,1,0.2)){
    
    dfnew <- df
    dfnew[(x < (meanShift + 0.25) & x > (meanShift - 0.25) & y < (meanShift +0.25) & y > (meanShift - 0.25)),"Class"] <- "2"
    dfnew$Class <- as.factor(dfnew$Class)
    
    openPDFEPS(paste("ImbalanceDS_mean_0",meanShift*10,sep=""))
    p <- ggplot(dfnew,aes(x,y,colour=Class)) + geom_point()
    print(p)
    dev.off()
    write.csv(dfnew,file=paste("ImbalanceDS_mean_0",meanShift*10,".csv",sep=""), row.names = FALSE)
    
  }
  
}





normals_imb <- function (n, cl = 2, r = sqrt(cl), sd = 0.5, prob = c(50,50), distance = 1) 
{
  e <- sample(0:(cl - 1), size = n, replace = TRUE, prob)
  m <- r * cbind(cos(pi/4 + e * 2 * pi/cl), sin(pi/4 + e * 2 * pi/cl))
  x <- matrix(rnorm(2 * n, sd = sd), ncol = 2) + m * distance
  retval <- data.frame( x=x[,1], y = x[,2], Class= factor(e))
  #retval <- list(x = x, classes = factor(e + 1))
  #class(retval) <- c("mlbench.2dnormals", "mlbench")
  retval
}




#totally different datasets

multi_normals_imb <- function(numDS, numEx, maxSD = 1, dist = 0.5)
{
  for (s in seq(0.5,maxSD,0.1)){
    
    for (d in seq(0,(numDS-1))){
      
      propC=95+d*(5/numDS)
      propc=5-d*(5/numDS)
      df<-normals_imb(numEx,sd = s, prob =c(propC,propc), distance = dist)
      
      openPDFEPS(paste("ImbDS2norm_",numEx,"_",s,"_(",propC,",",propc,")",sep=""))
      p <- ggplot(df,aes(x,y,colour=Class)) + geom_point()
      print(p)
      dev.off()
      write.csv(df,file=paste("ImbDS2norm_",numEx,"_",s,"_(",propC,",",propc,").csv",sep=""), row.names = FALSE)
      
      print(paste("ImbDS2norm_",numEx,"_",s,"_(",propC,",",propc,")",sep=""))
      
    }
  }
} 

#same dataset: change class

multi_normals_imb_2 <- function(numDS, numEx)
  {
    set.seed(288)
    propC=95
    propc=5
    numC=propC/100 * numEx
    numc= propc/100 * numEx
    df<-normals_imb(numEx,sd=1, prob =c(propC,propc))
    ggplot(df,aes(x,y,colour=Class)) + geom_point()
    numcReal <- (sum(df$Class == "2"))
    diff <- numc - numcReal
    subsetc <- which(df$Class == "2")
    subsetC <- which(df$Class == "1")
    if (diff < 0){ #sobran de la clase minoritaria
      df[sample(subsetc,abs(diff)),"Class"] <- "1"
    }else{ # sobran de la clase mayoritaria
      df[sample(subsetC,abs(diff)),"Class"] <- "2"
    }
   
    
        
    for (d in seq(1,(numDS))){
      
     
      propC=95+d*(5/numDS)
      propc=5-d*(5/numDS)
      numC=propC/100 * numEx
      numc= propc/100 * numEx
      toChangeClass = sum(df$Class == "2") - numc
      subsetc <- which(df$Class == "2")
      df[sample(subsetc, toChangeClass),"Class"] <- "1"
      
      openPDFEPS(paste("ImbDS2norm_",numEx,"_(",propC,",",propc,")",sep=""))
      p <- ggplot(df,aes(x,y,colour=Class)) + geom_point()
      print(p)
      dev.off()
      write.csv(df,file=paste("ImbDS2norm_",numEx,"_(",propC,",",propc,").csv",sep=""), row.names = FALSE)
      print(paste("ImbDS2norm_",numEx,"_(",propC,",",propc,")",sep=""))
      
      
    }
  
}

#same dataset: delete observations
multi_normals_imb_3 <- function(numDS=5, numEx=200, seed = 288, rep=1, propC = c(50,75,90,95,98), s = 1, dist = 1)
{
  set.seed(seed)
  propMaj = propC[1]
  propMin = 100 - propC[1]
  numC=propC[1]/100 * numEx
  numc= (100-propC[1])/100 * numEx
  df<-normals_imb(numEx,sd=1, prob =c(propMaj,propMin))
  numcReal <- (sum(df$Class == "1"))
  diff <- numc - numcReal
  
  while (diff != 0){
    df<-normals_imb(numEx,sd=s, prob =c(propMaj,propMin), distance = dist)
    numcReal <- (sum(df$Class == "1"))
    diff <- numc - numcReal
  }
  
  
  
  
  openPDFEPS(paste(rep,"_ImbDS2norm_e",numEx,"_sd",s,"_x",dist,"_(",propMaj,",",propMin,")",sep=""))
  p <- ggplot(df,aes(x,y,colour=Class)) + geom_point()
  print(p)
  dev.off()
  write.csv(df,file=paste(rep,"_ImbDS2norm_e",numEx,"_sd",s,"_x",dist,"_(",propMaj,",",propMin,").csv",sep=""), row.names = FALSE)
  print(paste(rep,"_ImbDS2norm_",numEx,"_(",propMaj,",",propMin,")",sep=""))
  
  orig.MinClassEx <- sum(df$Class == "1")
  
  for (d in seq(2,(numDS))){
    
    
    propMaj=propC[d]
    propMin.prev = 100 - propC[d-1]
    propMin=100-propC[d]
    dif = propMin.prev - propMin
    
    numC=propMaj/100 * nrow(df)
    numc= propMin/100 * nrow(df)
    toDeleteNumber= dif/100 * nrow(df)
    
    indices.Min <- which(df$Class == "1")
    indices.Maj <- which(df$Class == "0")
    
    dfTemp1 <- df[-c(sample(indices.Min, toDeleteNumber)),]
    dfTemp2 <- df[c(sample(indices.Maj, toDeleteNumber)),]
    df <- rbind(dfTemp1, dfTemp2)

    numEx = nrow(df)
    propMaj = length(which(df$Class == "0"))/nrow(df)*100
    propMin = length(which(df$Class == "1"))/nrow(df)*100
    
    
    
    openPDFEPS(paste(rep,"_ImbDS2norm_e",numEx,"_sd",s,"_x",dist,"_(",propMaj,",",propMin,")",sep=""))
    p <- ggplot(df,aes(x,y,colour=Class)) + geom_point()
    print(p)
    dev.off()
    write.csv(df,file=paste(rep,"_ImbDS2norm_e",numEx,"_sd",s,"_x",dist,"_(",propMaj,",",propMin,").csv",sep=""), row.names = FALSE)
    print(paste(rep,"_ImbDS2norm_",numEx,"_(",propMaj,",",propMin,")",sep=""))
    
    
  }
  
}

run <- function(){
  multi_normals_imb_3(s=1, dist=0.5, rep = 1)
  multi_normals_imb_3(s=0.5, dist=0.5, rep = 2)
  multi_normals_imb_3(s=0.25, dist=0.5, rep = 3)
  multi_normals_imb_3(s=1, dist=0.25, rep = 4)
  multi_normals_imb_3(s=0.5, dist=0.25, rep = 5)
  multi_normals_imb_3(s=0.25, dist=0.25, rep = 6)
  
  }
#same dataset: delete observations
#multi_normals_imb_3 <- function(numDS, numEx, seed)
 
repetitions_mni3 <- function(rep = 5){
  for (i in seq(1,rep,1)){
      multi_normals_imb_3(5, 200, i, i)
    }
  } 



reduceMinC_UCI <- function(ds="08_parkinsons.csv", numDS=5,seed=288){
  set.seed(seed)
  datos <- read.csv(ds)
  minC <- which.min(as.vector(table(datos$Class)))
  minC_class <- dimnames(table(datos$Class))[[1]][minC] 
  numC = nrow(datos[datos$Class==minC_class,])
  print(numC)
  numEx = nrow(datos)
  toDel = sample(row.names(datos[datos$Class==minC_class,]))
  numEx2Del = trunc(numC/numDS,0)
  
  for (i in 1:numDS){
   
    
    #numMin = nrow(datos[datos$Class==minC_class,])
    #print(numMin)
    toDelSegment = toDel[(((i-1)*numEx2Del)+1):(i*numEx2Del)]
    datos <- datos[!row.names(datos)%in%toDelSegment,]        
    
    #Data to print
    numExNew = nrow(datos)
    numClasses = length(unique(datos$Class))
    propClasses = paste(round(table(datos$Class)/numExNew,2))
    printPropClasses= paste(propClasses, collapse = " ")
    nameDS <- strsplit(ds,"[.]")[[1]][1] #keep just the name
    
    #PCA for plotting
    ir.pca <- prcomp(datos[,1:ncol(datos)-1],center = TRUE,scale. = TRUE) 
    do <- datos
    do$pc1 <- ir.pca$x[,'PC1']
    do$pc2 <- ir.pca$x[,'PC2']
  
    openPDFEPS(paste(nameDS,"_",numExNew,"_(",printPropClasses,")",sep=""))
    dis <- ggplot(do, aes(pc1,pc2, colour = factor(Class))) + geom_point() + theme_bw()       
    print(dis)
    dev.off()
    
    write.csv(datos,file=paste(nameDS,"_",numExNew,"_(",printPropClasses,").csv",sep=""), row.names = FALSE)
    print(paste(nameDS,"_",numExNew,"_(",printPropClasses,")",sep=""))
    
  }
}


  
balanceDS <- function(ds="1_hepatitis.csv", seed=288, number = 125){
  
  set.seed(seed)
  w <- read.csv(ds)
  minC <- which.min(as.vector(table(w$Class)))
  minC_class <- dimnames(table(w$Class))[[1]][minC] 
  numc = nrow(w[w$Class==minC_class,])
  
  majC <- which.max(as.vector(table(w$Class)))
  majC_class <- dimnames(table(w$Class))[[1]][majC] 
  numC = nrow(w[w$Class==majC_class,])
  
  if(numC == numc){
    majC_class = majC_class <- dimnames(table(w$Class))[[1]][majC+1] 
  }
  
  indices.Min <- which(w$Class == minC_class)
  indices.Maj <- which(w$Class == majC_class)
  
  if (numc<number){
    extraMin <- w[c(sample(indices.Min, number-numc)),]
    dfTemp1 <- rbind(w,extraMin) #oversampling 
  }else{
    dfTemp1 <- w[c(sample(indices.Min, number)),] #undersampling 
  }
  
  if (numC<number){
    extraMaj <- w[c(sample(indices.Maj, number-numC)),]
    #dfTemp2 <- rbind(w,extraMaj)
    dfTemp2 <- extraMaj
  }else{
    dfTemp2 <- w[c(sample(indices.Maj, number)),]
  }
  
  dfTemp3 <- rbind(dfTemp1, dfTemp2)
  dfTemp3 <- dfTemp3[sample(row.names(dfTemp3)),]
  row.names(dfTemp3) <- seq(1:nrow(dfTemp3))
  
  nameDS <- strsplit(ds,"[.]")[[1]][1] #keep just the name
  write.csv(dfTemp3,file=paste(nameDS,"_bal.csv",sep=""), row.names = FALSE)
  
}
  










