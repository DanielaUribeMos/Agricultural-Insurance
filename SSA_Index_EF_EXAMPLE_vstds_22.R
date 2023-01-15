###SSA
library(dplyr)
library(dbplyr)
library(tidyverse)
library(reshape)
#install.packages("Rssa")
library(Rssa)
setwd("/Users/DaniU/Documents/Semestre7/Econ/Ejercicios/Semana11")
set.seed(123)
#######################################################
##############FUNCIONES NO parte del paquete RSSA###
##Hankelizacisn:
UniHankel=function(Y,L){
  k<-length(Y)-L+1
  outer((1:L), (1:k), function(x,y) Y[(x+y-1)])
}

##Funcisn SVD
SVD<- function(Y,L){
  X<-UniHankel(Y,L)
  svd(X)
}

##Proceso de agrupamiento
Group<-function(Y,L,groups){
  I<-groups; p<-length(I)
  SVD<-SVD(Y,L)
  LambdaI<-matrix(diag(SVD$d) [I,I],p,p)
  SVD$u[,I]%*%LambdaI%*%t(SVD$v[,I])
}

###DiagonalAveraging
DiagAver<-function(X){
  L<-nrow(X); k<-ncol(X); N<-k+L-1
  D<-NULL
  for(j in 1:N){
    s1<-max(1, (j-N+L))
    s2<-min(L,j)
    place<-(s1:s2)+L*(((j+1-s1):(j+1-s2))-1)
    D[j]<-mean(X[place])
  }
  D
}

##Reconstruccisn
SSA.Rec<- function(Y,L, groups){
  N<-length(Y)
  I<-groups;p<-length(I)
  XI<-Group(Y,L, groups)
  Approx<-DiagAver(XI)
  Resid<-Y-Approx
  list(Approximation=Approx, Residual=Resid)
}

###W-correlation
W.corr<-function(Yt,L, groups){
  m<-length(groups); w.corr<-diag(m)
  N<-length(Yt)
  w<-((N+1)-abs((N+1)/2-L)-abs((N+1)/2-1:N)-
    abs(abs((N+1)/2-L)-abs((N+1)/2-1:N)))/2
wcorr<-function(i,j){
  Y1<-SSA.Rec(Yt,L,groups[[i]])$Approximation
  Y2<-SSA.Rec(Yt,L,groups[[j]])$Approximation
  sum(w*Y1*Y2)/sqrt(sum(w*Y1^2)*sum(w*Y2^2))}
for (i in 1:(m-1)){
  for (j in (i+1):m){
    w.corr[i,j]=w.corr[j,i]=wcorr(i,j)}}
rownames(w.corr)<-colnames(w.corr)<-groups
w.corr
}

###Plotting Images w.corr
Plt.Img<-function(x){
  min<-min(x)
  max<-max(x)
  yLabels<-rownames(x)
  xLabels<-colnames(x)
  if( is.null (xLabels)){
    xLabels <- c(1:ncol(x))
  }
  if (is.null(yLabels)){
    yLabels <- c(1:nrow(x))
  }
  layout(matrix(data=c(1,2), nrow=1, ncol=2),
         widths=c(4,1), heights=c(1,1))
  ColorRamp<-gray( seq(1,0, length=20))
  ColorLevels<-seq(min, max, length=length(ColorRamp))
                  par(mar=c(3,5,2.5,2))
                  image(1:length(xLabels), 1:length(yLabels),
                        t(x), col=colorRamp, xlab="",
                        ylab = "", axes=FALSE, zlim=c(min, max))
                  title(main=c("Image Plot"))
                  axis(BELOW<-1, at=1:length(xLabels),
                       labels=xLabels, cex.axis=0.7)
                  axis(LEFT<-2, at=1:length(yLabels),
                       labels=yLabels, las=HORIZONTAL<-1,
                       cex.axis=0.7)
                  box()
                  par(mar=c(3,2.5, 2.5,2))
                  image(1, ColorLevels,
                        matrix(data=ColorLevels, ncol=length(ColorLevels), nrow=1),
                        col=colorRamp,
                        xlab="", ylab="",
                        xaxt="n")
                  layout(1)
}

##Singular Values graph (7)
Sing.plt<-function(Y,L){
  lambda<-log(SVD(Y,L)$d)
  d<-length(lambda)
  win<-1:d
  plot.new()
  plot.window(xlim=range(win), ylim=range(lambda))
  usr=par("usr")
  rect(usr[1], usr[3], usr[2], usr[4])
  lines(win, lambda, lwd=2)
  points(win, lambda, pch=21, bg="gray")
  axis(1)
  axis(2)
  box()
  title(xlab="Number")
  title(ylab="Log.Singular Values")}

### Recurent forecasting SSA (9)
RSSA.Forecasting <- function(L,groups,h,Y){
  N <-length(Y)
  L <-min(L,(N-L+1))
  X <-UniHankel(Y,L)
  U <-matrix(svd(X)$u,L,L)
  pi <-array(U[L,groups], dim=length(groups))
  V2 <-sum(pi^2)
  m <-length(groups)
  Udelta=array(U[1:(L-1),groups], dim=c((L-1),m))
  A <-pi%*%t(Udelta)/(1-V2)
  yhat=array(0, dim=(N+h))
  yhat[1:N]<-SSA.Rec(Y,L,groups)$Approximation
  for(i in (N+1):(N+h))
  yhat[1]<-A%*%yhat[(1-L+1):(1-1)]
  yhat[(N+1):(N+h)]
}

###################################################
#################################################### SSA
## SSA 
cs2=read.csv("ptpm_col_example3.csv", header = TRUE, sep=";")
head(cs2)
library(psych)
describe(cs2$Station)
describe(cs2$Alt)
sapply(cs2, function(CL) length(unique(CL)))
#sapply(cs2, function(Stations) length(unique(Stations)))

head(cs2)
str(cs2)
cs2$Date<- as.Date(cs2$date, format='%d/%m/%Y')
head(cs2)
tail(cs2)
cs2_ssa=cs2[,c(-13)]
head(cs2_ssa)
str(cs2_ssa)
library(dplyr)
cs2_ssa$CL =as.character(cs2_ssa$CL)

#######################
indice <- function(x, lower_limit, lower_threshold, upper_limit, upper_threshold, max_indemnity){
  cond1 <- x <= lower_limit | x >= upper_limit
  cond2 <- x >= lower_threshold & x <= upper_threshold
  cond3 <- x >= lower_limit & x<= lower_threshold
  payment1 <- (lower_threshold-x)/(lower_threshold-lower_limit)*max_indemnity
  payment2 <- (x-upper_threshold)/(upper_limit-upper_threshold)*max_indemnity
  output <- ifelse(cond1, max_indemnity, 
                   ifelse(cond2, 0, 
                          ifelse(cond3, payment1, payment2)
                   )
  )
  return(output)
}



#########indice
###variables
max_indemnity=3200 ##cambiar
lower_limit = 35
lower_threshold = 60
upper_threshold = 140
upper_limit = 250
####################
##Promedio de lluvia por cluster
library(dbplyr)
library(plyr)
library(tidyverse)
library(magrittr)
library(psych)

cl_ptpm<-ddply(cs2_ssa, c("CL", "Date", "Crop_Season"), summarise,
               mean_pt = mean(PTPM)) %>% na.omit()
head(cl_ptpm)
tail(cl_ptpm)
str(cl_ptpm)

###
#CL1:
scl1_Wet <- dplyr::filter(cl_ptpm, CL == "1" & Crop_Season == 'Wet')
describe(scl1_Wet$mean_pt)
head(scl1_Wet)
scl1_Dry <- dplyr::filter(cl_ptpm, CL == "1" & Crop_Season == 'Dry')
describe(scl1_Dry$mean_pt)
head(scl1_Dry)
#CL2:
scl2_Wet <- dplyr::filter(cl_ptpm, CL == "2" & Crop_Season == 'Wet')
describe(scl2_Wet$mean_pt)
head(scl2_Wet)
scl2_Dry <- dplyr::filter(cl_ptpm, CL == "2" & Crop_Season == 'Dry')
describe(scl2_Dry$mean_pt)
head(scl2_Dry)
#CL3:
scl3_Wet <- dplyr::filter(cl_ptpm, CL == "3" & Crop_Season == 'Wet')
describe(scl3_Wet$mean_pt)
head(scl3_Wet)
scl3_Dry <- dplyr::filter(cl_ptpm, CL == "3" & Crop_Season == 'Dry')
describe(scl3_Dry$mean_pt)
head(scl3_Dry)

################################################
################## SSA


########################################## CLUSTER 1   WET ############################################
ssacl1_Wet<-ts(scl1_Wet$mean_pt, start=c(2010,04), frequency = 12)
ssacl1_Wet
plot(ssacl1_Wet)
library(Rssa)
###FASE 1. Descomposicisn
scl1_ssa_Wet=ssa(ssacl1_Wet)
summary(scl1_ssa_Wet)
##W-corrr matrix
wfort<-wcor(scl1_ssa_Wet, groups=1:30)
plot(wfort, grid=c(2,4,5,7))

##Grafica de valores singulares
Sing.plt(ssacl1_Wet, 24)

plot(scl1_ssa_Wet, type="vectors", idx=1:14)
plot(scl1_ssa_Wet, type= "paired", idx=2:11, plot.contrib = FALSE)  #los virtices son la frecuencia, determinan el periodo de las ondas seno.
#print(parestimate(s.fort, groups = list(2:3, 4:5), method = "pairs")) ##Muestra la frecuencia de los pares 2-3 y 4-5.

plot(wcor(scl1_ssa_Wet, groups = 1:30), scales=list(at=c(10, 20, 30)))

### Fase de Reconstruccisn:
r.cl1_Wet= reconstruct(scl1_ssa_Wet, groups= list(Trend=c(1,4), Seasonality=c(2:3, 5:6, 7:8 )))
#plot(r.cl1_Wet, add.residuals=TRUE, add.original=TRUE, plot.method="xyplot", main="Reconstruction Phase Cluster 1",
 #    superpose=TRUE, auto.key=list(columns=4))

####              Analisis Payoff por componentes Tendencia y Estacionalidad
r.cl1_S_Wet <- r.cl1_Wet$Seasonality+r.cl1_Wet$Trend

###Payouts:
#CL1_Wet
###kernel density
library(logKDE)
fitcl1_wet=logdensity(r.cl1_S_Wet, kernel = "logistic", n=1000)
inplot( fitcl1_wet, main = "Cluster 1 Wet Precipitation Density ",
      xlab="Precipitation (PTPM)", ylab="Density", 
      panel.first=grid(nx = NULL, ny = NULL, lwd=2, col = "lightgray", lty = "dotted"))

fitcl1_wet$x <- ifelse(fitcl1_wet$x < 0, .Machine$double.eps, fitcl1_wet$x) 
describe(fitcl1_wet$x)

datacl1_wet <- as.data.frame(x = fitcl1_wet$x)
##Payoff
payoffcl1_wet <- indice(datacl1_wet$`fitcl1_wet$x`, lower_limit, lower_threshold, upper_limit, upper_threshold, max_indemnity)
#comp3_wet <- as.data.frame(cbind(data3_wet, payoff3_wet))
#head(comp3_wet)
###Pago Indemnity payment due:
mean(payoffcl1_wet)
sd(payoffcl1_wet)

Percent_indemnizationcl1_wet=(mean(payoffcl1_wet)/max_indemnity)*100
Percent_indemnizationcl1_wet
#####################


#CL1 Dry SSA:
########################################## CLUSTER 1   Dry ############################################
ssacl1_Dry<-ts(scl1_Dry$mean_pt, start=c(2010,01), frequency = 12)
ssacl1_Dry
plot(ssacl1_Dry)

###FASE 1. Descomposicisn
scl1_ssa_Dry=ssa(ssacl1_Dry)
summary(scl1_ssa_Dry)
##W-corrr matrix
wfort<-wcor(scl1_ssa_Dry, groups=1:20)
plot(wfort, grid=c(2,4,5,7))

##Grafica de valores singulares
Sing.plt(ssacl1_Dry, 24)

plot(scl1_ssa_Dry, type="vectors", idx=1:14)
plot(scl1_ssa_Dry, type= "paired", idx=2:11, plot.contrib = FALSE)  #los virtices son la frecuencia, determinan el periodo de las ondas seno.
#print(parestimate(s.fort, groups = list(2:3, 4:5), method = "pairs")) ##Muestra la frecuencia de los pares 2-3 y 4-5.

plot(wcor(scl1_ssa_Dry, groups = 1:20), scales=list(at=c(10, 20, 30)))

### Fase de Reconstruccisn:
r.cl1_Dry= reconstruct(scl1_ssa_Dry, groups= list(Trend=c(1,4), Seasonality=c(2:3)))
#plot(r.cl1_Dry, add.residuals=TRUE, add.original=TRUE, plot.method="xyplot", main="Reconstruction Phase Cluster 1",
 #   superpose=TRUE, auto.key=list(columns=4))
######################################################################
####              Analisis Payoff por componentes Tendencia y Estacionalidad
r.cl1_S_Dry <- r.cl1_Dry$Seasonality+r.cl1_Dry$Trend

###Payouts:
#CL1_Dry
###kernel density
fitcl1_dry=logdensity(r.cl1_S_Dry, kernel = "logistic", n=1000)
plot( fitcl1_dry, main = "Cluster 1 Dry Precipitation Density ",
      xlab="Precipitation (PTPM)", ylab="Density", 
      panel.first=grid(nx = NULL, ny = NULL, lwd=2, col = "lightgray", lty = "dotted"))

fitcl1_dry$x <- ifelse(fitcl1_dry$x < 0, .Machine$double.eps, fitcl1_dry$x) 
describe(fitcl1_dry$x)

datacl1_dry <- as.data.frame(x = fitcl1_dry$x)
##Payoff
payoffcl1_dry <- indice(datacl1_dry$`fitcl1_dry$x`, lower_limit, lower_threshold, upper_limit, upper_threshold, max_indemnity)
#comp3_dry <- as.data.frame(cbind(data3_dry, payoff3_dry))
#head(comp3_dry)
###Pago Indemnity payment due:
mean(payoffcl1_dry)
sd(payoffcl1_dry)

Percent_indemnizationcl1_dry=(mean(payoffcl1_dry)/max_indemnity)*100
Percent_indemnizationcl1_dry
#######################################CONTINUE:
#--------------------------------------------------------------------

########################################## CLUSTER 2 WET ############################################
ssacl2_Wet<-ts(scl2_Wet$mean_pt, start=c(2010,04), frequency = 12)
ssacl2_Wet
plot(ssacl2_Wet)
library(Rssa)
###FASE 1. Descomposicisn
scl2_ssa_Wet=ssa(ssacl2_Wet)
summary(scl2_ssa_Wet)
##W-corrr matrix
wfort<-wcor(scl2_ssa_Wet, groups=1:30)
plot(wfort, grid=c(2,4,5,7))

##Grafica de valores singulares
Sing.plt(ssacl2_Wet, 24)

plot(scl2_ssa_Wet, type="vectors", idx=1:14)
plot(scl2_ssa_Wet, type= "paired", idx=2:11, plot.contrib = FALSE)  #los virtices son la frecuencia, determinan el periodo de las ondas seno.

#print(parestimate(s.fort, groups = list(2:3, 4:5), method = "pairs")) ##Muestra la frecuencia de los pares 2-3 y 4-5.

plot(wcor(scl2_ssa_Wet, groups = 1:30), scales=list(at=c(10, 20, 30)))

### Fase de Reconstruccisn:
r.cl2_Wet= reconstruct(scl2_ssa_Wet, groups= list(Trend=c(1,4,5,6,7,8), Seasonality=c(2:3)))
#plot(r.cl1_Wet, add.residuals=TRUE, add.original=TRUE, plot.method="xyplot", main="Reconstruction Phase Cluster 1",
#    superpose=TRUE, auto.key=list(columns=4))

######################################################################
####              Analisis Payoff por componentes Tendencia y Estacionalidad
r.cl2_S_Wet <- r.cl2_Wet$Seasonality+r.cl2_Wet$Trend

###Payouts:
#CL1_Wet
###kernel density
library(logKDE)
fitcl2_wet=logdensity(r.cl2_S_Wet, kernel = "logistic", n=1000)
inplot( fitcl2_wet, main = "Cluster 1 Wet Precipitation Density ",
        xlab="Precipitation (PTPM)", ylab="Density", 
        panel.first=grid(nx = NULL, ny = NULL, lwd=2, col = "lightgray", lty = "dotted"))

fitcl2_wet$x <- ifelse(fitcl2_wet$x < 0, .Machine$double.eps, fitcl2_wet$x) 
describe(fitcl2_wet$x)

datacl2_wet <- as.data.frame(x = fitcl2_wet$x)
##Payoff
payoffcl2_wet <- indice(datacl2_wet$`fitcl2_wet$x`, lower_limit, lower_threshold, upper_limit, upper_threshold, max_indemnity)
#comp3_wet <- as.data.frame(cbind(data3_wet, payoff3_wet))
#head(comp3_wet)
###Pago Indemnity payment due:
mean(payoffcl2_wet)
sd(payoffcl2_wet)

Percent_indemnizationcl2_wet=(mean(payoffcl2_wet)/max_indemnity)*100
Percent_indemnizationcl2_wet

########################################## CLUSTER 2 DRY ############################################
ssacl2_Dry<-ts(scl2_Dry$mean_pt, start=c(2010,04), frequency = 12)
ssacl2_Dry
plot(ssacl2_Dry)
library(Rssa)
###FASE 1. Descomposicisn
scl2_ssa_Dry=ssa(ssacl2_Dry)
summary(scl2_ssa_Dry)
##W-corrr matrix
wfort<-wcor(scl2_ssa_Dry, groups=1:30)
plot(wfort, grid=c(2,4,5,7))

##Grafica de valores singulares
Sing.plt(ssacl2_Dry, 24)

plot(scl2_ssa_Dry, type="vectors", idx=1:14)
plot(scl2_ssa_Dry, type= "paired", idx=2:11, plot.contrib = FALSE)  #los virtices son la frecuencia, determinan el periodo de las ondas seno.

#print(parestimate(s.fort, groups = list(2:3, 4:5), method = "pairs")) ##Muestra la frecuencia de los pares 2-3 y 4-5.

plot(wcor(scl2_ssa_Dry, groups = 1:30), scales=list(at=c(10, 20, 30)))

### Fase de Reconstruccisn:
r.cl2_Dry= reconstruct(scl2_ssa_Dry, groups= list(Trend=c(1,4,5,6,7,8), Seasonality=c(2:3)))
#plot(r.cl1_Wet, add.residuals=TRUE, add.original=TRUE, plot.method="xyplot", main="Reconstruction Phase Cluster 1",
#    superpose=TRUE, auto.key=list(columns=4))

######################################################################
####              Analisis Payoff por componentes Tendencia y Estacionalidad
r.cl2_S_Dry <- r.cl2_Dry$Seasonality+r.cl2_Dry$Trend

###Payouts:
#CL1_Wet
###kernel density
library(logKDE)
fitcl2_Dry=logdensity(r.cl2_S_Dry, kernel = "logistic", n=1000)
inplot( fitcl2_Dry, main = "Cluster 1 Dry Precipitation Density ",
        xlab="Precipitation (PTPM)", ylab="Density", 
        panel.first=grid(nx = NULL, ny = NULL, lwd=2, col = "lightgray", lty = "dotted"))

fitcl2_Dry$x <- ifelse(fitcl2_Dry$x < 0, .Machine$double.eps, fitcl2_Dry$x) 
describe(fitcl2_Dry$x)

datacl2_Dry <- as.data.frame(x = fitcl2_Dry$x)
##Payoff
payoffcl2_Dry <- indice(datacl2_Dry$`fitcl2_Dry$x`, lower_limit, lower_threshold, upper_limit, upper_threshold, max_indemnity)
#comp3_wet <- as.data.frame(cbind(data3_wet, payoff3_wet))
#head(comp3_wet)
###Pago Indemnity payment due:
mean(payoffcl2_Dry)
sd(payoffcl2_Dry)

Percent_indemnizationcl2_Dry=(mean(payoffcl2_Dry)/max_indemnity)*100
Percent_indemnizationcl2_Dry

########################################## CLUSTER 3 WET ############################################

ssacl3_Wet<-ts(scl3_Wet$mean_pt, start=c(2010,04), frequency = 12)
ssacl3_Wet
plot(ssacl3_Wet)
library(Rssa)
###FASE 1. Descomposicisn
scl3_ssa_Wet=ssa(ssacl3_Wet)
summary(scl3_ssa_Wet)
##W-corrr matrix
wfort<-wcor(scl3_ssa_Wet, groups=1:30)
plot(wfort, grid=c(2,4,5,7))

##Grafica de valores singulares
Sing.plt(ssacl3_Wet, 24)

plot(scl3_ssa_Wet, type="vectors", idx=1:14)
plot(scl3_ssa_Wet, type= "paired", idx=2:11, plot.contrib = FALSE)  #los virtices son la frecuencia, determinan el periodo de las ondas seno.

#print(parestimate(s.fort, groups = list(2:3, 4:5), method = "pairs")) ##Muestra la frecuencia de los pares 2-3 y 4-5.

plot(wcor(scl3_ssa_Wet, groups = 1:30), scales=list(at=c(10, 20, 30)))

### Fase de Reconstruccisn:
r.cl3_Wet= reconstruct(scl3_ssa_Wet, groups= list(Trend=c(1,4,5,6,7,8), Seasonality=c(2:3, 9:10)))
#plot(r.cl1_Wet, add.residuals=TRUE, add.original=TRUE, plot.method="xyplot", main="Reconstruction Phase Cluster 1",
#    superpose=TRUE, auto.key=list(columns=4))

r.cl3_S_Wet <- r.cl3_Wet$Seasonality+r.cl3_Wet$Trend

###Payouts:
#CL1_Wet
###kernel density
library(logKDE)
fitcl3_wet=logdensity(r.cl3_S_Wet, kernel = "logistic", n=1000)
inplot( fitcl3_wet, main = "Cluster 1 Wet Precipitation Density ",
        xlab="Precipitation (PTPM)", ylab="Density", 
        panel.first=grid(nx = NULL, ny = NULL, lwd=2, col = "lightgray", lty = "dotted"))

fitcl3_wet$x <- ifelse(fitcl3_wet$x < 0, .Machine$double.eps, fitcl3_wet$x) 
describe(fitcl3_wet$x)

datacl3_wet <- as.data.frame(x = fitcl3_wet$x)
##Payoff
payoffcl3_wet <- indice(datacl3_wet$`fitcl3_wet$x`, lower_limit, lower_threshold, upper_limit, upper_threshold, max_indemnity)
#comp3_wet <- as.data.frame(cbind(data3_wet, payoff3_wet))
#head(comp3_wet)
###Pago Indemnity payment due:
mean(payoffcl3_wet)
sd(payoffcl3_wet)

Percent_indemnizationcl3_wet=(mean(payoffcl3_wet)/max_indemnity)*100
Percent_indemnizationcl3_wet

########################################## CLUSTER 3 DRY ############################################

ssacl3_Dry<-ts(scl3_Dry$mean_pt, start=c(2010,04), frequency = 12)
ssacl3_Dry
plot(ssacl3_Dry)
library(Rssa)
###FASE 1. Descomposicisn
scl3_ssa_Dry=ssa(ssacl3_Dry)
summary(scl3_ssa_Dry)
##W-corrr matrix
wfort<-wcor(scl3_ssa_Dry, groups=1:30)
plot(wfort, grid=c(2,4,5,7))

##Grafica de valores singulares
Sing.plt(ssacl3_Dry, 24)

plot(scl3_ssa_Dry, type="vectors", idx=1:14)
plot(scl3_ssa_Dry, type= "paired", idx=2:11, plot.contrib = FALSE)  #los virtices son la frecuencia, determinan el periodo de las ondas seno.

#print(parestimate(s.fort, groups = list(2:3, 4:5), method = "pairs")) ##Muestra la frecuencia de los pares 2-3 y 4-5.

plot(wcor(scl3_ssa_Dry, groups = 1:30), scales=list(at=c(10, 20, 30)))

### Fase de Reconstruccisn:
r.cl3_Dry= reconstruct(scl3_ssa_Dry, groups= list(Trend=c(1,4,11, 12,15), Seasonality=c(2:3, 5:6, 7:10, 13:14)))
#plot(r.cl1_Wet, add.residuals=TRUE, add.original=TRUE, plot.method="xyplot", main="Reconstruction Phase Cluster 1",
#    superpose=TRUE, auto.key=list(columns=4))

r.cl3_S_Dry <- r.cl3_Dry$Seasonality+r.cl3_Dry$Trend

###Payouts:
#CL1_Wet
###kernel density
library(logKDE)
fitcl3_dry=logdensity(r.cl3_S_Dry, kernel = "logistic", n=1000)
inplot( fitcl3_dry, main = "Cluster 1 Wet Precipitation Density ",
        xlab="Precipitation (PTPM)", ylab="Density", 
        panel.first=grid(nx = NULL, ny = NULL, lwd=2, col = "lightgray", lty = "dotted"))

fitcl3_dry$x <- ifelse(fitcl3_dry$x < 0, .Machine$double.eps, fitcl3_dry$x) 
describe(fitcl3_wet$x)

datacl3_dry <- as.data.frame(x = fitcl3_dry$x)
##Payoff
payoffcl3_dry <- indice(datacl3_dry$`fitcl3_dry$x`, lower_limit, lower_threshold, upper_limit, upper_threshold, max_indemnity)
#comp3_wet <- as.data.frame(cbind(data3_wet, payoff3_wet))
#head(comp3_wet)
###Pago Indemnity payment due:
mean(payoffcl3_dry)
sd(payoffcl3_dry)

Percent_indemnizationcl3_dry=(mean(payoffcl3_dry)/max_indemnity)*100
Percent_indemnizationcl3_dry


###
#Resultados
# Clúster 1
# Wet: media= 1 192 295    indemnización= 39.7431
# Dry: media= 550 779.1    indemnización= 18.3593
# Clúster 2
# Wet: media= 1 874 234    indemnización= 62.4744 
# Dry: media= 1 523 873    indemnización= 50.7957
# Clúster 3
# Wet: media= 1 523 540    indemnización= 50.7846
# Dry: media= 872 507.8    indemnización= 29.0835
