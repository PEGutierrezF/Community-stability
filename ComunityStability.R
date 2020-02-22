
library(nlme)
laselva<-read.csv("D:/Curriculum/02_ Articulos/00 In progress/220 Community stability/Community-stability/laselva.csv",h=T) # read the accompanying csv file
head(laselva)
nrow(laselva)


## GRADUAL LINEAR DYNAMICS  #####

# where M is the diffusion constant 
# t denotes the time since beginning of the record, 
# and C is an error term representing the intercept. 

res_GLD<-list(data=NA)
aic_GLD<-list(data=NA)
GLD.Mod <- NA
print(GLD.Mod)

GLD.Mod <- try(gnls(dist ~ C + M*(samp_event), data = laselva, 
                        correlation=corAR1(), 
                        start=list(C=0,M=0)))

print(GLD.Mod)
print(class(GLD.Mod))
summary(GLD.Mod)

fitted(GLD.Mod) # ajuste de la curva
AIC(GLD.Mod) # AIC del modelo

options(digits=3)  # digit number
laselva$GLD <- unlist(fitted(GLD.Mod))
head(laselva)

plot(dist ~ samp_event, data=laselva)
curve(predict(GLD.Mod, newdata = data.frame(samp_event=x)), col = "black", add = TRUE)


### REVERSIBLE DYNAMICS  #####
## Defined by a double sigmoid mathematical fuction ##

## dist = distancia the bray curtis
## d_o = asymptotic height for onset
## d_r= asymptotic height for return
## teta_o = time of the onset, in which migration reaches one-half of its asymptotic height
## teta_r = time of the return, in which migration reaches one-half of its asymptotic height
## w_o = time elapsed between reaching one-half and three-quarters of the migration distance for onset
## w_r = time elapsed between reaching one-half and three-quarters of the migration distance for return
## t = time interval since beginning of the record.


res_Rev<-list(data=NA)
aic_Rev<-list(data=NA)

Rev.Mod <- try(gnls(dist ~ ((d_o)/(1+exp((teta_o-samp_event)/(w_o)))) - ((d_r)/(1+exp((teta_r-samp_event)/(w_r)))),
                    data = laselva, 
                    correlation=corAR1(),
                    start = list(d_o = 0.6 , d_r = 0.6, teta_o = 10,teta_r = 100, w_o = 1, w_r = 1),     
                    control = gnlsControl(nlsTol = 500)))

print(Rev.Mod )
print(class(Rev.Mod ))
summary(Rev.Mod )

fitted(Rev.Mod) # ajuste de la curva
AIC(Rev.Mod) # AIC del modelo

options(digits=3)  # digit number
laselva$Rev <- unlist(fitted(Rev.Mod))
head(laselva)
coef(Rev.Mod)

###### Reversible dynamics graph - Double sigmoidal #####

plot(dist ~ samp_event, data=laselva)
curve(predict(Rev.Mod , newdata = data.frame(samp_event=x)), col = "black", add = TRUE)

## STABLE BEHAVIOUR ####

res_stab<-list(data=NA)
aic_stab<-list(data=NA)
null.mod<-NA  
asym.HRmod<-NA
stab.Mod <- NA


### Step 1 Null Model
null.mod <- try(gnls(dist ~ A, data = laselva, 
                     correlation = corAR1(),
                     start = list(A = mean(laselva[,'dist'])), 
                     control=nlmeControl(maxIter=50, pnlsMaxIter=7, niterEM=25, returnObject=TRUE),verbose=TRUE))

## Step 2 Asymp model

asym.HRmod <- try(gnls(dist ~ Asym*(1-exp(lrc*year)), data = laselva, 
                       correlation = corAR1(),
                       start = list(Asym =summary(null.mod)$tTable[1],lrc=-0.059), 
                       control = gnlsControl(nlsTol = 50)))

## Step 3 Full stability model

stab.Mod <- try(gnls(dist ~ (Asym)*(1 - exp(lrc*(year))),data = laselva,
                     correlation = corAR1(),
                     start = c(Asym = summary(asym.HRmod)$tTable[1], lrc = summary(asym.HRmod)$tTable[2]),
                     control = gnlsControl(nlsTol = 50)))

print(class(stab.Mod))

fitted(stab.Mod)
AIC(stab.Mod)

laselva$stab <- unlist(fitted(stab.Mod))
head(laselva)


plot(dist ~ year, data=laselva)
curve(predict(stab.Mod, newdata = data.frame(year=x)), col = "black", add = TRUE)


## ABRUPT NONLINEAR BEHAVIOUR #### 


res_abrup<-list(data=NA)
aic_abrup<-list(data=NA)

  null.mod<-NA  
  asym.HRmod<-NA
  disp.Mod <- NA
  print(k)
  
  ## Step 1 NULL model
  null.mod <- try(gnls(dist ~ A, data = laselva, 
                       correlation = corAR1(),
                       start = c(A = mean(laselva[,'dist'])), 
                       control=nlmeControl(maxIter=1000, pnlsMaxIter=200, niterEM=400, returnObject=TRUE),verbose=F))
  
  ## Step 2 Asymp model
  asym.HRmod <- try(gnls(dist ~ Asym*(1-exp(lrc*year)), data = laselva, 
                         correlation = corAR1(),
                         start = c(Asym =summary(null.mod)$tTable[1],lrc=-0.059), 
                         control = gnlsControl(nlsTol = 100)))
 
   ## Step 3 Disp model 
  abrutp.Mod <- try(gnls(dist ~ (Asym)/(1 + exp((xmid-year)/scal)),data = laselva,
                       correlation = corAR1(), 
                       na.action = na.exclude,
                       start = c(Asym = summary(asym.HRmod)$tTable[1], xmid = 10, scal = 0.5),
                       control = gnlsControl(nlsTol = 100)))
  print(k)
  print(class(disp.Mod))
  
  fitted(abrutp.Mod)
  AIC(abrutp.Mod)
  
  laselva$abrupt <- unlist(fitted(abrutp.Mod))
  head(laselva)

  
  
  ##### Self-Starting Nls Logistic Model ####
  ## https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/SSlogis
  
  plot(dist ~ year, data=laselva)
  fit <- gnls(dist ~ SSlogis(year, Asym, xmid, scal), data=laselva)
  curve(predict( fit, newdata = data.frame(year=x)), add=TRUE)
  fitted(fit)
  AIC(fit)
  coef(fit)
  
  laselva$new <- unlist(fitted(fit))
  head(laselva)
  
  curve(predict(abrutp.Mod, newdata = data.frame(year=x)), add=TRUE)
  
  #####################################################
  ### Calculation of CC-score
  ## create a new object called "siteCC" based on object "site" above
  ## site is the updated dataframe from the gnls (nlme) fitted values 
  ## set x = siteCC$dist # this is observed distance
  ## set y_i = siteCC$fitted # this is model-fitted value
  ## CC-score will be caclulated from x and y 
  ## use y1 for linear, y2 for reversible, y3 for stable, y4 for abrupt
  
  
  laselva
  dat<-laselva
  head(dat)
  siteCC <- array(NA, dim=c(nrow(dat),7))
  colnames(siteCC)<-c("year", "x", "y1", "y2", "y3", "y4", "y5")
  siteCC<-data.frame(unlist(siteCC))
  head(siteCC)
  dim(siteCC)
  
  names(dat)
  siteCC$x<-dat$dist  

  #### select y from the columns, one at a time
  siteCC$y1<-dat$GLD 
  siteCC$y2<-dat$Rev
  siteCC$y3<-dat$stab
  siteCC$y4<-dat$abrupt 
  siteCC$y5<-dat$new
  
##CC score for gradual linear
  CC_GLD<-list(data=NA)
    num<-NA
    denom1<-NA
    denom2<-NA
    denom3<-NA
    num<-sum((siteCC$x-siteCC$y1)^2)
    denom1<-sum((siteCC$x-mean(siteCC$y1))^2)
    denom2<-sum((siteCC$y1-mean(siteCC$y1))^2)
    denom3<-length(siteCC$x)*(mean(siteCC$x)-mean(siteCC$y1))^2
    CC_GLD<- 1-((num)/(denom1+denom2+denom3))
    print(CC_GLD)

    ##CC score for reversible 
     
    CC_Rev<-list(data=NA)
      num<-NA
      denom1<-NA
      denom2<-NA
      denom3<-NA
      num<-sum((siteCC$x-siteCC$y2)^2)
      denom1<-sum((siteCC$x-mean(siteCC$y2))^2)
      denom2<-sum((siteCC$y2-mean(siteCC$y2))^2)
      denom3<-length(siteCC$x)*(mean(siteCC$x)-mean(siteCC$y2))^2
      CC_Rev<- 1-((num)/(denom1+denom2+denom3))
      print(CC_Rev)

      ##CC score for stable curve
      CC_stab<-list(data=NA)
        num<-NA
        denom1<-NA
        denom2<-NA
        denom3<-NA
        num<-sum((siteCC$x-siteCC$y3)^2)
        denom1<-sum((siteCC$x-mean(siteCC$y3))^2)
        denom2<-sum((siteCC$y3-mean(siteCC$y3))^2)
        denom3<-length(siteCC$x)*(mean(siteCC$x)-mean(siteCC$y3))^2
        CC_stab<- 1-((num)/(denom1+denom2+denom3))
        print(CC_stab)
        
        ##CC score for abrupt 
        
        CC_abrup<-list(data=NA)
          num<-NA
          denom1<-NA
          denom2<-NA
          denom3<-NA
          num<-sum((siteCC$x-siteCC$y4)^2)
          denom1<-sum((siteCC$x-mean(siteCC$y4))^2)
          denom2<-sum((siteCC$y4-mean(siteCC$y4))^2)
          denom3<-length(siteCC$x)*(mean(siteCC$x)-mean(siteCC$y4))^2
          CC_abrup<- 1-((num)/(denom1+denom2+denom3))
          print(CC_abrup)
          
          CC_fit<-list(data=NA)
          num<-NA
          denom1<-NA
          denom2<-NA
          denom3<-NA
          num<-sum((siteCC$x-siteCC$y5)^2)
          denom1<-sum((siteCC$x-mean(siteCC$y5))^2)
          denom2<-sum((siteCC$y4-mean(siteCC$y5))^2)
          denom3<-length(siteCC$x)*(mean(siteCC$x)-mean(siteCC$y5))^2
          CC_fit<- 1-((num)/(denom1+denom2+denom3))
          print(CC_fit)
        
          
  #### CC Sicrer #######
  ### Table 
          
          print(CC_GLD)
          print(CC_Rev)
          print(CC_stab)
          print(CC_abrup)
          print(CC_fit)
          
          AIC(GLD.Mod)
          AIC(Rev.Mod)
          AIC(stab.Mod)
          AIC(abrutp.Mod)
          AIC(fit)
       
          
plot(dist ~ year, data=laselva)
curve(predict(GLD.Mod, newdata = data.frame(year=x)), col = "black", add = TRUE)
curve(predict(Rev.Mod , newdata = data.frame(year=x)), col = "black", add = TRUE)
curve(predict(stab.Mod, newdata = data.frame(year=x)), col = "black", add = TRUE)
curve(predict(abrutp.Mod, newdata = data.frame(year=x)), add=TRUE)

plot(dist ~ year, data=laselva)
curve(predict(GLD.Mod, newdata = data.frame(year=x)), col = "black", add = TRUE)
curve(predict(Rev.Mod , newdata = data.frame(year=x)), col = "black", add = TRUE)
curve(predict(stab.Mod, newdata = data.frame(year=x)), col = "black", add = TRUE)
curve(predict(abrutp.Mod, newdata = data.frame(year=x)), add=TRUE)

Model <- c("Gradual", "Reversible", "Stable", "Abrutp", "Self_Starting")
CC <- c(print(CC_GLD), print(CC_Rev), print(CC_stab), print(CC_abrup), print(CC_fit))
AIC <- c(AIC(GLD.Mod),AIC(Rev.Mod),AIC(stab.Mod),AIC(abrutp.Mod),AIC(fit))

Results <- data.frame(Model,CC,AIC)
Results
