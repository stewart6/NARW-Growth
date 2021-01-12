# Supplemental code for Stewart et al., 2021
# Decreasing body sizes are a troubling trend in North Atlantic Right Whales
# This code will growth curve analyses for NARW photogrammetry data from 2000-2019.

# Required packages:
library(R2jags)

# Required Data Files:
# NARW Lengths.RData
# NARW_2PhaseGomp_Cov_PPred.jags

#######################
# 1) LOAD NARW_Data FILES
#######################

# Load the NARW Length and covariate data
# (make sure this .RData file is in your working directory)
load("NARW Lengths.RData")
# You should now have the NARW_Data object loaded in your R environment

# Notes: WhaleID corresponds to the North Atlantic Right Whale Consortium ID database.
# SurveyYear is the year that the photogrammetry image / measurement was taken.
# TL is the total length (in meters) of the measured whale.
# Type corresponds to the altitude measurement method (Radar, GPS, or Laser altimeter).
# BirthYear is the calendar year of birth, whereas BY is the birth year - 1980, for input as a covariate in the model.
# Births is the number of successful births by a measured whale prior to age 10 or prior to the measurement date.
# MomEntang is whether a measured whale's mother had evidence of new a severe entanglement injury while the measured whale was a dependent calf (1 = yes, 0 = no).
# MidGear is the midpoint between the minimum and maximum estimated durations of observed entanglements with attached gear, cumulative prior to age 10.
# See supplemental methods for details on covariates.


#############################
# 2) RUN THE MODEL
#############################


## Standardize the covariate data
Cov_Mult <- data.frame(MomEntang=NARW_Data$MomEntang, # Already formatted as min=0,  max=1, because this is modeled as a fixed effect
                       NBirths=NARW_Data$Births/max(NARW_Data$Births), # Divide observed number of births by maximum number of births so min=0, max=1
                       Gear=(NARW_Data$MidGear/max(NARW_Data$MidGear))) # Divide observed gear entanglement duration by max observed duration so min=1, max=1


## Make a list object including all input data for the model
jags.data <- list(Age = NARW_Data$Age, #Known age
                  TL = NARW_Data$TL, #Measured length
                  Nobs = length(NARW_Data$Age), #Number of length measurements
                  obs.types = max(as.numeric(as.factor(NARW_Data$Type))), #Number of different altitude measurement types
                  type = as.numeric(as.factor(NARW_Data$Type)), #Altimeter type turned into a factor variable (1,2,3)
                  Phase = NARW_Data$Phase, #Growth phase that each measurement belongs to (age <1 or >1)
                  Cov = Cov_Mult, #Standardized covariates (excluding birth year)
                  Ncov = dim(Cov_Mult)[2], #Number of covariates in the Cov data frame
                  BY = NARW_Data$BY, #Birth year
                  n.years = max(NARW_Data$BY)) #Number of birth years


## Specify MCMC settings
ni <- 100000
nt <- 50
nb <- 50000
nc <- 3

## Specify which parameters to track
parameters <- c("Linf","C","K","sigma","loglik","TL.mean","L.slope","L.BY.slope","Post.Pred")

## Run the model
NARWGom2Phase_PP <- jags(jags.data, inits=NULL, parameters, "NARW_2PhaseGomp_Cov_PPred.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


#############################
# 3) VISUALIZE RESULTS
#############################

# Attach model output to R environment
attach.jags(NARWGom2Phase_PP)



### Growth Curves (Figure 2)

## Custom function for plotting predicted two-phase Gompertz growth curves based on posterior distributions of L, K and C, plus covariate effects on L (birth year and/or entanglement)
plot2PhaseGOMCI<-function(xvals,Linf,K,COL="#C1CDCD50",Lslope=0,L.BY.slope=0,Entang=0,BY=0,DENS=NA,WEIGHT=1){
  range1 <- seq(xvals[1],xvals[2],by=0.05)
  range2 <- seq(xvals[2],xvals[3],by=0.05)
  
  ymat1 <- array(dim=c(dim(Linf)[1],length(range1)))
  ymat2 <- array(dim=c(dim(Linf)[1],length(range2)))
  
  if(length(L.BY.slope)==1) L.BY.slope=rep(0,length(Linf))
  if(length(Lslope)==1) Lslope=rep(0,length(Linf))
  
  
  Linf.mu <- array(dim=dim(Linf))
  for(i in 1:dim(Linf)[1]){
    for(p in 1:dim(Linf)[2]){
      Linf.mu[i,p] <- Linf[i,p] + Lslope[i]*Entang + L.BY.slope[i]*BY
    }#p
    ymat1[i,] <- Linf.mu[i,1]*exp(-C[i,1]*exp(-K[i,1]*range1))
    ymat2[i,] <- Linf.mu[i,2]*exp(-C[i,2]*exp(-K[i,2]*range2))
    
  }#i
  
  for(p in 1:dim(Linf)[2]){
    assign(paste0("UCI",p), apply(get(paste0("ymat",p)),2,quantile,probs=c(0.975)) )
    assign(paste0("LCI",p), apply(get(paste0("ymat",p)),2,quantile,probs=c(0.025)) )
    polygon(x=c(get(paste0("range",p)),get(paste0("range",p))[length(get(paste0("range",p))):1]), y=c(get(paste0("LCI",p)),get(paste0("UCI",p))[length(get(paste0("UCI",p))):1]),col=COL,border = NA, density=DENS,lwd=WEIGHT, lend=2)
    lines(x=get(paste0("range",p)),y=apply(get(paste0("ymat",p)),2,median))
  }
}


## Add some columns to the data frame for plotting:

# The AgeCEX column sets custom point sizes for length measurements based on their birth year
# Note there are three birth year groups for visual clarity, although the model estimated effects for every birth year
NARW_Data$AgeCEX <- 0.75
NARW_Data$AgeCEX[which(NARW_Data$BirthYear<1994)] <- 0.3
NARW_Data$AgeCEX[which(NARW_Data$BirthYear>2006)] <- 1.25

# Create a Gear column that corresponds to the standardized model input data for attached gear entanglement duration
NARW_Data$Gear<-(NARW_Data$MidGear/max(NARW_Data$MidGear))


## Plot figure:

# Set plotting parameters so there are 4 panels
par(mfrow=c(2,2))

#Panel 1, Birth Year Effect, No Entanglement
par(mai=c(0.2,0.8,0.6,0))
plot(x=NARW_Data$Age,y=NARW_Data$TL,pch=19,xlab="",ylim=c(4,15), cex=NARW_Data$AgeCEX, main="", xaxt='n',ylab="Total Length (m)",cex.axis=1.5,cex.lab=1.5)
axis(1, labels=FALSE)
plot2PhaseGOMCI(xvals = c(0,1,38),Linf=Linf,K=K,COL="#888B95") # Baseline growth curve with no birth year or entanglement effects
plot2PhaseGOMCI(xvals = c(0,1,38),Linf=Linf,K=K,L.BY.slope=L.BY.slope,BY=40,COL="#FF7F2490") # Growth curve with birth year effects
points(x=NARW_Data$Age,y=NARW_Data$TL,pch=19,cex=NARW_Data$AgeCEX,col="black")
legend(x=10,y=12,legend=c("1981-1993","1994-2006","2007-2019"),pch=19,pt.cex=c(0.3,0.75,1.25),title="Birth Year",bty='n',x.intersp=0.5,cex=1.5)
text(x=0.5,y=14.7,labels="(A)",cex=1.5)

#Panel 2, No Birth Year Effect, Maximum Entanglement Duration
par(mai=c(0.2,0.2,0.6,0.6))
plot(x=NARW_Data$Age,y=NARW_Data$TL,pch=19,ylab="",xlab="",ylim=c(4,15), cex=NARW_Data$AgeCEX, main="",col="gray70",xaxt='n',yaxt='n',cex.axis=1.5,cex.lab=1.5)
axis(1,labels=FALSE)
axis(2,labels=FALSE)
plot2PhaseGOMCI(xvals = c(0,1,38),Linf=Linf,K=K,COL="#888B95") # Baseline growth curve with no birth year or entanglement effects
plot2PhaseGOMCI(xvals = c(0,1,38),Linf=Linf,K=K,Lslope=L.slope[,3],Entang=1,COL="#00308C90") # Growth curve with entanglement effects
points(x=NARW_Data$Age,y=NARW_Data$TL,pch=19,cex=NARW_Data$AgeCEX,col="gray70")
points(x=NARW_Data$Age[which(NARW_Data$Gear>0)],y=NARW_Data$TL[which(NARW_Data$Gear>0)],cex=NARW_Data$AgeCEX[which(NARW_Data$Gear>0)],col="black",pch=19)
legend(x=10,y=10,legend=c("Whale Not Entangled","Whale Entangled"),pch=19,pt.cex=1,col=c("gray70","black"),bty='n',x.intersp=0.5,cex=1.5)
text(x=0.5,y=14.7,labels="(B)",cex=1.5)

#Panel 3, No Birth Year Effect, Maternal Entanglement
par(mai=c(0.8,0.8,0,0))
plot(x=NARW_Data$Age,y=NARW_Data$TL,pch=19,xlab="Age",ylim=c(4,15), cex=NARW_Data$AgeCEX, main="",col="gray70",ylab="Total Length (m)",cex.axis=1.5,cex.lab=1.5)
plot2PhaseGOMCI(xvals = c(0,1,38),Linf=Linf,K=K,COL="#888B95") # Baseline growth curve with no birth year or entanglement effects
plot2PhaseGOMCI(xvals = c(0,1,38),Linf=Linf,K=K,Lslope=L.slope[,1],Entang=1,COL="#6EADF390") # Growth curve with maternal entanglement effects
points(x=NARW_Data$Age,y=NARW_Data$TL,pch=19,cex=NARW_Data$AgeCEX,col="gray70")
points(x=NARW_Data$Age[which(NARW_Data$MomEntang==1)],y=NARW_Data$TL[which(NARW_Data$MomEntang==1)],cex=NARW_Data$AgeCEX[which(NARW_Data$MomEntang==1)],col="black",pch=19)
legend(x=10,y=10,legend=c("Mother Not Entangled","Mother Entangled"),pch=19,pt.cex=1,col=c("gray70","black"),bty='n',x.intersp=0.5,cex=1.5)
text(x=0.5,y=14.7,labels="(C)",cex=1.5)

#Panel 4, Birth Year Effect Plus Maximum Entanglement Duration
par(mai=c(0.8,0.2,0,0.6))
plot(x=NARW_Data$Age,y=NARW_Data$TL,pch=19,ylab="",xlab="Age",ylim=c(4,15), cex=NARW_Data$AgeCEX,col="gray70", main="",yaxt='n',cex.axis=1.5,cex.lab=1.5)
axis(2,labels=F)
plot2PhaseGOMCI(xvals = c(0,1,38),Linf=Linf,K=K,Lslope=L.slope,COL="#888B95") # Baseline growth curve with no birth year or entanglement effects
plot2PhaseGOMCI(xvals = c(0,1,38),Linf=Linf,K=K,Lslope=L.slope[,3],Entang=1,L.BY.slope=L.BY.slope,BY=40,COL="#FF7F2490") # Growth curve with birth year and entanglement effects
points(x=NARW_Data$Age,y=NARW_Data$TL,pch=19,cex=NARW_Data$AgeCEX,col="gray70")
points(x=NARW_Data$Age[which(NARW_Data$Gear>0)],y=NARW_Data$TL[which(NARW_Data$Gear>0)],cex=NARW_Data$AgeCEX[which(NARW_Data$Gear>0)],col="black",pch=19)
points(x=NARW_Data$Age[which(NARW_Data$MomEntang==1)],y=NARW_Data$TL[which(NARW_Data$MomEntang==1)],cex=NARW_Data$AgeCEX[which(NARW_Data$MomEntang==1)],col="black",pch=19)
legend(x=10,y=10,legend=c("Whale or Mother Not Entangled","Whale or Mother Entangled"),pch=19,pt.cex=1,col=c("gray70","black"),bty='n',x.intersp=0.5,cex=1.5)
text(x=0.5,y=14.7,labels="(D)",cex=1.5)

# Note: Figure 2D in the manuscript has a blue/orange striped curve; the blue stripes were added manually.





### Covariate Effects Violin Plot (Figure 3)

## Combine all covariate effects into a single data frame
Cov_Effects <- as.data.frame(cbind(L.slope,L.BY.slope*39)) # Birth Year slope is multiplied by 39 so that all covariates are plotted with their maximum effect sizes
colnames(Cov_Effects) <- c("Mother Entangled","Number of Lactations","Gear Entanglement Duration","Birth Year")

## Format for ggplot
Covs_AS <- gather(Cov_Effects,key="Covariate",value="Slope")

## Create violin plot ggplot object
CovViolinPlot <- ggplot(data=transform(Covs_AS,Covariate=factor(Covariate,levels=c("Birth Year","Gear Entanglement Duration","Mother Entangled","Number of Lactations"))),aes(x=Covariate,y=Slope,fill=Covariate)) +
  geom_violin(color=NA,trim=TRUE,scale="width") + 
  scale_fill_manual(values=c("#FF7F2490","#00308C90","#6EADF390","#1D696390")) +
  labs(x="Covariate",y="Effect on Asymptotic Length (m)") +
  geom_boxplot(width=0.15, fill='white',outlier.shape = NA)+
  ylim(-3,2) +
  theme_minimal() + theme(axis.line = element_line(colour = "black"),
                          text = element_text(size = 16),
                          panel.grid.minor.y = element_blank(),
                          legend.position = "none") +
  geom_hline(yintercept = 0,size=1)

## Plot violin plots
CovViolinPlot





### Posterior Predictive Checks (Figure S2)

## Randomly select 20 observed lengths
PPSelect <- sample(1:202,size=20)

## Set plot parameters
par(mfrow=c(5,4))
par(mar=rep(2,4))

## Plot predicted and observed lengths
for(i in 1:length(PPSelect)){
  
  SelectedObs <- PPSelect[i]
  hist(Post.Pred[,PPSelect[i]],main="",xlab="Predicted Total Length") # Plot posterior prediction distribution
  abline(v=NARW_Data$TL[SelectedObs],col="red",lwd=2) # Plot observed length
  abline(v=quantile(Post.Pred[,PPSelect[i]],c(0.025,0.975)),lty=2) # Plot 95% prediction intervals
  
}








