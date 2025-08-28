### Design
# rm(list=ls())
set.seed(123)

## Illustrating the "general" approach for a time-to-events data
library(DoseFinding, lib.loc = "~/RStudio/packages")


# source("powMCT1.R")

dose = c(0, 5, 25, 50, 100)
ngrp=length(dose)
surv0=0.5   ## Median survival for control
HR=0.6 #hazard ratio
surv1=surv0/HR  ## Median survival for Test arm
lambda0=log(2)/surv0
lambda1=log(2)/surv1


emax1 <- guesst(d=dose[4],p=0.5,model="emax")
emax2 <- guesst(d=dose[3],p=0.8,model="emax")
exp <- guesst(d=dose[4],p=0.1,model="exponential",Maxd=dose[ngrp])
logit <- guesst(d=c(dose[3],dose[4]),p=c(0.1,0.8),"logistic",Maxd=dose[ngrp])
betam1 <- guesst(d=dose[2],p=0.3,"betaMod", scal=120, dMax=50, Maxd=dose[ngrp])


models1 <- Mods(emax=c(emax1, emax2),
                 linear=NULL,
                 exponential=c(exp),
                 logistic=c(logit),
                 betaMod=c(betam1),
                 doses=dose,direction=c("decreasing"),
		     placEff=log(lambda0), maxEff=(log(lambda1)-log(lambda0))) 

plot(models1)

## Get the three model set ups and plot them
y0=getResp(models1, doses=dose)
y=y0[1:dim(y0)[1],]
y

## plot the dose response curve using median survival time
survtall=log(2)/exp(y0[1:dim(y0)[1],])
survtall

nshape=dim(y0)[2]
plot(dose, survtall[,1], type='l', ylab='median survival (years)')
for (i in 2:nshape) lines(dose,survtall[,i])


#########################################################################################
############ Use S for calculating contrast and power, use S0 to derive critical value ##
#########################################################################################

etotal=242
## allocation ratio nk/n0, each arm against PBO, the first number is always 1
eta<-c(1, 1, 1, 1, 1) 
#eta<-c(1, 1/2, 1/2, 1/2, 1) 

beta <-y-log(lambda0)

# 0.5 is added as correction per the paper
p0=1/apply(exp(0.5*beta)*eta,2, sum) # sum of columns
pk=(exp(0.5*beta)*eta)%*%diag(p0)

## S matrix
S=array(0,c(ngrp-1,ngrp-1,nshape)) # fill with 0s, 3-dimensional as given in c(...)
for (i in 1:nshape) {
  S[,,i]=1/p0[i]
  diag(S[,,i])=1/pk[-1,i]+1/p0[i]
}
S=S/etotal

## s0 matrix
S0=diag(sum(eta)+sum(eta)/eta[-1])
S0[upper.tri(S0)]=S0[lower.tri(S0)]=sum(eta)
S0=S0/etotal
dosPlac <- dose[-1]

model1 <- Mods(emax=emax1, doses=dose, placEff=log(lambda0), maxEff=(log(lambda1)-log(lambda0))) 
model2 <- Mods(emax=emax2, doses=dose, placEff=log(lambda0), maxEff=(log(lambda1)-log(lambda0))) 
model3 <- Mods(linear=NULL, doses=dose, placEff=log(lambda0), maxEff=(log(lambda1)-log(lambda0))) 
model4 <- Mods(exponential=exp, doses=dose, placEff=log(lambda0), maxEff=(log(lambda1)-log(lambda0))) 
model5 <- Mods(logistic=logit, doses=dose, placEff=log(lambda0), maxEff=(log(lambda1)-log(lambda0))) 
model6 <- Mods(betaMod=betam1, doses=dose, placEff=log(lambda0), maxEff=(log(lambda1)-log(lambda0))) # ist doch oben schon alles definiert???

contrast1 <- optContr(model1, doses=dosPlac, S=S[,,1], placAdj=T)$contMat
contrast2 <- optContr(model2, doses=dosPlac, S=S[,,2], placAdj=T)$contMat
contrast3 <- optContr(model3, doses=dosPlac, S=S[,,3], placAdj=T)$contMat
contrast4 <- optContr(model4, doses=dosPlac, S=S[,,4], placAdj=T)$contMat
contrast5 <- optContr(model5, doses=dosPlac, S=S[,,5], placAdj=T)$contMat
contrast6 <- optContr(model6, doses=dosPlac, S=S[,,6], placAdj=T)$contMat

contrast <- cbind(contrast1, contrast2, contrast3, contrast4, contrast5, contrast6)

### powMCT is doseFinding package is modified to take different covariance matrix 
### for different dose response shapes. df is set to Inf for normal approximation
mcpmod.power1=powMCT1(contrast, altModels = model1, alpha = 0.05, S=S[,,1], S0=S0, placAdj=T, df=Inf)
mcpmod.power2=powMCT1(contrast, altModels = model2, alpha = 0.05, S=S[,,2], S0=S0, placAdj=T, df=Inf)
mcpmod.power3=powMCT1(contrast, altModels = model3, alpha = 0.05, S=S[,,3], S0=S0, placAdj=T, df=Inf)
mcpmod.power4=powMCT1(contrast, altModels = model4, alpha = 0.05, S=S[,,4], S0=S0, placAdj=T, df=Inf)
mcpmod.power5=powMCT1(contrast, altModels = model5, alpha = 0.05, S=S[,,5], S0=S0, placAdj=T, df=Inf)
mcpmod.power6=powMCT1(contrast, altModels = model6, alpha = 0.05, S=S[,,6], S0=S0, placAdj=T, df=Inf)


mcpmod.power=c(mcpmod.power1, mcpmod.power2, mcpmod.power3, mcpmod.power4, mcpmod.power5, mcpmod.power6)
round(mcpmod.power,3)
final1 <- round(mean(mcpmod.power),3)

###################################################################

mcpmod.power12=powMCT2(contrast, altModels = model1, alpha = 0.05, S=S[,,1], S0=S0, placAdj=T, df=Inf)
mcpmod.power22=powMCT2(contrast, altModels = model2, alpha = 0.05, S=S[,,2], S0=S0, placAdj=T, df=Inf)
mcpmod.power32=powMCT2(contrast, altModels = model3, alpha = 0.05, S=S[,,3], S0=S0, placAdj=T, df=Inf)
mcpmod.power42=powMCT2(contrast, altModels = model4, alpha = 0.05, S=S[,,4], S0=S0, placAdj=T, df=Inf)
mcpmod.power52=powMCT2(contrast, altModels = model5, alpha = 0.05, S=S[,,5], S0=S0, placAdj=T, df=Inf)
mcpmod.power62=powMCT2(contrast, altModels = model6, alpha = 0.05, S=S[,,6], S0=S0, placAdj=T, df=Inf)


mcpmod.power2=c(mcpmod.power12, mcpmod.power22, mcpmod.power32, mcpmod.power42, mcpmod.power52, mcpmod.power62)
round(mcpmod.power2,3)
final2 <- round(mean(mcpmod.power2),3)


final1
final2