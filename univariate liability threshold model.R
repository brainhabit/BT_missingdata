#==========
# univariate twin analysis 
# used to test liability threshold models for experiment-level missing
#
# author: Giorgia Bussu
# project: BT missing data
# version: June 2024
#==========
rm(list=ls())


require(OpenMx)
source('C:/Users/myfolder/twin_modelling/tutorial/miFunctions.R')

###############################################################################################################
### prepare data:
### import and check data

data <- read.csv(file='twin_bin_sensitivity_data.csv',header=T,
                 sep=',')
names(data);dim(data)

### select variable for analysis:

Vars <- 'fomo_bin'
nv <- 1
ntv <- nv*2
nth <- 1
(selVars <- paste(Vars,c(rep(1,nv),rep(2,nv)),sep=''))

### select subsets for analysis:

mz <- subset(data,zygosity=='MZ',selVars)
dz <- subset(data,zygosity=='DZ',selVars)

# convert the variable to a factor:

mz <- mxFactor(x=mz[,selVars],levels=c(0:nth))
dz <- mxFactor(x=dz[,selVars],levels=c(0:nth))

##### assumptions testing and twin correlations #####

### fit a fully saturated model

# matrix for the mean, which is fixed to 0:

expMean <- mxMatrix(type='Zero',nrow=1,ncol=ntv,name='ExpMean')

# matrices for the thresholds:

expThreshMZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=1,
                        labels=labFull('th_mz',1,ntv),
                        name='ExpThreshMZ')
expThreshDZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=1,
                        labels=labFull('th_dz',1,ntv),
                        name='ExpThreshDZ')

# twin correlations:

expCorMZ <- mxMatrix(type='Stand',nrow=ntv,ncol=ntv,free=T,values=.5,
                     label='cor_mz',lbound=-.9999,ubound=.9999,
                     name='ExpCorMZ')
expCorDZ <- mxMatrix(type='Stand',nrow=ntv,ncol=ntv,free=T,values=.5,
                     label='cor_dz',lbound=-.9999,ubound=.9999,
                     name='ExpCorDZ')

# specify observed data:

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# specify objectives:

objMZ <- mxExpectationNormal(covariance='ExpCorMZ',
                             means='ExpMean',dimnames=selVars,
                             thresholds='ExpThreshMZ')
objDZ <- mxExpectationNormal(covariance='ExpCorDZ',
                             means='ExpMean',dimnames=selVars,
                             thresholds='ExpThreshDZ')

# specify estimation method:

funcML <- mxFitFunctionML()

# specify the data groups:

modelMZ <- mxModel('MZ',expMean,expThreshMZ,expCorMZ,dataMZ,objMZ,
                   funcML)
modelDZ <- mxModel('DZ',expMean,expThreshDZ,expCorDZ,dataDZ,objDZ,
                   funcML)

# combine data groups:

multi <- mxFitFunctionMultigroup(c('MZ','DZ'))

# confidence intervals:

ci <- mxCI(c('MZ.ExpCorMZ','DZ.ExpCorDZ'))

# combine all model objects:

SatModel <- mxModel('Sat',modelMZ,modelDZ,ci,multi)

# fit the model:

SatFit <- mxTryHardOrdinal(SatModel,intervals=T)
summary(SatFit)

### test assumptions

# assumption 1: equal thresholds within twin pairs:

SatModel2 <- mxModel(SatFit,name='Sat2')
SatModel2 <- omxSetParameters(SatModel2,labels=labFull('th_mz',nth,ntv),
                              free=T,values=1,newlabel='th_mz')
SatModel2 <- omxSetParameters(SatModel2,labels=labFull('th_dz',nth,ntv),
                              free=T,values=1,newlabel='th_dz')
SatFit2 <- mxTryHardOrdinal(SatModel2,intervals=T)
summary(SatFit2)

# assumption 2: 

SatModel3 <- mxModel(SatFit2,name='Sat3')
SatModel3 <- omxSetParameters(SatModel3,labels=c('th_mz','th_dz'),free=T,
                              values=1,newlabel='th')
SatFit3 <- mxTryHardOrdinal(SatModel3,intervals=T)
summary(SatFit3)

# compare the fit of the models:

mxCompare(SatFit,c(SatFit2,SatFit3))

# significance correlations
# assumption 2: 

SatModel_corr <- mxModel(SatFit3,name='SatCorr')
SatModel_corr <- omxSetParameters(SatModel_corr,labels=c('cor_dz'),free=F,
                              values=0)
SatFit_corr <- mxTryHardOrdinal(SatModel_corr,intervals=T)
summary(SatFit_corr)

mxCompare(SatFit3,SatFit_corr)


##### ACE model #####

### full ACE model

# path coefficients:

pathA <- mxMatrix(type='Full',nrow=nv,ncol=nv,free=T,values=1,label='a_1_1',
                  name='a')
pathC <- mxMatrix(type='Full',nrow=nv,ncol=nv,free=T,values=1,label='c_1_1',
                  name='c')
pathE <- mxMatrix(type='Full',nrow=nv,ncol=nv,free=T,values=1,label='e_1_1',
                  name='e')

# calculate variance components:

varA <- mxAlgebra(a^2,name='A')
varC <- mxAlgebra(c^2,name='C')
varE <- mxAlgebra(e^2,name='E')

# calculate the total variance:

varP <- mxAlgebra(A+C+E,name='V')

# constrain the total variance to equal 1:

matU <- mxMatrix(type='Unit',nrow=nv,ncol=nv,name='U')
varConst <- mxConstraint(V==U,name='VarConst')

# mean:

expMean <- mxMatrix(type='Zero',nrow=1,ncol=ntv,name='ExpMean')

# threshold:

expThresh <- mxMatrix(type='Full',nrow=nth,ncol=ntv,free=T,values=1,label='th',
                      name='ExpThresh')

# variance/covariance:

expCovMZ <- mxAlgebra(rbind(cbind(V,A+C),
                            cbind(A+C,V)),name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V,0.5%x%A+C),
                            cbind(0.5%x%A+C,V)),name='ExpCovDZ')

# convert the variance components to proportions:

estVC <- mxAlgebra(cbind(A/V,C/V,E/V),name='EstVC')

# observed data:

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

#  objectives:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMean',dimnames=selVars,
                             thresholds='ExpThresh')
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMean',dimnames=selVars,
                             thresholds='ExpThresh')
funcML <- mxFitFunctionML()

# specify data groups:

pars <- list(pathA,pathC,pathE,varA,varC,varE,varP,matU,expMean,expThresh)

modelMZ <- mxModel('MZ',pars,varConst,expCovMZ,estVC,dataMZ,objMZ,funcML)
modelDZ <- mxModel('DZ',pars,varConst,expCovDZ,estVC,dataDZ,objDZ,funcML)

# combine data groups:

multi <- mxFitFunctionMultigroup(c('MZ','DZ'))

# confidence intervals:

ci <- mxCI(c('MZ.EstVC'))

# combine all model objects:

ModelACE <- mxModel('ACE',modelMZ,modelDZ,multi,ci)

# fit the model:

FitACE <- mxTryHardOrdinal(ModelACE,intervals=T)
summary(FitACE)

# compare fit to the saturated model:

mxCompare(SatFit,FitACE)

### nested models:

# AE model:

ModelAE <- mxModel(FitACE,name='AE')
ModelAE <- omxSetParameters(ModelAE,label='c_1_1',free=F,values=0)
FitAE <- mxTryHardOrdinal(ModelAE,intervals=T)
aefit<-summary(FitAE)

# CE model:

ModelCE <- mxModel(FitACE,name='CE')
ModelCE <- omxSetParameters(ModelCE,label='a_1_1',free=F,values=0)
FitCE <- mxTryHardOrdinal(ModelCE,intervals=T)
cefit<-summary(FitCE)

# E model:

ModelE <- mxModel(FitACE,name='E')
ModelE <- omxSetParameters(ModelE,labels=c('a_1_1','c_1_1'),free=F,values=0)
FitE <- mxTryHardOrdinal(ModelE,intervals=T)
summary(FitE)

# compare fit with the ACE model:

mxCompare(FitACE,c(FitAE,FitCE,FitE))

# AC for p purposes
ModelAC <- mxModel(FitACE,name='AC')
ModelAC <- omxSetParameters(ModelAC,label='e_1_1',free=F,values=0)
FitAC <- mxTryHardOrdinal(ModelAC,intervals=T)
summary(FitAC)
































