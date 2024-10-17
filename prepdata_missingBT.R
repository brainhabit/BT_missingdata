#=========================================================================================================
# Main for project: Genetic influences on missing data across experimental measures in infancy 
#
# used to prepare data for subsequent twin analysis
#
# author: Giorgia Bussu
# project: BT missing data
# version: June 2024
#==========================================================================================================

rm(list=ls())

### set a working directory (otherwise create a new R project):

setwd('C:/Users/myfolder/missingdata_BATSS')
list.files() # show files in the working directory

### import and check data:

data_fomo <- read.csv('batss_fomo.csv')
names(data_fomo);dim(data_fomo)

data_emi <- readxl::read_xlsx('Missingdata_EMI.xlsx')
names(data_emi);dim(data_emi)

data_plr <- read.csv('BT_PLR.csv')
data_plr<-data_plr[-1]
names(data_plr);dim(data_plr)

data_demo <- readxl::read_xlsx('C:/Users/myfolder/Demographics_5m.xlsx')
names(data_demo);dim(data_demo)

data_exclusion<-readxl::read_xlsx('C:/Users/myfolder/BabyTwins background summary and exclusion 20201102.xlsx')
names(data_exclusion); dim(data_exclusion)

data_background<- readxl::read_xlsx('C:/Users/myfolder/Background_5m.xlsx')

data_pregnancy<- readxl::read_xlsx('C:/Users/myfolder/pregnancy_time_data.xlsx')

### matching participants across files
indx<-match(data_emi$filename,data_demo$`Participant EASE Code`)
data_demo<-data_demo[indx,]

indx_ex<-match(data_emi$filename,data_exclusion$Code)
data_exclusion<-data_exclusion[indx_ex,]

indx_back<-match(data_emi$filename,data_background$Kod)
data_background<-data_background[indx_back,]

indx_preg<-match(data_emi$filename,data_pregnancy$ID)
data_pregnancy<-data_pregnancy[indx_preg,]

indx_fomo<-match(data_emi$filename,data_fomo$subject)

indx_plr<-match(data_emi$filename,data_plr$id)

### creating separate data frame for analysis, collating binary/continuous missing data info for each of the 3 experiments

# ET EMI
my_data<-data_emi

# EEG GMP
my_data$fomo_cont <- rep(NA,length(my_data$valid_trials))

for(i in 1:length(indx_fomo)){if(!is.na(indx_fomo[i])){my_data$fomo_cont[i]=data_fomo$motion2Hz_nrVEPs[indx_fomo[i]]}}

my_data$fomo_bin <- rep(NA,length(my_data$valid_trials))

for(i in 1:length(indx_fomo)){if(!is.na(indx_fomo[i])){my_data$fomo_bin[i]=data_fomo$pred_ok[indx_fomo[i]]}}

names(my_data)<-c('id','emi_cont','emi_bin','fomo_cont','fomo_bin')

# ET PLR
my_data$plr_cont <- data_plr$nrValidTrials[indx_plr]

my_data$plr_bin <- data_plr$binary[indx_plr]


### exclusion based on general criteria

excl_indx<-which(data_exclusion$`Exclusion version A`==1)

my_data<-my_data[-excl_indx,]

# after removal based on general exclusion criteria, I have data from every remaining participant for the EMI, while there is 34 missing data for GMP and 75 missing data for the GAP
# now set those to 0 since those have been excluded anyways

my_data$fomo_bin[is.na(my_data$fomo_bin)]<-0
my_data$plr_bin[is.na(my_data$plr_bin)]<-0

# re-coding to 1= missing and 0 = not missing
my_data$fomo_bin<-!my_data$fomo_bin
my_data$plr_bin<-!my_data$plr_bin
my_data$emi_bin<-!my_data$emi_bin

# composite missing data variable for univariate twin modelling
my_data$composite_missing<-my_data$emi_bin+my_data$fomo_bin+my_data$plr_bin

# overall 348 infants provided no missing data, 185 missing in 1 experiment only, 54 missing in 2 experiments, 7 missing in all experiments
# based on the distribution, we decided to go for a binary measure 0 = not missing vs 1 = missing in any 1 or more experiments

my_data$composite_bin <- ifelse(my_data$composite_missing>0,1,my_data$composite_missing)


# clean data files from excluded participants
data_demo_clean<-data_demo[-excl_indx,]
data_background_clean<-data_background[-excl_indx,]
data_pregnancy_clean<-data_pregnancy[-excl_indx,]

# check pregnancy
preg_day<-data_pregnancy_clean$Gdays
preg_day[which(is.na(preg_day))]<-0

pregnancy_term<-data_pregnancy_clean$Gweeks*7+preg_day

### dataset included in our analyses

data<-data.frame(cbind(my_data,data_demo_clean$Gender,data_demo_clean$`Age at Date of Assessment (days)`,pregnancy_term,data_demo_clean$`Bio Mum Age`,data_demo_clean$`Bio Dad Age`,data_demo_clean$`A. Highest level of education`,data_demo_clean$`B. Highest level of education`,data_background_clean$`F16. Ungefär hur hög är familjens* gemensamma månadsinkomst innan skatt (lön och andra ersättningar/bidrag)?`,data_demo_clean$TWAB,data_demo_clean$`Twin pair no.`))
names(data)<-c(names(my_data),'sex','age','term_age','Mum_age','Dad_age','A_edu','B_edu','family_income','TWAB','Twinpair')

# parental education level to MAX within family
data$A_edu_num<-as.numeric(as.factor(data$A_edu))
data$A_edu_num[which(data$A_edu_num==4)]<-6
data$A_edu_num[which(data$A_edu_num==5)]<-4
data$A_edu_num[which(data$A_edu_num==6)]<-5

data$B_edu_num<-as.numeric(as.factor(data$B_edu))
data$B_edu_num[which(data$B_edu_num==4)]<-6
data$B_edu_num[which(data$B_edu_num==5)]<-4
data$B_edu_num[which(data$B_edu_num==6)]<-5

#data$edu_max<-pmax(data$A_edu_num,data$B_edu_num)
data$edu_mean<-rowMeans(cbind(data$A_edu_num,data$B_edu_num),na.rm = T)


# parental age to mean across mum and dad (NA.RM=TRUE for 3 missing Dad age)
data$Mum_age<-as.numeric(data$Mum_age)
data$Dad_age<-as.numeric(data$Dad_age)
data$parental_age<-rowMeans(cbind(data$Mum_age,data$Dad_age),na.rm = T)

## add zygosity
data_exclusion_clean<-data_exclusion[-excl_indx,]
data$zygosity<-data_exclusion_clean$Zygosity

names(data);dim(data)

table(data$zygosity)

## table demo and descriptives ##

library(tidyverse)
library(finalfit)

data$sex_factor<-as.factor(data$sex)
data$family_income[which(data$family_income=='vet ej')]<-'unknown'
data$income_factor<-as.factor(data$family_income)
data$zyg_factor<-as.factor(data$zygosity)
explanatory <- c('age','sex_factor','term_age','edu_factor','income_factor','parental_age')
data$edu_factor<-as.factor(round(data$edu_mean))
data <- data %>%
  mutate(
    income_factor = ff_label(income_factor, "Family income")
  )
data <- data %>%
  mutate(
    edu_factor = ff_label(edu_factor, "Mean parental education")
  )
data <- data %>%
  mutate(
    term_age = ff_label(term_age, "Gestation age (days)")
  )
data <- data %>%
  mutate(
    age = ff_label(age, "Age (days)")
  )
data <- data %>%
  mutate(
    sex_factor = ff_label(sex_factor, "Sex")
  )
table1 <- data %>%
  summary_factorlist('zyg_factor', explanatory,
                     p=TRUE, na_include=TRUE,
                     add_dependent_label=TRUE,
                     dependent_label_prefix = "Zygosity: "
  )
write.table(table1, file = "table1.txt", sep = ",", quote = FALSE, row.names = F)

# check how many pairs are incomplete by making a variable that counts 
#   the frequency of their pair ID:

data <- merge(data,data.frame(table(Twinpair=data$Twinpair)),by='Twinpair')
table(data$Freq) # check how many pair IDs appear only once 

data$Twinpair[which(data$Freq==1)]
data_test<-rbind(data[1:97,],data[97:188,],data[188:293,],data[293:532,],data[532:length(data$id),])
data_test <- merge(data_test,data.frame(table(id=data_test$id)),by='id')
naindx<-matrix(which(data_test$Freq.y==2),nrow=2)[1,]


data_test$emi_cont[naindx]<-NA
data_test$emi_bin[naindx]<-NA
data_test$fomo_cont[naindx]<-NA
data_test$fomo_bin[naindx]<-NA
data_test$plr_cont[naindx]<-NA
data_test$plr_bin[naindx]<-NA

data_test$composite_missing[naindx]<-NA
data_test$composite_bin[naindx]<-NA

data_test$TWAB[naindx[c(1)]]<-2
data_test$TWAB[naindx[c(2,3,4)]]<-1

data<-data_test

dim(data) 

# remove frequency variables
data <- data[-c(25,26)]

##### transform data for twin analysis #####

# binary sex: 0=Females; 1=Males.
data$sex[which(data$sex=='Male')]<-1
data$sex[which(data$sex=='Female')]<-0
data$sex<-as.numeric(data$sex)

# ordinal discrete income
data$income<-as.numeric(as.factor(data$family_income))
data$income[which(data$income==2)]<-12
data$income[which(data$income==1)]<-2
data$income[which(data$income==11)]<-1
data$income[which(data$income==12)]<-11

# remove char income
data<-data[,-c(18)]

# remove unused variables
data<-data[,-c(14:17,19:20)]

# numeric variables
vars <- colnames(data)[c(2:16,18)] 
data[vars]<-lapply(data[vars],as.numeric)

### check the distribution of the missing data scales:

vars <- colnames(data)[c(3:10)]
lapply(data[vars],psych::describe)

## composite binary variable skewness = 0.44, continuous ones around -1

## no need to scale the composite missing binary variable

# transforming continuous variables to meet criteria of skewness around 0.3

data$emi_t <- I(data$emi_cont)^2
data$fomo_t <- I(data$fomo_cont)^2
data$plr_t <- I(data$plr_cont)^2

####### check potential covariate effects

library(drgee)

summary(gee( formula = composite_bin~age,data = data, clusterid = "Twinpair", cond = F))
# no significant age effect 

summary(gee( formula = composite_bin~sex,data = data, clusterid = "Twinpair", cond = F))
# no significant effect of sex

summary(gee( formula = composite_bin~parental_age,data = data, clusterid = "Twinpair", cond = F))
# no significant effect of parental age

summary(gee( formula = composite_bin~term_age,data = data, clusterid = "Twinpair", cond = F))
# significant effect of term age (b=-.007, se=.003, p=.022)

summary(gee( formula = composite_bin~edu_mean,data = data, clusterid = "Twinpair", cond = F))
# no significant effect of education

summary(gee( formula = composite_bin~income,data = data, clusterid = "Twinpair", cond = F))
# no significant effect of income 

## scale key variables
data$emi_ts <- scale(data$emi_t)
data$fomo_ts <- scale(data$fomo_t)
data$plr_ts <- scale(data$plr_t)

############################ Randomized split of twin pairs for twin format ################################################

### give each person a random number (only use two possible numbers, and never 0 and 1):

set.seed(2023)
npairs <- nrow(data)/2
rand <- c(rep(4,npairs),rep(5,npairs))
data$rand <- sample(rand) # randomly choose 1 number for each person

### divide the twins into two subsets: twin 1 and twin 2

twin1 <- subset(data,TWAB==1)
twin2 <- subset(data,TWAB==2)

### now we need to flip twin 2's random number so that it is the opposite of twin 1:

randVar <- data.frame(cbind(twin1$Twinpair,twin1$rand))
colnames(randVar) <- c('Twinpair','rand')

# now take away twin 2's random number and replace it with the twin 1 random number:

twin2 <- twin2[-c(25)] # before working on continuous measures was 21
twin2 <- merge(twin2,randVar,by='Twinpair')

# now we check that the frequency of each number is the same in both files (i.e.
#   twin 1 and twin 2 should both have the same random number):

table(twin1$rand);table(twin2$rand)

# for twin 1, convert the random numbers to 0s and 1s:

twin1$rand[twin1$rand==4] <- 0
twin1$rand[twin1$rand==5] <- 1

# and then recode twin 2, but in the opposite direction to twin 1:

twin2$rand[twin2$rand==4] <- 1
twin2$rand[twin2$rand==5] <- 0

# check that it worked:

table(twin1$rand);table(twin2$rand) # twin 1 should have as many 0s as twin 2 has 1s

# now join twin 1 and twin 2 back together:

data <- data.frame(rbind(twin1,twin2))
names(data);dim(data)

# sort the data by pairnr and check the first few rows:

data <- data[order(data[,2]),]
head(data)

# split the data by random number:

twin1 <- subset(data,rand==1)
twin2 <- subset(data,rand==0)

# check that the dimensions are the same:

dim(twin1);dim(twin2)

# add '1' to the end of the variables for twin 1 and '2' for twin 2:

twin1labs <- paste(colnames(twin1),1,sep='')
twin2labs <- paste(colnames(twin2),2,sep='')

names(twin1) <- twin1labs
names(twin2) <- twin2labs

##########################################################################################################################

# combine data so that 1 pair each line
dataD <- data.frame(cbind(twin1,twin2))

# remove unused variables:

dataD <- dataD[-c(11:13,15,16,18:21,25,27,36:46,50)] # before was dataD[-c(4,6,8,11:13,15,16,18:23,25,27,29,32:42)]

# relabel a few variables:

names(dataD)[names(dataD)=='Twinpair1'] <- 'Twinpair'
names(dataD)[names(dataD)=='id1'] <- 'id'
names(dataD)[names(dataD)=='zygosity1'] <- 'zygosity'
names(dataD)[names(dataD)=='TWAB1'] <- 'TWAB'


#### done. Save two versions of the data: 

write.csv(data,'long_sensitivity_data.csv',row.names = F)
write.csv(dataD,'twin_bin_sensitivity_data.csv',row.names = F)

