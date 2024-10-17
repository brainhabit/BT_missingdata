#==========
# GEE analysis
# used to test phenotypic associations
#
# author: Giorgia Bussu
# project: BT missing data
# version: June 2024
#==========

rm(list=ls())

data <- read.csv(file='C:/Users/myfolder/long_wcontinuous_data.csv',header=T,
                 sep=',')
names(data);dim(data)

quest<-read.csv(file='C:/Users/myfolder/imputed_data.csv',header=T,
                sep=',')
names(quest);dim(quest)


# match files
indx<-match(quest$id,data$id)

# matching data
my_data<-data[indx,]

data<-quest
data$missing<-my_data$emi_ts

####### check phenotypic associations

library(drgee)

p<-rep(0,17)

test<-summary(gee( formula = missing~sex,data = data, clusterid = "Twinpair", cond = F))
p[1]<-test$coefficients[2,4]


test<-summary(gee( formula = missing~age_s,data = data, clusterid = "Twinpair", cond = F))
p[2]<-test$coefficients[2,4]


test<-summary(gee( formula = missing~page_s,data = data, clusterid = "Twinpair", cond = F))
p[3]<-test$coefficients[2,4]


test<-summary(gee( formula = missing~tage_s,data = data, clusterid = "Twinpair", cond = F))
p[4]<-test$coefficients[2,4]


test<-summary(gee( formula = missing~edu_s,data = data, clusterid = "Twinpair", cond = F))
p[5]<-test$coefficients[2,4]


test<-summary(gee( formula = missing~income_s,data = data, clusterid = "Twinpair", cond = F))
p[6]<-test$coefficients[2,4]


test<-summary(gee( formula = missing~qchat,data = data, clusterid = "Twinpair", cond = F))
p[7]<-test$coefficients[2,4]

test<-summary(gee( formula = missing~vabs_sc,data = data, clusterid = "Twinpair", cond = F))
p[8]<-test$coefficients[2,4]

test<-summary(gee( formula = missing~vabs_m,data = data, clusterid = "Twinpair", cond = F))
p[9]<-test$coefficients[2,4]

test<-summary(gee( formula = missing~surgency,data = data, clusterid = "Twinpair", cond = F))
p[10]<-test$coefficients[2,4]

test<-summary(gee( formula = missing~negaff,data = data, clusterid = "Twinpair", cond = F))
p[11]<-test$coefficients[2,4]

test<-summary(gee( formula = missing~eff,data = data, clusterid = "Twinpair", cond = F))
p[12]<-test$coefficients[2,4]

test<-summary(gee( formula = missing~sens,data = data, clusterid = "Twinpair", cond = F))
p[13]<-test$coefficients[2,4]

test<-summary(gee( formula = missing~avoid,data = data, clusterid = "Twinpair", cond = F))
p[14]<-test$coefficients[2,4]

test<-summary(gee( formula = missing~seek,data = data, clusterid = "Twinpair", cond = F))
p[15]<-test$coefficients[2,4]

test<-summary(gee( formula = missing~lowreg,data = data, clusterid = "Twinpair", cond = F))
p[16]<-test$coefficients[2,4]

test<-summary(gee( formula = missing~pgs+V1+V2+V3+V4+V5+V6+V7+V8+V9+V10,data = data, clusterid = "Twinpair", cond = F))
p[17]<-test$coefficients[2,4]

p.adjust(p,method='holm')
