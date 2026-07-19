library(tidyverse)

file <- "https://hbiostat.org/data/repo/rhc.csv"

rhc <- read_csv(file, col_select = -1)



rhc_df <-
  rhc %>%
  dplyr::select(age = age,
                sex = sex,
                race = race,
                edu = edu,
                income = income,
                ninsclas = ninsclas,
                cat1 = cat1,
                ca = ca,
                wtkilo1 = wtkilo1,
                temp1 = temp1,
                meanbp1 = meanbp1,
                resp1 = resp1,
                hrt1 = hrt1,
                pafi1 = pafi1,
                paco21 = paco21,
                ph1 = ph1,
                cardiohx,
                chfhx,
                "dementhx", "psychhx", "chrpulhx",
                "renalhx", "liverhx", "gibledhx", "malighx", "immunhx", "transhx",
                "amihx",
                "surv2md1", "das2d3pc",
                "aps1", "scoma1", "wblc1",
                "alb1", "hema1", "bili1", "crea1", "sod1",
                "pot1", "dnr1",
                "resp", "card", "neuro", "gastr", "renal", "meta", "hema", "seps",
                "trauma", "ortho", "adld3p", "urin1",
                death30 = dth30,
                death180 = death,
                time = t3d30,
                RHC = swang1) ### this is the treatment



cnames <- c("cat1", "cat2", "ca", "sadmdte", "dschdte", "dthdte", "lstctdte",
            "death", "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx",
            "renalhx", "liverhx", "gibledhx", "malighx", "immunhx", "transhx",
            "amihx", "age", "sex", "edu", "surv2md1", "das2d3pc", "t3d30",
            "dth30", "aps1", "scoma1", "meanbp1", "wblc1", "hrt1", "resp1",
            "temp1", "pafi1", "alb1", "hema1", "bili1", "crea1", "sod1",
            "pot1", "paco21", "ph1", "swang1", "wtkilo1", "dnr1", "ninsclas",
            "resp", "card", "neuro", "gastr", "renal", "meta", "hema", "seps",
            "trauma", "ortho", "adld3p", "urin1", "race", "income", "ptid")


clabels <- c("Primary disease category", "Secondary disease category",
             "Cancer",
             "Study Admission Date", "Hospital Discharge Date", "Date of Death",
             "Date of Last Contact",
             "Death at any time up to 180 Days",
             "Acute MI, Peripheral Vascular Disease, Severe Cardiovascular Symptoms (NYHA-Class III), Very Severe Cardiovascular Symptoms (NYHA-Class IV)",
             "Congestive Heart Failure",
             "Dementia, Stroke or Cerebral Infarct, Parkinson's Disease",
             "Psychiatric History, Active Psychosis or Severe Depression",
             "Chronic Pulmonary Disease, Severe Pulmonary Disease, Very Severe Pulmonary Disease",
             "Chronic Renal Disease, Chronic Hemodialysis or Peritoneal Dialysis",
             "Cirrhosis, Hepatic Failure",
             "Upper GI Bleeding",
             "Solid Tumor, Metastatic Disease, Chronic Leukemia/Myeloma, Acute Leukemia, Lymphoma",
             "Immunosupperssion, Organ Transplant, HIV Positivity, Diabetes Mellitus Without End Organ Damage, Diabetes Mellitus With End Organ Damage, Connective Tissue Disease",
             "Transfer (> 24 Hours) from Another Hospital",
             "Definite Myocardial Infarction",
             "Age", "Sex", "Years of education", "Prob Surv > 2yrs",  #"Support model estimate of the prob. of surviving 2 months",
             "Duke Activity Status Index", "death indicator", "time to death within 30", "APACHE score", "Glasgow Coma Score",
             "Mean blood pressure", "White Blood Count)", "Heart rate", "Respiratory rate", "Temperature",
             "PaO2/FIO2 ratio", "Albumin", "Hematocrit", "Bilirubin", "Creatinine", "Sodium", "Potassium",
             "PaCo2", "PH", "Right Heart Catheterization (RHC)", "Weight (kg)", "DNR status on day1", "Medical insurance",
             "Respiratory Diagnosis", "Cardiovascular Diagnosis", "Neurological Diagnosis", "Gastrointestinal Diagnosis",
             "Renal Diagnosis", "Metabolic Diagnosis", "Hematologic Diagnosis", "Sepsis Diagnosis", "Trauma Diagnosis",
             "Orthopedic Diagnosis", "ADL", "Urine output", "Race", "Income", "Patient ID")




### this is the outcome -- death within 30 days

rhc_df$death <- as.numeric(as.factor(rhc_df$death30)) -1

### this is the treatment -- make it 1 for treated (RHC) and 0 for untreated (control)
rhc_df$RHC <- as.numeric(as.factor(rhc_df$RHC)) - 1



## these are the baseline covariates to control for/balance
vars_2_balance <- c("age",
                    "sex",
                    "race",
                    "edu",
                    "income",
                    "ninsclas",
                    "cat1",
                    "ca",
                    "wtkilo1",
                    "temp1",
                    "meanbp1",
                    "resp1",
                    "hrt1",
                    "pafi1",
                    "paco21",
                    "ph1",
                    "cardiohx",
                    "chfhx",
                    "dementhx", "psychhx", "chrpulhx",
                    "renalhx", "liverhx", "gibledhx", "malighx", "immunhx", "transhx",
                    "amihx",
                    "surv2md1", "das2d3pc",
                    "aps1", "scoma1", "wblc1",
                    "alb1", "hema1", "bili1", "crea1", "sod1",
                    "pot1", "dnr1",
                    "resp", "card", "neuro", "gastr", "renal", "meta", "hema", "seps",
                    "trauma", "ortho")

Xmain=rhc_df[,vars_2_balance]
amain=rhc_df$RHC
ymain=rhc_df$death
Xmain=lapply(Xmain, function(x) if(is.character(x)) factor(x) else x)
Xmain=data.frame(Xmain)
form=formula(paste("~",paste(colnames(Xmain),collapse="+")))
des=model.matrix(form,data = Xmain)
des=des[,-1]#the numerical design matrix
s=c(1:5735)
x=as.matrix(des[s,])
y=ymain[s]
a=amain[s]
library(cobalt)
fit_real=forest_balance(x,a,y,cross.fitting = T,num.folds = 2)
weights=fit_real$weights
w=weights

dfxa  <- data.frame(a = a, x = x)
dfxay <- data.frame(a = a, x = x, y = y)

## model for treatment
formula_a <- as.formula(paste0("a ~ ", paste(colnames(dfxa)[-1], collapse = " + ")))

m=love.plot(formula_a, data = dfxa, weights = list(ebw_grf_multi = w),
            int = F,var.order = "adjusted", poly = 1,abs=T,
            un = TRUE)
tab= bal.tab(formula_a, data = dfxa, weights = list(ebw_grf_multi = w),
             int = F,var.order = "adjusted", poly = 1,abs=T,
             un = TRUE)
table=tab$Balance
change=table$Diff.Un-table$Diff.Adj
ref_change=change/(table$Diff.Un)
ref_change=(ref_change+abs(ref_change))/2
balance1=table$Diff.Adj
k=order(change)[46:65]
feat=colnames(des)[k]
newdata=des[,feat]
dfxa  <- data.frame(a = a,  newdata)
formula_a <- as.formula(paste0("a ~ ", paste(colnames(dfxa)[-1], collapse = " + ")))
m=love.plot(formula_a, data = dfxa, weights = list(ebw_grf_multi = w),
            int = F,var.order = "adjusted", poly = 1,abs=T,
            un = TRUE)
pval=matrix(0,ncol=2,nrow=20)
for(i in 1:20){
  pval[i,1]=t.test(des[,k[i]]~a,alternate=two.sided)$p.value
  pval[i,2]=t.test(des[,k[i]]~y,alternate=two.sided)$p.value
  
}
pval=data.frame(pval)
pval$feat=feat


pvalentire=matrix(0,ncol=2,nrow=65)
for(i in 1:65){
  if(table$Type[i]=="Contin."){
    pvalentire[i,1]=t.test(des[,i]~a,alternate=two.sided)$p.value
    pvalentire[i,2]=t.test(des[,i]~y,alternate=two.sided)$p.value
    
  }else{
    table1=table(des[,i],a)
    table2=table(des[,i],y)
    pvalentire[i,1]=chisq.test(table1)$p.value
    pvalentire[i,2]=chisq.test(table2)$p.value
    
    
  }
  
}

plot(change,pvalentire[,1])
plot(change,pvalentire[,2])

min=order(change)[1:20]
feat_min=colnames(des)[min]
pval_min=matrix(0,ncol=2,nrow=20)
for(i in 1:20){
  pval_min[i,]=pvalentire[min[i],]
}
pval_min=data.frame(pval_min)
pval_min$feat=feat_min


data=data.frame(Balance=change,TreatmentPval=pvalentire[,1],OutcomePval=pvalentire[,2])
library(ggplot2)
g1=ggplot(data,aes(x=Balance,y=log10(-log10(TreatmentPval))))+geom_point()+theme_minimal()
g2=ggplot(data,aes(x=Balance,y=log10(-log10((OutcomePval)))))+geom_point()+theme_minimal()
library(patchwork)
g1+g2

ggplot(data,aes(x=log10(-log10(TreatmentPval)),y=log10(-log10((OutcomePval))),size=Balance))+
  geom_point(alpha=.5)+theme_minimal()+labs(x="log of negative log of Treatment p-value", y="log of negative log of Outcome p-value")


test_stat=matrix(0,ncol=2,nrow=65)
for(i in 1:65){
  if(table$Type[i]=="Contin."){
    test_stat[i,1]=abs(t.test(des[,i]~a,alternate=two.sided)$statistic)
    test_stat[i,2]=abs(t.test(des[,i]~y,alternate=two.sided)$statistic)
    
  }else{
    table1=table(des[,i],a)
    table2=table(des[,i],y)
    test_stat[i,1]=sqrt(chisq.test(table1)$statistic)
    test_stat[i,2]=sqrt(chisq.test(table2)$statistic)
    
    
  }
  
}
data_test=data.frame(Change_in_balance=change,Treatment_teststatistic=test_stat[,1],Outcome_teststatistic=test_stat[,2])
library(ggplot2)
g1=ggplot(data_test,aes(x=Change_in_balance,y=Treatment_teststatistic))+geom_point()+labs(
  x = "Change in Balance ",         # x-axis label
  y = "Treatment Test Statistic"         # y-axis label
  # legend label (based on the 'fill' aesthetic)
)+theme_minimal()
g2=ggplot(data_test,aes(x=Change_in_balance,y=Outcome_teststatistic))+geom_point()+labs(
  x = "Change in Balance ",         # x-axis label
  y = "Outcome Test Statistic"         # y-axis label
  # legend label (based on the 'fill' aesthetic)
)+theme_minimal()
library(patchwork)
g1+g2

ggplot(data_test,aes(x=Treatment_teststatistic,y=Outcome_teststatistic,size=Change_in_balance))+
  geom_point(alpha=.5)+theme_minimal()+labs(x="Treatment Test Statistic", y="Outcome Test Statistic", size=" Change in Balance")



