library(osqp)
library(MASS)
library(gbm)
library(ranger)
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
simrep=function(s){
  x=as.matrix(des[s,])
  y=ymain[s]
  a=amain[s]
  dfxa  <- data.frame(a = a, x = x)
  formula_a <- as.formula(paste0("a ~ ", paste(colnames(dfxa)[-1], collapse = " + ")))
  rf_fit_a <- ranger(formula_a, data = dfxa,
                     num.trees = 1000)
  prop=rf_fit_a$predictions
  prop[a==0]=1-prop[a==0]
  invprop=1/prop
  return(weighted.mean(y[a==1],w=invprop[a==1])-weighted.mean(y[a==0],w=invprop[a==0]))
  
}
set.seed(124)
resmatrdpropinv=numeric(40)
for(i in 1:40){
  s=sample(c(1:5735),5735,replace = T)
  resmatrdpropinv[i]=simrep(s)
  
}
write.table(resmatrdpropinv,"resmatrdpropinv2",row.names = F,quote = F,col.names = TRUE)