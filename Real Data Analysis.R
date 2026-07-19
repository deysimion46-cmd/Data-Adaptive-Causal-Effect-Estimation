library(cobalt)
weights=read.table(file.choose(),header = F)
w=weights[,1]
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
g1=ggplot(data,aes(x=Balance,y=TreatmentPval))+geom_point()+theme_minimal()+scale_y_log10()
g2=ggplot(data,aes(x=Balance,y=OutcomePval))+geom_point()+theme_minimal()+scale_y_log10()
library(patchwork)
g1+g2

ggplot(data,aes(x=log10(-log10(TreatmentPval)),y=log10(-log10((OutcomePval))),size=Balance))+
  geom_point(alpha=.5)+theme_minimal()+labs(x="Treatment p-value", y="Outcome p-value")


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


