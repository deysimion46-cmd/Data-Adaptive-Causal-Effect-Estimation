setwd("/Users/simionde/Downloads/results")
library(ggplot2)
library(tidyverse)
library(patchwork)
res_mat=read.table("resmat141",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv141",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
g1=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + theme_bw(base_size = 18)+labs(x=NULL)+
  theme(axis.text.x = element_blank())+ggtitle("n=500,p=10")
res_mat=read.table("resmat143",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv143",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
#2nd plot
g3=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + labs(x=NULL)+ theme_bw(base_size = 18)+
  theme(axis.text.x = element_blank())+labs(y=NULL)+ggtitle("n=500,p=100")

res_mat=read.table("resmat147",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv147",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
g7=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ggtitle("n=2000,p=10")

res_mat=read.table("resmat149",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv149",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
g9=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1) )+labs(y=NULL)+ggtitle("n=2000,p=100")

rplot14=g1+g3+g7+g9
ggsave("Model143.pdf",plot = rplot14,device = "pdf",height = 9,width = 11)


library(ggplot2)
library(tidyverse)
library(patchwork)
res_mat=read.table("resmat21",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv21",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
g1=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + theme_bw(base_size = 18)+labs(x=NULL)+
  theme(axis.text.x = element_blank())+ggtitle("n=500,p=10")
res_mat=read.table("resmat23",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv23",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
#2nd plot
g3=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + labs(x=NULL)+ theme_bw(base_size = 18)+
  theme(axis.text.x = element_blank())+labs(y=NULL)+ggtitle("n=500,p=100")

res_mat=read.table("resmat27",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv27",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
g7=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ggtitle("n=2000,p=10")

res_mat=read.table("resmat29",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv29",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
g9=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1) )+labs(y=NULL)+ggtitle("n=2000,p=100")

rplot2=g1+g3+g7+g9
ggsave("Model23.pdf",plot = rplot2,device = "pdf",height = 9,width = 11)




library(ggplot2)
library(tidyverse)
library(patchwork)
res_mat=read.table("resmat81",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv81",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
g1=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + theme_bw(base_size = 18)+labs(x=NULL)+
  theme(axis.text.x = element_blank())+ggtitle("n=500,p=10")
res_mat=read.table("resmat83",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv83",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
#2nd plot
g3=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + labs(x=NULL)+ theme_bw(base_size = 18)+
  theme(axis.text.x = element_blank())+labs(y=NULL)+ggtitle("n=500,p=100")

res_mat=read.table("resmat87",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv87",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
g7=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ggtitle("n=2000,p=10")

res_mat=read.table("resmat89",header=T)
res_mat=res_mat[,c(2,7,10,16)]
colnames(res_mat)[1]="Gaussian Kernel MMD"
colnames(res_mat)[2]="RF kernel MMD"
colnames(res_mat)[3]="Logistic IPW"
resprop=read.table("respropinv89",header=T)
colnames(resprop)=c("RF IPW")
res_mat=cbind.data.frame(res_mat,resprop)
g9=res_mat %>% as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "Method",
               values_to = "Error") %>%
  ggplot(aes(x = Method, y = Error)) +
  geom_boxplot() + 
  labs(x=NULL) + theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust=1) )+labs(y=NULL)+ggtitle("n=2000,p=100")

rplot8=g1+g3+g7+g9
ggsave("Model83.pdf",plot = rplot2,device = "pdf",height = 9,width = 11)









