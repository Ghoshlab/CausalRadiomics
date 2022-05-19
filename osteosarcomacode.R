# R code for latent causal with radiomics on osteosarcoma dataset

osteo <- read.csv("osteosarcoma.csv",header=T)
y <- ifelse(osteo[,8]=="effective group",1,0)

#recode to stage II versus not
stagedic <- ifelse(osteo[,9] == unique(osteo[,9])[2],1,0)

# This requires the PLS package
library(pls)
pls1 <- mvr(radiomics~stagedic,scale=T,method="kernelpls")
ce <- lm(y~pls1$scores)

# bootstrap this
beta <- NULL
B <- 1000
for (i in 1:B) {
  set.seed(3*i)
 tmp <- sample(1:102,size=102,replace=T) 
 tmppls1 <- mvr(radiomics[tmp]~stagedic[tmp],scale=T,method="kernelpls")
 tmpce <- lm(y[tmp]~tmppls1$scores)
 beta <- c(beta,as.numeric(coef(tmpce)[2]))
 cat(i,"\n")
}


# Other analyses

pca1 <- prcomp(radiomics,scale.=T)

ce2 <- lm(y~pca1$x[,1]+stagedic)
