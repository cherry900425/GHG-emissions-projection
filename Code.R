#things to run before any session
if(!require(R2jags)) install.packages("R2jags")
library(msm)
library(readxl)
#import data
data_gdp <- read_excel("E:/Data1_update/data_original.xlsx", sheet="gdp")
data_gdp <- data.frame(data_gdp)
data_fh <- read_excel("E:/Data1_update/data_original.xlsx", sheet="fh")
data_fh <- data.frame(data_fh)
data_fv <- read_excel("E:/Data1_update/data_original.xlsx", sheet="fv")
data_fv <- data.frame(data_fv)
data_fg <- read_excel("E:/Data1_update/data_original.xlsx", sheet="fg")
data_fg <- data.frame(data_fg)
data_s <- read_excel("E:/Data1_update/data_original.xlsx", sheet="s")
data_s <- data.frame(data_s)
data_v <- read_excel("E:/Data1_update/data_original.xlsx", sheet="v")
data_v <- data.frame(data_v)
data_g <- read_excel("E:/Data1_update/data_original.xlsx", sheet="g")
data_g <- data.frame(data_g)
data_y <- read_excel("E:/Data1_update/data_original.xlsx", sheet="y")
data_y <- data.frame(data_y)
data_p <- read_excel("E:/Data1_update/data_original.xlsx", sheet="pop")
data_p <- data.frame(data_p)
data_r <- read_excel("E:/Data1_update/data_R.xlsx", sheet="ratio")
data_r <- data.frame(data_r)

n_2050=56
n_2018=24
n_2015=21
n_region=40

gdp=(data_gdp [1:24,2:41])
x=log(data_gdp [1:24,2:41])
x_2023=log(data_gdp [1:29,2:41])
fh=(data_fh[,3:42])
fv=(data_fv[,3:42])
fg=(data_fg[,3:42])
so=(data_s[,3:42])
vo=(data_v[,3:42])
go=(data_g[,3:42])
yo=(data_y[,3:42])
ro=(data_r[,2:41])
#adjust factors with zero value to non-zero for LMDI decomposition analysis
for (j in 1:40){
  for (i in 1:nrow(fh)){
    if (fh[i,j] > 0){fh[i,j]=fh[i,j]}
    else {fh[i,j]<-10^(-50)}
    if (so[i,j] > 0){so[i,j]=so[i,j]}
    else {so[i,j]<-10^(-50)}
  }
}
for (j in 1:40){
  for (i in 1:nrow(fv)){  
    if (fv[i,j] > 0){fv[i,j]=fv[i,j]}
    else {fv[i,j]<-10^(-50)}
    if (vo[i,j] > 0){vo[i,j]=vo[i,j]}
    else {vo[i,j]<-10^(-50)}
  }
}
for (j in 1:40){
  for (i in 1:nrow(fg)){  
    if (fg[i,j] > 0){fg[i,j]=fg[i,j]}
    else {fg[i,j]<-10^(-50)}
    if (go[i,j] > 0){go[i,j]=go[i,j]}
    else {go[i,j]<-10^(-50)}
  }
}
for (j in 1:40){
  for (i in 1:nrow(yo)){
    if (yo[i,j] >- 0){yo[i,j]=yo[i,j]}
    else {yo[i,j]<-10^(-50)}
  }
}

eih=log(fh)
eiv=log(fv)
eig=log(fg)
s=log(so)
v=log(vo)
g=log(go)
y=log(yo)
p=data_p[,2:41]

n.iter = 100000
n.burnin = 3000
n.thin = 10
n.chains = 7

#parameters for gdp
x_pre=log(data_gdp[1:56,2:41])
d<- matrix(NA, nrow =n_2050-1,ncol=n_region)
for  (t in 1:n_2050-1){
  for  (j in 1:n_region){
    d[t,j]<-x_pre[t+1,j]-x_pre[t,j]
  }
}

library(R2jags)

jagsdata_s1 <- with(data_gdp, list(gdp=gdp,x=x,n_2018=24,n_region=40))

lm1_jags <- function(){
  # Likelihood:
  for (j in 1:n_region){
    for (t in 2:n_2018){
      gdp[t,j] ~ dnorm(mu_gdp[t,j], tau_gdp[j]) # tau is precision (1 / variance)
      mu_gdp[t,j] <- exp(alpha_gdp[j]+x[t-1,j])
    }
  }

  # Priors:
  for (j in 1:n_region){ 
    alpha_gdp[j] ~ dnorm(mu_alpha_gdp[j], tau_alpha_gdp[j]) # intercept
    mu_alpha_gdp[j]~ dunif(0,1) 
    sigma_alpha_gdp[j] ~ dunif(0, 0.1) # standard deviation
    tau_alpha_gdp[j] <- 1 / (sigma_alpha_gdp[j] * sigma_alpha_gdp[j]) # sigma^2 doesn't work in JAGS
    tau_gdp[j] <- 1 / (sigma_gdp[j] * sigma_gdp[j]) # sigma^2 doesn't work in JAGS
    sigma_gdp[j] ~ dunif(0,5) # standard deviation
  }
}

params <- c("alpha_gdp[1]","sigma_gdp[1]","alpha_gdp[2]","sigma_gdp[2]","alpha_gdp[3]","sigma_gdp[3]",
            "alpha_gdp[4]","sigma_gdp[4]","alpha_gdp[5]","sigma_gdp[5]","alpha_gdp[6]","sigma_gdp[6]", 
            "alpha_gdp[7]","sigma_gdp[7]","alpha_gdp[8]","sigma_gdp[8]","alpha_gdp[9]","sigma_gdp[9]",
            "alpha_gdp[10]","sigma_gdp[10]","alpha_gdp[11]", "sigma_gdp[11]","alpha_gdp[12]","sigma_gdp[12]",
            "alpha_gdp[13]","sigma_gdp[13]","alpha_gdp[14]","sigma_gdp[14]","alpha_gdp[15]","sigma_gdp[15]",
            "alpha_gdp[16]","sigma_gdp[16]","alpha_gdp[17]","sigma_gdp[17]","alpha_gdp[18]","sigma_gdp[18]",
            "alpha_gdp[19]","sigma_gdp[19]","alpha_gdp[20]","sigma_gdp[20]","alpha_gdp[21]","sigma_gdp[21]", 
            "alpha_gdp[22]","sigma_gdp[22]","alpha_gdp[23]","sigma_gdp[23]","alpha_gdp[24]","sigma_gdp[24]",
            "alpha_gdp[25]","sigma_gdp[25]","alpha_gdp[26]", "sigma_gdp[26]","alpha_gdp[27]","sigma_gdp[27]",
            "alpha_gdp[28]","sigma_gdp[28]","alpha_gdp[29]","sigma_gdp[29]","alpha_gdp[30]","sigma_gdp[30]",
            "alpha_gdp[31]","sigma_gdp[31]","alpha_gdp[32]","sigma_gdp[32]","alpha_gdp[33]","sigma_gdp[33]",
            "alpha_gdp[34]","sigma_gdp[34]","alpha_gdp[35]","sigma_gdp[35]","alpha_gdp[36]","sigma_gdp[36]", 
            "alpha_gdp[37]","sigma_gdp[37]","alpha_gdp[38]","sigma_gdp[38]","alpha_gdp[39]","sigma_gdp[39]",
            "alpha_gdp[40]","sigma_gdp[40]","sigma_alpha_gdp[1]","sigma_alpha_gdp[2]","sigma_alpha_gdp[3]",
            "sigma_alpha_gdp[4]","sigma_alpha_gdp[5]","sigma_alpha_gdp[6]","sigma_alpha_gdp[7]","sigma_alpha_gdp[8]",
            "sigma_alpha_gdp[9]","sigma_alpha_gdp[10]","sigma_alpha_gdp[11]","sigma_alpha_gdp[12]","sigma_alpha_gdp[13]",
            "sigma_alpha_gdp[14]","sigma_alpha_gdp[15]","sigma_alpha_gdp[16]","sigma_alpha_gdp[17]","sigma_alpha_gdp[18]",
            "sigma_alpha_gdp[19]","sigma_alpha_gdp[20]","sigma_alpha_gdp[21]","sigma_alpha_gdp[22]","sigma_alpha_gdp[23]",
            "sigma_alpha_gdp[24]","sigma_alpha_gdp[25]","sigma_alpha_gdp[26]","sigma_alpha_gdp[27]","sigma_alpha_gdp[28]",
            "sigma_alpha_gdp[29]","sigma_alpha_gdp[30]","sigma_alpha_gdp[31]","sigma_alpha_gdp[32]","sigma_alpha_gdp[33]",
            "sigma_alpha_gdp[34]","sigma_alpha_gdp[35]","sigma_alpha_gdp[36]","sigma_alpha_gdp[37]","sigma_alpha_gdp[38]","sigma_alpha_gdp[39]","sigma_alpha_gdp[40]")

fit_lm1 <- jags(data = jagsdata_s1,  parameters.to.save = params, model.file = lm1_jags,
                n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC = F)

fit_lm1

lm1_mcmc <- as.mcmc(fit_lm1)
lm1_mcmc_combi <- as.mcmc(rbind(lm1_mcmc[[1]], lm1_mcmc[[2]], lm1_mcmc[[3]],lm1_mcmc[[4]], lm1_mcmc[[5]], lm1_mcmc[[6]],lm1_mcmc[[7]]))
lm1_mcmc_gdp <- do.call(cbind, rep(list("alpha_gdp[1]"=lm1_mcmc_combi[,"alpha_gdp[1]"],"alpha_gdp[2]"=lm1_mcmc_combi[,"alpha_gdp[2]"],
                                        "alpha_gdp[3]"=lm1_mcmc_combi[,"alpha_gdp[3]"],"alpha_gdp[4]"=lm1_mcmc_combi[,"alpha_gdp[4]"],
                                        "alpha_gdp[5]"=lm1_mcmc_combi[,"alpha_gdp[5]"],"alpha_gdp[6]"=lm1_mcmc_combi[,"alpha_gdp[6]"],
                                        "alpha_gdp[7]"=lm1_mcmc_combi[,"alpha_gdp[7]"],"alpha_gdp[8]"=lm1_mcmc_combi[,"alpha_gdp[8]"],
                                        "alpha_gdp[9]"=lm1_mcmc_combi[,"alpha_gdp[9]"],"alpha_gdp[10]"=lm1_mcmc_combi[,"alpha_gdp[10]"],
                                        "alpha_gdp[11]"=lm1_mcmc_combi[,"alpha_gdp[11]"],"alpha_gdp[12]"=lm1_mcmc_combi[,"alpha_gdp[12]"],
                                        "alpha_gdp[13]"=lm1_mcmc_combi[,"alpha_gdp[13]"],"alpha_gdp[14]"=lm1_mcmc_combi[,"alpha_gdp[14]"],
                                        "alpha_gdp[15]"=lm1_mcmc_combi[,"alpha_gdp[15]"],"alpha_gdp[16]"=lm1_mcmc_combi[,"alpha_gdp[16]"],
                                        "alpha_gdp[17]"=lm1_mcmc_combi[,"alpha_gdp[17]"],"alpha_gdp[18]"=lm1_mcmc_combi[,"alpha_gdp[18]"],
                                        "alpha_gdp[19]"=lm1_mcmc_combi[,"alpha_gdp[19]"],"alpha_gdp[20]"=lm1_mcmc_combi[,"alpha_gdp[20]"],
                                        "alpha_gdp[21]"=lm1_mcmc_combi[,"alpha_gdp[21]"],"alpha_gdp[22]"=lm1_mcmc_combi[,"alpha_gdp[22]"],
                                        "alpha_gdp[23]"=lm1_mcmc_combi[,"alpha_gdp[23]"],"alpha_gdp[24]"=lm1_mcmc_combi[,"alpha_gdp[24]"],
                                        "alpha_gdp[25]"=lm1_mcmc_combi[,"alpha_gdp[25]"],"alpha_gdp[26]"=lm1_mcmc_combi[,"alpha_gdp[26]"],
                                        "alpha_gdp[27]"=lm1_mcmc_combi[,"alpha_gdp[27]"],"alpha_gdp[28]"=lm1_mcmc_combi[,"alpha_gdp[28]"],
                                        "alpha_gdp[29]"=lm1_mcmc_combi[,"alpha_gdp[29]"],"alpha_gdp[30]"=lm1_mcmc_combi[,"alpha_gdp[30]"],
                                        "alpha_gdp[31]"=lm1_mcmc_combi[,"alpha_gdp[31]"],"alpha_gdp[42]"=lm1_mcmc_combi[,"alpha_gdp[32]"],
                                        "alpha_gdp[33]"=lm1_mcmc_combi[,"alpha_gdp[33]"],"alpha_gdp[44]"=lm1_mcmc_combi[,"alpha_gdp[34]"],
                                        "alpha_gdp[35]"=lm1_mcmc_combi[,"alpha_gdp[35]"],"alpha_gdp[46]"=lm1_mcmc_combi[,"alpha_gdp[36]"],
                                        "alpha_gdp[37]"=lm1_mcmc_combi[,"alpha_gdp[37]"],"alpha_gdp[48]"=lm1_mcmc_combi[,"alpha_gdp[38]"],
                                        "alpha_gdp[39]"=lm1_mcmc_combi[,"alpha_gdp[39]"],"alpha_gdp[40]"=lm1_mcmc_combi[,"alpha_gdp[40]"],
                                        "sigma_gdp[1]"=lm1_mcmc_combi[,"sigma_gdp[1]"],"sigma_gdp[2]"=lm1_mcmc_combi[,"sigma_gdp[2]"],
                                        "sigma_gdp[3]"=lm1_mcmc_combi[,"sigma_gdp[3]"],"sigma_gdp[4]"=lm1_mcmc_combi[,"sigma_gdp[4]"],
                                        "sigma_gdp[5]"=lm1_mcmc_combi[,"sigma_gdp[5]"],"sigma_gdp[6]"=lm1_mcmc_combi[,"sigma_gdp[6]"],
                                        "sigma_gdp[7]"=lm1_mcmc_combi[,"sigma_gdp[7]"],"sigma_gdp[8]"=lm1_mcmc_combi[,"sigma_gdp[8]"],
                                        "sigma_gdp[9]"=lm1_mcmc_combi[,"sigma_gdp[9]"],"sigma_gdp[10]"=lm1_mcmc_combi[,"sigma_gdp[10]"],
                                        "sigma_gdp[11]"=lm1_mcmc_combi[,"sigma_gdp[11]"],"sigma_gdp[12]"=lm1_mcmc_combi[,"sigma_gdp[12]"],
                                        "sigma_gdp[13]"=lm1_mcmc_combi[,"sigma_gdp[13]"],"sigma_gdp[14]"=lm1_mcmc_combi[,"sigma_gdp[14]"],
                                        "sigma_gdp[15]"=lm1_mcmc_combi[,"sigma_gdp[15]"],"sigma_gdp[16]"=lm1_mcmc_combi[,"sigma_gdp[16]"],
                                        "sigma_gdp[17]"=lm1_mcmc_combi[,"sigma_gdp[17]"],"sigma_gdp[18]"=lm1_mcmc_combi[,"sigma_gdp[18]"],
                                        "sigma_gdp[19]"=lm1_mcmc_combi[,"sigma_gdp[19]"],"sigma_gdp[20]"=lm1_mcmc_combi[,"sigma_gdp[20]"],
                                        "sigma_gdp[21]"=lm1_mcmc_combi[,"sigma_gdp[21]"],"sigma_gdp[22]"=lm1_mcmc_combi[,"sigma_gdp[22]"],
                                        "sigma_gdp[23]"=lm1_mcmc_combi[,"sigma_gdp[23]"],"sigma_gdp[24]"=lm1_mcmc_combi[,"sigma_gdp[24]"],
                                        "sigma_gdp[25]"=lm1_mcmc_combi[,"sigma_gdp[25]"],"sigma_gdp[26]"=lm1_mcmc_combi[,"sigma_gdp[26]"],
                                        "sigma_gdp[27]"=lm1_mcmc_combi[,"sigma_gdp[27]"],"sigma_gdp[28]"=lm1_mcmc_combi[,"sigma_gdp[28]"],
                                        "sigma_gdp[29]"=lm1_mcmc_combi[,"sigma_gdp[29]"],"sigma_gdp[30]"=lm1_mcmc_combi[,"sigma_gdp[30]"],
                                        "sigma_gdp[31]"=lm1_mcmc_combi[,"sigma_gdp[31]"],"sigma_gdp[32]"=lm1_mcmc_combi[,"sigma_gdp[32]"],
                                        "sigma_gdp[33]"=lm1_mcmc_combi[,"sigma_gdp[33]"],"sigma_gdp[34]"=lm1_mcmc_combi[,"sigma_gdp[34]"],
                                        "sigma_gdp[35]"=lm1_mcmc_combi[,"sigma_gdp[35]"],"sigma_gdp[36]"=lm1_mcmc_combi[,"sigma_gdp[26]"],
                                        "sigma_gdp[37]"=lm1_mcmc_combi[,"sigma_gdp[37]"],"sigma_gdp[38]"=lm1_mcmc_combi[,"sigma_gdp[38]"],
                                        "sigma_gdp[39]"=lm1_mcmc_combi[,"sigma_gdp[39]"],"sigma_gdp[40]"=lm1_mcmc_combi[,"sigma_gdp[40]"],
                                        "sigma_alpha_gdp[1]"=lm1_mcmc_combi[,"sigma_alpha_gdp[1]"],"sigma_alpha_gdp[2]"=lm1_mcmc_combi[,"sigma_alpha_gdp[2]"],
                                        "sigma_alpha_gdp[3]"=lm1_mcmc_combi[,"sigma_alpha_gdp[3]"],"sigma_alpha_gdp[4]"=lm1_mcmc_combi[,"sigma_alpha_gdp[4]"],
                                        "sigma_alpha_gdp[5]"=lm1_mcmc_combi[,"sigma_alpha_gdp[5]"],"sigma_alpha_gdp[6]"=lm1_mcmc_combi[,"sigma_alpha_gdp[6]"],
                                        "sigma_alpha_gdp[7]"=lm1_mcmc_combi[,"sigma_alpha_gdp[7]"],"sigma_alpha_gdp[8]"=lm1_mcmc_combi[,"sigma_alpha_gdp[8]"],
                                        "sigma_alpha_gdp[9]"=lm1_mcmc_combi[,"sigma_alpha_gdp[9]"],"sigma_alpha_gdp[10]"=lm1_mcmc_combi[,"sigma_alpha_gdp[10]"],
                                        "sigma_alpha_gdp[11]"=lm1_mcmc_combi[,"sigma_alpha_gdp[11]"],"sigma_alpha_gdp[12]"=lm1_mcmc_combi[,"sigma_alpha_gdp[12]"],
                                        "sigma_alpha_gdp[13]"=lm1_mcmc_combi[,"sigma_alpha_gdp[13]"],"sigma_alpha_gdp[14]"=lm1_mcmc_combi[,"sigma_alpha_gdp[14]"],
                                        "sigma_alpha_gdp[15]"=lm1_mcmc_combi[,"sigma_alpha_gdp[15]"],"sigma_alpha_gdp[16]"=lm1_mcmc_combi[,"sigma_alpha_gdp[16]"],
                                        "sigma_alpha_gdp[17]"=lm1_mcmc_combi[,"sigma_alpha_gdp[17]"],"sigma_alpha_gdp[18]"=lm1_mcmc_combi[,"sigma_alpha_gdp[18]"],
                                        "sigma_alpha_gdp[19]"=lm1_mcmc_combi[,"sigma_alpha_gdp[19]"],"sigma_alpha_gdp[20]"=lm1_mcmc_combi[,"sigma_alpha_gdp[20]"],
                                        "sigma_alpha_gdp[21]"=lm1_mcmc_combi[,"sigma_alpha_gdp[21]"],"sigma_alpha_gdp[22]"=lm1_mcmc_combi[,"sigma_alpha_gdp[22]"],
                                        "sigma_alpha_gdp[23]"=lm1_mcmc_combi[,"sigma_alpha_gdp[23]"],"sigma_alpha_gdp[24]"=lm1_mcmc_combi[,"sigma_alpha_gdp[24]"],
                                        "sigma_alpha_gdp[25]"=lm1_mcmc_combi[,"sigma_alpha_gdp[25]"],"sigma_alpha_gdp[26]"=lm1_mcmc_combi[,"sigma_alpha_gdp[26]"],
                                        "sigma_alpha_gdp[27]"=lm1_mcmc_combi[,"sigma_alpha_gdp[27]"],"sigma_alpha_gdp[28]"=lm1_mcmc_combi[,"sigma_alpha_gdp[28]"],
                                        "sigma_alpha_gdp[29]"=lm1_mcmc_combi[,"sigma_alpha_gdp[29]"],"sigma_alpha_gdp[30]"=lm1_mcmc_combi[,"sigma_alpha_gdp[30]"],
                                        "sigma_alpha_gdp[31]"=lm1_mcmc_combi[,"sigma_alpha_gdp[31]"],"sigma_alpha_gdp[32]"=lm1_mcmc_combi[,"sigma_alpha_gdp[32]"],
                                        "sigma_alpha_gdp[33]"=lm1_mcmc_combi[,"sigma_alpha_gdp[33]"],"sigma_alpha_gdp[34]"=lm1_mcmc_combi[,"sigma_alpha_gdp[34]"],
                                        "sigma_alpha_gdp[35]"=lm1_mcmc_combi[,"sigma_alpha_gdp[35]"],"sigma_alpha_gdp[36]"=lm1_mcmc_combi[,"sigma_alpha_gdp[26]"],
                                        "sigma_alpha_gdp[37]"=lm1_mcmc_combi[,"sigma_alpha_gdp[37]"],"sigma_alpha_gdp[38]"=lm1_mcmc_combi[,"sigma_alpha_gdp[38]"],
                                        "sigma_alpha_gdp[39]"=lm1_mcmc_combi[,"sigma_alpha_gdp[39]"],"sigma_alpha_gdp[40]"=lm1_mcmc_combi[,"sigma_alpha_gdp[40]"]), 1)) 
                                    
#parametesrs for population 
x1_pre=log(data_p[1:56,2:41])
d1<- matrix(NA, nrow =n_2050-1,ncol=n_region)
for  (t in 1:n_2050-1){
  for  (j in 1:n_region){
    d1[t,j]<-x1_pre[t+1,j]-x1_pre[t,j]
  }
}
x1=log(data_p [1:26,2:41])


library(R2jags)

jagsdata_s1 <- with(data_p, list(p=p, x1=x1,n_region=40))

lm1_jags <- function(){
  # Likelihood:
  for (j in 1:n_region){
    for (t in 2:26){
      p[t,j] ~ dnorm(mu_p[t,j], tau_p[j]) # tau is precision (1 / variance)
      mu_p[t,j] <- exp(alpha_p[j]+x1[t-1,j])
    }
  }
  
  # Priors:
  for (j in 1:n_region){ 
    alpha_p[j] ~ dnorm(mu_alpha_p[j], tau_alpha_p[j]) # intercept
    mu_alpha_p[j]~ dunif(0,1) 
    sigma_alpha_p[j] ~ dunif(0, 0.03) # standard deviation
    tau_alpha_p[j] <- 1 / (sigma_alpha_p[j] * sigma_alpha_p[j]) # sigma^2 doesn't work in JAGS
    tau_p[j] <- 1 / (sigma_p[j] * sigma_p[j]) # sigma^2 doesn't work in JAGS
    sigma_p[j] ~ dunif(0,500) # standard deviation
  }
}

params <- c("alpha_p[1]","sigma_p[1]","alpha_p[2]","sigma_p[2]","alpha_p[3]","sigma_p[3]",
            "alpha_p[4]","sigma_p[4]","alpha_p[5]","sigma_p[5]","alpha_p[6]","sigma_p[6]", 
            "alpha_p[7]","sigma_p[7]","alpha_p[8]","sigma_p[8]","alpha_p[9]","sigma_p[9]",
            "alpha_p[10]","sigma_p[10]","alpha_p[11]", "sigma_p[11]","alpha_p[12]","sigma_p[12]",
            "alpha_p[13]","sigma_p[13]","alpha_p[14]","sigma_p[14]","alpha_p[15]","sigma_p[15]",
            "alpha_p[16]","sigma_p[16]","alpha_p[17]","sigma_p[17]","alpha_p[18]","sigma_p[18]",
            "alpha_p[19]","sigma_p[19]","alpha_p[20]","sigma_p[20]","alpha_p[21]","sigma_p[21]", 
            "alpha_p[22]","sigma_p[22]","alpha_p[23]","sigma_p[23]","alpha_p[24]","sigma_p[24]",
            "alpha_p[25]","sigma_p[25]","alpha_p[26]","sigma_p[26]","alpha_p[27]","sigma_p[27]",
            "alpha_p[28]","sigma_p[28]","alpha_p[29]","sigma_p[29]","alpha_p[30]","sigma_p[30]",
            "alpha_p[31]","sigma_p[31]","alpha_p[32]","sigma_p[32]","alpha_p[33]","sigma_p[33]",
            "alpha_p[34]","sigma_p[34]","alpha_p[35]","sigma_p[35]","alpha_p[36]","sigma_p[36]", 
            "alpha_p[37]","sigma_p[37]","alpha_p[38]","sigma_p[38]","alpha_p[39]","sigma_p[39]",
            "alpha_p[40]","sigma_p[40]","sigma_alpha_p[1]","sigma_alpha_p[2]","sigma_alpha_p[3]",
            "sigma_alpha_p[4]","sigma_alpha_p[5]","sigma_alpha_p[6]","sigma_alpha_p[7]","sigma_alpha_p[8]",
            "sigma_alpha_p[9]","sigma_alpha_p[10]","sigma_alpha_p[11]","sigma_alpha_p[12]","sigma_alpha_p[13]",
            "sigma_alpha_p[14]","sigma_alpha_p[15]","sigma_alpha_p[16]","sigma_alpha_p[17]","sigma_alpha_p[18]",
            "sigma_alpha_p[19]","sigma_alpha_p[20]","sigma_alpha_p[21]","sigma_alpha_p[22]","sigma_alpha_p[23]",
            "sigma_alpha_p[24]","sigma_alpha_p[25]","sigma_alpha_p[26]","sigma_alpha_p[27]","sigma_alpha_p[28]",
            "sigma_alpha_p[29]","sigma_alpha_p[30]","sigma_alpha_p[31]","sigma_alpha_p[32]","sigma_alpha_p[33]",
            "sigma_alpha_p[34]","sigma_alpha_p[35]","sigma_alpha_p[36]","sigma_alpha_p[37]","sigma_alpha_p[38]","sigma_alpha_p[39]","sigma_alpha_p[40]")

fit_lm1 <- jags(data = jagsdata_s1,  parameters.to.save = params, model.file = lm1_jags,
                n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC = F)

fit_lm1

lm1_mcmc <- as.mcmc(fit_lm1)
lm1_mcmc_combi <- as.mcmc(rbind(lm1_mcmc[[1]], lm1_mcmc[[2]], lm1_mcmc[[3]],lm1_mcmc[[4]], lm1_mcmc[[5]], lm1_mcmc[[6]],lm1_mcmc[[7]]))
lm1_mcmc_p <- do.call(cbind, rep(list("alpha_p[1]"=lm1_mcmc_combi[,"alpha_p[1]"],"alpha_p[2]"=lm1_mcmc_combi[,"alpha_p[2]"],
                                      "alpha_p[3]"=lm1_mcmc_combi[,"alpha_p[3]"],"alpha_p[4]"=lm1_mcmc_combi[,"alpha_p[4]"],
                                      "alpha_p[5]"=lm1_mcmc_combi[,"alpha_p[5]"],"alpha_p[6]"=lm1_mcmc_combi[,"alpha_p[6]"],
                                      "alpha_p[7]"=lm1_mcmc_combi[,"alpha_p[7]"],"alpha_p[8]"=lm1_mcmc_combi[,"alpha_p[8]"],
                                      "alpha_p[9]"=lm1_mcmc_combi[,"alpha_p[9]"],"alpha_p[10]"=lm1_mcmc_combi[,"alpha_p[10]"],
                                      "alpha_p[11]"=lm1_mcmc_combi[,"alpha_p[11]"],"alpha_p[12]"=lm1_mcmc_combi[,"alpha_p[12]"],
                                      "alpha_p[13]"=lm1_mcmc_combi[,"alpha_p[13]"],"alpha_p[14]"=lm1_mcmc_combi[,"alpha_p[14]"],
                                      "alpha_p[15]"=lm1_mcmc_combi[,"alpha_p[15]"],"alpha_p[16]"=lm1_mcmc_combi[,"alpha_p[16]"],
                                      "alpha_p[17]"=lm1_mcmc_combi[,"alpha_p[17]"],"alpha_p[18]"=lm1_mcmc_combi[,"alpha_p[18]"],
                                      "alpha_p[19]"=lm1_mcmc_combi[,"alpha_p[19]"],"alpha_p[20]"=lm1_mcmc_combi[,"alpha_p[20]"],
                                      "alpha_p[21]"=lm1_mcmc_combi[,"alpha_p[21]"],"alpha_p[22]"=lm1_mcmc_combi[,"alpha_p[22]"],
                                      "alpha_p[23]"=lm1_mcmc_combi[,"alpha_p[23]"],"alpha_p[24]"=lm1_mcmc_combi[,"alpha_p[24]"],
                                      "alpha_p[25]"=lm1_mcmc_combi[,"alpha_p[25]"],"alpha_p[26]"=lm1_mcmc_combi[,"alpha_p[26]"],
                                      "alpha_p[27]"=lm1_mcmc_combi[,"alpha_p[27]"],"alpha_p[28]"=lm1_mcmc_combi[,"alpha_p[28]"],
                                      "alpha_p[29]"=lm1_mcmc_combi[,"alpha_p[29]"],"alpha_p[30]"=lm1_mcmc_combi[,"alpha_p[30]"],
                                      "alpha_p[31]"=lm1_mcmc_combi[,"alpha_p[31]"],"alpha_p[32]"=lm1_mcmc_combi[,"alpha_p[32]"],
                                      "alpha_p[33]"=lm1_mcmc_combi[,"alpha_p[33]"],"alpha_p[34]"=lm1_mcmc_combi[,"alpha_p[34]"],
                                      "alpha_p[35]"=lm1_mcmc_combi[,"alpha_p[35]"],"alpha_p[36]"=lm1_mcmc_combi[,"alpha_p[36]"],
                                      "alpha_p[37]"=lm1_mcmc_combi[,"alpha_p[37]"],"alpha_p[38]"=lm1_mcmc_combi[,"alpha_p[38]"],
                                      "alpha_p[39]"=lm1_mcmc_combi[,"alpha_p[39]"],"alpha_p[40]"=lm1_mcmc_combi[,"alpha_p[40]"],
                                      "sigma_p[1]"=lm1_mcmc_combi[,"sigma_p[1]"],"sigma_p[2]"=lm1_mcmc_combi[,"sigma_p[2]"],
                                      "sigma_p[3]"=lm1_mcmc_combi[,"sigma_p[3]"],"sigma_p[4]"=lm1_mcmc_combi[,"sigma_p[4]"],
                                      "sigma_p[5]"=lm1_mcmc_combi[,"sigma_p[5]"],"sigma_p[6]"=lm1_mcmc_combi[,"sigma_p[6]"],
                                      "sigma_p[7]"=lm1_mcmc_combi[,"sigma_p[7]"],"sigma_p[8]"=lm1_mcmc_combi[,"sigma_p[8]"],
                                      "sigma_p[9]"=lm1_mcmc_combi[,"sigma_p[9]"],"sigma_p[10]"=lm1_mcmc_combi[,"sigma_p[10]"],
                                      "sigma_p[11]"=lm1_mcmc_combi[,"sigma_p[11]"],"sigma_p[12]"=lm1_mcmc_combi[,"sigma_p[12]"],
                                      "sigma_p[13]"=lm1_mcmc_combi[,"sigma_p[13]"],"sigma_p[14]"=lm1_mcmc_combi[,"sigma_p[14]"],
                                      "sigma_p[15]"=lm1_mcmc_combi[,"sigma_p[15]"],"sigma_p[16]"=lm1_mcmc_combi[,"sigma_p[16]"],
                                      "alpha_p[17]"=lm1_mcmc_combi[,"sigma_p[17]"],"sigma_p[18]"=lm1_mcmc_combi[,"sigma_p[18]"],
                                      "alpha_p[19]"=lm1_mcmc_combi[,"sigma_p[19]"],"sigma_p[20]"=lm1_mcmc_combi[,"sigma_p[20]"],
                                      "sigma_p[21]"=lm1_mcmc_combi[,"sigma_p[21]"],"sigma_p[22]"=lm1_mcmc_combi[,"sigma_p[22]"],
                                      "sigma_p[23]"=lm1_mcmc_combi[,"sigma_p[23]"],"sigma_p[24]"=lm1_mcmc_combi[,"sigma_p[24]"],
                                      "sigma_p[25]"=lm1_mcmc_combi[,"sigma_p[25]"],"sigma_p[26]"=lm1_mcmc_combi[,"sigma_p[26]"],
                                      "sigma_p[27]"=lm1_mcmc_combi[,"sigma_p[27]"],"sigma_p[28]"=lm1_mcmc_combi[,"sigma_p[28]"],
                                      "sigma_p[29]"=lm1_mcmc_combi[,"sigma_p[29]"],"sigma_p[30]"=lm1_mcmc_combi[,"sigma_p[30]"],
                                      "sigma_p[31]"=lm1_mcmc_combi[,"sigma_p[31]"],"sigma_p[32]"=lm1_mcmc_combi[,"sigma_p[32]"],
                                      "sigma_p[33]"=lm1_mcmc_combi[,"sigma_p[33]"],"sigma_p[34]"=lm1_mcmc_combi[,"sigma_p[34]"],
                                      "sigma_p[35]"=lm1_mcmc_combi[,"sigma_p[35]"],"sigma_p[36]"=lm1_mcmc_combi[,"sigma_p[36]"],
                                      "sigma_p[37]"=lm1_mcmc_combi[,"sigma_p[37]"],"sigma_p[38]"=lm1_mcmc_combi[,"sigma_p[38]"],
                                      "sigma_p[39]"=lm1_mcmc_combi[,"sigma_p[39]"],"sigma_p[40]"=lm1_mcmc_combi[,"sigma_p[40]"],
                                      "sigma_alpha_p[1]"=lm1_mcmc_combi[,"sigma_alpha_p[1]"],"sigma_alpha_p[2]"=lm1_mcmc_combi[,"sigma_alpha_p[2]"],
                                      "sigma_alpha_p[3]"=lm1_mcmc_combi[,"sigma_alpha_p[3]"],"sigma_alpha_p[4]"=lm1_mcmc_combi[,"sigma_alpha_p[4]"],
                                      "sigma_alpha_p[5]"=lm1_mcmc_combi[,"sigma_alpha_p[5]"],"sigma_alpha_p[6]"=lm1_mcmc_combi[,"sigma_alpha_p[6]"],
                                      "sigma_alpha_p[7]"=lm1_mcmc_combi[,"sigma_alpha_p[7]"],"sigma_alpha_p[8]"=lm1_mcmc_combi[,"sigma_alpha_p[8]"],
                                      "sigma_alpha_p[9]"=lm1_mcmc_combi[,"sigma_alpha_p[9]"],"sigma_alpha_p[10]"=lm1_mcmc_combi[,"sigma_alpha_p[10]"],
                                      "sigma_alpha_p[11]"=lm1_mcmc_combi[,"sigma_alpha_p[11]"],"sigma_alpha_p[12]"=lm1_mcmc_combi[,"sigma_alpha_p[12]"],
                                      "sigma_alpha_p[13]"=lm1_mcmc_combi[,"sigma_alpha_p[13]"],"sigma_alpha_p[14]"=lm1_mcmc_combi[,"sigma_alpha_p[14]"],
                                      "sigma_alpha_p[15]"=lm1_mcmc_combi[,"sigma_alpha_p[15]"],"sigma_alpha_p[16]"=lm1_mcmc_combi[,"sigma_alpha_p[16]"],
                                      "sigma_alpha_p[17]"=lm1_mcmc_combi[,"sigma_alpha_p[17]"],"sigma_alpha_p[18]"=lm1_mcmc_combi[,"sigma_alpha_p[18]"],
                                      "sigma_alpha_p[19]"=lm1_mcmc_combi[,"sigma_alpha_p[19]"],"sigma_alpha_p[20]"=lm1_mcmc_combi[,"sigma_alpha_p[20]"],
                                      "sigma_alpha_p[21]"=lm1_mcmc_combi[,"sigma_alpha_p[21]"],"sigma_alpha_p[22]"=lm1_mcmc_combi[,"sigma_alpha_p[22]"],
                                      "sigma_alpha_p[23]"=lm1_mcmc_combi[,"sigma_alpha_p[23]"],"sigma_alpha_p[24]"=lm1_mcmc_combi[,"sigma_alpha_p[24]"],
                                      "sigma_alpha_p[25]"=lm1_mcmc_combi[,"sigma_alpha_p[25]"],"sigma_alpha_p[26]"=lm1_mcmc_combi[,"sigma_alpha_p[26]"],
                                      "sigma_alpha_p[27]"=lm1_mcmc_combi[,"sigma_alpha_p[27]"],"sigma_alpha_p[28]"=lm1_mcmc_combi[,"sigma_alpha_p[28]"],
                                      "sigma_alpha_p[29]"=lm1_mcmc_combi[,"sigma_alpha_p[29]"],"sigma_alpha_p[30]"=lm1_mcmc_combi[,"sigma_alpha_p[30]"],
                                      "sigma_alpha_p[31]"=lm1_mcmc_combi[,"sigma_alpha_p[31]"],"sigma_alpha_p[32]"=lm1_mcmc_combi[,"sigma_alpha_p[32]"],
                                      "sigma_alpha_p[33]"=lm1_mcmc_combi[,"sigma_alpha_p[33]"],"sigma_alpha_p[34]"=lm1_mcmc_combi[,"sigma_alpha_p[34]"],
                                      "sigma_alpha_p[35]"=lm1_mcmc_combi[,"sigma_alpha_p[35]"],"sigma_alpha_p[36]"=lm1_mcmc_combi[,"sigma_alpha_p[26]"],
                                      "sigma_alpha_p[37]"=lm1_mcmc_combi[,"sigma_alpha_p[37]"],"sigma_alpha_p[38]"=lm1_mcmc_combi[,"sigma_alpha_p[38]"],
                                      "sigma_alpha_p[39]"=lm1_mcmc_combi[,"sigma_alpha_p[39]"],"sigma_alpha_p[40]"=lm1_mcmc_combi[,"sigma_alpha_p[40]"]), 1)) 
#parameters for final-demand-to-gdp ratio
x2=(data_r [1:21,2:41])
x2_pre=(data_r[1:56,2:41])
d2<- matrix(NA, nrow =n_2050-1,ncol=n_region)
for  (t in 1:n_2050-1){
  for  (j in 1:n_region){
    d2[t,j]<-x2_pre[t+1,j]-x2_pre[t,j]
  }
}

library(R2jags)

jagsdata_s1 <- with(data_r, list(x2=x2,n_2015=21,n_region=40))

lm1_jags <- function(){
  # Likelihood:
  for (j in 1:n_region){
    for (t in 2:n_2015){
      r[t,j] ~ dnorm(mu_r[t,j], tau_r[j]) # tau is precision (1 / variance)
      mu_r[t,j] <- (alpha_r[j]+x2[t-1,j])
    }
  }
  
  # Priors:
  for (j in 1:n_region){ 
    alpha_r[j] ~ dnorm(mu_alpha_r[j], tau_alpha_r[j]) # intercept
    mu_alpha_r[j]~ dunif(-1,1) 
    sigma_alpha_r[j] ~ dunif(0,0.05) # standard deviation
    tau_alpha_r[j] <- 1 / (sigma_alpha_r[j] * sigma_alpha_r[j]) # sigma^2 doesn't work in JAGS
    tau_r[j] <- 1 / (sigma_r[j] * sigma_r[j]) # sigma^2 doesn't work in JAGS
    sigma_r[j] ~ dunif(0,0.05) # standard deviation
  }
}

params <- c("alpha_r[1]","sigma_r[1]","alpha_r[2]","sigma_r[2]","alpha_r[3]","sigma_r[3]",
            "alpha_r[4]","sigma_r[4]","alpha_r[5]","sigma_r[5]","alpha_r[6]","sigma_r[6]", 
            "alpha_r[7]","sigma_r[7]","alpha_r[8]","sigma_r[8]","alpha_r[9]","sigma_r[9]",
            "alpha_r[10]","sigma_r[10]","alpha_r[11]", "sigma_r[11]","alpha_r[12]","sigma_r[12]",
            "alpha_r[13]","sigma_r[13]","alpha_r[14]","sigma_r[14]","alpha_r[15]","sigma_r[15]",
            "alpha_r[16]","sigma_r[16]","alpha_r[17]","sigma_r[17]","alpha_r[18]","sigma_r[18]",
            "alpha_r[19]","sigma_r[19]","alpha_r[20]","sigma_r[20]","alpha_r[21]","sigma_r[21]", 
            "alpha_r[22]","sigma_r[22]","alpha_r[23]","sigma_r[23]","alpha_r[24]","sigma_r[24]",
            "alpha_r[25]","sigma_r[25]","alpha_r[26]","sigma_r[26]","alpha_r[27]","sigma_r[27]",
            "alpha_r[28]","sigma_r[28]","alpha_r[29]","sigma_r[29]","alpha_r[30]","sigma_r[30]",
            "alpha_r[31]","sigma_r[31]","alpha_r[32]","sigma_r[32]","alpha_r[33]","sigma_r[33]",
            "alpha_r[34]","sigma_r[34]","alpha_r[35]","sigma_r[35]","alpha_r[36]","sigma_r[36]", 
            "alpha_r[37]","sigma_r[37]","alpha_r[38]","sigma_r[38]","alpha_r[39]","sigma_r[39]",
            "alpha_r[40]","sigma_r[40]","sigma_alpha_r[1]","sigma_alpha_r[2]","sigma_alpha_r[3]",
            "sigma_alpha_r[4]","sigma_alpha_r[5]","sigma_alpha_r[6]","sigma_alpha_r[7]","sigma_alpha_r[8]",
            "sigma_alpha_r[9]","sigma_alpha_r[10]","sigma_alpha_r[11]","sigma_alpha_r[12]","sigma_alpha_r[13]",
            "sigma_alpha_r[14]","sigma_alpha_r[15]","sigma_alpha_r[16]","sigma_alpha_r[17]","sigma_alpha_r[18]",
            "sigma_alpha_r[19]","sigma_alpha_r[20]","sigma_alpha_r[21]","sigma_alpha_r[22]","sigma_alpha_r[23]",
            "sigma_alpha_r[24]","sigma_alpha_r[25]","sigma_alpha_r[26]","sigma_alpha_r[27]","sigma_alpha_r[28]",
            "sigma_alpha_r[29]","sigma_alpha_r[30]","sigma_alpha_r[31]","sigma_alpha_r[32]","sigma_alpha_r[33]",
            "sigma_alpha_r[34]","sigma_alpha_r[35]","sigma_alpha_r[36]","sigma_alpha_r[37]","sigma_alpha_r[38]","sigma_alpha_r[39]","sigma_alpha_r[40]")

fit_lm1 <- jags(data = jagsdata_s1,  parameters.to.save = params, model.file = lm1_jags,
                n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC = F)

fit_lm1

lm1_mcmc <- as.mcmc(fit_lm1)
lm1_mcmc_combi <- as.mcmc(rbind(lm1_mcmc[[1]], lm1_mcmc[[2]], lm1_mcmc[[3]],lm1_mcmc[[4]], lm1_mcmc[[5]], lm1_mcmc[[6]],lm1_mcmc[[7]]))
lm1_mcmc_r <- do.call(cbind, rep(list("alpha_r[1]"=lm1_mcmc_combi[,"alpha_r[1]"],"alpha_r[2]"=lm1_mcmc_combi[,"alpha_r[2]"],
                                      "alpha_r[3]"=lm1_mcmc_combi[,"alpha_r[3]"],"alpha_r[4]"=lm1_mcmc_combi[,"alpha_r[4]"],
                                      "alpha_r[5]"=lm1_mcmc_combi[,"alpha_r[5]"],"alpha_r[6]"=lm1_mcmc_combi[,"alpha_r[6]"],
                                      "alpha_r[7]"=lm1_mcmc_combi[,"alpha_r[7]"],"alpha_r[8]"=lm1_mcmc_combi[,"alpha_r[8]"],
                                      "alpha_r[9]"=lm1_mcmc_combi[,"alpha_r[9]"],"alpha_r[10]"=lm1_mcmc_combi[,"alpha_r[10]"],
                                      "alpha_r[11]"=lm1_mcmc_combi[,"alpha_r[11]"],"alpha_r[12]"=lm1_mcmc_combi[,"alpha_r[12]"],
                                      "alpha_r[13]"=lm1_mcmc_combi[,"alpha_r[13]"],"alpha_r[14]"=lm1_mcmc_combi[,"alpha_r[14]"],
                                      "alpha_r[15]"=lm1_mcmc_combi[,"alpha_r[15]"],"alpha_r[16]"=lm1_mcmc_combi[,"alpha_r[16]"],
                                      "alpha_r[17]"=lm1_mcmc_combi[,"alpha_r[17]"],"alpha_r[18]"=lm1_mcmc_combi[,"alpha_r[18]"],
                                      "alpha_r[19]"=lm1_mcmc_combi[,"alpha_r[19]"],"alpha_r[20]"=lm1_mcmc_combi[,"alpha_r[20]"],
                                      "alpha_r[21]"=lm1_mcmc_combi[,"alpha_r[21]"],"alpha_r[22]"=lm1_mcmc_combi[,"alpha_r[22]"],
                                      "alpha_r[23]"=lm1_mcmc_combi[,"alpha_r[23]"],"alpha_r[24]"=lm1_mcmc_combi[,"alpha_r[24]"],
                                      "alpha_r[25]"=lm1_mcmc_combi[,"alpha_r[25]"],"alpha_r[26]"=lm1_mcmc_combi[,"alpha_r[26]"],
                                      "alpha_r[27]"=lm1_mcmc_combi[,"alpha_r[27]"],"alpha_r[28]"=lm1_mcmc_combi[,"alpha_r[28]"],
                                      "alpha_r[29]"=lm1_mcmc_combi[,"alpha_r[29]"],"alpha_r[30]"=lm1_mcmc_combi[,"alpha_r[30]"],
                                      "alpha_r[31]"=lm1_mcmc_combi[,"alpha_r[31]"],"alpha_r[32]"=lm1_mcmc_combi[,"alpha_r[32]"],
                                      "alpha_r[33]"=lm1_mcmc_combi[,"alpha_r[33]"],"alpha_r[34]"=lm1_mcmc_combi[,"alpha_r[34]"],
                                      "alpha_r[35]"=lm1_mcmc_combi[,"alpha_r[35]"],"alpha_r[36]"=lm1_mcmc_combi[,"alpha_r[36]"],
                                      "alpha_r[37]"=lm1_mcmc_combi[,"alpha_r[37]"],"alpha_r[38]"=lm1_mcmc_combi[,"alpha_r[38]"],
                                      "alpha_r[39]"=lm1_mcmc_combi[,"alpha_r[39]"],"alpha_r[40]"=lm1_mcmc_combi[,"alpha_r[40]"],
                                      "sigma_r[1]"=lm1_mcmc_combi[,"sigma_r[1]"],"sigma_r[2]"=lm1_mcmc_combi[,"sigma_r[2]"],
                                      "sigma_r[3]"=lm1_mcmc_combi[,"sigma_r[3]"],"sigma_r[4]"=lm1_mcmc_combi[,"sigma_r[4]"],
                                      "sigma_r[5]"=lm1_mcmc_combi[,"sigma_r[5]"],"sigma_r[6]"=lm1_mcmc_combi[,"sigma_r[6]"],
                                      "sigma_r[7]"=lm1_mcmc_combi[,"sigma_r[7]"],"sigma_r[8]"=lm1_mcmc_combi[,"sigma_r[8]"],
                                      "sigma_r[9]"=lm1_mcmc_combi[,"sigma_r[9]"],"sigma_r[10]"=lm1_mcmc_combi[,"sigma_r[10]"],
                                      "sigma_r[11]"=lm1_mcmc_combi[,"sigma_r[11]"],"sigma_r[12]"=lm1_mcmc_combi[,"sigma_r[12]"],
                                      "sigma_r[13]"=lm1_mcmc_combi[,"sigma_r[13]"],"sigma_r[14]"=lm1_mcmc_combi[,"sigma_r[14]"],
                                      "sigma_r[15]"=lm1_mcmc_combi[,"sigma_r[15]"],"sigma_r[16]"=lm1_mcmc_combi[,"sigma_r[16]"],
                                      "alpha_r[17]"=lm1_mcmc_combi[,"sigma_r[17]"],"sigma_r[18]"=lm1_mcmc_combi[,"sigma_r[18]"],
                                      "alpha_r[19]"=lm1_mcmc_combi[,"sigma_r[19]"],"sigma_r[20]"=lm1_mcmc_combi[,"sigma_r[20]"],
                                      "sigma_r[21]"=lm1_mcmc_combi[,"sigma_r[21]"],"sigma_r[22]"=lm1_mcmc_combi[,"sigma_r[22]"],
                                      "sigma_r[23]"=lm1_mcmc_combi[,"sigma_r[23]"],"sigma_r[24]"=lm1_mcmc_combi[,"sigma_r[24]"],
                                      "sigma_r[25]"=lm1_mcmc_combi[,"sigma_r[25]"],"sigma_r[26]"=lm1_mcmc_combi[,"sigma_r[26]"],
                                      "sigma_r[27]"=lm1_mcmc_combi[,"sigma_r[27]"],"sigma_r[28]"=lm1_mcmc_combi[,"sigma_r[28]"],
                                      "sigma_r[29]"=lm1_mcmc_combi[,"sigma_r[29]"],"sigma_r[30]"=lm1_mcmc_combi[,"sigma_r[30]"],
                                      "sigma_r[31]"=lm1_mcmc_combi[,"sigma_r[31]"],"sigma_r[32]"=lm1_mcmc_combi[,"sigma_r[32]"],
                                      "sigma_r[33]"=lm1_mcmc_combi[,"sigma_r[33]"],"sigma_r[34]"=lm1_mcmc_combi[,"sigma_r[34]"],
                                      "sigma_r[35]"=lm1_mcmc_combi[,"sigma_r[35]"],"sigma_r[36]"=lm1_mcmc_combi[,"sigma_r[36]"],
                                      "sigma_r[37]"=lm1_mcmc_combi[,"sigma_r[37]"],"sigma_r[38]"=lm1_mcmc_combi[,"sigma_r[38]"],
                                      "sigma_r[39]"=lm1_mcmc_combi[,"sigma_r[39]"],"sigma_r[40]"=lm1_mcmc_combi[,"sigma_r[40]"],
                                      "sigma_alpha_r[1]"=lm1_mcmc_combi[,"sigma_alpha_r[1]"],"sigma_alpha_r[2]"=lm1_mcmc_combi[,"sigma_alpha_r[2]"],
                                      "sigma_alpha_r[3]"=lm1_mcmc_combi[,"sigma_alpha_r[3]"],"sigma_alpha_r[4]"=lm1_mcmc_combi[,"sigma_alpha_r[4]"],
                                      "sigma_alpha_r[5]"=lm1_mcmc_combi[,"sigma_alpha_r[5]"],"sigma_alpha_r[6]"=lm1_mcmc_combi[,"sigma_alpha_r[6]"],
                                      "sigma_alpha_r[7]"=lm1_mcmc_combi[,"sigma_alpha_r[7]"],"sigma_alpha_r[8]"=lm1_mcmc_combi[,"sigma_alpha_r[8]"],
                                      "sigma_alpha_r[9]"=lm1_mcmc_combi[,"sigma_alpha_r[9]"],"sigma_alpha_r[10]"=lm1_mcmc_combi[,"sigma_alpha_r[10]"],
                                      "sigma_alpha_r[11]"=lm1_mcmc_combi[,"sigma_alpha_r[11]"],"sigma_alpha_r[12]"=lm1_mcmc_combi[,"sigma_alpha_r[12]"],
                                      "sigma_alpha_r[13]"=lm1_mcmc_combi[,"sigma_alpha_r[13]"],"sigma_alpha_r[14]"=lm1_mcmc_combi[,"sigma_alpha_r[14]"],
                                      "sigma_alpha_r[15]"=lm1_mcmc_combi[,"sigma_alpha_r[15]"],"sigma_alpha_r[16]"=lm1_mcmc_combi[,"sigma_alpha_r[16]"],
                                      "sigma_alpha_r[17]"=lm1_mcmc_combi[,"sigma_alpha_r[17]"],"sigma_alpha_r[18]"=lm1_mcmc_combi[,"sigma_alpha_r[18]"],
                                      "sigma_alpha_r[19]"=lm1_mcmc_combi[,"sigma_alpha_r[19]"],"sigma_alpha_r[20]"=lm1_mcmc_combi[,"sigma_alpha_r[20]"],
                                      "sigma_alpha_r[21]"=lm1_mcmc_combi[,"sigma_alpha_r[21]"],"sigma_alpha_r[22]"=lm1_mcmc_combi[,"sigma_alpha_r[22]"],
                                      "sigma_alpha_r[23]"=lm1_mcmc_combi[,"sigma_alpha_r[23]"],"sigma_alpha_r[24]"=lm1_mcmc_combi[,"sigma_alpha_r[24]"],
                                      "sigma_alpha_r[25]"=lm1_mcmc_combi[,"sigma_alpha_r[25]"],"sigma_alpha_r[26]"=lm1_mcmc_combi[,"sigma_alpha_r[26]"],
                                      "sigma_alpha_r[27]"=lm1_mcmc_combi[,"sigma_alpha_r[27]"],"sigma_alpha_r[28]"=lm1_mcmc_combi[,"sigma_alpha_r[28]"],
                                      "sigma_alpha_r[29]"=lm1_mcmc_combi[,"sigma_alpha_r[29]"],"sigma_alpha_r[30]"=lm1_mcmc_combi[,"sigma_alpha_r[30]"],
                                      "sigma_alpha_r[31]"=lm1_mcmc_combi[,"sigma_alpha_r[31]"],"sigma_alpha_r[32]"=lm1_mcmc_combi[,"sigma_alpha_r[32]"],
                                      "sigma_alpha_r[33]"=lm1_mcmc_combi[,"sigma_alpha_r[33]"],"sigma_alpha_r[34]"=lm1_mcmc_combi[,"sigma_alpha_r[34]"],
                                      "sigma_alpha_r[35]"=lm1_mcmc_combi[,"sigma_alpha_r[35]"],"sigma_alpha_r[36]"=lm1_mcmc_combi[,"sigma_alpha_r[26]"],
                                      "sigma_alpha_r[37]"=lm1_mcmc_combi[,"sigma_alpha_r[37]"],"sigma_alpha_r[38]"=lm1_mcmc_combi[,"sigma_alpha_r[38]"],
                                      "sigma_alpha_r[39]"=lm1_mcmc_combi[,"sigma_alpha_r[39]"],"sigma_alpha_r[40]"=lm1_mcmc_combi[,"sigma_alpha_r[40]"]), 1)) 
#parameters for benchmarks of emission intensity and consumption patterns factors
jagsdata_s1 <- with(data_gdp, list(x=x,x_2023=x_2023,
                                   eih=eih,eiv=eiv,eig=eig,s=s,v=v,y=yo,y_ref=yo,
                                   n_2015=21,n_region=40))

lm1_jags <- function(){
  # Likelihood:
  
  for (t in 1:n_2015){
    eih[t,33] ~ dnorm(mu_fh1[t], tau_fh1) # tau is precision (1 / variance)
    mu_fh1[t] <- (alpha_fh1 + beta_fh1* (x[t,33]))
    eih[t+21,34] ~ dnorm(mu_fh2[t], tau_fh2) # tau is precision (1 / variance)
    mu_fh2[t] <- (alpha_fh2 + beta_fh2* (x[t,34]))
    eih[t+21*2,33] ~ dnorm(mu_fh3[t], tau_fh3) # tau is precision (1 / variance)
    mu_fh3[t] <- (alpha_fh3 + beta_fh3* (x[t,33]))
    eih[t+21*3,20] ~ dnorm(mu_fh4[t], tau_fh4) # tau is precision (1 / variance)
    mu_fh4[t] <- (alpha_fh4 + beta_fh4* (x[t,20]))
    eih[t+21*4,34] ~ dnorm(mu_fh5[t], tau_fh5) # tau is precision (1 / variance)
    mu_fh5[t] <- (alpha_fh5 + beta_fh5* (x[t,34]))
    eih[t+21*5,33] ~ dnorm(mu_fh6[t], tau_fh6) # tau is precision (1 / variance)
    mu_fh6[t] <- (alpha_fh6 + beta_fh6* (x[t,33]))
    
    eiv[t,35] ~ dnorm(mu_fv1[t], tau_fv1) # tau is precision (1 / variance)
    mu_fv1[t] <- (alpha_fv1 + beta_fv1* (x[t,35]))
    eiv[t+21,33] ~ dnorm(mu_fv2[t], tau_fv2) # tau is precision (1 / variance)
    mu_fv2[t] <- (alpha_fv2 + beta_fv2* (x[t,33]))
    eiv[t+21*2,33] ~ dnorm(mu_fv3[t], tau_fv3) # tau is precision (1 / variance)
    mu_fv3[t] <- (alpha_fv3 + beta_fv3* (x[t,33]))
    eiv[t+21*3,33] ~ dnorm(mu_fv4[t], tau_fv4) # tau is precision (1 / variance)
    mu_fv4[t] <- (alpha_fv4 + beta_fv4* (x[t,33]))
    eiv[t+21*4,20] ~ dnorm(mu_fv5[t], tau_fv5) # tau is precision (1 / variance)
    mu_fv5[t] <- (alpha_fv5 + beta_fv5* (x[t,20]))
    
    eig[t,20] ~ dnorm(mu_fg1[t], tau_fg1) # tau is precision (1 / variance)
    mu_fg1[t] <- (alpha_fg1 + beta_fg1* (x[t,20]))
  }
  for (j in 1:n_region){
    for (t in 1:n_2015){
      s[t,j] ~ dnorm(mu_s1[t,j], tau_s1) # tau is precision (1 / variance)
      mu_s1[t,j] <- (alpha_s1 + beta_s1* (x[t,j]))
      s[t+21,j] ~ dnorm(mu_s2[t,j], tau_s2) # tau is precision (1 / variance)
      mu_s2[t,j] <- (alpha_s2 + beta_s2* (x[t,j]))
      s[t+21*2,j] ~ dnorm(mu_s3[t,j], tau_s3) # tau is precision (1 / variance)
      mu_s3[t,j] <- (alpha_s3 + beta_s3* (x[t,j]))
      s[t+21*3,j] ~ dnorm(mu_s4[t,j], tau_s4) # tau is precision (1 / variance)
      mu_s4[t,j] <- (alpha_s4 + beta_s4* (x[t,j]))
      s[t+21*4,j] ~ dnorm(mu_s5[t,j], tau_s5) # tau is precision (1 / variance)
      mu_s5[t,j] <- (alpha_s5 + beta_s5* (x[t,j]))
      s[t+21*5,j] ~ dnorm(mu_s6[t,j], tau_s6) # tau is precision (1 / variance)
      mu_s6[t,j] <- (alpha_s6 + beta_s6* (x[t,j]))
      
      v[t,j] ~ dnorm(mu_v1[t,j], tau_v1) # tau is precision (1 / variance)
      mu_v1[t,j] <- (alpha_v1 + beta_v1* (x[t,j]))
      v[t+21,j] ~ dnorm(mu_v2[t,j], tau_v2) # tau is precision (1 / variance)
      mu_v2[t,j] <- (alpha_v2 + beta_v2* (x[t,j]))
      v[t+21*2,j] ~ dnorm(mu_v3[t,j], tau_v3) # tau is precision (1 / variance)
      mu_v3[t,j] <- (alpha_v3 + beta_v3* (x[t,j]))
      v[t+21*3,j] ~ dnorm(mu_v4[t,j], tau_v4) # tau is precision (1 / variance)
      mu_v4[t,j] <- (alpha_v4 + beta_v4* (x[t,j]))
      v[t+21*4,j] ~ dnorm(mu_v5[t,j], tau_v5) # tau is precision (1 / variance)
      mu_v5[t,j] <- (alpha_v5 + beta_v5* (x[t,j]))
    }
  }
  
  for (j in 1:36){
    for (t in 1:29){
      y[t,j] ~ dnorm(mu_y1[t,j], tau_y1) # tau is precision (1 / variance)
      mu_y1[t,j] <- (alpha_y1 + beta_y1* (x_2023[t,j]))
      y[t+29,j] ~ dnorm(mu_y2[t,j], tau_y2) # tau is precision (1 / variance)
      mu_y2[t,j] <- (alpha_y2 + beta_y2* (x_2023[t,j]))
      y[t+29*2,j] ~ dnorm(mu_y3[t,j], tau_y3) # tau is precision (1 / variance)
      mu_y3[t,j] <- (alpha_y3 + beta_y3* (x_2023[t,j]))
    }
  }
  
  for (t in 1:29){
    y_ref[t,37] ~ dnorm(mu_y1_ref[t], tau_y1_ref) # tau is precision (1 / variance)
    mu_y1_ref[t] <- (alpha_y1_ref + beta2_y1_ref* (x_2023[t,37]* x_2023[t,37])+ beta1_y1_ref* (x_2023[t,37]))
    y_ref[t+29,37] ~ dnorm(mu_y2_ref[t], tau_y2_ref) # tau is precision (1 / variance)
    mu_y2_ref[t] <- (alpha_y2_ref + beta2_y2_ref* (x_2023[t,37]* x_2023[t,37])+ beta1_y2_ref* (x_2023[t,37]))
    y_ref[t+29*2,37] ~ dnorm(mu_y3_ref[t], tau_y3_ref) # tau is precision (1 / variance)
    mu_y3_ref[t] <- (alpha_y3_ref + beta2_y3_ref* (x_2023[t,37]* x_2023[t,37])+ beta1_y3_ref* (x_2023[t,37]))
  }
  
  # Priors£º
  alpha_fh1 ~ dnorm(mu_alpha_fh1, tau_alpha_fh1) # intercept
  beta_fh1 ~ dnorm(mu_beta_fh1, tau_beta_fh1) # slope
  mu_alpha_fh1~ dunif(-30,30) 
  mu_beta_fh1~ dunif(-5, 5) 
  sigma_alpha_fh1 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh1 <- 1 / (sigma_alpha_fh1 * sigma_alpha_fh1) # sigma^2 doesn't work in JAGS
  sigma_beta_fh1~ dunif(0, 1) # standard deviation
  tau_beta_fh1<- 1 / (sigma_beta_fh1 * sigma_beta_fh1) # sigma^2 doesn't work in JAGS
  sigma_fh1 ~ dunif(0, 0.05) # standard deviation
  tau_fh1 <- 1 / (sigma_fh1 * sigma_fh1) # sigma^2 doesn't work in JAGS
  
  alpha_fh2 ~ dnorm(mu_alpha_fh2, tau_alpha_fh2) # intercept
  beta_fh2 ~ dnorm(mu_beta_fh2, tau_beta_fh2) # slope
  mu_alpha_fh2~ dunif(-30,30) 
  mu_beta_fh2~ dunif(-5,5) 
  sigma_alpha_fh2 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh2 <- 1 / (sigma_alpha_fh2 * sigma_alpha_fh2) # sigma^2 doesn't work in JAGS
  sigma_beta_fh2~ dunif(0, 1) # standard deviation
  tau_beta_fh2<- 1 / (sigma_beta_fh2 * sigma_beta_fh2) # sigma^2 doesn't work in JAGS
  sigma_fh2 ~ dunif(0, 0.05) # standard deviation
  tau_fh2 <- 1 / (sigma_fh2 * sigma_fh2) # sigma^2 doesn't work in JAGS
  
  alpha_fh3 ~ dnorm(mu_alpha_fh3, tau_alpha_fh3) # intercept
  beta_fh3 ~ dnorm(mu_beta_fh3, tau_beta_fh3) # slope
  mu_alpha_fh3~ dunif(-30,30) 
  mu_beta_fh3~ dunif(-5, 5) 
  sigma_alpha_fh3 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh3 <- 1 / (sigma_alpha_fh3 * sigma_alpha_fh3) # sigma^2 doesn't work in JAGS
  sigma_beta_fh3~ dunif(0,1) # standard deviation
  tau_beta_fh3<- 1 / (sigma_beta_fh3 * sigma_beta_fh3) # sigma^2 doesn't work in JAGS
  sigma_fh3 ~ dunif(0, 0.05)# standard deviation
  tau_fh3 <- 1 / (sigma_fh3 * sigma_fh3) # sigma^2 doesn't work in JAGS
  
  alpha_fh4 ~ dnorm(mu_alpha_fh4, tau_alpha_fh4) # intercept
  beta_fh4 ~ dnorm(mu_beta_fh4, tau_beta_fh4) # slope
  mu_alpha_fh4~ dunif(-30,30) 
  mu_beta_fh4~ dunif(-5, 5) 
  sigma_alpha_fh4 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh4 <- 1 / (sigma_alpha_fh4 * sigma_alpha_fh4) # sigma^2 doesn't work in JAGS
  sigma_beta_fh4~ dunif(0, 1) # standard deviation
  tau_beta_fh4<- 1 / (sigma_beta_fh4 * sigma_beta_fh4) # sigma^2 doesn't work in JAGS
  sigma_fh4 ~ dunif(0, 0.05) # standard deviation
  tau_fh4 <- 1 / (sigma_fh4 * sigma_fh4) # sigma^2 doesn't work in JAGS
  
  alpha_fh5 ~ dnorm(mu_alpha_fh5, tau_alpha_fh5) # intercept
  beta_fh5 ~ dnorm(mu_beta_fh5, tau_beta_fh5) # slope
  mu_alpha_fh5~ dunif(-30,30) 
  mu_beta_fh5~ dunif(-5, 5) 
  sigma_alpha_fh5 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh5 <- 1 / (sigma_alpha_fh5 * sigma_alpha_fh5) # sigma^2 doesn't work in JAGS
  sigma_beta_fh5~ dunif(0, 1) # standard deviation
  tau_beta_fh5<- 1 / (sigma_beta_fh5* sigma_beta_fh5) # sigma^2 doesn't work in JAGS
  sigma_fh5 ~ dunif(0, 0.05) # standard deviation
  tau_fh5 <- 1 / (sigma_fh5 * sigma_fh5) # sigma^2 doesn't work in JAGS
  
  alpha_fh6 ~ dnorm(mu_alpha_fh6, tau_alpha_fh6) # intercept
  beta_fh6 ~ dnorm(mu_beta_fh6, tau_beta_fh6) # slope
  mu_alpha_fh6~ dunif(-30,30) 
  mu_beta_fh6~ dunif(-5, 5) 
  sigma_alpha_fh6 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh6 <- 1 / (sigma_alpha_fh6 * sigma_alpha_fh6) # sigma^2 doesn't work in JAGS
  sigma_beta_fh6~ dunif(0, 1) # standard deviation
  tau_beta_fh6<- 1 / (sigma_beta_fh6* sigma_beta_fh6) # sigma^2 doesn't work in JAGS
  sigma_fh6 ~ dunif(0, 0.05) # standard deviation
  tau_fh6 <- 1 / (sigma_fh6 * sigma_fh6) # sigma^2 doesn't work in JAGS
  
  alpha_fv1 ~ dnorm(mu_alpha_fv1, tau_alpha_fv1) # intercept
  beta_fv1 ~ dnorm(mu_beta_fv1, tau_beta_fv1) # slope
  mu_alpha_fv1~ dunif(-30,30) 
  mu_beta_fv1~ dunif(-5, 5) 
  sigma_alpha_fv1 ~ dunif(0, 1) # standard deviation
  tau_alpha_fv1 <- 1 / (sigma_alpha_fv1 * sigma_alpha_fv1) # sigma^2 doesn't work in JAGS
  sigma_beta_fv1~ dunif(0,1) # standard deviation
  tau_beta_fv1<- 1 / (sigma_beta_fv1 * sigma_beta_fv1) # sigma^2 doesn't work in JAGS
  sigma_fv1 ~ dunif(0, 0.05) # standard deviation
  tau_fv1 <- 1 / (sigma_fv1 * sigma_fv1) # sigma^2 doesn't work in JAGS
  
  alpha_fv2 ~ dnorm(mu_alpha_fv2, tau_alpha_fv2) # intercept
  beta_fv2 ~ dnorm(mu_beta_fv2, tau_beta_fv2) # slope
  mu_alpha_fv2~ dunif(-30,30) 
  mu_beta_fv2~ dunif(-5, 5) 
  sigma_alpha_fv2 ~ dunif(0, 1) # standard deviation
  tau_alpha_fv2 <- 1 / (sigma_alpha_fv2 * sigma_alpha_fv2) # sigma^2 doesn't work in JAGS
  sigma_beta_fv2~ dunif(0,1) # standard deviation
  tau_beta_fv2<- 1 / (sigma_beta_fv2 * sigma_beta_fv2) # sigma^2 doesn't work in JAGS
  sigma_fv2 ~ dunif(0, 0.05) # standard deviation
  tau_fv2 <- 1 / (sigma_fv2 * sigma_fv2) # sigma^2 doesn't work in JAGS
  
  alpha_fv3 ~ dnorm(mu_alpha_fv3, tau_alpha_fv3) # intercept
  beta_fv3 ~ dnorm(mu_beta_fv3, tau_beta_fv3) # slope
  mu_alpha_fv3~ dunif(-30,30) 
  mu_beta_fv3~ dunif(-5, 5) 
  sigma_alpha_fv3 ~ dunif(0, 1) # standard deviation
  tau_alpha_fv3 <- 1 / (sigma_alpha_fv3 * sigma_alpha_fv3) # sigma^2 doesn't work in JAGS
  sigma_beta_fv3~ dunif(0,1) # standard deviation
  tau_beta_fv3<- 1 / (sigma_beta_fv3 * sigma_beta_fv3) # sigma^2 doesn't work in JAGS
  sigma_fv3 ~ dunif(0, 0.05) # standard deviation
  tau_fv3 <- 1 / (sigma_fv3 * sigma_fv3) # sigma^2 doesn't work in JAGS
  
  alpha_fv4 ~ dnorm(mu_alpha_fv4, tau_alpha_fv4) # intercept
  beta_fv4 ~ dnorm(mu_beta_fv4, tau_beta_fv4) # slope
  mu_alpha_fv4~ dunif(-30,30) 
  mu_beta_fv4~ dunif(-5, 5) 
  sigma_alpha_fv4 ~ dunif(0, 1) # standard deviation
  tau_alpha_fv4 <- 1 / (sigma_alpha_fv4 * sigma_alpha_fv4) # sigma^2 doesn't work in JAGS
  sigma_beta_fv4~ dunif(0,1) # standard deviation
  tau_beta_fv4<- 1 / (sigma_beta_fv4 * sigma_beta_fv4) # sigma^2 doesn't work in JAGS
  sigma_fv4 ~ dunif(0, 0.05) # standard deviation
  tau_fv4 <- 1 / (sigma_fv4 * sigma_fv4) # sigma^2 doesn't work in JAGS
  
  alpha_fv5 ~ dnorm(mu_alpha_fv5, tau_alpha_fv5) # intercept
  beta_fv5 ~ dnorm(mu_beta_fv5, tau_beta_fv5) # slope
  mu_alpha_fv5~ dunif(-30,30) 
  mu_beta_fv5~ dunif(-5, 5) 
  sigma_alpha_fv5 ~ dunif(0, 1) # standard deviation
  tau_alpha_fv5 <- 1 / (sigma_alpha_fv5 * sigma_alpha_fv5) # sigma^2 doesn't work in JAGS
  sigma_beta_fv5~ dunif(0,1) # standard deviation
  tau_beta_fv5<- 1 / (sigma_beta_fv5 * sigma_beta_fv5) # sigma^2 doesn't work in JAGS
  sigma_fv5 ~ dunif(0, 0.05) # standard deviation
  tau_fv5 <- 1 / (sigma_fv5 * sigma_fv5) # sigma^2 doesn't work in JAGS
  
  alpha_fg1 ~ dnorm(mu_alpha_fg1, tau_alpha_fg1) # intercept
  beta_fg1 ~ dnorm(mu_beta_fg1, tau_beta_fg1) # slope
  mu_alpha_fg1~ dunif(-30,30) 
  mu_beta_fg1~ dunif(-5, 5) 
  sigma_alpha_fg1 ~ dunif(0, 1) # standard deviation
  tau_alpha_fg1 <- 1 / (sigma_alpha_fg1 * sigma_alpha_fg1) # sigma^2 doesn't work in JAGS
  sigma_beta_fg1~ dunif(0,1) # standard deviation
  tau_beta_fg1<- 1 / (sigma_beta_fg1 * sigma_beta_fg1) # sigma^2 doesn't work in JAGS
  sigma_fg1 ~ dunif(0, 0.05) # standard deviation
  tau_fg1 <- 1 / (sigma_fg1 * sigma_fg1) # sigma^2 doesn't work in JAGS
  
  alpha_s1 ~ dnorm(mu_alpha_s1, tau_alpha_s1) # intercept
  beta_s1 ~ dnorm(mu_beta_s1, tau_beta_s1) # slope
  mu_alpha_s1~ dunif(-1, 1) 
  mu_beta_s1~ dunif(-1, 1) 
  sigma_alpha_s1 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s1 <- 1 / (sigma_alpha_s1 * sigma_alpha_s1) # sigma^2 doesn't work in JAGS
  sigma_beta_s1~ dunif(0, 0.1) # standard deviation
  tau_beta_s1<- 1 / (sigma_beta_s1 * sigma_beta_s1) # sigma^2 doesn't work in JAGS
  sigma_s1 ~ dunif(0,0.02)# standard deviation
  tau_s1 <- 1 / (sigma_s1 * sigma_s1) # sigma^2 doesn't work in JAGS
  
  alpha_s2 ~ dnorm(mu_alpha_s2, tau_alpha_s2) # intercept
  beta_s2 ~ dnorm(mu_beta_s2, tau_beta_s2) # slope
  mu_alpha_s2~ dunif(-1, 1) 
  mu_beta_s2~ dunif(-1, 1) 
  sigma_alpha_s2 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s2 <- 1 / (sigma_alpha_s2 * sigma_alpha_s2) # sigma^2 doesn't work in JAGS
  sigma_beta_s2~ dunif(0, 0.1) # standard deviation
  tau_beta_s2<- 1 / (sigma_beta_s2 * sigma_beta_s2) # sigma^2 doesn't work in JAGS
  sigma_s2 ~ dunif(0, 0.02) # standard deviation
  tau_s2 <- 1 / (sigma_s2 * sigma_s2) # sigma^2 doesn't work in JAGS
  
  alpha_s3 ~ dnorm(mu_alpha_s3, tau_alpha_s3) # intercept
  beta_s3 ~ dnorm(mu_beta_s3, tau_beta_s3) # slope
  mu_alpha_s3~ dunif(-1, 1) 
  mu_beta_s3~ dunif(-1, 1)  
  sigma_alpha_s3 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s3 <- 1 / (sigma_alpha_s3 * sigma_alpha_s3) # sigma^2 doesn't work in JAGS
  sigma_beta_s3~ dunif(0, 0.1) # standard deviation
  tau_beta_s3<- 1 / (sigma_beta_s3 * sigma_beta_s3) # sigma^2 doesn't work in JAGS
  sigma_s3 ~ dunif(0, 0.02) # standard deviation
  tau_s3 <- 1 / (sigma_s3 * sigma_s3) # sigma^2 doesn't work in JAGS
  
  alpha_s4 ~ dnorm(mu_alpha_s4, tau_alpha_s4) # intercept
  beta_s4 ~ dnorm(mu_beta_s4, tau_beta_s4) # slope
  mu_alpha_s4~ dunif(-1, 1) 
  mu_beta_s4~ dunif(-1, 1) 
  sigma_alpha_s4 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s4 <- 1 / (sigma_alpha_s4 * sigma_alpha_s4) # sigma^2 doesn't work in JAGS
  sigma_beta_s4~ dunif(0, 0.1) # standard deviation
  tau_beta_s4<- 1 / (sigma_beta_s4 * sigma_beta_s4) # sigma^2 doesn't work in JAGS
  sigma_s4 ~ dunif(0, 0.02) # standard deviation
  tau_s4 <- 1 / (sigma_s4 * sigma_s4) # sigma^2 doesn't work in JAGS
  
  alpha_s5 ~ dnorm(mu_alpha_s5, tau_alpha_s5) # intercept
  beta_s5 ~ dnorm(mu_beta_s5, tau_beta_s5) # slope
  mu_alpha_s5~ dunif(-1, 1) 
  mu_beta_s5~ dunif(-1, 1) 
  sigma_alpha_s5 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s5 <- 1 / (sigma_alpha_s5 * sigma_alpha_s5) # sigma^2 doesn't work in JAGS
  sigma_beta_s5~ dunif(0, 0.1) # standard deviation
  tau_beta_s5<- 1 / (sigma_beta_s5* sigma_beta_s5) # sigma^2 doesn't work in JAGS
  sigma_s5 ~dunif(0, 0.02) # standard deviation
  tau_s5 <- 1 / (sigma_s5 * sigma_s5) # sigma^2 doesn't work in JAGS
  
  alpha_s6 ~ dnorm(mu_alpha_s6, tau_alpha_s6) # intercept
  beta_s6 ~ dnorm(mu_beta_s6, tau_beta_s6) # slope
  mu_alpha_s6~ dunif(-1, 1) 
  mu_beta_s6~ dunif(-1, 1) 
  sigma_alpha_s6 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s6 <- 1 / (sigma_alpha_s6 * sigma_alpha_s6) # sigma^2 doesn't work in JAGS
  sigma_beta_s6~ dunif(0, 0.1) # standard deviation
  tau_beta_s6<- 1 / (sigma_beta_s6* sigma_beta_s6) # sigma^2 doesn't work in JAGS
  sigma_s6 ~ dunif(0, 0.02) # standard deviation
  tau_s6 <- 1 / (sigma_s6 * sigma_s6) # sigma^2 doesn't work in JAGS
  
  alpha_v1 ~ dnorm(mu_alpha_v1, tau_alpha_v1) # intercept
  beta_v1 ~ dnorm(mu_beta_v1, tau_beta_v1) # slope
  mu_alpha_v1~ dunif(-1, 1) 
  mu_beta_v1~ dunif(-1, 1) 
  sigma_alpha_v1 ~ dunif(0, 1) # standard deviation
  tau_alpha_v1 <- 1 / (sigma_alpha_v1 * sigma_alpha_v1) # sigma^2 doesn't work in JAGS
  sigma_beta_v1~ dunif(0, 1) # standard deviation
  tau_beta_v1<- 1 / (sigma_beta_v1 * sigma_beta_v1) # sigma^2 doesn't work in JAGS
  sigma_v1 ~ dunif(0, 0.02) # standard deviation
  tau_v1 <- 1 / (sigma_v1 * sigma_v1) # sigma^2 doesn't work in JAGS
  
  alpha_v2 ~ dnorm(mu_alpha_v2, tau_alpha_v2) # intercept
  beta_v2 ~ dnorm(mu_beta_v2, tau_beta_v2) # slope
  mu_alpha_v2~ dunif(-1, 1) 
  mu_beta_v2~ dunif(-1, 1) 
  sigma_alpha_v2 ~ dunif(0, 1) # standard deviation
  tau_alpha_v2 <- 1 / (sigma_alpha_v2 * sigma_alpha_v2) # sigma^2 doesn't work in JAGS
  sigma_beta_v2~ dunif(0, 1) # standard deviation
  tau_beta_v2<- 1 / (sigma_beta_v2 * sigma_beta_v2) # sigma^2 doesn't work in JAGS
  sigma_v2 ~ dunif(0, 0.02) # standard deviation
  tau_v2 <- 1 / (sigma_v2 * sigma_v2) # sigma^2 doesn't work in JAGS
  
  alpha_v3 ~ dnorm(mu_alpha_v3, tau_alpha_v3) # intercept
  beta_v3 ~ dnorm(mu_beta_v3, tau_beta_v3) # slope
  mu_alpha_v3~ dunif(-1, 1) 
  mu_beta_v3~ dunif(-1, 1) 
  sigma_alpha_v3 ~ dunif(0, 1) # standard deviation
  tau_alpha_v3 <- 1 / (sigma_alpha_v3 * sigma_alpha_v3) # sigma^2 doesn't work in JAGS
  sigma_beta_v3~ dunif(0, 1) # standard deviation
  tau_beta_v3<- 1 / (sigma_beta_v3 * sigma_beta_v3) # sigma^2 doesn't work in JAGS
  sigma_v3 ~ dunif(0, 0.02) # standard deviation
  tau_v3 <- 1 / (sigma_v3 * sigma_v3) # sigma^2 doesn't work in JAGS
  
  alpha_v4 ~ dnorm(mu_alpha_v4, tau_alpha_v4) # intercept
  beta_v4 ~ dnorm(mu_beta_v4, tau_beta_v4) # slope
  mu_alpha_v4~ dunif(-1, 1) 
  mu_beta_v4~ dunif(-1, 1) 
  sigma_alpha_v4 ~ dunif(0, 1) # standard deviation
  tau_alpha_v4 <- 1 / (sigma_alpha_v4 * sigma_alpha_v4) # sigma^2 doesn't work in JAGS
  sigma_beta_v4~ dunif(0, 1) # standard deviation
  tau_beta_v4<- 1 / (sigma_beta_v4 * sigma_beta_v4) # sigma^2 doesn't work in JAGS
  sigma_v4 ~ dunif(0, 0.02) # standard deviation
  tau_v4 <- 1 / (sigma_v4 * sigma_v4) # sigma^2 doesn't work in JAGS
  
  alpha_v5 ~ dnorm(mu_alpha_v5, tau_alpha_v5) # intercept
  beta_v5 ~ dnorm(mu_beta_v5, tau_beta_v5) # slope
  mu_alpha_v5~ dunif(-1, 1) 
  mu_beta_v5~ dunif(-1, 1) 
  sigma_alpha_v5 ~ dunif(0, 1) # standard deviation
  tau_alpha_v5 <- 1 / (sigma_alpha_v5 * sigma_alpha_v5) # sigma^2 doesn't work in JAGS
  sigma_beta_v5~ dunif(0, 1) # standard deviation
  tau_beta_v5<- 1 / (sigma_beta_v5* sigma_beta_v5) # sigma^2 doesn't work in JAGS
  sigma_v5 ~ dunif(0, 0.02) # standard deviation
  tau_v5 <- 1 / (sigma_v5 * sigma_v5) # sigma^2 doesn't work in JAGS
  
  alpha_y1 ~ dnorm(mu_alpha_y1, tau_alpha_y1) # intercept
  beta_y1 ~ dnorm(mu_beta_y1, tau_beta_y1) # slope
  mu_alpha_y1~ dunif(-1,1) 
  mu_beta_y1~ dunif(-1,1) 
  sigma_alpha_y1 ~ dunif(0, 1) # standard deviation
  tau_alpha_y1 <- 1 / (sigma_alpha_y1 * sigma_alpha_y1) # sigma^2 doesn't work in JAGS
  sigma_beta_y1~ dunif(0, 1) # standard deviation
  tau_beta_y1<- 1 / (sigma_beta_y1 * sigma_beta_y1) # sigma^2 doesn't work in JAGS
  sigma_y1 ~ dunif(0, 0.02) # standard deviation
  tau_y1 <- 1 / (sigma_y1 * sigma_y1) # sigma^2 doesn't work in JAGS
  
  alpha_y2 ~ dnorm(mu_alpha_y2, tau_alpha_y2) # intercept
  beta_y2 ~ dnorm(mu_beta_y2, tau_beta_y2) # slope
  mu_alpha_y2~ dunif(-1,1) 
  mu_beta_y2~ dunif(-1,1) 
  sigma_alpha_y2 ~ dunif(0, 1) # standard deviation
  tau_alpha_y2 <- 1 / (sigma_alpha_y2 * sigma_alpha_y2) # sigma^2 doesn't work in JAGS
  sigma_beta_y2~ dunif(0, 1) # standard deviation
  tau_beta_y2<- 1 / (sigma_beta_y2 * sigma_beta_y2) # sigma^2 doesn't work in JAGS
  sigma_y2 ~ dunif(0, 0.02) # standard deviation
  tau_y2 <- 1 / (sigma_y2* sigma_y2) # sigma^2 doesn't work in JAGS
  
  alpha_y3 ~ dnorm(mu_alpha_y3, tau_alpha_y3) # intercept
  beta_y3 ~ dnorm(mu_beta_y3, tau_beta_y3) # slope
  mu_alpha_y3~ dunif(-1,1) 
  mu_beta_y3~ dunif(-1,1) 
  sigma_alpha_y3 ~ dunif(0, 1) # standard deviation
  tau_alpha_y3 <- 1 / (sigma_alpha_y3 * sigma_alpha_y3) # sigma^2 doesn't work in JAGS
  sigma_beta_y3~ dunif(0, 1) # standard deviation
  tau_beta_y3<- 1 / (sigma_beta_y3 * sigma_beta_y3) # sigma^2 doesn't work in JAGS
  sigma_y3 ~ dunif(0, 0.02) # standard deviation
  tau_y3 <- 1 / (sigma_y3 * sigma_y3) # sigma^2 doesn't work in JAGS
  
  alpha_y1_ref ~ dnorm(mu_alpha_y1_ref, tau_alpha_y1_ref) # intercept
  beta1_y1_ref ~ dnorm(mu_beta1_y1_ref, tau_beta1_y1_ref) # slope
  beta2_y1_ref ~ dnorm(mu_beta2_y1_ref, tau_beta2_y1_ref) # slope
  mu_alpha_y1_ref~ dunif(-1,1) 
  mu_beta1_y1_ref~ dunif(-1,1) 
  mu_beta2_y1_ref~ dunif(-1,1) 
  sigma_alpha_y1_ref ~ dunif(0, 1) # standard deviation
  tau_alpha_y1_ref <- 1 / (sigma_alpha_y1_ref * sigma_alpha_y1_ref) # sigma^2 doesn't work in JAGS
  sigma_beta1_y1_ref~ dunif(0, 1) # standard deviation
  sigma_beta2_y1_ref~ dunif(0, 1) # standard deviation
  tau_beta1_y1_ref<- 1 / (sigma_beta1_y1_ref * sigma_beta1_y1_ref) # sigma^2 doesn't work in JAGS
  tau_beta2_y1_ref<- 1 / (sigma_beta2_y1_ref * sigma_beta2_y1_ref) # sigma^2 doesn't work in JAGS
  sigma_y1_ref ~ dunif(0, 0.02) # standard deviation
  tau_y1_ref <- 1 / (sigma_y1_ref * sigma_y1_ref) # sigma^2 doesn't work in JAGS
  
  alpha_y2_ref ~ dnorm(mu_alpha_y2_ref, tau_alpha_y2_ref) # intercept
  beta1_y2_ref ~ dnorm(mu_beta1_y2_ref, tau_beta1_y2_ref) # slope
  beta2_y2_ref ~ dnorm(mu_beta2_y2_ref, tau_beta2_y2_ref) # slope
  mu_alpha_y2_ref~ dunif(-1,1) 
  mu_beta1_y2_ref~ dunif(-1,1) 
  mu_beta2_y2_ref~ dunif(-1,1) 
  sigma_alpha_y2_ref ~ dunif(0, 1) # standard deviation
  tau_alpha_y2_ref <- 1 / (sigma_alpha_y2_ref * sigma_alpha_y2_ref) # sigma^2 doesn't work in JAGS
  sigma_beta1_y2_ref~ dunif(0, 1) # standard deviation
  sigma_beta2_y2_ref~ dunif(0, 1) # standard deviation
  tau_beta1_y2_ref<- 1 / (sigma_beta1_y2_ref * sigma_beta1_y2_ref) # sigma^2 doesn't work in JAGS
  tau_beta2_y2_ref<- 1 / (sigma_beta2_y2_ref * sigma_beta2_y2_ref) # sigma^2 doesn't work in JAGS
  sigma_y2_ref ~ dunif(0, 0.02) # standard deviation
  tau_y2_ref <- 1 / (sigma_y2_ref * sigma_y2_ref) # sigma^2 doesn't work in JAGS
  
  alpha_y3_ref ~ dnorm(mu_alpha_y3_ref, tau_alpha_y3_ref) # intercept
  beta1_y3_ref ~ dnorm(mu_beta1_y3_ref, tau_beta1_y3_ref) # slope
  beta2_y3_ref ~ dnorm(mu_beta2_y3_ref, tau_beta2_y3_ref) # slope
  mu_alpha_y3_ref~ dunif(-1,1) 
  mu_beta1_y3_ref~ dunif(-1,1) 
  mu_beta2_y3_ref~ dunif(-1,1) 
  sigma_alpha_y3_ref ~ dunif(0, 1) # standard deviation
  tau_alpha_y3_ref <- 1 / (sigma_alpha_y3_ref * sigma_alpha_y3_ref) # sigma^2 doesn't work in JAGS
  sigma_beta1_y3_ref~ dunif(0, 1) # standard deviation
  sigma_beta2_y3_ref~ dunif(0, 1) # standard deviation
  tau_beta1_y3_ref<- 1 / (sigma_beta1_y3_ref * sigma_beta1_y3_ref) # sigma^2 doesn't work in JAGS
  tau_beta2_y3_ref<- 1 / (sigma_beta2_y3_ref * sigma_beta2_y3_ref) # sigma^2 doesn't work in JAGS
  sigma_y3_ref ~ dunif(0, 0.02) # standard deviation
  tau_y3_ref <- 1 / (sigma_y3_ref * sigma_y3_ref) # sigma^2 doesn't work in JAGS
}

params <- c("alpha_fh1", "beta_fh1", "sigma_fh1","alpha_fh2", "beta_fh2", "sigma_fh2","alpha_fh3", "beta_fh3", "sigma_fh3","alpha_fh4", "beta_fh4", "sigma_fh4",
            "alpha_fh5", "beta_fh5", "sigma_fh5","alpha_fh6", "beta_fh6", "sigma_fh6",
            "alpha_fv1", "beta_fv1", "sigma_fv1","alpha_fv2", "beta_fv2", "sigma_fv2","alpha_fv3", "beta_fv3", "sigma_fv3","alpha_fv4", "beta_fv4", "sigma_fv4",
            "alpha_fv5", "beta_fv5", "sigma_fv5",
            "alpha_fg1", "beta_fg1", "sigma_fg1",
            "alpha_s1", "beta_s1", "sigma_s1",
            "alpha_s2", "beta_s2", "sigma_s2",
            "alpha_s3", "beta_s3", "sigma_s3",
            "alpha_s4", "beta_s4", "sigma_s4",
            "alpha_s5", "beta_s5", "sigma_s5",
            "alpha_s6", "beta_s6", "sigma_s6",
            "alpha_v1", "beta_v1", "sigma_v1",
            "alpha_v2", "beta_v2", "sigma_v2",
            "alpha_v3", "beta_v3", "sigma_v3",
            "alpha_v4", "beta_v4", "sigma_v4",
            "alpha_v5", "beta_v5", "sigma_v5",
            "alpha_y1", "beta_y1", "sigma_y1",
            "alpha_y2", "beta_y2", "sigma_y2",
            "alpha_y3", "beta_y3", "sigma_y3",
            "alpha_y1_ref", "beta1_y1_ref", "beta2_y1_ref","sigma_y1_ref",
            "alpha_y2_ref", "beta1_y2_ref", "beta2_y2_ref","sigma_y2_ref",
            "alpha_y3_ref", "beta1_y3_ref", "beta2_y3_ref","sigma_y3_ref")
fit_lm1 <- jags(data = jagsdata_s1,  parameters.to.save = params,model.file = lm1_jags,
                n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC = F)

lm1_mcmc <- as.mcmc(fit_lm1)
lm1_mcmc_combi <- as.mcmc(rbind(lm1_mcmc[[1]], lm1_mcmc[[2]], lm1_mcmc[[3]],lm1_mcmc[[4]],lm1_mcmc[[5]],lm1_mcmc[[6]],lm1_mcmc[[7]]))
lm1_mcmc_combi_rep <- do.call(rbind, rep(list(lm1_mcmc_combi), 1)) # replication

lm1_mcmc_ref <- do.call(cbind, rep(list("alpha_fh1"=lm1_mcmc_combi[,"alpha_fh1"],"alpha_fh2"=lm1_mcmc_combi[,"alpha_fh2"],"alpha_fh3"=lm1_mcmc_combi[,"alpha_fh3"],"alpha_fh4"=lm1_mcmc_combi[,"alpha_fh4"],"alpha_fh5"=lm1_mcmc_combi[,"alpha_fh5"],"alpha_fh6"=lm1_mcmc_combi[,"alpha_fh6"],
                                        "beta_fh1"=lm1_mcmc_combi[,"beta_fh1"],"beta_fh2"=lm1_mcmc_combi[,"beta_fh2"],"beta_fh3"=lm1_mcmc_combi[,"beta_fh3"],"beta_fh4"=lm1_mcmc_combi[,"beta_fh4"], "beta_fh5"=lm1_mcmc_combi[,"beta_fh5"],"beta_fh6"=lm1_mcmc_combi[,"beta_fh6"], 
                                        "sigma_fh1"=lm1_mcmc_combi[,"sigma_fh1"],"sigma_fh2"=lm1_mcmc_combi[,"sigma_fh2"],"sigma_fh3"=lm1_mcmc_combi[,"sigma_fh3"],"sigma_fh4"=lm1_mcmc_combi[,"sigma_fh4"], "sigma_fh5"=lm1_mcmc_combi[,"sigma_fh5"],"sigma_fh6"=lm1_mcmc_combi[,"sigma_fh6"],
                                        "alpha_fv1"=lm1_mcmc_combi[,"alpha_fv1"],"alpha_fv2"=lm1_mcmc_combi[,"alpha_fv2"],"alpha_fv3"=lm1_mcmc_combi[,"alpha_fv3"],"alpha_fv4"=lm1_mcmc_combi[,"alpha_fv4"],"alpha_fv5"=lm1_mcmc_combi[,"alpha_fv5"],
                                        "beta_fv1"=lm1_mcmc_combi[,"beta_fv1"],"beta_fv2"=lm1_mcmc_combi[,"beta_fv2"],"beta_fv3"=lm1_mcmc_combi[,"beta_fv3"],"beta_fv4"=lm1_mcmc_combi[,"beta_fv4"], "beta_fv5"=lm1_mcmc_combi[,"beta_fv5"],
                                        "sigma_fv1"=lm1_mcmc_combi[,"sigma_fv1"],"sigma_fv2"=lm1_mcmc_combi[,"sigma_fv2"],"sigma_fv3"=lm1_mcmc_combi[,"sigma_fv3"],"sigma_fv4"=lm1_mcmc_combi[,"sigma_fv4"], "sigma_fv5"=lm1_mcmc_combi[,"sigma_fv5"],
                                        "alpha_fg1"=lm1_mcmc_combi[,"alpha_fg1"],"beta_fg1"=lm1_mcmc_combi[,"beta_fg1"],"sigma_fg1"=lm1_mcmc_combi[,"sigma_fg1"],
                                        "alpha_s1"=lm1_mcmc_combi[,"alpha_s1"],"alpha_s2"=lm1_mcmc_combi[,"alpha_s2"],"alpha_s3"=lm1_mcmc_combi[,"alpha_s3"],"alpha_s4"=lm1_mcmc_combi[,"alpha_s4"],
                                        "alpha_s5"=lm1_mcmc_combi[,"alpha_s5"],"alpha_s6"=lm1_mcmc_combi[,"alpha_s6"],
                                        "beta_s1"=lm1_mcmc_combi[,"beta_s1"],"beta_s2"=lm1_mcmc_combi[,"beta_s2"],"beta_s3"=lm1_mcmc_combi[,"beta_s3"],"beta_s4"=lm1_mcmc_combi[,"beta_s4"],
                                        "beta_s5"=lm1_mcmc_combi[,"beta_s5"],"beta_s6"=lm1_mcmc_combi[,"beta_s6"],
                                        "sigma_s1"=lm1_mcmc_combi[,"sigma_s1"],"sigma_s2"=lm1_mcmc_combi[,"sigma_s2"],"sigma_s3"=lm1_mcmc_combi[,"sigma_s3"],"sigma_s4"=lm1_mcmc_combi[,"sigma_s4"],
                                        "sigma_s5"=lm1_mcmc_combi[,"sigma_s5"],"sigma_s6"=lm1_mcmc_combi[,"sigma_s6"],
                                        "alpha_v1"=lm1_mcmc_combi[,"alpha_v1"],"alpha_v2"=lm1_mcmc_combi[,"alpha_v2"],"alpha_v3"=lm1_mcmc_combi[,"alpha_v3"],"alpha_v4"=lm1_mcmc_combi[,"alpha_v4"],
                                        "alpha_v5"=lm1_mcmc_combi[,"alpha_v5"],
                                        "beta_v1"=lm1_mcmc_combi[,"beta_v1"],"beta_v2"=lm1_mcmc_combi[,"beta_v2"],"beta_v3"=lm1_mcmc_combi[,"beta_v3"],"beta_v4"=lm1_mcmc_combi[,"beta_v4"],
                                        "beta_v5"=lm1_mcmc_combi[,"beta_v5"],
                                        "sigma_v1"=lm1_mcmc_combi[,"sigma_v1"],"sigma_v2"=lm1_mcmc_combi[,"sigma_v2"],"sigma_v3"=lm1_mcmc_combi[,"sigma_v3"],"sigma_v4"=lm1_mcmc_combi[,"sigma_v4"],
                                        "sigma_v5"=lm1_mcmc_combi[,"sigma_v5"],
                                        "alpha_y1"=lm1_mcmc_combi[,"alpha_y1"],"alpha_y2"=lm1_mcmc_combi[,"alpha_y2"],"alpha_y3"=lm1_mcmc_combi[,"alpha_y3"],
                                        "beta_y1"=lm1_mcmc_combi[,"beta_y1"],"beta_y2"=lm1_mcmc_combi[,"beta_y2"],"beta_y3"=lm1_mcmc_combi[,"beta_y3"],
                                        "sigma_y1"=lm1_mcmc_combi[,"sigma_y1"],"sigma_y2"=lm1_mcmc_combi[,"sigma_y2"],"sigma_y3"=lm1_mcmc_combi[,"sigma_y3"]), 1))
lm1_mcmc_ref1 <- do.call(cbind, rep(list("alpha_y1_ref"=lm1_mcmc_combi[,"alpha_y1_ref"],"alpha_y2_ref"=lm1_mcmc_combi[,"alpha_y2_ref"],"alpha_y3_ref"=lm1_mcmc_combi[,"alpha_y3_ref"],
                                         "beta1_y1_ref"=lm1_mcmc_combi[,"beta1_y1_ref"],"beta1_y2_ref"=lm1_mcmc_combi[,"beta1_y2_ref"],"beta1_y3_ref"=lm1_mcmc_combi[,"beta1_y3_ref"],
                                         "beta2_y1_ref"=lm1_mcmc_combi[,"beta2_y1_ref"],"beta2_y2_ref"=lm1_mcmc_combi[,"beta2_y2_ref"],"beta2_y3_ref"=lm1_mcmc_combi[,"beta2_y3_ref"],
                                         "sigma_y1_ref"=lm1_mcmc_combi[,"sigma_y1_ref"],"sigma_y2_ref"=lm1_mcmc_combi[,"sigma_y2_ref"],"sigma_y3_ref"=lm1_mcmc_combi[,"sigma_y3_ref"]), 1)) 

#Parameters for emission intensity and consumption patterns
tau_gdp_pre1<- matrix(NA, nrow =1,ncol=n_region)
for (j in 1:n_region){
  tau_gdp_pre1[,j]=1/(sd(lm1_mcmc_gdp[,j])^2)
}

tau_gdp_pre2<- matrix(NA, nrow =1,ncol=n_region)
tau_gdp_pre3<- matrix(NA, nrow =1,ncol=n_region)
for (j in 1:n_region){
  tau_gdp_pre2[,j]=(mean(lm1_mcmc_gdp[,j+40]))
  tau_gdp_pre3[,j]=1/(sd(lm1_mcmc_gdp[,j+40])^2)
}

library(R2jags)

jagsdata_s1 <- with(data_gdp, list(x=x,x_2023=x_2023,d=d,tau_gdp_pre1=tau_gdp_pre1,tau_gdp_pre2=tau_gdp_pre2,tau_gdp_pre3=tau_gdp_pre3,
                                   eih=eih,eiv=eiv,eig=eig,
                                   s=s,v=v,y=yo,
                                   eih1_j=eih,eih2_j=eih,eih3_j=eih,eih4_j=eih,eih5_j=eih,eih6_j=eih,
                                   eiv1_j=eiv,eiv2_j=eiv,eiv3_j=eiv,eiv4_j=eiv,eiv5_j=eiv,
                                   eig1_j=eig,
                                   s1_j=s,s2_j=s,s3_j=s,s4_j=s,s5_j=s,s6_j=s,
                                   v1_j=v,v2_j=v,v3_j=v,v4_j=v,v5_j=v,
                                   y1_j=yo,y2_j=yo,y3_j=yo,y_ref=yo,
                                   n_2015=21,n_region=40))

lm1_jags <- function(){
  # Likelihood:
  for (j in 1:n_region){
    for (t in 1:24){x_new[t,j]<-(x[t,j])}
    for (t in 25:25){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error1[j]}
    for (t in 26:26){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error2[j]}
    for (t in 27:27){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error3[j]}
    for (t in 28:28){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error4[j]}
    for (t in 29:29){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error5[j]}
    for (t in 30:30){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error6[j]}
    for (t in 31:31){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error7[j]}
    for (t in 32:32){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error8[j]}
    for (t in 33:33){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error9[j]}
    for (t in 34:34){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error10[j]}
    for (t in 35:35){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error11[j]}
    for (t in 36:36){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error12[j]}
    for (t in 37:37){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error13[j]}
    for (t in 38:38){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error14[j]}
    for (t in 39:39){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error15[j]}
    for (t in 40:40){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error16[j]}
    for (t in 41:41){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error17[j]}
    for (t in 42:42){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error18[j]}
    for (t in 43:43){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error19[j]}
    for (t in 44:44){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error20[j]}
    for (t in 45:45){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error21[j]}
    for (t in 46:46){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error22[j]}
    for (t in 47:47){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error23[j]}
    for (t in 48:48){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error24[j]}
    for (t in 49:49){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error25[j]}
    for (t in 50:50){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error26[j]}
    for (t in 51:51){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error27[j]}
    for (t in 52:52){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error28[j]}
    for (t in 53:53){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error29[j]}
    for (t in 54:54){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error30[j]}
    for (t in 55:55){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error31[j]}
    for (t in 56:56){x_new[t,j] <- x_new[t-1,j]+d[t-1,j]+alpha_gdp_pre[j]+error32[j]}
  }
  
  for (t in 1:n_2015){
    eih[t,33] ~ dnorm(mu_fh1[t], tau_fh1) # tau is precision (1 / variance)
    mu_fh1[t] <- (alpha_fh1 + beta_fh1* (x[t,33]))
    eih[t+21,34] ~ dnorm(mu_fh2[t], tau_fh2) # tau is precision (1 / variance)
    mu_fh2[t] <- (alpha_fh2 + beta_fh2* (x[t,34]))
    eih[t+21*2,33] ~ dnorm(mu_fh3[t], tau_fh3) # tau is precision (1 / variance)
    mu_fh3[t] <- (alpha_fh3 + beta_fh3* (x[t,33]))
    eih[t+21*3,20] ~ dnorm(mu_fh4[t], tau_fh4) # tau is precision (1 / variance)
    mu_fh4[t] <- (alpha_fh4 + beta_fh4* (x[t,20]))
    eih[t+21*4,34] ~ dnorm(mu_fh5[t], tau_fh5) # tau is precision (1 / variance)
    mu_fh5[t] <- (alpha_fh5 + beta_fh5* (x[t,34]))
    eih[t+21*5,33] ~ dnorm(mu_fh6[t], tau_fh6) # tau is precision (1 / variance)
    mu_fh6[t] <- (alpha_fh6 + beta_fh6* (x[t,33]))
    
    eiv[t,35] ~ dnorm(mu_fv1[t], tau_fv1) # tau is precision (1 / variance)
    mu_fv1[t] <- (alpha_fv1 + beta_fv1* (x[t,35]))
    eiv[t+21,33] ~ dnorm(mu_fv2[t], tau_fv2) # tau is precision (1 / variance)
    mu_fv2[t] <- (alpha_fv2 + beta_fv2* (x[t,33]))
    eiv[t+21*2,33] ~ dnorm(mu_fv3[t], tau_fv3) # tau is precision (1 / variance)
    mu_fv3[t] <- (alpha_fv3 + beta_fv3* (x[t,33]))
    eiv[t+21*3,33] ~ dnorm(mu_fv4[t], tau_fv4) # tau is precision (1 / variance)
    mu_fv4[t] <- (alpha_fv4 + beta_fv4* (x[t,33]))
    eiv[t+21*4,20] ~ dnorm(mu_fv5[t], tau_fv5) # tau is precision (1 / variance)
    mu_fv5[t] <- (alpha_fv5 + beta_fv5* (x[t,20]))
    
    eig[t,20] ~ dnorm(mu_fg1[t], tau_fg1) # tau is precision (1 / variance)
    mu_fg1[t] <- (alpha_fg1 + beta_fg1* (x[t,20]))
  }
  
  for (j in 1:n_region){
    for (t in 1:n_2015){
      s[t,j] ~ dnorm(mu_s1[t,j], tau_s1) # tau is precision (1 / variance)
      mu_s1[t,j] <- (alpha_s1 + beta_s1* (x[t,j]))
      s[t+21,j] ~ dnorm(mu_s2[t,j], tau_s2) # tau is precision (1 / variance)
      mu_s2[t,j] <- (alpha_s2 + beta_s2* (x[t,j]))
      s[t+21*2,j] ~ dnorm(mu_s3[t,j], tau_s3) # tau is precision (1 / variance)
      mu_s3[t,j] <- (alpha_s3 + beta_s3* (x[t,j]))
      s[t+21*3,j] ~ dnorm(mu_s4[t,j], tau_s4) # tau is precision (1 / variance)
      mu_s4[t,j] <- (alpha_s4 + beta_s4* (x[t,j]))
      s[t+21*4,j] ~ dnorm(mu_s5[t,j], tau_s5) # tau is precision (1 / variance)
      mu_s5[t,j] <- (alpha_s5 + beta_s5* (x[t,j]))
      s[t+21*5,j] ~ dnorm(mu_s6[t,j], tau_s6) # tau is precision (1 / variance)
      mu_s6[t,j] <- (alpha_s6 + beta_s6* (x[t,j]))
      
      v[t,j] ~ dnorm(mu_v1[t,j], tau_v1) # tau is precision (1 / variance)
      mu_v1[t,j] <- (alpha_v1 + beta_v1* (x[t,j]))
      v[t+21,j] ~ dnorm(mu_v2[t,j], tau_v2) # tau is precision (1 / variance)
      mu_v2[t,j] <- (alpha_v2 + beta_v2* (x[t,j]))
      v[t+21*2,j] ~ dnorm(mu_v3[t,j], tau_v3) # tau is precision (1 / variance)
      mu_v3[t,j] <- (alpha_v3 + beta_v3* (x[t,j]))
      v[t+21*3,j] ~ dnorm(mu_v4[t,j], tau_v4) # tau is precision (1 / variance)
      mu_v4[t,j] <- (alpha_v4 + beta_v4* (x[t,j]))
      v[t+21*4,j] ~ dnorm(mu_v5[t,j], tau_v5) # tau is precision (1 / variance)
      mu_v5[t,j] <- (alpha_v5 + beta_v5* (x[t,j]))
    }
  }
  
  for (j in 1:36){
    for (t in 1:29){
      y[t,j] ~ dnorm(mu_y1[t,j], tau_y1) # tau is precision (1 / variance)
      mu_y1[t,j] <- (alpha_y1 + beta_y1* (x_new[t,j]))
      y[t+29,j] ~ dnorm(mu_y2[t,j], tau_y2) # tau is precision (1 / variance)
      mu_y2[t,j] <- (alpha_y2 + beta_y2* (x_new[t,j]))
      y[t+29*2,j] ~ dnorm(mu_y3[t,j], tau_y3) # tau is precision (1 / variance)
      mu_y3[t,j] <- (alpha_y3 + beta_y3* (x_new[t,j]))
    }
  }
  
  for (t in 1:n_2015){ 
    eihr1[t] <- (alpha_fh1 + beta_fh1 *(x[t,33]))
    eihr2[t] <- (alpha_fh2 + beta_fh2 *(x[t,34]))
    eihr3[t] <- (alpha_fh3 + beta_fh3 *(x[t,33]))
    eihr4[t] <- (alpha_fh4 + beta_fh4 *(x[t,20]))
    eihr5[t] <- (alpha_fh5 + beta_fh5 *(x[t,34]))
    eihr6[t] <- (alpha_fh6 + beta_fh6 *(x[t,33]))
    
    eivr1[t] <- (alpha_fv1 + beta_fv1 *(x[t,35]))
    eivr2[t] <- (alpha_fv2 + beta_fv2 *(x[t,33]))
    eivr3[t] <- (alpha_fv3 + beta_fv3 *(x[t,33]))
    eivr4[t] <- (alpha_fv4 + beta_fv4 *(x[t,33]))
    eivr5[t] <- (alpha_fv5 + beta_fv5 *(x[t,20]))
    
    eigr1[t] <- (alpha_fg1 + beta_fg1 *(x[t,20]))
  }
  
  for (j in 1:n_region){
    for (t in 1:n_2015){  
      sr1[t,j] <- (alpha_s1 + beta_s1 *(x[t,j]))
      sr2[t,j] <- (alpha_s2 + beta_s2 *(x[t,j]))
      sr3[t,j] <- (alpha_s3 + beta_s3 *(x[t,j]))
      sr4[t,j] <- (alpha_s4 + beta_s4 *(x[t,j]))
      sr5[t,j] <- (alpha_s5 + beta_s5 *(x[t,j]))
      sr6[t,j] <- (alpha_s6 + beta_s6 *(x[t,j]))
      
      vr1[t,j] <- (alpha_v1 + beta_v1 *(x[t,j]))
      vr2[t,j] <- (alpha_v2 + beta_v2 *(x[t,j]))
      vr3[t,j] <- (alpha_v3 + beta_v3 *(x[t,j]))
      vr4[t,j] <- (alpha_v4 + beta_v4 *(x[t,j]))
      vr5[t,j] <- (alpha_v5 + beta_v5 *(x[t,j]))
    }
  }
  
  for (j in 1:37){
    for (t in 1:56){
      yr1[t,j] <- (alpha_y1 + beta_y1 *(x_new[t,j]))
      yr2[t,j] <- (alpha_y2 + beta_y2 *(x_new[t,j]))
      yr3[t,j] <- (alpha_y3 + beta_y3 *(x_new[t,j]))
    }
  }
  
  
  for (j in 1:n_region){
    for (t in 2:n_2015){
      eih1_j[t,j] ~ dnorm(mu_j_fh1[t,j], tau_j_fh1[j]) # tau is precision (1 / variance)
      mu_j_fh1[t,j] <- eihr1[t]+ar_fh1[j]*(eih1_j[t-1,j]-eihr1[t-1])
      eih2_j[t+21,j] ~ dnorm(mu_j_fh2[t,j], tau_j_fh2[j]) # tau is precision (1 / variance)
      mu_j_fh2[t,j] <- eihr2[t]+ar_fh2[j]*(eih2_j[t-1+21,j]-eihr2[t-1])
      eih3_j[t+21*2,j] ~ dnorm(mu_j_fh3[t,j], tau_j_fh3[j]) # tau is precision (1 / variance)
      mu_j_fh3[t,j] <- eihr3[t]+ar_fh3[j]*(eih3_j[t-1+21*2,j]-eihr3[t-1])
      eih4_j[t+21*3,j] ~ dnorm(mu_j_fh4[t,j], tau_j_fh4[j]) # tau is precision (1 / variance)
      mu_j_fh4[t,j] <- eihr4[t]+ar_fh4[j]*(eih4_j[t-1+21*3,j]-eihr4[t-1])
      eih5_j[t+21*4,j] ~ dnorm(mu_j_fh5[t,j], tau_j_fh5[j]) # tau is precision (1 / variance)
      mu_j_fh5[t,j] <- eihr5[t]+ar_fh5[j]*(eih5_j[t-1+21*4,j]-eihr5[t-1])
      eih6_j[t+21*5,j] ~ dnorm(mu_j_fh6[t,j], tau_j_fh6[j]) # tau is precision (1 / variance)
      mu_j_fh6[t,j] <- eihr6[t]+ar_fh6[j]*(eih6_j[t-1+21*5,j]-eihr6[t-1])
      
      eiv1_j[t,j] ~ dnorm(mu_j_fv1[t,j], tau_j_fv1[j]) # tau is precision (1 / variance)
      mu_j_fv1[t,j] <- eivr1[t]+ar_fv1[j]*(eiv1_j[t-1,j]-eivr1[t-1])
      eiv2_j[t+21,j] ~ dnorm(mu_j_fv2[t,j], tau_j_fv2[j]) # tau is precision (1 / variance)
      mu_j_fv2[t,j] <- eivr2[t]+ar_fv2[j]*(eiv2_j[t-1+21,j]-eivr2[t-1])
      eiv3_j[t+21*2,j] ~ dnorm(mu_j_fv3[t,j], tau_j_fv3[j]) # tau is precision (1 / variance)
      mu_j_fv3[t,j] <- eivr3[t]+ar_fv3[j]*(eiv3_j[t-1+21*2,j]-eivr3[t-1])
      eiv4_j[t+21*3,j] ~ dnorm(mu_j_fv4[t,j], tau_j_fv4[j]) # tau is precision (1 / variance)
      mu_j_fv4[t,j] <- eivr4[t]+ar_fv4[j]*(eiv4_j[t-1+21*3,j]-eivr4[t-1])
      eiv5_j[t+21*4,j] ~ dnorm(mu_j_fv5[t,j], tau_j_fv5[j]) # tau is precision (1 / variance)
      mu_j_fv5[t,j] <- eivr5[t]+ar_fv5[j]*(eiv5_j[t-1+21*4,j]-eivr5[t-1])
      
      eig1_j[t,j] ~ dnorm(mu_j_fg1[t,j], tau_j_fg1[j]) # tau is precision (1 / variance)
      mu_j_fg1[t,j] <- eigr1[t]+ar_fg1[j]*(eig1_j[t-1,j]-eigr1[t-1])
      
      s1_j[t,j] ~ dnorm(mu_j_s1[t,j], tau_j_s1[j]) # tau is precision (1 / variance)
      mu_j_s1[t,j] <- sr1[t,j]+ar_s1[j]*(s1_j[t-1,j]-sr1[t-1,j])
      s2_j[t+21,j] ~ dnorm(mu_j_s2[t,j], tau_j_s2[j]) # tau is precision (1 / variance)
      mu_j_s2[t,j] <- sr2[t,j]+ar_s2[j]*(s2_j[t-1+21,j]-sr2[t-1,j])
      s3_j[t+21*2,j] ~ dnorm(mu_j_s3[t,j], tau_j_s3[j]) # tau is precision (1 / variance)
      mu_j_s3[t,j] <- sr3[t,j]+ar_s3[j]*(s3_j[t-1+21*2,j]-sr3[t-1,j])
      s4_j[t+21*3,j] ~ dnorm(mu_j_s4[t,j], tau_j_s4[j]) # tau is precision (1 / variance)
      mu_j_s4[t,j] <- sr4[t,j]+ar_s4[j]*(s4_j[t-1+21*3,j]-sr4[t-1,j])
      s5_j[t+21*4,j] ~ dnorm(mu_j_s5[t,j], tau_j_s5[j]) # tau is precision (1 / variance)
      mu_j_s5[t,j] <- sr5[t,j]+ar_s5[j]*(s5_j[t-1+21*4,j]-sr5[t-1,j])
      s6_j[t+21*5,j] ~ dnorm(mu_j_s6[t,j], tau_j_s6[j]) # tau is precision (1 / variance)
      mu_j_s6[t,j] <- sr6[t,j]+ar_s6[j]*(s6_j[t-1+21*5,j]-sr6[t-1,j])
      
      v1_j[t,j] ~ dnorm(mu_j_v1[t,j], tau_j_v1[j]) # tau is precision (1 / variance)
      mu_j_v1[t,j] <- vr1[t,j]+ar_v1[j]*(v1_j[t-1,j]-vr1[t-1,j])
      v2_j[t+21,j] ~ dnorm(mu_j_v2[t,j], tau_j_v2[j]) # tau is precision (1 / variance)
      mu_j_v2[t,j] <- vr2[t,j]+ar_v2[j]*(v2_j[t-1+21,j]-vr2[t-1,j])
      v3_j[t+21*2,j] ~ dnorm(mu_j_v3[t,j], tau_j_v3[j]) # tau is precision (1 / variance)
      mu_j_v3[t,j] <- vr3[t,j]+ar_v3[j]*(v3_j[t-1+21*2,j]-vr3[t-1,j])
      v4_j[t+21*3,j] ~ dnorm(mu_j_v4[t,j], tau_j_v4[j]) # tau is precision (1 / variance)
      mu_j_v4[t,j] <- vr4[t,j]+ar_v4[j]*(v4_j[t-1+21*3,j]-vr4[t-1,j])
      v5_j[t+21*4,j] ~ dnorm(mu_j_v5[t,j], tau_j_v5[j]) # tau is precision (1 / variance)
      mu_j_v5[t,j] <- vr5[t,j]+ar_v5[j]*(v5_j[t-1+21*4,j]-vr5[t-1,j])
    }
  }
  
  for (j in 1:36){
    for (t in 2:29){     
      y1_j[t,j] ~ dnorm(mu_j_y1[t,j], tau_j_y1[j]) # tau is precision (1 / variance)
      mu_j_y1[t,j] <- yr1[t,j]+ar_y1[j]*(y1_j[t-1,j]-yr1[t-1,j])
      y2_j[t+29,j] ~ dnorm(mu_j_y2[t,j], tau_j_y2[j]) # tau is precision (1 / variance)
      mu_j_y2[t,j] <- yr2[t,j]+ar_y2[j]*(y2_j[t-1+29,j]-yr2[t-1,j])
      y3_j[t+29*2,j] ~ dnorm(mu_j_y3[t,j], tau_j_y3[j]) # tau is precision (1 / variance)
      mu_j_y3[t,j] <- yr3[t,j]+ar_y3[j]*(y3_j[t-1+29*2,j]-yr3[t-1,j])
    }
  }
  
  for (j in 37:37){
    for (t in 17:29){     
      y1_j[t,j] ~ dnorm(mu_j_y1[t,j], tau_j_y1[j]) # tau is precision (1 / variance)
      mu_j_y1[t,j] <- yr1[t,j]+ar_y1[j]*(y1_j[t-1,j]-yr1[t-1,j])
      y2_j[t+29,j] ~ dnorm(mu_j_y2[t,j], tau_j_y2[j]) # tau is precision (1 / variance)
      mu_j_y2[t,j] <- yr2[t,j]+ar_y2[j]*(y2_j[t-1+29,j]-yr2[t-1,j])
      y3_j[t+29*2,j] ~ dnorm(mu_j_y3[t,j], tau_j_y3[j]) # tau is precision (1 / variance)
      mu_j_y3[t,j] <- yr3[t,j]+ar_y3[j]*(y3_j[t-1+29*2,j]-yr3[t-1,j])
    }
  }
  
  for (t in 1:29){
    y_ref[t,37] ~ dnorm(mu_y1_ref[t], tau_y1_ref) # tau is precision (1 / variance)
    mu_y1_ref[t] <- (alpha_y1_ref + beta2_y1_ref* (x_2023[t,37]* x_2023[t,37])+ beta1_y1_ref* (x_2023[t,37]))
    y_ref[t+29,37] ~ dnorm(mu_y2_ref[t], tau_y2_ref) # tau is precision (1 / variance)
    mu_y2_ref[t] <- (alpha_y2_ref + beta2_y2_ref* (x_2023[t,37]* x_2023[t,37])+ beta1_y2_ref* (x_2023[t,37]))
    y_ref[t+29*2,37] ~ dnorm(mu_y3_ref[t], tau_y3_ref) # tau is precision (1 / variance)
    mu_y3_ref[t] <- (alpha_y3_ref + beta2_y3_ref* (x_2023[t,37]* x_2023[t,37])+ beta1_y3_ref* (x_2023[t,37]))
  }
  
  for (j in 38:40){
    for (t in 1:29){
      yref1[t,j] <- (alpha_y1_ref + beta2_y1_ref* (x_2023[t,j]* x_2023[t,j])+ beta1_y1_ref* (x_2023[t,j]))
      yref2[t,j] <- (alpha_y2_ref + beta2_y2_ref* (x_2023[t,j]* x_2023[t,j])+ beta1_y2_ref* (x_2023[t,j]))
      yref3[t,j] <- (alpha_y3_ref + beta2_y3_ref* (x_2023[t,j]* x_2023[t,j])+ beta1_y3_ref* (x_2023[t,j]))
    }
  }
  
  for (j in 38:40){
    for (t in 2:29){     
      y1_j[t,j] ~ dnorm(mu_j_y1[t,j], tau_j_y1[j]) # tau is precision (1 / variance)
      mu_j_y1[t,j] <- yref1[t,j]+ar_y1[j]*(y1_j[t-1,j]-yref1[t-1,j])
      y2_j[t+29,j] ~ dnorm(mu_j_y2[t,j], tau_j_y2[j]) # tau is precision (1 / variance)
      mu_j_y2[t,j] <- yref2[t,j]+ar_y2[j]*(y2_j[t-1+29,j]-yref1[t-1,j])
      y3_j[t+29*2,j] ~ dnorm(mu_j_y3[t,j], tau_j_y3[j]) # tau is precision (1 / variance)
      mu_j_y3[t,j] <- yref3[t,j]+ar_y3[j]*(y3_j[t-1+29*2,j]-yref1[t-1,j])
    }
  }
  
  # Priors:
  for (j in 1:n_region){ 
    alpha_gdp_pre[j]~ dnorm(0, tau_gdp_pre1[,j]) 
    error1[j]~dnorm(0, tau_error[j]) 
    error2[j]~dnorm(0, tau_error[j]) 
    error3[j]~dnorm(0, tau_error[j]) 
    error4[j]~dnorm(0, tau_error[j]) 
    error5[j]~dnorm(0, tau_error[j]) 
    error6[j]~dnorm(0, tau_error[j]) 
    error7[j]~dnorm(0, tau_error[j]) 
    error8[j]~dnorm(0, tau_error[j]) 
    error9[j]~dnorm(0, tau_error[j]) 
    error10[j]~dnorm(0, tau_error[j]) 
    error11[j]~dnorm(0, tau_error[j]) 
    error12[j]~dnorm(0, tau_error[j]) 
    error13[j]~dnorm(0, tau_error[j]) 
    error14[j]~dnorm(0, tau_error[j]) 
    error15[j]~dnorm(0, tau_error[j]) 
    error16[j]~dnorm(0, tau_error[j]) 
    error17[j]~dnorm(0, tau_error[j]) 
    error18[j]~dnorm(0, tau_error[j]) 
    error19[j]~dnorm(0, tau_error[j]) 
    error20[j]~dnorm(0, tau_error[j]) 
    error21[j]~dnorm(0, tau_error[j]) 
    error22[j]~dnorm(0, tau_error[j]) 
    error23[j]~dnorm(0, tau_error[j]) 
    error24[j]~dnorm(0, tau_error[j]) 
    error25[j]~dnorm(0, tau_error[j]) 
    error26[j]~dnorm(0, tau_error[j]) 
    error27[j]~dnorm(0, tau_error[j]) 
    error28[j]~dnorm(0, tau_error[j]) 
    error29[j]~dnorm(0, tau_error[j]) 
    error30[j]~dnorm(0, tau_error[j]) 
    error31[j]~dnorm(0, tau_error[j]) 
    error32[j]~dnorm(0, tau_error[j]) 
    tau_error[j]<- 1 / (sigma_error[j] * sigma_error[j]) 
    sigma_error[j] ~ dnorm(tau_gdp_pre2[,j], tau_gdp_pre3[,j])
  }
  
  alpha_fh1 ~ dnorm(mu_alpha_fh1, tau_alpha_fh1) # intercept
  beta_fh1 ~ dnorm(mu_beta_fh1, tau_beta_fh1) # slope
  mu_alpha_fh1~ dunif(-30,30) 
  mu_beta_fh1~ dunif(-5, 5) 
  sigma_alpha_fh1 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh1 <- 1 / (sigma_alpha_fh1 * sigma_alpha_fh1) # sigma^2 doesn't work in JAGS
  sigma_beta_fh1~ dunif(0, 1) # standard deviation
  tau_beta_fh1<- 1 / (sigma_beta_fh1 * sigma_beta_fh1) # sigma^2 doesn't work in JAGS
  sigma_fh1 ~ dunif(0, 0.05) # standard deviation
  tau_fh1 <- 1 / (sigma_fh1 * sigma_fh1) # sigma^2 doesn't work in JAGS
  
  alpha_fh2 ~ dnorm(mu_alpha_fh2, tau_alpha_fh2) # intercept
  beta_fh2 ~ dnorm(mu_beta_fh2, tau_beta_fh2) # slope
  mu_alpha_fh2~ dunif(-30,30) 
  mu_beta_fh2~ dunif(-5,5) 
  sigma_alpha_fh2 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh2 <- 1 / (sigma_alpha_fh2 * sigma_alpha_fh2) # sigma^2 doesn't work in JAGS
  sigma_beta_fh2~ dunif(0, 1) # standard deviation
  tau_beta_fh2<- 1 / (sigma_beta_fh2 * sigma_beta_fh2) # sigma^2 doesn't work in JAGS
  sigma_fh2 ~ dunif(0, 0.05) # standard deviation
  tau_fh2 <- 1 / (sigma_fh2 * sigma_fh2) # sigma^2 doesn't work in JAGS
  
  alpha_fh3 ~ dnorm(mu_alpha_fh3, tau_alpha_fh3) # intercept
  beta_fh3 ~ dnorm(mu_beta_fh3, tau_beta_fh3) # slope
  mu_alpha_fh3~ dunif(-30,30) 
  mu_beta_fh3~ dunif(-5, 5) 
  sigma_alpha_fh3 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh3 <- 1 / (sigma_alpha_fh3 * sigma_alpha_fh3) # sigma^2 doesn't work in JAGS
  sigma_beta_fh3~ dunif(0,1) # standard deviation
  tau_beta_fh3<- 1 / (sigma_beta_fh3 * sigma_beta_fh3) # sigma^2 doesn't work in JAGS
  sigma_fh3 ~ dunif(0, 0.05)# standard deviation
  tau_fh3 <- 1 / (sigma_fh3 * sigma_fh3) # sigma^2 doesn't work in JAGS
  
  alpha_fh4 ~ dnorm(mu_alpha_fh4, tau_alpha_fh4) # intercept
  beta_fh4 ~ dnorm(mu_beta_fh4, tau_beta_fh4) # slope
  mu_alpha_fh4~ dunif(-30,30) 
  mu_beta_fh4~ dunif(-5, 5) 
  sigma_alpha_fh4 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh4 <- 1 / (sigma_alpha_fh4 * sigma_alpha_fh4) # sigma^2 doesn't work in JAGS
  sigma_beta_fh4~ dunif(0, 1) # standard deviation
  tau_beta_fh4<- 1 / (sigma_beta_fh4 * sigma_beta_fh4) # sigma^2 doesn't work in JAGS
  sigma_fh4 ~ dunif(0, 0.05) # standard deviation
  tau_fh4 <- 1 / (sigma_fh4 * sigma_fh4) # sigma^2 doesn't work in JAGS
  
  alpha_fh5 ~ dnorm(mu_alpha_fh5, tau_alpha_fh5) # intercept
  beta_fh5 ~ dnorm(mu_beta_fh5, tau_beta_fh5) # slope
  mu_alpha_fh5~ dunif(-30,30) 
  mu_beta_fh5~ dunif(-5, 5) 
  sigma_alpha_fh5 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh5 <- 1 / (sigma_alpha_fh5 * sigma_alpha_fh5) # sigma^2 doesn't work in JAGS
  sigma_beta_fh5~ dunif(0, 1) # standard deviation
  tau_beta_fh5<- 1 / (sigma_beta_fh5* sigma_beta_fh5) # sigma^2 doesn't work in JAGS
  sigma_fh5 ~ dunif(0, 0.05) # standard deviation
  tau_fh5 <- 1 / (sigma_fh5 * sigma_fh5) # sigma^2 doesn't work in JAGS
  
  alpha_fh6 ~ dnorm(mu_alpha_fh6, tau_alpha_fh6) # intercept
  beta_fh6 ~ dnorm(mu_beta_fh6, tau_beta_fh6) # slope
  mu_alpha_fh6~ dunif(-30,30) 
  mu_beta_fh6~ dunif(-5, 5) 
  sigma_alpha_fh6 ~ dunif(0, 1) # standard deviation
  tau_alpha_fh6 <- 1 / (sigma_alpha_fh6 * sigma_alpha_fh6) # sigma^2 doesn't work in JAGS
  sigma_beta_fh6~ dunif(0, 1) # standard deviation
  tau_beta_fh6<- 1 / (sigma_beta_fh6* sigma_beta_fh6) # sigma^2 doesn't work in JAGS
  sigma_fh6 ~ dunif(0, 0.05) # standard deviation
  tau_fh6 <- 1 / (sigma_fh6 * sigma_fh6) # sigma^2 doesn't work in JAGS
  
  alpha_fv1 ~ dnorm(mu_alpha_fv1, tau_alpha_fv1) # intercept
  beta_fv1 ~ dnorm(mu_beta_fv1, tau_beta_fv1) # slope
  mu_alpha_fv1~ dunif(-30,30) 
  mu_beta_fv1~ dunif(-5, 5) 
  sigma_alpha_fv1 ~ dunif(0, 1) # standard deviation
  tau_alpha_fv1 <- 1 / (sigma_alpha_fv1 * sigma_alpha_fv1) # sigma^2 doesn't work in JAGS
  sigma_beta_fv1~ dunif(0,1) # standard deviation
  tau_beta_fv1<- 1 / (sigma_beta_fv1 * sigma_beta_fv1) # sigma^2 doesn't work in JAGS
  sigma_fv1 ~ dunif(0, 0.05) # standard deviation
  tau_fv1 <- 1 / (sigma_fv1 * sigma_fv1) # sigma^2 doesn't work in JAGS
  
  alpha_fv2 ~ dnorm(mu_alpha_fv2, tau_alpha_fv2) # intercept
  beta_fv2 ~ dnorm(mu_beta_fv2, tau_beta_fv2) # slope
  mu_alpha_fv2~ dunif(-30,30) 
  mu_beta_fv2~ dunif(-5, 5) 
  sigma_alpha_fv2 ~ dunif(0, 1) # standard deviation
  tau_alpha_fv2 <- 1 / (sigma_alpha_fv2 * sigma_alpha_fv2) # sigma^2 doesn't work in JAGS
  sigma_beta_fv2~ dunif(0,1) # standard deviation
  tau_beta_fv2<- 1 / (sigma_beta_fv2 * sigma_beta_fv2) # sigma^2 doesn't work in JAGS
  sigma_fv2 ~ dunif(0, 0.05) # standard deviation
  tau_fv2 <- 1 / (sigma_fv2 * sigma_fv2) # sigma^2 doesn't work in JAGS
  
  alpha_fv3 ~ dnorm(mu_alpha_fv3, tau_alpha_fv3) # intercept
  beta_fv3 ~ dnorm(mu_beta_fv3, tau_beta_fv3) # slope
  mu_alpha_fv3~ dunif(-30,30) 
  mu_beta_fv3~ dunif(-5, 5) 
  sigma_alpha_fv3 ~ dunif(0, 1) # standard deviation
  tau_alpha_fv3 <- 1 / (sigma_alpha_fv3 * sigma_alpha_fv3) # sigma^2 doesn't work in JAGS
  sigma_beta_fv3~ dunif(0,1) # standard deviation
  tau_beta_fv3<- 1 / (sigma_beta_fv3 * sigma_beta_fv3) # sigma^2 doesn't work in JAGS
  sigma_fv3 ~ dunif(0, 0.05) # standard deviation
  tau_fv3 <- 1 / (sigma_fv3 * sigma_fv3) # sigma^2 doesn't work in JAGS
  
  alpha_fv4 ~ dnorm(mu_alpha_fv4, tau_alpha_fv4) # intercept
  beta_fv4 ~ dnorm(mu_beta_fv4, tau_beta_fv4) # slope
  mu_alpha_fv4~ dunif(-30,30) 
  mu_beta_fv4~ dunif(-5, 5) 
  sigma_alpha_fv4 ~ dunif(0, 1) # standard deviation
  tau_alpha_fv4 <- 1 / (sigma_alpha_fv4 * sigma_alpha_fv4) # sigma^2 doesn't work in JAGS
  sigma_beta_fv4~ dunif(0,1) # standard deviation
  tau_beta_fv4<- 1 / (sigma_beta_fv4 * sigma_beta_fv4) # sigma^2 doesn't work in JAGS
  sigma_fv4 ~ dunif(0, 0.05) # standard deviation
  tau_fv4 <- 1 / (sigma_fv4 * sigma_fv4) # sigma^2 doesn't work in JAGS
  
  alpha_fv5 ~ dnorm(mu_alpha_fv5, tau_alpha_fv5) # intercept
  beta_fv5 ~ dnorm(mu_beta_fv5, tau_beta_fv5) # slope
  mu_alpha_fv5~ dunif(-30,30) 
  mu_beta_fv5~ dunif(-5, 5) 
  sigma_alpha_fv5 ~ dunif(0, 1) # standard deviation
  tau_alpha_fv5 <- 1 / (sigma_alpha_fv5 * sigma_alpha_fv5) # sigma^2 doesn't work in JAGS
  sigma_beta_fv5~ dunif(0,1) # standard deviation
  tau_beta_fv5<- 1 / (sigma_beta_fv5 * sigma_beta_fv5) # sigma^2 doesn't work in JAGS
  sigma_fv5 ~ dunif(0, 0.05) # standard deviation
  tau_fv5 <- 1 / (sigma_fv5 * sigma_fv5) # sigma^2 doesn't work in JAGS
  
  alpha_fg1 ~ dnorm(mu_alpha_fg1, tau_alpha_fg1) # intercept
  beta_fg1 ~ dnorm(mu_beta_fg1, tau_beta_fg1) # slope
  mu_alpha_fg1~ dunif(-30,30) 
  mu_beta_fg1~ dunif(-5, 5) 
  sigma_alpha_fg1 ~ dunif(0, 1) # standard deviation
  tau_alpha_fg1 <- 1 / (sigma_alpha_fg1 * sigma_alpha_fg1) # sigma^2 doesn't work in JAGS
  sigma_beta_fg1~ dunif(0,1) # standard deviation
  tau_beta_fg1<- 1 / (sigma_beta_fg1 * sigma_beta_fg1) # sigma^2 doesn't work in JAGS
  sigma_fg1 ~ dunif(0, 0.05) # standard deviation
  tau_fg1 <- 1 / (sigma_fg1 * sigma_fg1) # sigma^2 doesn't work in JAGS
  
  alpha_s1 ~ dnorm(mu_alpha_s1, tau_alpha_s1) # intercept
  beta_s1 ~ dnorm(mu_beta_s1, tau_beta_s1) # slope
  mu_alpha_s1~ dunif(-1, 1) 
  mu_beta_s1~ dunif(-1, 1) 
  sigma_alpha_s1 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s1 <- 1 / (sigma_alpha_s1 * sigma_alpha_s1) # sigma^2 doesn't work in JAGS
  sigma_beta_s1~ dunif(0, 0.1) # standard deviation
  tau_beta_s1<- 1 / (sigma_beta_s1 * sigma_beta_s1) # sigma^2 doesn't work in JAGS
  sigma_s1 ~ dunif(0,0.02)# standard deviation
  tau_s1 <- 1 / (sigma_s1 * sigma_s1) # sigma^2 doesn't work in JAGS
  
  alpha_s2 ~ dnorm(mu_alpha_s2, tau_alpha_s2) # intercept
  beta_s2 ~ dnorm(mu_beta_s2, tau_beta_s2) # slope
  mu_alpha_s2~ dunif(-1, 1) 
  mu_beta_s2~ dunif(-1, 1) 
  sigma_alpha_s2 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s2 <- 1 / (sigma_alpha_s2 * sigma_alpha_s2) # sigma^2 doesn't work in JAGS
  sigma_beta_s2~ dunif(0, 0.1) # standard deviation
  tau_beta_s2<- 1 / (sigma_beta_s2 * sigma_beta_s2) # sigma^2 doesn't work in JAGS
  sigma_s2 ~ dunif(0, 0.02) # standard deviation
  tau_s2 <- 1 / (sigma_s2 * sigma_s2) # sigma^2 doesn't work in JAGS
  
  alpha_s3 ~ dnorm(mu_alpha_s3, tau_alpha_s3) # intercept
  beta_s3 ~ dnorm(mu_beta_s3, tau_beta_s3) # slope
  mu_alpha_s3~ dunif(-1, 1) 
  mu_beta_s3~ dunif(-1, 1)  
  sigma_alpha_s3 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s3 <- 1 / (sigma_alpha_s3 * sigma_alpha_s3) # sigma^2 doesn't work in JAGS
  sigma_beta_s3~ dunif(0, 0.1) # standard deviation
  tau_beta_s3<- 1 / (sigma_beta_s3 * sigma_beta_s3) # sigma^2 doesn't work in JAGS
  sigma_s3 ~ dunif(0, 0.02) # standard deviation
  tau_s3 <- 1 / (sigma_s3 * sigma_s3) # sigma^2 doesn't work in JAGS
  
  alpha_s4 ~ dnorm(mu_alpha_s4, tau_alpha_s4) # intercept
  beta_s4 ~ dnorm(mu_beta_s4, tau_beta_s4) # slope
  mu_alpha_s4~ dunif(-1, 1) 
  mu_beta_s4~ dunif(-1, 1) 
  sigma_alpha_s4 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s4 <- 1 / (sigma_alpha_s4 * sigma_alpha_s4) # sigma^2 doesn't work in JAGS
  sigma_beta_s4~ dunif(0, 0.1) # standard deviation
  tau_beta_s4<- 1 / (sigma_beta_s4 * sigma_beta_s4) # sigma^2 doesn't work in JAGS
  sigma_s4 ~ dunif(0, 0.02) # standard deviation
  tau_s4 <- 1 / (sigma_s4 * sigma_s4) # sigma^2 doesn't work in JAGS
  
  alpha_s5 ~ dnorm(mu_alpha_s5, tau_alpha_s5) # intercept
  beta_s5 ~ dnorm(mu_beta_s5, tau_beta_s5) # slope
  mu_alpha_s5~ dunif(-1, 1) 
  mu_beta_s5~ dunif(-1, 1) 
  sigma_alpha_s5 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s5 <- 1 / (sigma_alpha_s5 * sigma_alpha_s5) # sigma^2 doesn't work in JAGS
  sigma_beta_s5~ dunif(0, 0.1) # standard deviation
  tau_beta_s5<- 1 / (sigma_beta_s5* sigma_beta_s5) # sigma^2 doesn't work in JAGS
  sigma_s5 ~dunif(0, 0.02) # standard deviation
  tau_s5 <- 1 / (sigma_s5 * sigma_s5) # sigma^2 doesn't work in JAGS
  
  alpha_s6 ~ dnorm(mu_alpha_s6, tau_alpha_s6) # intercept
  beta_s6 ~ dnorm(mu_beta_s6, tau_beta_s6) # slope
  mu_alpha_s6~ dunif(-1, 1) 
  mu_beta_s6~ dunif(-1, 1) 
  sigma_alpha_s6 ~ dunif(0, 0.1) # standard deviation
  tau_alpha_s6 <- 1 / (sigma_alpha_s6 * sigma_alpha_s6) # sigma^2 doesn't work in JAGS
  sigma_beta_s6~ dunif(0, 0.1) # standard deviation
  tau_beta_s6<- 1 / (sigma_beta_s6* sigma_beta_s6) # sigma^2 doesn't work in JAGS
  sigma_s6 ~ dunif(0, 0.02) # standard deviation
  tau_s6 <- 1 / (sigma_s6 * sigma_s6) # sigma^2 doesn't work in JAGS
  
  alpha_v1 ~ dnorm(mu_alpha_v1, tau_alpha_v1) # intercept
  beta_v1 ~ dnorm(mu_beta_v1, tau_beta_v1) # slope
  mu_alpha_v1~ dunif(-1, 1) 
  mu_beta_v1~ dunif(-1, 1) 
  sigma_alpha_v1 ~ dunif(0, 1) # standard deviation
  tau_alpha_v1 <- 1 / (sigma_alpha_v1 * sigma_alpha_v1) # sigma^2 doesn't work in JAGS
  sigma_beta_v1~ dunif(0, 1) # standard deviation
  tau_beta_v1<- 1 / (sigma_beta_v1 * sigma_beta_v1) # sigma^2 doesn't work in JAGS
  sigma_v1 ~ dunif(0, 0.02) # standard deviation
  tau_v1 <- 1 / (sigma_v1 * sigma_v1) # sigma^2 doesn't work in JAGS
  
  alpha_v2 ~ dnorm(mu_alpha_v2, tau_alpha_v2) # intercept
  beta_v2 ~ dnorm(mu_beta_v2, tau_beta_v2) # slope
  mu_alpha_v2~ dunif(-1, 1) 
  mu_beta_v2~ dunif(-1, 1) 
  sigma_alpha_v2 ~ dunif(0, 1) # standard deviation
  tau_alpha_v2 <- 1 / (sigma_alpha_v2 * sigma_alpha_v2) # sigma^2 doesn't work in JAGS
  sigma_beta_v2~ dunif(0, 1) # standard deviation
  tau_beta_v2<- 1 / (sigma_beta_v2 * sigma_beta_v2) # sigma^2 doesn't work in JAGS
  sigma_v2 ~ dunif(0, 0.02) # standard deviation
  tau_v2 <- 1 / (sigma_v2 * sigma_v2) # sigma^2 doesn't work in JAGS
  
  alpha_v3 ~ dnorm(mu_alpha_v3, tau_alpha_v3) # intercept
  beta_v3 ~ dnorm(mu_beta_v3, tau_beta_v3) # slope
  mu_alpha_v3~ dunif(-1, 1) 
  mu_beta_v3~ dunif(-1, 1) 
  sigma_alpha_v3 ~ dunif(0, 1) # standard deviation
  tau_alpha_v3 <- 1 / (sigma_alpha_v3 * sigma_alpha_v3) # sigma^2 doesn't work in JAGS
  sigma_beta_v3~ dunif(0, 1) # standard deviation
  tau_beta_v3<- 1 / (sigma_beta_v3 * sigma_beta_v3) # sigma^2 doesn't work in JAGS
  sigma_v3 ~ dunif(0, 0.02) # standard deviation
  tau_v3 <- 1 / (sigma_v3 * sigma_v3) # sigma^2 doesn't work in JAGS
  
  alpha_v4 ~ dnorm(mu_alpha_v4, tau_alpha_v4) # intercept
  beta_v4 ~ dnorm(mu_beta_v4, tau_beta_v4) # slope
  mu_alpha_v4~ dunif(-1, 1) 
  mu_beta_v4~ dunif(-1, 1) 
  sigma_alpha_v4 ~ dunif(0, 1) # standard deviation
  tau_alpha_v4 <- 1 / (sigma_alpha_v4 * sigma_alpha_v4) # sigma^2 doesn't work in JAGS
  sigma_beta_v4~ dunif(0, 1) # standard deviation
  tau_beta_v4<- 1 / (sigma_beta_v4 * sigma_beta_v4) # sigma^2 doesn't work in JAGS
  sigma_v4 ~ dunif(0, 0.02) # standard deviation
  tau_v4 <- 1 / (sigma_v4 * sigma_v4) # sigma^2 doesn't work in JAGS
  
  alpha_v5 ~ dnorm(mu_alpha_v5, tau_alpha_v5) # intercept
  beta_v5 ~ dnorm(mu_beta_v5, tau_beta_v5) # slope
  mu_alpha_v5~ dunif(-1, 1) 
  mu_beta_v5~ dunif(-1, 1) 
  sigma_alpha_v5 ~ dunif(0, 1) # standard deviation
  tau_alpha_v5 <- 1 / (sigma_alpha_v5 * sigma_alpha_v5) # sigma^2 doesn't work in JAGS
  sigma_beta_v5~ dunif(0, 1) # standard deviation
  tau_beta_v5<- 1 / (sigma_beta_v5* sigma_beta_v5) # sigma^2 doesn't work in JAGS
  sigma_v5 ~ dunif(0, 0.02) # standard deviation
  tau_v5 <- 1 / (sigma_v5 * sigma_v5) # sigma^2 doesn't work in JAGS
  
  alpha_y1 ~ dnorm(mu_alpha_y1, tau_alpha_y1) # intercept
  beta_y1 ~ dnorm(mu_beta_y1, tau_beta_y1) # slope
  mu_alpha_y1~ dunif(-1,1) 
  mu_beta_y1~ dunif(-1,1) 
  sigma_alpha_y1 ~ dunif(0, 1) # standard deviation
  tau_alpha_y1 <- 1 / (sigma_alpha_y1 * sigma_alpha_y1) # sigma^2 doesn't work in JAGS
  sigma_beta_y1~ dunif(0, 1) # standard deviation
  tau_beta_y1<- 1 / (sigma_beta_y1 * sigma_beta_y1) # sigma^2 doesn't work in JAGS
  sigma_y1 ~ dunif(0, 0.02) # standard deviation
  tau_y1 <- 1 / (sigma_y1 * sigma_y1) # sigma^2 doesn't work in JAGS
  
  alpha_y2 ~ dnorm(mu_alpha_y2, tau_alpha_y2) # intercept
  beta_y2 ~ dnorm(mu_beta_y2, tau_beta_y2) # slope
  mu_alpha_y2~ dunif(-1,1) 
  mu_beta_y2~ dunif(-1,1) 
  sigma_alpha_y2 ~ dunif(0, 1) # standard deviation
  tau_alpha_y2 <- 1 / (sigma_alpha_y2 * sigma_alpha_y2) # sigma^2 doesn't work in JAGS
  sigma_beta_y2~ dunif(0, 1) # standard deviation
  tau_beta_y2<- 1 / (sigma_beta_y2 * sigma_beta_y2) # sigma^2 doesn't work in JAGS
  sigma_y2 ~ dunif(0, 0.02) # standard deviation
  tau_y2 <- 1 / (sigma_y2* sigma_y2) # sigma^2 doesn't work in JAGS
  
  alpha_y3 ~ dnorm(mu_alpha_y3, tau_alpha_y3) # intercept
  beta_y3 ~ dnorm(mu_beta_y3, tau_beta_y3) # slope
  mu_alpha_y3~ dunif(-1,1) 
  mu_beta_y3~ dunif(-1,1) 
  sigma_alpha_y3 ~ dunif(0, 1) # standard deviation
  tau_alpha_y3 <- 1 / (sigma_alpha_y3 * sigma_alpha_y3) # sigma^2 doesn't work in JAGS
  sigma_beta_y3~ dunif(0, 1) # standard deviation
  tau_beta_y3<- 1 / (sigma_beta_y3 * sigma_beta_y3) # sigma^2 doesn't work in JAGS
  sigma_y3 ~ dunif(0, 0.02) # standard deviation
  tau_y3 <- 1 / (sigma_y3 * sigma_y3) # sigma^2 doesn't work in JAGS
  
  alpha_y1_ref ~ dnorm(mu_alpha_y1_ref, tau_alpha_y1_ref) # intercept
  beta1_y1_ref ~ dnorm(mu_beta1_y1_ref, tau_beta1_y1_ref) # slope
  beta2_y1_ref ~ dnorm(mu_beta2_y1_ref, tau_beta2_y1_ref) # slope
  mu_alpha_y1_ref~ dunif(-1,1) 
  mu_beta1_y1_ref~ dunif(-1,1) 
  mu_beta2_y1_ref~ dunif(-1,1) 
  sigma_alpha_y1_ref ~ dunif(0, 1) # standard deviation
  tau_alpha_y1_ref <- 1 / (sigma_alpha_y1_ref * sigma_alpha_y1_ref) # sigma^2 doesn't work in JAGS
  sigma_beta1_y1_ref~ dunif(0, 1) # standard deviation
  sigma_beta2_y1_ref~ dunif(0, 1) # standard deviation
  tau_beta1_y1_ref<- 1 / (sigma_beta1_y1_ref * sigma_beta1_y1_ref) # sigma^2 doesn't work in JAGS
  tau_beta2_y1_ref<- 1 / (sigma_beta2_y1_ref * sigma_beta2_y1_ref) # sigma^2 doesn't work in JAGS
  sigma_y1_ref ~ dunif(0, 0.02) # standard deviation
  tau_y1_ref <- 1 / (sigma_y1_ref * sigma_y1_ref) # sigma^2 doesn't work in JAGS
  
  alpha_y2_ref ~ dnorm(mu_alpha_y2_ref, tau_alpha_y2_ref) # intercept
  beta1_y2_ref ~ dnorm(mu_beta1_y2_ref, tau_beta1_y2_ref) # slope
  beta2_y2_ref ~ dnorm(mu_beta2_y2_ref, tau_beta2_y2_ref) # slope
  mu_alpha_y2_ref~ dunif(-1,1) 
  mu_beta1_y2_ref~ dunif(-1,1) 
  mu_beta2_y2_ref~ dunif(-1,1) 
  sigma_alpha_y2_ref ~ dunif(0, 1) # standard deviation
  tau_alpha_y2_ref <- 1 / (sigma_alpha_y2_ref * sigma_alpha_y2_ref) # sigma^2 doesn't work in JAGS
  sigma_beta1_y2_ref~ dunif(0, 1) # standard deviation
  sigma_beta2_y2_ref~ dunif(0, 1) # standard deviation
  tau_beta1_y2_ref<- 1 / (sigma_beta1_y2_ref * sigma_beta1_y2_ref) # sigma^2 doesn't work in JAGS
  tau_beta2_y2_ref<- 1 / (sigma_beta2_y2_ref * sigma_beta2_y2_ref) # sigma^2 doesn't work in JAGS
  sigma_y2_ref ~ dunif(0, 0.02) # standard deviation
  tau_y2_ref <- 1 / (sigma_y2_ref * sigma_y2_ref) # sigma^2 doesn't work in JAGS
  
  alpha_y3_ref ~ dnorm(mu_alpha_y3_ref, tau_alpha_y3_ref) # intercept
  beta1_y3_ref ~ dnorm(mu_beta1_y3_ref, tau_beta1_y3_ref) # slope
  beta2_y3_ref ~ dnorm(mu_beta2_y3_ref, tau_beta2_y3_ref) # slope
  mu_alpha_y3_ref~ dunif(-1,1) 
  mu_beta1_y3_ref~ dunif(-1,1) 
  mu_beta2_y3_ref~ dunif(-1,1) 
  sigma_alpha_y3_ref ~ dunif(0, 1) # standard deviation
  tau_alpha_y3_ref <- 1 / (sigma_alpha_y3_ref * sigma_alpha_y3_ref) # sigma^2 doesn't work in JAGS
  sigma_beta1_y3_ref~ dunif(0, 1) # standard deviation
  sigma_beta2_y3_ref~ dunif(0, 1) # standard deviation
  tau_beta1_y3_ref<- 1 / (sigma_beta1_y3_ref * sigma_beta1_y3_ref) # sigma^2 doesn't work in JAGS
  tau_beta2_y3_ref<- 1 / (sigma_beta2_y3_ref * sigma_beta2_y3_ref) # sigma^2 doesn't work in JAGS
  sigma_y3_ref ~ dunif(0, 0.02) # standard deviation
  tau_y3_ref <- 1 / (sigma_y3_ref * sigma_y3_ref) # sigma^2 doesn't work in JAGS
  
  for (j in 1:n_region){ 
    ar_fh1[j]~dnorm(phi_fh1[j],tau_ar_fh1[j]);T(0,1) # parameter
    tau_ar_fh1[j] <- 1/ (sigma_ar_fh1[j] * sigma_ar_fh1[j]) 
    sigma_ar_fh1[j]~ dunif(0, 1) # standard deviation
    phi_fh1[j]~ dunif(0, 1) 
    tau_j_fh1[j] <- 1/ (sigma_j_fh1[j] * sigma_j_fh1[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fh1[j] ~ dunif(0,0.03) # standard deviation
    
    ar_fh2[j]~dnorm(phi_fh2[j],tau_ar_fh2[j]);T(0,1) # parameter
    tau_ar_fh2[j] <- 1/ (sigma_ar_fh2[j] * sigma_ar_fh2[j]) 
    sigma_ar_fh2[j]~ dunif(0, 1) # standard deviation
    phi_fh2[j]~ dunif(0, 1) 
    tau_j_fh2[j] <- 1/ (sigma_j_fh2[j] * sigma_j_fh2[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fh2[j] ~ dunif(0,0.03) # standard deviation
    
    ar_fh3[j]~dnorm(phi_fh3[j],tau_ar_fh3[j]);T(0,1) # parameter
    tau_ar_fh3[j] <- 1/ (sigma_ar_fh3[j] * sigma_ar_fh3[j]) 
    sigma_ar_fh3[j]~ dunif(0, 1) # standard deviation
    phi_fh3[j]~ dunif(0, 1) 
    tau_j_fh3[j] <- 1/ (sigma_j_fh3[j] * sigma_j_fh3[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fh3[j] ~ dunif(0,0.03) # standard deviation
    
    ar_fh4[j]~dnorm(phi_fh4[j],tau_ar_fh4[j]);T(0,1) # parameter
    tau_ar_fh4[j] <- 1/ (sigma_ar_fh4[j] * sigma_ar_fh4[j]) 
    sigma_ar_fh4[j]~ dunif(0, 1) # standard deviation
    phi_fh4[j]~ dunif(0, 1) 
    tau_j_fh4[j] <- 1/ (sigma_j_fh4[j] * sigma_j_fh4[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fh4[j] ~ dunif(0,0.03) # standard deviation
    
    ar_fh5[j]~dnorm(phi_fh5[j],tau_ar_fh5[j]);T(0,1) # parameter
    tau_ar_fh5[j] <- 1/ (sigma_ar_fh5[j] * sigma_ar_fh5[j]) 
    sigma_ar_fh5[j]~ dunif(0, 1) # standard deviation
    phi_fh5[j]~ dunif(0, 1) 
    tau_j_fh5[j] <- 1/ (sigma_j_fh5[j] * sigma_j_fh5[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fh5[j] ~ dunif(0,0.03) # standard deviation
    
    ar_fh6[j]~dnorm(phi_fh6[j],tau_ar_fh6[j]);T(0,1) # parameter
    tau_ar_fh6[j] <- 1/ (sigma_ar_fh6[j] * sigma_ar_fh6[j]) 
    sigma_ar_fh6[j]~ dunif(0, 1) # standard deviation
    phi_fh6[j]~ dunif(0, 1) 
    tau_j_fh6[j] <- 1/ (sigma_j_fh6[j] * sigma_j_fh6[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fh6[j] ~ dunif(0,0.03) # standard deviation
    
    ar_fv1[j]~dnorm(phi_fv1[j],tau_ar_fv1[j]);T(0,1) # parameter
    tau_ar_fv1[j] <- 1/ (sigma_ar_fv1[j] * sigma_ar_fv1[j]) 
    sigma_ar_fv1[j]~ dunif(0, 1) # standard deviation
    phi_fv1[j]~ dunif(0, 1) 
    tau_j_fv1[j] <- 1/ (sigma_j_fv1[j] * sigma_j_fv1[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fv1[j] ~ dunif(0,0.03) # standard deviation
    
    ar_fv2[j]~dnorm(phi_fv2[j],tau_ar_fv2[j]);T(0,1) # parameter
    tau_ar_fv2[j] <- 1/ (sigma_ar_fv2[j] * sigma_ar_fv2[j]) 
    sigma_ar_fv2[j]~ dunif(0, 1) # standard deviation
    phi_fv2[j]~ dunif(0, 1) 
    tau_j_fv2[j] <- 1/ (sigma_j_fv2[j] * sigma_j_fv2[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fv2[j] ~ dunif(0,0.03) # standard deviation
    
    ar_fv3[j]~dnorm(phi_fv3[j],tau_ar_fv3[j]);T(0,1) # parameter
    tau_ar_fv3[j] <- 1/ (sigma_ar_fv3[j] * sigma_ar_fv3[j]) 
    sigma_ar_fv3[j]~ dunif(0, 1) # standard deviation
    phi_fv3[j]~ dunif(0, 1) 
    tau_j_fv3[j] <- 1/ (sigma_j_fv3[j] * sigma_j_fv3[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fv3[j] ~ dunif(0,0.03) # standard deviation
    
    ar_fv4[j]~dnorm(phi_fv4[j],tau_ar_fv4[j]);T(0,1) # parameter
    tau_ar_fv4[j] <- 1/ (sigma_ar_fv4[j] * sigma_ar_fv4[j]) 
    sigma_ar_fv4[j]~ dunif(0, 1) # standard deviation
    phi_fv4[j]~ dunif(0, 1) 
    tau_j_fv4[j] <- 1/ (sigma_j_fv4[j] * sigma_j_fv4[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fv4[j] ~ dunif(0,0.03) # standard deviation
    
    ar_fv5[j]~dnorm(phi_fv5[j],tau_ar_fv5[j]);T(0,1) # parameter
    tau_ar_fv5[j] <- 1/ (sigma_ar_fv5[j] * sigma_ar_fv5[j]) 
    sigma_ar_fv5[j]~ dunif(0, 1) # standard deviation
    phi_fv5[j]~ dunif(0, 1) 
    tau_j_fv5[j] <- 1/ (sigma_j_fv5[j] * sigma_j_fv5[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fv5[j] ~ dunif(0,0.03) # standard deviation
    
    ar_fg1[j]~dnorm(phi_fg1[j],tau_ar_fg1[j]);T(0,1) # parameter
    tau_ar_fg1[j] <- 1/ (sigma_ar_fg1[j] * sigma_ar_fg1[j]) 
    sigma_ar_fg1[j]~ dunif(0, 1) # standard deviation
    phi_fg1[j]~ dunif(0, 1) 
    tau_j_fg1[j] <- 1/ (sigma_j_fg1[j] * sigma_j_fg1[j]) # sigma^2 doesn't work in JAGS
    sigma_j_fg1[j] ~ dunif(0,0.03) # standard deviation
    
    ar_s1[j]~dnorm(phi_s1[j],tau_ar_s1[j]);T(0,1) # parameter
    tau_ar_s1[j] <- 1/ (sigma_ar_s1[j] * sigma_ar_s1[j]) 
    sigma_ar_s1[j]~ dunif(0, 1) # standard deviation
    phi_s1[j]~ dunif(0, 1) 
    tau_j_s1[j] <- 1/ (sigma_j_s1[j] * sigma_j_s1[j]) # sigma^2 doesn't work in JAGS
    sigma_j_s1[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_s2[j]~dnorm(phi_s2[j],tau_ar_s2[j]);T(0,1) # parameter
    tau_ar_s2[j] <- 1/ (sigma_ar_s2[j] * sigma_ar_s2[j]) 
    sigma_ar_s2[j]~ dunif(0, 1) # standard deviation
    phi_s2[j]~ dunif(0, 1) 
    tau_j_s2[j] <- 1/ (sigma_j_s2[j] * sigma_j_s2[j]) # sigma^2 doesn't work in JAGS
    sigma_j_s2[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_s3[j]~dnorm(phi_s3[j],tau_ar_s3[j]);T(0,1) # parameter
    tau_ar_s3[j] <- 1/ (sigma_ar_s3[j] * sigma_ar_s3[j]) 
    sigma_ar_s3[j]~ dunif(0, 1) # standard deviation
    phi_s3[j]~ dunif(0, 1) 
    tau_j_s3[j] <- 1/ (sigma_j_s3[j] * sigma_j_s3[j]) # sigma^2 doesn't work in JAGS
    sigma_j_s3[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_s4[j]~dnorm(phi_s4[j],tau_ar_s4[j]);T(0,1) # parameter
    tau_ar_s4[j] <- 1/ (sigma_ar_s4[j] * sigma_ar_s4[j]) 
    sigma_ar_s4[j]~ dunif(0, 1) # standard deviation
    phi_s4[j]~ dunif(0, 1) 
    tau_j_s4[j] <- 1/ (sigma_j_s4[j] * sigma_j_s4[j]) # sigma^2 doesn't work in JAGS
    sigma_j_s4[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_s5[j]~dnorm(phi_s5[j],tau_ar_s5[j]);T(0,1) # parameter
    tau_ar_s5[j] <- 1/ (sigma_ar_s5[j] * sigma_ar_s5[j]) 
    sigma_ar_s5[j]~ dunif(0, 1) # standard deviation
    phi_s5[j]~ dunif(0, 1) 
    tau_j_s5[j] <- 1/ (sigma_j_s5[j] * sigma_j_s5[j]) # sigma^2 doesn't work in JAGS
    sigma_j_s5[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_s6[j]~dnorm(phi_s6[j],tau_ar_s6[j]);T(0,1) # parameter
    tau_ar_s6[j] <- 1/ (sigma_ar_s6[j] * sigma_ar_s6[j]) 
    sigma_ar_s6[j]~ dunif(0, 1) # standard deviation
    phi_s6[j]~ dunif(0, 1) 
    tau_j_s6[j] <- 1/ (sigma_j_s6[j] * sigma_j_s6[j]) # sigma^2 doesn't work in JAGS
    sigma_j_s6[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_v1[j]~dnorm(phi_v1[j],tau_ar_v1[j]);T(0,1) # parameter
    tau_ar_v1[j] <- 1/ (sigma_ar_v1[j] * sigma_ar_v1[j]) 
    sigma_ar_v1[j]~ dunif(0, 1) # standard deviation
    phi_v1[j]~ dunif(0, 1) 
    tau_j_v1[j] <- 1/ (sigma_j_v1[j] * sigma_j_v1[j]) # sigma^2 doesn't work in JAGS
    sigma_j_v1[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_v2[j]~dnorm(phi_v2[j],tau_ar_v2[j]);T(0,1) # parameter
    tau_ar_v2[j] <- 1/ (sigma_ar_v2[j] * sigma_ar_v2[j]) 
    sigma_ar_v2[j]~ dunif(0, 1) # standard deviation
    phi_v2[j]~ dunif(0, 1) 
    tau_j_v2[j] <- 1/ (sigma_j_v2[j] * sigma_j_v2[j]) # sigma^2 doesn't work in JAGS
    sigma_j_v2[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_v3[j]~dnorm(phi_v3[j],tau_ar_v3[j]);T(0,1) # parameter
    tau_ar_v3[j] <- 1/ (sigma_ar_v3[j] * sigma_ar_v3[j]) 
    sigma_ar_v3[j]~ dunif(0, 1) # standard deviation
    phi_v3[j]~ dunif(0, 1) 
    tau_j_v3[j] <- 1/ (sigma_j_v3[j] * sigma_j_v3[j]) # sigma^2 doesn't work in JAGS
    sigma_j_v3[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_v4[j]~dnorm(phi_v4[j],tau_ar_v4[j]);T(0,1) # parameter
    tau_ar_v4[j] <- 1/ (sigma_ar_v4[j] * sigma_ar_v4[j]) 
    sigma_ar_v4[j]~ dunif(0, 1) # standard deviation
    phi_v4[j]~ dunif(0, 1) 
    tau_j_v4[j] <- 1/ (sigma_j_v4[j] * sigma_j_v4[j]) # sigma^2 doesn't work in JAGS
    sigma_j_v4[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_v5[j]~dnorm(phi_v5[j],tau_ar_v5[j]);T(0,1) # parameter
    tau_ar_v5[j] <- 1/ (sigma_ar_v5[j] * sigma_ar_v5[j]) 
    sigma_ar_v5[j]~ dunif(0, 1) # standard deviation
    phi_v5[j]~ dunif(0, 1) 
    tau_j_v5[j] <- 1/ (sigma_j_v5[j] * sigma_j_v5[j]) # sigma^2 doesn't work in JAGS
    sigma_j_v5[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_y1[j]~dnorm(phi_y1[j],tau_ar_y1[j]);T(0,1) # parameter
    tau_ar_y1[j] <- 1/ (sigma_ar_y1[j] * sigma_ar_y1[j]) 
    sigma_ar_y1[j]~ dunif(0, 1) # standard deviation
    phi_y1[j]~ dunif(0, 1) 
    tau_j_y1[j] <- 1/ (sigma_j_y1[j] * sigma_j_y1[j]) # sigma^2 doesn't work in JAGS
    sigma_j_y1[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_y2[j]~dnorm(phi_y2[j],tau_ar_y2[j]);T(0,1) # parameter
    tau_ar_y2[j] <- 1/ (sigma_ar_y2[j] * sigma_ar_y2[j]) 
    sigma_ar_y2[j]~ dunif(0, 1) # standard deviation
    phi_y2[j]~ dunif(0, 1) 
    tau_j_y2[j] <- 1/ (sigma_j_y2[j] * sigma_j_y2[j]) # sigma^2 doesn't work in JAGS
    sigma_j_y2[j] ~ dunif(0, 0.02) # standard deviation
    
    ar_y3[j]~dnorm(phi_y3[j],tau_ar_y3[j]);T(0,1) # parameter
    tau_ar_y3[j] <- 1/ (sigma_ar_y3[j] * sigma_ar_y3[j]) 
    sigma_ar_y3[j]~ dunif(0, 1) # standard deviation
    phi_y3[j]~ dunif(0, 1) 
    tau_j_y3[j] <- 1/ (sigma_j_y3[j] * sigma_j_y3[j]) # sigma^2 doesn't work in JAGS
    sigma_j_y3[j] ~ dunif(0, 0.02) # standard deviation
  }
}

params <- c("alpha_gdp_pre[1]","alpha_gdp_pre[2]","alpha_gdp_pre[3]","alpha_gdp_pre[4]","alpha_gdp_pre[5]",
            "alpha_gdp_pre[6]","alpha_gdp_pre[7]","alpha_gdp_pre[8]","alpha_gdp_pre[9]","alpha_gdp_pre[10]",
            "alpha_gdp_pre[11]","alpha_gdp_pre[12]","alpha_gdp_pre[13]","alpha_gdp_pre[14]","alpha_gdp_pre[15]",
            "alpha_gdp_pre[16]","alpha_gdp_pre[17]","alpha_gdp_pre[18]","alpha_gdp_pre[19]","alpha_gdp_pre[20]",
            "alpha_gdp_pre[21]","alpha_gdp_pre[22]","alpha_gdp_pre[23]","alpha_gdp_pre[24]","alpha_gdp_pre[25]",
            "alpha_gdp_pre[26]","alpha_gdp_pre[27]","alpha_gdp_pre[28]","alpha_gdp_pre[29]","alpha_gdp_pre[30]",
            "alpha_gdp_pre[31]","alpha_gdp_pre[32]","alpha_gdp_pre[33]","alpha_gdp_pre[34]","alpha_gdp_pre[35]",
            "alpha_gdp_pre[36]","alpha_gdp_pre[37]","alpha_gdp_pre[38]","alpha_gdp_pre[39]","alpha_gdp_pre[40]",
            "error1[1]","error1[2]","error1[3]","error1[4]","error1[5]","error1[6]","error1[7]","error1[8]","error1[9]","error1[10]",
            "error1[11]","error1[12]","error1[13]","error1[14]","error1[15]", "error1[16]","error1[17]","error1[18]","error1[19]","error1[20]",
            "error1[21]","error1[22]","error1[23]","error1[24]","error1[25]", "error1[26]","error1[27]","error1[28]","error1[29]","error1[30]",
            "error1[31]","error1[32]","error1[33]","error1[34]","error1[35]", "error1[36]","error1[37]","error1[38]","error1[39]","error1[40]",
            "error2[1]","error2[2]","error2[3]","error2[4]","error2[5]","error2[6]","error2[7]","error2[8]","error2[9]","error2[10]",
            "error2[11]","error2[12]","error2[13]","error2[14]","error2[15]","error2[16]","error2[17]","error2[18]","error2[19]","error2[20]",
            "error2[21]","error2[22]","error2[23]","error2[24]","error2[25]", "error2[26]","error2[27]","error2[28]","error2[29]","error2[30]",
            "error2[31]","error2[32]","error2[33]","error2[34]","error2[35]", "error2[36]","error2[37]","error2[38]","error2[39]","error2[40]",
            "error3[1]","error3[2]","error3[3]","error3[4]","error3[5]","error3[6]","error3[7]","error3[8]","error3[9]","error3[10]",
            "error3[11]","error3[12]","error3[13]","error3[14]","error3[15]", "error3[16]","error3[17]","error3[18]","error3[19]","error3[20]",
            "error3[21]","error3[22]","error3[23]","error3[24]","error3[25]", "error3[26]","error3[27]","error3[28]","error3[29]","error3[30]",
            "error3[31]","error3[32]","error3[33]","error3[34]","error3[35]", "error3[36]","error3[37]","error3[38]","error3[39]","error3[40]",           
            "error4[1]","error4[2]","error4[3]","error4[4]","error4[5]","error4[6]","error4[7]","error4[8]","error4[9]","error4[10]",
            "error4[11]","error4[12]","error4[13]","error4[14]","error4[15]", "error4[16]","error4[17]","error4[18]","error4[19]","error4[20]",
            "error4[21]","error4[22]","error4[23]","error4[24]","error4[25]", "error4[26]","error4[27]","error4[28]","error4[29]","error4[30]",
            "error4[31]","error4[32]","error4[33]","error4[34]","error4[35]", "error4[36]","error4[37]","error4[38]","error4[39]","error4[40]",
            "error5[1]","error5[2]","error5[3]","error5[4]","error5[5]","error5[6]","error5[7]","error5[8]","error5[9]","error5[10]",
            "error5[11]","error5[12]","error5[13]","error5[14]","error5[15]", "error5[16]","error5[17]","error5[18]","error5[19]","error5[20]",
            "error5[21]","error5[22]","error5[23]","error5[24]","error5[25]", "error5[26]","error5[27]","error5[28]","error5[29]","error5[30]",
            "error5[31]","error5[32]","error5[33]","error5[34]","error5[35]", "error5[36]","error5[37]","error5[38]","error5[39]","error5[40]",
            "error6[1]","error6[2]","error6[3]","error6[4]","error6[5]","error6[6]","error6[7]","error6[8]","error6[9]","error6[10]",
            "error6[11]","error6[12]","error6[13]","error6[14]","error6[15]","error6[16]","error6[17]","error6[18]","error6[19]","error6[20]",
            "error6[21]","error6[22]","error6[23]","error6[24]","error6[25]", "error6[26]","error6[27]","error6[28]","error6[29]","error6[30]",
            "error6[31]","error6[32]","error6[33]","error6[34]","error6[35]", "error6[36]","error6[37]","error6[38]","error6[39]","error6[40]",
            "error7[1]","error7[2]","error7[3]","error7[4]","error7[5]","error7[6]","error7[7]","error7[8]","error7[9]","error7[10]",
            "error7[11]","error7[12]","error7[13]","error7[14]","error7[15]", "error7[16]","error7[17]","error7[18]","error7[19]","error7[20]",
            "error7[21]","error7[22]","error7[23]","error7[24]","error7[25]", "error7[26]","error7[27]","error7[28]","error7[29]","error7[30]",
            "error7[31]","error7[32]","error7[33]","error7[34]","error7[35]", "error7[36]","error7[37]","error7[38]","error7[39]","error7[40]",
            "error8[1]","error8[2]","error8[3]","error8[4]","error8[5]","error8[6]","error8[7]","error8[8]","error8[9]","error8[10]",
            "error8[11]","error8[12]","error8[13]","error8[14]","error8[15]","error8[16]","error8[17]","error8[18]","error8[19]","error8[20]",
            "error8[21]","error8[22]","error8[23]","error8[24]","error8[25]", "error8[26]","error8[27]","error8[28]","error8[29]","error8[30]",
            "error8[31]","error8[32]","error8[33]","error8[34]","error8[35]", "error8[36]","error8[37]","error8[38]","error8[39]","error8[40]",
            "error9[1]","error9[2]","error9[3]","error9[4]","error9[5]","error9[6]","error9[7]","error9[8]","error9[9]","error9[10]",
            "error9[11]","error9[12]","error9[13]","error9[14]","error9[15]", "error9[16]","error9[17]","error9[18]","error9[19]","error9[20]",
            "error9[21]","error9[22]","error9[23]","error9[24]","error9[25]", "error9[26]","error9[27]","error9[28]","error9[29]","error9[30]",
            "error9[31]","error9[32]","error9[33]","error9[34]","error9[35]", "error9[36]","error9[37]","error9[38]","error9[39]","error9[40]",
            "error10[1]","error10[2]","error10[3]","error10[4]","error10[5]","error10[6]","error10[7]","error10[8]","error10[9]","error10[10]",
            "error10[11]","error10[12]","error10[13]","error10[14]","error10[15]", "error10[16]","error10[17]","error10[18]","error10[19]","error10[20]",
            "error10[21]","error10[22]","error10[23]","error10[24]","error10[25]", "error10[26]","error10[27]","error10[28]","error10[29]","error10[30]",
            "error10[31]","error10[32]","error10[33]","error10[34]","error10[35]", "error10[36]","error10[37]","error10[38]","error10[39]","error10[40]",
            "error11[1]","error11[2]","error11[3]","error11[4]","error11[5]","error11[6]","error11[7]","error11[8]","error11[9]","error11[10]",
            "error11[11]","error11[12]","error11[13]","error11[14]","error11[15]", "error11[16]","error11[17]","error11[18]","error11[19]","error11[20]",
            "error11[21]","error11[22]","error11[23]","error11[24]","error11[25]", "error11[26]","error11[27]","error11[28]","error11[29]","error11[30]",
            "error11[31]","error11[32]","error11[33]","error11[34]","error11[35]", "error11[36]","error11[37]","error11[38]","error11[39]","error11[40]",
            "error12[1]","error12[2]","error12[3]","error12[4]","error12[5]","error12[6]","error12[7]","error12[8]","error12[9]","error12[10]",
            "error12[11]","error12[12]","error12[13]","error12[14]","error12[15]", "error12[16]","error12[17]","error12[18]","error12[19]","error12[20]",
            "error12[21]","error12[22]","error12[23]","error12[24]","error12[25]", "error12[26]","error12[27]","error12[28]","error12[29]","error12[30]",
            "error12[31]","error12[32]","error12[33]","error12[34]","error12[35]", "error12[36]","error12[37]","error12[38]","error12[39]","error12[40]",
            "error13[1]","error13[2]","error13[3]","error13[4]","error13[5]","error13[6]","error13[7]","error13[8]","error13[9]","error13[10]",
            "error13[11]","error13[12]","error13[13]","error13[14]","error13[15]", "error13[16]","error13[17]","error13[18]","error13[19]","error13[20]",
            "error13[21]","error13[22]","error13[23]","error13[24]","error13[25]", "error13[26]","error13[27]","error13[28]","error13[29]","error13[30]",
            "error13[31]","error13[32]","error13[33]","error13[34]","error13[35]", "error13[36]","error13[37]","error13[38]","error13[39]","error13[40]",
            "error14[1]","error14[2]","error14[3]","error14[4]","error14[5]","error14[6]","error14[7]","error14[8]","error14[9]","error14[10]",
            "error14[11]","error14[12]","error14[13]","error14[14]","error14[15]", "error14[16]","error14[17]","error14[18]","error14[19]","error14[20]",
            "error14[21]","error14[22]","error14[23]","error14[24]","error14[25]", "error14[26]","error14[27]","error14[28]","error14[29]","error14[30]",
            "error14[31]","error14[32]","error14[33]","error14[34]","error14[35]", "error14[36]","error14[37]","error14[38]","error14[39]","error14[40]",
            "error15[1]","error15[2]","error15[3]","error15[4]","error15[5]","error15[6]","error15[7]","error15[8]","error15[9]","error15[10]",
            "error15[11]","error15[12]","error15[13]","error15[14]","error15[15]", "error15[16]","error15[17]","error15[18]","error15[19]","error15[20]",
            "error15[21]","error15[22]","error15[23]","error15[24]","error15[25]", "error15[26]","error15[27]","error15[28]","error15[29]","error15[30]",
            "error15[31]","error15[32]","error15[33]","error15[34]","error15[35]", "error15[36]","error15[37]","error15[38]","error15[39]","error15[40]",
            "error16[1]","error16[2]","error16[3]","error16[4]","error16[5]","error16[6]","error16[7]","error16[8]","error16[9]","error16[10]",
            "error16[11]","error16[12]","error16[13]","error16[14]","error16[15]", "error16[16]","error16[17]","error16[18]","error16[19]","error16[20]",
            "error16[21]","error16[22]","error16[23]","error16[24]","error16[25]", "error16[26]","error16[27]","error16[28]","error16[29]","error16[30]",
            "error16[31]","error16[32]","error16[33]","error16[34]","error16[35]", "error16[36]","error16[37]","error16[38]","error16[39]","error16[40]",
            "error17[1]","error17[2]","error17[3]","error17[4]","error17[5]","error17[6]","error17[7]","error17[8]","error17[9]","error17[10]",
            "error17[11]","error17[12]","error17[13]","error17[14]","error17[15]", "error17[16]","error17[17]","error17[18]","error17[19]","error17[20]",
            "error17[21]","error17[22]","error17[23]","error17[24]","error17[25]", "error17[26]","error17[27]","error17[28]","error17[29]","error17[30]",
            "error17[31]","error17[32]","error17[33]","error17[34]","error17[35]", "error17[36]","error17[37]","error17[38]","error17[39]","error17[40]",
            "error18[1]","error18[2]","error18[3]","error18[4]","error18[5]","error18[6]","error18[7]","error18[8]","error18[9]","error18[10]",
            "error18[11]","error18[12]","error18[13]","error18[14]","error18[15]","error18[16]","error18[17]","error18[18]","error18[19]","error18[20]",
            "error18[21]","error18[22]","error18[23]","error18[24]","error18[25]", "error18[26]","error18[27]","error18[28]","error18[29]","error18[30]",
            "error18[31]","error18[32]","error18[33]","error18[34]","error18[35]", "error18[36]","error18[37]","error18[38]","error18[39]","error18[40]",
            "error19[1]","error19[2]","error19[3]","error19[4]","error19[5]","error19[6]","error19[7]","error19[8]","error19[9]","error19[10]",
            "error19[11]","error19[12]","error19[13]","error19[14]","error19[15]", "error19[16]","error19[17]","error19[18]","error19[19]","error19[20]",
            "error19[21]","error19[22]","error19[23]","error19[24]","error19[25]", "error19[26]","error19[27]","error19[28]","error19[29]","error19[30]",
            "error19[31]","error19[32]","error19[33]","error19[34]","error19[35]", "error19[36]","error19[37]","error19[38]","error19[39]","error19[40]",
            "error20[1]","error20[2]","error20[3]","error20[4]","error20[5]","error20[6]","error20[7]","error20[8]","error20[9]","error20[10]",
            "error20[11]","error20[12]","error20[13]","error20[14]","error20[15]", "error20[16]","error20[17]","error20[18]","error20[19]","error20[20]",
            "error20[21]","error20[22]","error20[23]","error20[24]","error20[25]", "error20[26]","error20[27]","error20[28]","error20[29]","error20[30]",
            "error20[31]","error20[32]","error20[33]","error20[34]","error20[35]", "error20[36]","error20[37]","error20[38]","error20[39]","error20[40]",
            "error21[1]","error21[2]","error21[3]","error21[4]","error21[5]","error21[6]","error21[7]","error21[8]","error21[9]","error21[10]",
            "error21[11]","error21[12]","error21[13]","error21[14]","error21[15]", "error21[16]","error21[17]","error21[18]","error21[19]","error21[20]",
            "error21[21]","error21[22]","error21[23]","error21[24]","error21[25]", "error21[26]","error21[27]","error21[28]","error21[29]","error21[30]",
            "error21[31]","error21[32]","error21[33]","error21[34]","error21[35]", "error21[36]","error21[37]","error21[38]","error21[39]","error21[40]",
            "error22[1]","error22[2]","error22[3]","error22[4]","error22[5]","error22[6]","error22[7]","error22[8]","error22[9]","error22[10]",
            "error22[11]","error22[12]","error22[13]","error22[14]","error22[15]", "error22[16]","error22[17]","error22[18]","error22[19]","error22[20]",
            "error22[21]","error22[22]","error22[23]","error22[24]","error22[25]", "error22[26]","error22[27]","error22[28]","error22[29]","error22[30]",
            "error22[31]","error22[32]","error22[33]","error22[34]","error22[35]", "error22[36]","error22[37]","error22[38]","error22[39]","error22[40]",
            "error23[1]","error23[2]","error23[3]","error23[4]","error23[5]","error23[6]","error23[7]","error23[8]","error23[9]","error23[10]",
            "error23[11]","error23[12]","error23[13]","error23[14]","error23[15]", "error23[16]","error23[17]","error23[18]","error23[19]","error23[20]",
            "error23[21]","error23[22]","error23[23]","error23[24]","error23[25]", "error23[26]","error23[27]","error23[28]","error23[29]","error23[30]",
            "error23[31]","error23[32]","error23[33]","error23[34]","error23[35]", "error23[36]","error23[37]","error23[38]","error23[39]","error23[40]",
            "error24[1]","error24[2]","error24[3]","error24[4]","error24[5]","error24[6]","error24[7]","error24[8]","error24[9]","error24[10]",
            "error24[11]","error24[12]","error24[13]","error24[14]","error24[15]", "error24[16]","error24[17]","error24[18]","error24[19]","error24[20]",
            "error24[21]","error24[22]","error24[23]","error24[24]","error24[25]", "error24[26]","error24[27]","error24[28]","error24[29]","error24[30]",
            "error24[31]","error24[32]","error24[33]","error24[34]","error24[35]", "error24[36]","error24[37]","error24[38]","error24[39]","error24[40]",
            "error25[1]","error25[2]","error25[3]","error25[4]","error25[5]","error25[6]","error25[7]","error25[8]","error25[9]","error25[10]",
            "error25[11]","error25[12]","error25[13]","error25[14]","error25[15]", "error25[16]","error25[17]","error25[18]","error25[19]","error25[20]",
            "error25[21]","error25[22]","error25[23]","error25[24]","error25[25]", "error25[26]","error25[27]","error25[28]","error25[29]","error25[30]",
            "error25[31]","error25[32]","error25[33]","error25[34]","error25[35]", "error25[36]","error25[37]","error25[38]","error25[39]","error25[40]",
            "error26[1]","error26[2]","error26[3]","error26[4]","error26[5]","error26[6]","error26[7]","error26[8]","error26[9]","error26[10]",
            "error26[11]","error26[12]","error26[13]","error26[14]","error26[15]", "error26[16]","error26[17]","error26[18]","error26[19]","error26[20]",
            "error26[21]","error26[22]","error26[23]","error26[24]","error26[25]", "error26[26]","error26[27]","error26[28]","error26[29]","error26[30]",
            "error26[31]","error26[32]","error26[33]","error26[34]","error26[35]", "error26[36]","error26[37]","error26[38]","error26[39]","error26[40]",
            "error27[1]","error27[2]","error27[3]","error27[4]","error27[5]","error27[6]","error27[7]","error27[8]","error27[9]","error27[10]",
            "error27[11]","error27[12]","error27[13]","error27[14]","error27[15]", "error27[16]","error27[17]","error27[18]","error27[19]","error27[20]",
            "error27[21]","error27[22]","error27[23]","error27[24]","error27[25]", "error27[26]","error27[27]","error27[28]","error27[29]","error27[30]",
            "error27[31]","error27[32]","error27[33]","error27[34]","error27[35]", "error27[36]","error27[37]","error27[38]","error27[39]","error27[40]",
            "error28[1]","error28[2]","error28[3]","error28[4]","error28[5]","error28[6]","error28[7]","error28[8]","error28[9]","error28[10]",
            "error28[11]","error28[12]","error28[13]","error28[14]","error28[15]", "error28[16]","error28[17]","error28[18]","error28[19]","error28[20]",
            "error28[21]","error28[22]","error28[23]","error28[24]","error28[25]", "error28[26]","error28[27]","error28[28]","error28[29]","error28[30]",
            "error28[31]","error28[32]","error28[33]","error28[34]","error28[35]", "error28[36]","error28[37]","error28[38]","error28[39]","error28[40]",
            "error29[1]","error29[2]","error29[3]","error29[4]","error29[5]","error29[6]","error29[7]","error29[8]","error29[9]","error29[10]",
            "error29[11]","error29[12]","error29[13]","error29[14]","error29[15]", "error29[16]","error29[17]","error29[18]","error29[19]","error29[20]",
            "error29[21]","error29[22]","error29[23]","error29[24]","error29[25]", "error29[26]","error29[27]","error29[28]","error29[29]","error29[30]",
            "error29[31]","error29[32]","error29[33]","error29[34]","error29[35]", "error29[36]","error29[37]","error29[38]","error29[39]","error29[40]",
            "error30[1]","error30[2]","error30[3]","error30[4]","error30[5]","error30[6]","error30[7]","error30[8]","error30[9]","error30[10]",
            "error30[11]","error30[12]","error30[13]","error30[14]","error30[15]", "error30[16]","error30[17]","error30[18]","error30[19]","error30[20]",
            "error30[21]","error30[22]","error30[23]","error30[24]","error30[25]", "error30[26]","error30[27]","error30[28]","error30[29]","error30[30]",
            "error30[31]","error30[32]","error30[33]","error30[34]","error30[35]", "error30[36]","error30[37]","error30[38]","error30[39]","error30[40]",
            "error31[1]","error31[2]","error31[3]","error31[4]","error31[5]","error31[6]","error31[7]","error31[8]","error31[9]","error31[10]",
            "error31[11]","error31[12]","error31[13]","error31[14]","error31[15]", "error31[16]","error31[17]","error31[18]","error31[19]","error31[20]",
            "error31[21]","error31[22]","error31[23]","error31[24]","error31[25]", "error31[26]","error31[27]","error31[28]","error31[29]","error31[30]",
            "error31[31]","error31[32]","error31[33]","error31[34]","error31[35]", "error31[36]","error31[37]","error31[38]","error31[39]","error31[40]",          
            "error32[1]","error32[2]","error32[3]","error32[4]","error32[5]","error32[6]","error32[7]","error32[8]","error32[9]","error32[10]",
            "error32[11]","error32[12]","error32[13]","error32[14]","error32[15]", "error32[16]","error32[17]","error32[18]","error32[19]","error32[20]",
            "error32[21]","error32[22]","error32[23]","error32[24]","error32[25]", "error32[26]","error32[27]","error32[28]","error32[29]","error32[30]",
            "error32[31]","error32[32]","error32[33]","error32[34]","error32[35]", "error32[36]","error32[37]","error32[38]","error32[39]","error32[40]",
            "alpha_fh1", "beta_fh1", "sigma_fh1","alpha_fh2", "beta_fh2", "sigma_fh2","alpha_fh3", "beta_fh3", "sigma_fh3","alpha_fh4", "beta_fh4", "sigma_fh4",
            "alpha_fh5", "beta_fh5", "sigma_fh5","alpha_fh6", "beta_fh6", "sigma_fh6",
            "alpha_fv1", "beta_fv1", "sigma_fv1","alpha_fv2", "beta_fv2", "sigma_fv2","alpha_fv3", "beta_fv3", "sigma_fv3","alpha_fv4", "beta_fv4", "sigma_fv4",
            "alpha_fv5", "beta_fv5", "sigma_fv5",
            "alpha_fg1", "beta_fg1", "sigma_fg1",
            "ar_fh1[1]","sigma_j_fh1[1]","ar_fh1[2]", "sigma_j_fh1[2]", "ar_fh1[3]", "sigma_j_fh1[3]", "ar_fh1[4]", "sigma_j_fh1[4]",
            "ar_fh1[5]","sigma_j_fh1[5]","ar_fh1[6]", "sigma_j_fh1[6]", "ar_fh1[7]", "sigma_j_fh1[7]", "ar_fh1[8]", "sigma_j_fh1[8]",
            "ar_fh1[9]","sigma_j_fh1[9]","ar_fh1[10]","sigma_j_fh1[10]","ar_fh1[11]","sigma_j_fh1[11]","ar_fh1[12]","sigma_j_fh1[12]",
            "ar_fh1[13]", "sigma_j_fh1[13]","ar_fh1[14]", "sigma_j_fh1[14]","ar_fh1[15]", "sigma_j_fh1[15]","ar_fh1[16]","sigma_j_fh1[16]",
            "ar_fh1[17]", "sigma_j_fh1[17]","ar_fh1[18]", "sigma_j_fh1[18]","ar_fh1[19]", "sigma_j_fh1[19]","ar_fh1[20]", "sigma_j_fh1[20]",
            "ar_fh1[21]", "sigma_j_fh1[21]","ar_fh1[22]", "sigma_j_fh1[22]","ar_fh1[23]", "sigma_j_fh1[23]","ar_fh1[24]", "sigma_j_fh1[24]",
            "ar_fh1[25]", "sigma_j_fh1[25]","ar_fh1[26]", "sigma_j_fh1[26]","ar_fh1[27]", "sigma_j_fh1[27]","ar_fh1[28]", "sigma_j_fh1[28]",
            "ar_fh1[29]", "sigma_j_fh1[29]","ar_fh1[30]", "sigma_j_fh1[30]","ar_fh1[31]", "sigma_j_fh1[31]","ar_fh1[32]", "sigma_j_fh1[32]",
            "ar_fh1[33]", "sigma_j_fh1[33]","ar_fh1[34]", "sigma_j_fh1[34]","ar_fh1[35]", "sigma_j_fh1[35]","ar_fh1[36]", "sigma_j_fh1[36]",
            "ar_fh1[37]", "sigma_j_fh1[37]","ar_fh1[38]", "sigma_j_fh1[38]","ar_fh1[39]", "sigma_j_fh1[39]","ar_fh1[40]", "sigma_j_fh1[40]",
            "ar_fh2[1]","sigma_j_fh2[1]","ar_fh2[2]", "sigma_j_fh2[2]", "ar_fh2[3]", "sigma_j_fh2[3]", "ar_fh2[4]",  "sigma_j_fh2[4]",
            "ar_fh2[5]","sigma_j_fh2[5]","ar_fh2[6]", "sigma_j_fh2[6]", "ar_fh2[7]", "sigma_j_fh2[7]", "ar_fh2[8]",  "sigma_j_fh2[8]",
            "ar_fh2[9]","sigma_j_fh2[9]","ar_fh2[10]","sigma_j_fh2[10]","ar_fh2[11]","sigma_j_fh2[11]","ar_fh2[12]", "sigma_j_fh2[12]",
            "ar_fh2[13]", "sigma_j_fh2[13]","ar_fh2[14]", "sigma_j_fh2[14]","ar_fh2[15]","sigma_j_fh2[15]","ar_fh2[16]",  "sigma_j_fh2[16]",
            "ar_fh2[17]", "sigma_j_fh2[17]","ar_fh2[18]", "sigma_j_fh2[18]","ar_fh2[19]", "sigma_j_fh2[19]","ar_fh2[20]", "sigma_j_fh2[20]",
            "ar_fh2[21]", "sigma_j_fh2[21]","ar_fh2[22]", "sigma_j_fh2[22]","ar_fh2[23]", "sigma_j_fh2[23]","ar_fh2[24]", "sigma_j_fh2[24]",
            "ar_fh2[25]", "sigma_j_fh2[25]","ar_fh2[26]", "sigma_j_fh2[26]","ar_fh2[27]", "sigma_j_fh2[27]","ar_fh2[28]", "sigma_j_fh2[28]",
            "ar_fh2[29]", "sigma_j_fh2[29]","ar_fh2[30]", "sigma_j_fh2[30]","ar_fh2[31]", "sigma_j_fh2[31]","ar_fh2[32]", "sigma_j_fh2[32]",
            "ar_fh2[33]", "sigma_j_fh2[33]","ar_fh2[34]", "sigma_j_fh2[34]","ar_fh2[35]", "sigma_j_fh2[35]","ar_fh2[36]", "sigma_j_fh2[36]",
            "ar_fh2[37]", "sigma_j_fh2[37]","ar_fh2[38]", "sigma_j_fh2[38]","ar_fh2[39]", "sigma_j_fh2[39]","ar_fh2[40]", "sigma_j_fh2[40]",
            "ar_fh3[1]","sigma_j_fh3[1]","ar_fh3[2]", "sigma_j_fh3[2]", "ar_fh3[3]",  "sigma_j_fh3[3]", "ar_fh3[4]",  "sigma_j_fh3[4]",
            "ar_fh3[5]","sigma_j_fh3[5]","ar_fh3[6]", "sigma_j_fh3[6]", "ar_fh3[7]",  "sigma_j_fh3[7]", "ar_fh3[8]",  "sigma_j_fh3[8]",
            "ar_fh3[9]","sigma_j_fh3[9]","ar_fh3[10]","sigma_j_fh3[10]","ar_fh3[11]", "sigma_j_fh3[11]","ar_fh3[12]", "sigma_j_fh3[12]",
            "ar_fh3[13]", "sigma_j_fh3[13]","ar_fh3[14]", "sigma_j_fh3[14]","ar_fh3[15]", "sigma_j_fh3[15]","ar_fh3[16]","sigma_j_fh3[16]",
            "ar_fh3[17]", "sigma_j_fh3[17]","ar_fh3[18]", "sigma_j_fh3[18]","ar_fh3[19]", "sigma_j_fh3[19]","ar_fh3[20]","sigma_j_fh3[20]",
            "ar_fh3[21]", "sigma_j_fh3[21]","ar_fh3[22]", "sigma_j_fh3[22]","ar_fh3[23]", "sigma_j_fh3[23]","ar_fh3[24]","sigma_j_fh3[24]",
            "ar_fh3[25]", "sigma_j_fh3[25]","ar_fh3[26]", "sigma_j_fh3[26]","ar_fh3[27]", "sigma_j_fh3[27]","ar_fh3[28]","sigma_j_fh3[28]",
            "ar_fh3[29]", "sigma_j_fh3[29]","ar_fh3[30]", "sigma_j_fh3[30]","ar_fh3[31]", "sigma_j_fh3[31]","ar_fh3[32]","sigma_j_fh3[32]",
            "ar_fh3[33]", "sigma_j_fh3[33]","ar_fh3[34]", "sigma_j_fh3[34]","ar_fh3[35]", "sigma_j_fh3[35]","ar_fh3[36]","sigma_j_fh3[36]",
            "ar_fh3[37]", "sigma_j_fh3[37]","ar_fh3[38]", "sigma_j_fh3[38]","ar_fh3[39]", "sigma_j_fh3[39]","ar_fh3[40]","sigma_j_fh3[40]",
            "ar_fh4[1]","sigma_j_fh4[1]","ar_fh4[2]",  "sigma_j_fh4[2]", "ar_fh4[3]",  "sigma_j_fh4[3]", "ar_fh4[4]","sigma_j_fh4[4]",
            "ar_fh4[5]","sigma_j_fh4[5]","ar_fh4[6]",  "sigma_j_fh4[6]", "ar_fh4[7]",  "sigma_j_fh4[7]", "ar_fh4[8]", "sigma_j_fh4[8]",
            "ar_fh4[9]","sigma_j_fh4[9]","ar_fh4[10]", "sigma_j_fh4[10]","ar_fh4[11]", "sigma_j_fh4[11]","ar_fh4[12]", "sigma_j_fh4[12]",
            "ar_fh4[13]", "sigma_j_fh4[13]","ar_fh4[14]", "sigma_j_fh4[14]","ar_fh4[15]", "sigma_j_fh4[15]","ar_fh4[16]","sigma_j_fh4[16]",
            "ar_fh4[17]", "sigma_j_fh4[17]","ar_fh4[18]", "sigma_j_fh4[18]","ar_fh4[19]", "sigma_j_fh4[19]","ar_fh4[20]","sigma_j_fh4[20]",
            "ar_fh4[21]", "sigma_j_fh4[21]","ar_fh4[22]", "sigma_j_fh4[22]","ar_fh4[23]", "sigma_j_fh4[23]","ar_fh4[24]","sigma_j_fh4[24]",
            "ar_fh4[25]", "sigma_j_fh4[25]","ar_fh4[26]", "sigma_j_fh4[26]","ar_fh4[27]", "sigma_j_fh4[27]","ar_fh4[28]","sigma_j_fh4[28]",
            "ar_fh4[29]", "sigma_j_fh4[29]","ar_fh4[30]", "sigma_j_fh4[30]","ar_fh4[31]", "sigma_j_fh4[31]","ar_fh4[32]","sigma_j_fh4[32]",
            "ar_fh4[33]", "sigma_j_fh4[33]","ar_fh4[34]", "sigma_j_fh4[34]","ar_fh4[35]", "sigma_j_fh4[35]","ar_fh4[36]","sigma_j_fh4[36]",
            "ar_fh4[37]", "sigma_j_fh4[37]","ar_fh4[38]", "sigma_j_fh4[38]","ar_fh4[39]", "sigma_j_fh4[39]","ar_fh4[40]","sigma_j_fh4[40]",
            "ar_fh5[1]","sigma_j_fh5[1]","ar_fh5[2]", "sigma_j_fh5[2]","ar_fh5[3]", "sigma_j_fh5[3]","ar_fh5[4]","sigma_j_fh5[4]",
            "ar_fh5[5]","sigma_j_fh5[5]","ar_fh5[6]", "sigma_j_fh5[6]","ar_fh5[7]", "sigma_j_fh5[7]","ar_fh5[8]", "sigma_j_fh5[8]",
            "ar_fh5[9]","sigma_j_fh5[9]","ar_fh5[10]", "sigma_j_fh5[10]","ar_fh5[11]", "sigma_j_fh5[11]","ar_fh5[12]", "sigma_j_fh5[12]",
            "ar_fh5[13]", "sigma_j_fh5[13]","ar_fh5[14]", "sigma_j_fh5[14]","ar_fh5[15]", "sigma_j_fh5[15]","ar_fh5[16]","sigma_j_fh5[16]",
            "ar_fh5[17]", "sigma_j_fh5[17]","ar_fh5[18]", "sigma_j_fh5[18]","ar_fh5[19]", "sigma_j_fh5[19]","ar_fh5[20]","sigma_j_fh5[20]",
            "ar_fh5[21]", "sigma_j_fh5[21]","ar_fh5[22]", "sigma_j_fh5[22]","ar_fh5[23]", "sigma_j_fh5[23]","ar_fh5[24]","sigma_j_fh5[24]",
            "ar_fh5[25]", "sigma_j_fh5[25]","ar_fh5[26]", "sigma_j_fh5[26]","ar_fh5[27]", "sigma_j_fh5[27]","ar_fh5[28]","sigma_j_fh5[28]",
            "ar_fh5[29]", "sigma_j_fh5[29]","ar_fh5[30]", "sigma_j_fh5[30]","ar_fh5[31]", "sigma_j_fh5[31]","ar_fh5[32]","sigma_j_fh5[32]",
            "ar_fh5[33]", "sigma_j_fh5[33]","ar_fh5[34]", "sigma_j_fh5[34]","ar_fh5[35]", "sigma_j_fh5[35]","ar_fh5[36]","sigma_j_fh5[36]",
            "ar_fh5[37]", "sigma_j_fh5[37]","ar_fh5[38]", "sigma_j_fh5[38]","ar_fh5[39]", "sigma_j_fh5[39]","ar_fh5[40]","sigma_j_fh5[40]",
            "ar_fh6[1]","sigma_j_fh6[1]","ar_fh6[2]", "sigma_j_fh6[2]", "ar_fh6[3]", "sigma_j_fh6[3]","ar_fh6[4]","sigma_j_fh6[4]",
            "ar_fh6[5]","sigma_j_fh6[5]","ar_fh6[6]", "sigma_j_fh6[6]", "ar_fh6[7]", "sigma_j_fh6[7]","ar_fh6[8]", "sigma_j_fh6[8]",
            "ar_fh6[9]","sigma_j_fh6[9]","ar_fh6[10]","sigma_j_fh6[10]","ar_fh6[11]", "sigma_j_fh6[11]","ar_fh6[12]", "sigma_j_fh6[12]",
            "ar_fh6[13]", "sigma_j_fh6[13]","ar_fh6[14]", "sigma_j_fh6[14]","ar_fh6[15]", "sigma_j_fh6[15]","ar_fh6[16]","sigma_j_fh6[16]",
            "ar_fh6[17]", "sigma_j_fh6[17]","ar_fh6[18]", "sigma_j_fh6[18]","ar_fh6[19]", "sigma_j_fh6[19]","ar_fh6[20]","sigma_j_fh6[20]",
            "ar_fh6[21]", "sigma_j_fh6[21]","ar_fh6[22]", "sigma_j_fh6[22]","ar_fh6[23]", "sigma_j_fh6[23]","ar_fh6[24]","sigma_j_fh6[24]",
            "ar_fh6[25]", "sigma_j_fh6[25]","ar_fh6[26]", "sigma_j_fh6[26]","ar_fh6[27]", "sigma_j_fh6[27]","ar_fh6[28]","sigma_j_fh6[28]",
            "ar_fh6[29]", "sigma_j_fh6[29]","ar_fh6[30]", "sigma_j_fh6[30]","ar_fh6[31]", "sigma_j_fh6[31]","ar_fh6[32]","sigma_j_fh6[32]",
            "ar_fh6[33]", "sigma_j_fh6[33]","ar_fh6[34]", "sigma_j_fh6[34]","ar_fh6[35]", "sigma_j_fh6[35]","ar_fh6[36]","sigma_j_fh6[36]",
            "ar_fh6[37]", "sigma_j_fh6[37]","ar_fh6[38]", "sigma_j_fh6[38]","ar_fh6[39]", "sigma_j_fh6[39]","ar_fh6[40]","sigma_j_fh6[40]",
            "ar_fv1[1]","sigma_j_fv1[1]","ar_fv1[2]", "sigma_j_fv1[2]","ar_fv1[3]", "sigma_j_fv1[3]","ar_fv1[4]","sigma_j_fv1[4]",
            "ar_fv1[5]", "sigma_j_fv1[5]","ar_fv1[6]", "sigma_j_fv1[6]","ar_fv1[7]", "sigma_j_fv1[7]","ar_fv1[8]", "sigma_j_fv1[8]",
            "ar_fv1[9]", "sigma_j_fv1[9]","ar_fv1[10]", "sigma_j_fv1[10]","ar_fv1[11]", "sigma_j_fv1[11]","ar_fv1[12]", "sigma_j_fv1[12]",
            "ar_fv1[13]", "sigma_j_fv1[13]","ar_fv1[14]", "sigma_j_fv1[14]","ar_fv1[15]", "sigma_j_fv1[15]","ar_fv1[16]","sigma_j_fv1[16]",
            "ar_fv1[17]", "sigma_j_fv1[17]","ar_fv1[18]", "sigma_j_fv1[18]","ar_fv1[19]", "sigma_j_fv1[19]","ar_fv1[20]", "sigma_j_fv1[20]",
            "ar_fv1[21]", "sigma_j_fv1[21]","ar_fv1[22]", "sigma_j_fv1[22]","ar_fv1[23]", "sigma_j_fv1[23]","ar_fv1[24]", "sigma_j_fv1[24]",
            "ar_fv1[25]", "sigma_j_fv1[25]","ar_fv1[26]", "sigma_j_fv1[26]","ar_fv1[27]", "sigma_j_fv1[27]","ar_fv1[28]", "sigma_j_fv1[28]",
            "ar_fv1[29]", "sigma_j_fv1[29]","ar_fv1[30]", "sigma_j_fv1[30]","ar_fv1[31]", "sigma_j_fv1[31]","ar_fv1[32]", "sigma_j_fv1[32]",
            "ar_fv1[33]", "sigma_j_fv1[33]","ar_fv1[34]", "sigma_j_fv1[34]","ar_fv1[35]", "sigma_j_fv1[35]","ar_fv1[36]", "sigma_j_fv1[36]",
            "ar_fv1[37]", "sigma_j_fv1[37]","ar_fv1[38]", "sigma_j_fv1[38]","ar_fv1[39]", "sigma_j_fv1[39]","ar_fv1[40]", "sigma_j_fv1[40]",
            "ar_fv2[1]","sigma_j_fv2[1]","ar_fv2[2]", "sigma_j_fv2[2]","ar_fv2[3]", "sigma_j_fv2[3]","ar_fv2[4]","sigma_j_fv2[4]",
            "ar_fv2[5]", "sigma_j_fv2[5]","ar_fv2[6]", "sigma_j_fv2[6]","ar_fv2[7]", "sigma_j_fv2[7]","ar_fv2[8]", "sigma_j_fv2[8]",
            "ar_fv2[9]", "sigma_j_fv2[9]","ar_fv2[10]", "sigma_j_fv2[10]","ar_fv2[11]", "sigma_j_fv2[11]","ar_fv2[12]", "sigma_j_fv2[12]",
            "ar_fv2[13]", "sigma_j_fv2[13]","ar_fv2[14]", "sigma_j_fv2[14]", "ar_fv2[15]","sigma_j_fv2[15]","ar_fv2[16]","sigma_j_fv2[16]",
            "ar_fv2[17]", "sigma_j_fv2[17]","ar_fv2[18]", "sigma_j_fv2[18]","ar_fv2[19]", "sigma_j_fv2[19]","ar_fv2[20]", "sigma_j_fv2[20]",
            "ar_fv2[21]", "sigma_j_fv2[21]","ar_fv2[22]", "sigma_j_fv2[22]","ar_fv2[23]", "sigma_j_fv2[23]","ar_fv2[24]", "sigma_j_fv2[24]",
            "ar_fv2[25]", "sigma_j_fv2[25]","ar_fv2[26]", "sigma_j_fv2[26]","ar_fv2[27]", "sigma_j_fv2[27]","ar_fv2[28]", "sigma_j_fv2[28]",
            "ar_fv2[29]", "sigma_j_fv2[29]","ar_fv2[30]", "sigma_j_fv2[30]","ar_fv2[31]", "sigma_j_fv2[31]","ar_fv2[32]", "sigma_j_fv2[32]",
            "ar_fv2[33]", "sigma_j_fv2[33]","ar_fv2[34]", "sigma_j_fv2[34]","ar_fv2[35]", "sigma_j_fv2[35]","ar_fv2[36]", "sigma_j_fv2[36]",
            "ar_fv2[37]", "sigma_j_fv2[37]","ar_fv2[38]", "sigma_j_fv2[38]","ar_fv2[39]", "sigma_j_fv2[39]","ar_fv2[40]", "sigma_j_fv2[40]",
            "ar_fv3[1]","sigma_j_fv3[1]","ar_fv3[2]", "sigma_j_fv3[2]","ar_fv3[3]", "sigma_j_fv3[3]","ar_fv3[4]","sigma_j_fv3[4]",
            "ar_fv3[5]", "sigma_j_fv3[5]","ar_fv3[6]", "sigma_j_fv3[6]","ar_fv3[7]", "sigma_j_fv3[7]","ar_fv3[8]", "sigma_j_fv3[8]",
            "ar_fv3[9]", "sigma_j_fv3[9]","ar_fv3[10]","sigma_j_fv3[10]", "ar_fv3[11]", "sigma_j_fv3[11]","ar_fv3[12]","sigma_j_fv3[12]",
            "ar_fv3[13]", "sigma_j_fv3[13]","ar_fv3[14]", "sigma_j_fv3[14]","ar_fv3[15]", "sigma_j_fv3[15]","ar_fv3[16]","sigma_j_fv3[16]",
            "ar_fv3[17]", "sigma_j_fv3[17]","ar_fv3[18]", "sigma_j_fv3[18]","ar_fv3[19]", "sigma_j_fv3[19]","ar_fv3[20]","sigma_j_fv3[20]",
            "ar_fv3[21]", "sigma_j_fv3[21]","ar_fv3[22]", "sigma_j_fv3[22]","ar_fv3[23]", "sigma_j_fv3[23]","ar_fv3[24]","sigma_j_fv3[24]",
            "ar_fv3[25]", "sigma_j_fv3[25]","ar_fv3[26]", "sigma_j_fv3[26]","ar_fv3[27]", "sigma_j_fv3[27]","ar_fv3[28]","sigma_j_fv3[28]",
            "ar_fv3[29]", "sigma_j_fv3[29]","ar_fv3[30]", "sigma_j_fv3[30]","ar_fv3[31]", "sigma_j_fv3[31]","ar_fv3[32]","sigma_j_fv3[32]",
            "ar_fv3[33]", "sigma_j_fv3[33]","ar_fv3[34]", "sigma_j_fv3[34]","ar_fv3[35]", "sigma_j_fv3[35]","ar_fv3[36]","sigma_j_fv3[36]",
            "ar_fv3[37]", "sigma_j_fv3[37]","ar_fv3[38]", "sigma_j_fv3[38]","ar_fv3[39]", "sigma_j_fv3[39]","ar_fv3[40]","sigma_j_fv3[40]",
            "ar_fv4[1]","sigma_j_fv4[1]","ar_fv4[2]", "sigma_j_fv4[2]","ar_fv4[3]", "sigma_j_fv4[3]","ar_fv4[4]","sigma_j_fv4[4]",
            "ar_fv4[5]", "sigma_j_fv4[5]","ar_fv4[6]", "sigma_j_fv4[6]","ar_fv4[7]", "sigma_j_fv4[7]","ar_fv4[8]", "sigma_j_fv4[8]",
            "ar_fv4[9]", "sigma_j_fv4[9]","ar_fv4[10]", "sigma_j_fv4[10]","ar_fv4[11]", "sigma_j_fv4[11]","ar_fv4[12]", "sigma_j_fv4[12]",
            "ar_fv4[13]", "sigma_j_fv4[13]","ar_fv4[14]", "sigma_j_fv4[14]", "ar_fv4[15]","sigma_j_fv4[15]","ar_fv4[16]","sigma_j_fv4[16]",
            "ar_fv4[17]", "sigma_j_fv4[17]","ar_fv4[18]", "sigma_j_fv4[18]","ar_fv4[19]", "sigma_j_fv4[19]","ar_fv4[20]","sigma_j_fv4[20]",
            "ar_fv4[21]", "sigma_j_fv4[21]","ar_fv4[22]", "sigma_j_fv4[22]","ar_fv4[23]", "sigma_j_fv4[23]","ar_fv4[24]","sigma_j_fv4[24]",
            "ar_fv4[25]", "sigma_j_fv4[25]","ar_fv4[26]", "sigma_j_fv4[26]","ar_fv4[27]", "sigma_j_fv4[27]","ar_fv4[28]","sigma_j_fv4[28]",
            "ar_fv4[29]", "sigma_j_fv4[29]","ar_fv4[30]", "sigma_j_fv4[30]","ar_fv4[31]", "sigma_j_fv4[31]","ar_fv4[32]","sigma_j_fv4[32]",
            "ar_fv4[33]", "sigma_j_fv4[33]","ar_fv4[34]", "sigma_j_fv4[34]","ar_fv4[35]", "sigma_j_fv4[35]","ar_fv4[36]","sigma_j_fv4[36]",
            "ar_fv4[37]", "sigma_j_fv4[37]","ar_fv4[38]", "sigma_j_fv4[38]","ar_fv4[39]", "sigma_j_fv4[39]","ar_fv4[40]","sigma_j_fv4[40]",
            "ar_fv5[1]","sigma_j_fv5[1]","ar_fv5[2]", "sigma_j_fv5[2]","ar_fv5[3]", "sigma_j_fv5[3]","ar_fv5[4]","sigma_j_fv5[4]",
            "ar_fv5[5]", "sigma_j_fv5[5]","ar_fv5[6]", "sigma_j_fv5[6]","ar_fv5[7]", "sigma_j_fv5[7]","ar_fv5[8]", "sigma_j_fv5[8]",
            "ar_fv5[9]", "sigma_j_fv5[9]","ar_fv5[10]", "sigma_j_fv5[10]","ar_fv5[11]", "sigma_j_fv5[11]","ar_fv5[12]", "sigma_j_fv5[12]",
            "ar_fv5[13]", "sigma_j_fv5[13]","ar_fv5[14]", "sigma_j_fv5[14]","ar_fv5[15]", "sigma_j_fv5[15]","ar_fv5[16]","sigma_j_fv5[16]",
            "ar_fv5[17]", "sigma_j_fv5[17]","ar_fv5[18]", "sigma_j_fv5[18]","ar_fv5[19]", "sigma_j_fv5[19]","ar_fv5[20]","sigma_j_fv5[20]",
            "ar_fv5[21]", "sigma_j_fv5[21]","ar_fv5[22]", "sigma_j_fv5[22]","ar_fv5[23]", "sigma_j_fv5[23]","ar_fv5[24]","sigma_j_fv5[24]",
            "ar_fv5[25]", "sigma_j_fv5[25]","ar_fv5[26]", "sigma_j_fv5[26]","ar_fv5[27]", "sigma_j_fv5[27]","ar_fv5[28]","sigma_j_fv5[28]",
            "ar_fv5[29]", "sigma_j_fv5[29]","ar_fv5[30]", "sigma_j_fv5[30]","ar_fv5[31]", "sigma_j_fv5[31]","ar_fv5[32]","sigma_j_fv5[32]",
            "ar_fv5[33]", "sigma_j_fv5[33]","ar_fv5[34]", "sigma_j_fv5[34]","ar_fv5[35]", "sigma_j_fv5[35]","ar_fv5[36]","sigma_j_fv5[36]",
            "ar_fv5[37]", "sigma_j_fv5[37]","ar_fv5[38]", "sigma_j_fv5[38]","ar_fv5[39]", "sigma_j_fv5[39]","ar_fv5[40]","sigma_j_fv5[40]",
            "ar_fg1[1]","sigma_j_fg1[1]","ar_fg1[2]", "sigma_j_fg1[2]","ar_fg1[3]", "sigma_j_fg1[3]","ar_fg1[4]","sigma_j_fg1[4]",
            "ar_fg1[5]", "sigma_j_fg1[5]","ar_fg1[6]", "sigma_j_fg1[6]","ar_fg1[7]", "sigma_j_fg1[7]","ar_fg1[8]", "sigma_j_fg1[8]",
            "ar_fg1[9]", "sigma_j_fg1[9]","ar_fg1[10]", "sigma_j_fg1[10]","ar_fg1[11]", "sigma_j_fg1[11]","ar_fg1[12]", "sigma_j_fg1[12]",
            "ar_fg1[13]", "sigma_j_fg1[13]","ar_fg1[14]", "sigma_j_fg1[14]","ar_fg1[15]", "sigma_j_fg1[15]","ar_fg1[16]","sigma_j_fg1[16]",
            "ar_fg1[17]", "sigma_j_fg1[17]","ar_fg1[18]", "sigma_j_fg1[18]","ar_fg1[19]", "sigma_j_fg1[19]","ar_fg1[20]", "sigma_j_fg1[20]",
            "ar_fg1[21]", "sigma_j_fg1[21]","ar_fg1[22]", "sigma_j_fg1[22]","ar_fg1[23]", "sigma_j_fg1[23]","ar_fg1[24]", "sigma_j_fg1[24]",
            "ar_fg1[25]", "sigma_j_fg1[25]","ar_fg1[26]", "sigma_j_fg1[26]","ar_fg1[27]", "sigma_j_fg1[27]","ar_fg1[28]", "sigma_j_fg1[28]",
            "ar_fg1[29]", "sigma_j_fg1[29]","ar_fg1[30]", "sigma_j_fg1[30]","ar_fg1[31]", "sigma_j_fg1[31]","ar_fg1[32]", "sigma_j_fg1[32]",
            "ar_fg1[33]", "sigma_j_fg1[33]","ar_fg1[34]", "sigma_j_fg1[34]","ar_fg1[35]", "sigma_j_fg1[35]","ar_fg1[36]", "sigma_j_fg1[36]",
            "ar_fg1[37]", "sigma_j_fg1[37]","ar_fg1[38]", "sigma_j_fg1[38]","ar_fg1[39]", "sigma_j_fg1[39]","ar_fg1[40]", "sigma_j_fg1[40]",
            "alpha_s1", "beta_s1", "sigma_s1",
            "ar_s1[1]","sigma_j_s1[1]","ar_s1[2]", "sigma_j_s1[2]","ar_s1[3]", "sigma_j_s1[3]","ar_s1[4]","sigma_j_s1[4]",
            "ar_s1[5]", "sigma_j_s1[5]","ar_s1[6]", "sigma_j_s1[6]","ar_s1[7]", "sigma_j_s1[7]","ar_s1[8]", "sigma_j_s1[8]",
            "ar_s1[9]", "sigma_j_s1[9]","ar_s1[10]", "sigma_j_s1[10]","ar_s1[11]", "sigma_j_s1[11]","ar_s1[12]", "sigma_j_s1[12]",
            "ar_s1[13]", "sigma_j_s1[13]","ar_s1[14]", "sigma_j_s1[14]", "ar_s1[15]", "sigma_j_s1[15]","ar_s1[16]","sigma_j_s1[16]",
            "ar_s1[17]", "sigma_j_s1[17]","ar_s1[18]", "sigma_j_s1[18]","ar_s1[19]", "sigma_j_s1[19]","ar_s1[20]", "sigma_j_s1[20]",
            "ar_s1[21]", "sigma_j_s1[21]","ar_s1[22]", "sigma_j_s1[22]","ar_s1[23]", "sigma_j_s1[23]","ar_s1[24]", "sigma_j_s1[24]",
            "ar_s1[25]", "sigma_j_s1[25]","ar_s1[26]", "sigma_j_s1[26]","ar_s1[27]", "sigma_j_s1[27]","ar_s1[28]", "sigma_j_s1[28]",
            "ar_s1[29]", "sigma_j_s1[29]","ar_s1[30]", "sigma_j_s1[30]","ar_s1[31]", "sigma_j_s1[31]","ar_s1[32]", "sigma_j_s1[32]",
            "ar_s1[33]", "sigma_j_s1[33]","ar_s1[34]", "sigma_j_s1[34]","ar_s1[35]", "sigma_j_s1[35]","ar_s1[36]", "sigma_j_s1[36]",
            "ar_s1[37]", "sigma_j_s1[37]","ar_s1[38]", "sigma_j_s1[38]","ar_s1[39]", "sigma_j_s1[39]","ar_s1[40]", "sigma_j_s1[40]",
            "alpha_s2", "beta_s2", "sigma_s2",
            "ar_s2[1]","sigma_j_s2[1]","ar_s2[2]", "sigma_j_s2[2]","ar_s2[3]", "sigma_j_s2[3]","ar_s2[4]","sigma_j_s2[4]",
            "ar_s2[5]", "sigma_j_s2[5]","ar_s2[6]", "sigma_j_s2[6]","ar_s2[7]", "sigma_j_s2[7]","ar_s2[8]", "sigma_j_s2[8]",
            "ar_s2[9]", "sigma_j_s2[9]","ar_s2[10]", "sigma_j_s2[10]","ar_s2[11]", "sigma_j_s2[11]","ar_s2[12]", "sigma_j_s2[12]",
            "ar_s2[13]", "sigma_j_s2[13]","ar_s2[14]", "sigma_j_s2[14]", "ar_s2[15]", "sigma_j_s2[15]","ar_s2[16]","sigma_j_s2[16]",
            "ar_s2[17]", "sigma_j_s2[17]","ar_s2[18]", "sigma_j_s2[18]","ar_s2[19]", "sigma_j_s2[19]","ar_s2[20]", "sigma_j_s2[20]",
            "ar_s2[21]", "sigma_j_s2[21]","ar_s2[22]", "sigma_j_s2[22]","ar_s2[23]", "sigma_j_s2[23]","ar_s2[24]", "sigma_j_s2[24]",
            "ar_s2[25]", "sigma_j_s2[25]","ar_s2[26]", "sigma_j_s2[26]","ar_s2[27]", "sigma_j_s2[27]","ar_s2[28]", "sigma_j_s2[28]",
            "ar_s2[29]", "sigma_j_s2[29]","ar_s2[30]", "sigma_j_s2[30]","ar_s2[31]", "sigma_j_s2[31]","ar_s2[32]", "sigma_j_s2[32]",
            "ar_s2[33]", "sigma_j_s2[33]","ar_s2[34]", "sigma_j_s2[34]","ar_s2[35]", "sigma_j_s2[35]","ar_s2[36]", "sigma_j_s2[36]",
            "ar_s2[37]", "sigma_j_s2[37]","ar_s2[38]", "sigma_j_s2[38]","ar_s2[39]", "sigma_j_s2[39]","ar_s2[40]", "sigma_j_s2[40]",
            "alpha_s3", "beta_s3", "sigma_s3",
            "ar_s3[1]","sigma_j_s3[1]","ar_s3[2]", "sigma_j_s3[2]","ar_s3[3]", "sigma_j_s3[3]","ar_s3[4]","sigma_j_s3[4]",
            "ar_s3[5]", "sigma_j_s3[5]","ar_s3[6]", "sigma_j_s3[6]","ar_s3[7]", "sigma_j_s3[7]","ar_s3[8]", "sigma_j_s3[8]",
            "ar_s3[9]", "sigma_j_s3[9]","ar_s3[10]", "sigma_j_s3[10]","ar_s3[11]", "sigma_j_s3[11]","ar_s3[12]", "sigma_j_s3[12]",
            "ar_s3[13]", "sigma_j_s3[13]","ar_s3[14]", "sigma_j_s3[14]", "ar_s3[15]", "sigma_j_s3[15]","ar_s3[16]","sigma_j_s3[16]",
            "ar_s3[17]", "sigma_j_s3[17]","ar_s3[18]", "sigma_j_s3[18]","ar_s3[19]", "sigma_j_s3[19]","ar_s3[20]", "sigma_j_s3[20]",
            "ar_s3[21]", "sigma_j_s3[21]","ar_s3[22]", "sigma_j_s3[22]","ar_s3[23]", "sigma_j_s3[23]","ar_s3[24]", "sigma_j_s3[24]",
            "ar_s3[25]", "sigma_j_s3[25]","ar_s3[26]", "sigma_j_s3[26]","ar_s3[27]", "sigma_j_s3[27]","ar_s3[28]", "sigma_j_s3[28]",
            "ar_s3[29]", "sigma_j_s3[29]","ar_s3[30]", "sigma_j_s3[30]","ar_s3[31]", "sigma_j_s3[31]","ar_s3[32]", "sigma_j_s3[32]",
            "ar_s3[33]", "sigma_j_s3[33]","ar_s3[34]", "sigma_j_s3[34]","ar_s3[35]", "sigma_j_s3[35]","ar_s3[36]", "sigma_j_s3[36]",
            "ar_s3[37]", "sigma_j_s3[37]","ar_s3[38]", "sigma_j_s3[38]","ar_s3[39]", "sigma_j_s3[39]","ar_s3[40]", "sigma_j_s3[40]",
            "alpha_s4", "beta_s4", "sigma_s4",
            "ar_s4[1]","sigma_j_s4[1]","ar_s4[2]", "sigma_j_s4[2]","ar_s4[3]", "sigma_j_s4[3]","ar_s4[4]","sigma_j_s4[4]",
            "ar_s4[5]", "sigma_j_s4[5]","ar_s4[6]", "sigma_j_s4[6]","ar_s4[7]", "sigma_j_s4[7]","ar_s4[8]", "sigma_j_s4[8]",
            "ar_s4[9]", "sigma_j_s4[9]","ar_s4[10]", "sigma_j_s4[10]","ar_s4[11]", "sigma_j_s4[11]","ar_s4[12]", "sigma_j_s4[12]",
            "ar_s4[13]", "sigma_j_s4[13]","ar_s4[14]", "sigma_j_s4[14]", "ar_s4[15]", "sigma_j_s4[15]","ar_s4[16]","sigma_j_s4[16]",
            "ar_s4[17]", "sigma_j_s4[17]","ar_s4[18]", "sigma_j_s4[18]","ar_s4[19]", "sigma_j_s4[19]","ar_s4[20]", "sigma_j_s4[20]",
            "ar_s4[21]", "sigma_j_s4[21]","ar_s4[22]", "sigma_j_s4[22]","ar_s4[23]", "sigma_j_s4[23]","ar_s4[24]", "sigma_j_s4[24]",
            "ar_s4[25]", "sigma_j_s4[25]","ar_s4[26]", "sigma_j_s4[26]","ar_s4[27]", "sigma_j_s4[27]","ar_s4[28]", "sigma_j_s4[28]",
            "ar_s4[29]", "sigma_j_s4[29]","ar_s4[30]", "sigma_j_s4[30]","ar_s4[31]", "sigma_j_s4[31]","ar_s4[32]", "sigma_j_s4[32]",
            "ar_s4[33]", "sigma_j_s4[33]","ar_s4[34]", "sigma_j_s4[34]","ar_s4[35]", "sigma_j_s4[35]","ar_s4[36]", "sigma_j_s4[36]",
            "ar_s4[37]", "sigma_j_s4[37]","ar_s4[38]", "sigma_j_s4[38]","ar_s4[39]", "sigma_j_s4[39]","ar_s4[40]", "sigma_j_s4[40]",
            "alpha_s5", "beta_s5", "sigma_s5",
            "ar_s5[1]","sigma_j_s5[1]","ar_s5[2]", "sigma_j_s5[2]","ar_s5[3]", "sigma_j_s5[3]","ar_s5[4]","sigma_j_s5[4]",
            "ar_s5[5]", "sigma_j_s5[5]","ar_s5[6]", "sigma_j_s5[6]","ar_s5[7]", "sigma_j_s5[7]","ar_s5[8]", "sigma_j_s5[8]",
            "ar_s5[9]", "sigma_j_s5[9]","ar_s5[10]", "sigma_j_s5[10]","ar_s5[11]", "sigma_j_s5[11]","ar_s5[12]", "sigma_j_s5[12]",
            "ar_s5[13]", "sigma_j_s5[13]","ar_s5[14]", "sigma_j_s5[14]", "ar_s5[15]", "sigma_j_s5[15]","ar_s5[16]","sigma_j_s5[16]",
            "ar_s5[17]", "sigma_j_s5[17]","ar_s5[18]", "sigma_j_s5[18]","ar_s5[19]", "sigma_j_s5[19]","ar_s5[20]", "sigma_j_s5[20]",
            "ar_s5[21]", "sigma_j_s5[21]","ar_s5[22]", "sigma_j_s5[22]","ar_s5[23]", "sigma_j_s5[23]","ar_s5[24]", "sigma_j_s5[24]",
            "ar_s5[25]", "sigma_j_s5[25]","ar_s5[26]", "sigma_j_s5[26]","ar_s5[27]", "sigma_j_s5[27]","ar_s5[28]", "sigma_j_s5[28]",
            "ar_s5[29]", "sigma_j_s5[29]","ar_s5[30]", "sigma_j_s5[30]","ar_s5[31]", "sigma_j_s5[31]","ar_s5[32]", "sigma_j_s5[32]",
            "ar_s5[33]", "sigma_j_s5[33]","ar_s5[34]", "sigma_j_s5[34]","ar_s5[35]", "sigma_j_s5[35]","ar_s5[36]", "sigma_j_s5[36]",
            "ar_s5[37]", "sigma_j_s5[37]","ar_s5[38]", "sigma_j_s5[38]","ar_s5[39]", "sigma_j_s5[39]","ar_s5[40]", "sigma_j_s5[40]",
            "alpha_s6", "beta_s6", "sigma_s6",
            "ar_s6[1]","sigma_j_s6[1]","ar_s6[2]", "sigma_j_s6[2]","ar_s6[3]", "sigma_j_s6[3]","ar_s6[4]","sigma_j_s6[4]",
            "ar_s6[5]", "sigma_j_s6[5]","ar_s6[6]", "sigma_j_s6[6]","ar_s6[7]", "sigma_j_s6[7]","ar_s6[8]", "sigma_j_s6[8]",
            "ar_s6[9]", "sigma_j_s6[9]","ar_s6[10]", "sigma_j_s6[10]","ar_s6[11]", "sigma_j_s6[11]","ar_s6[12]", "sigma_j_s6[12]",
            "ar_s6[13]", "sigma_j_s6[13]","ar_s6[14]", "sigma_j_s6[14]", "ar_s6[15]", "sigma_j_s6[15]","ar_s6[16]","sigma_j_s6[16]",
            "ar_s6[17]", "sigma_j_s6[17]","ar_s6[18]", "sigma_j_s6[18]","ar_s6[19]", "sigma_j_s6[19]","ar_s6[20]", "sigma_j_s6[20]",
            "ar_s6[21]", "sigma_j_s6[21]","ar_s6[22]", "sigma_j_s6[22]","ar_s6[23]", "sigma_j_s6[23]","ar_s6[24]", "sigma_j_s6[24]",
            "ar_s6[25]", "sigma_j_s6[25]","ar_s6[26]", "sigma_j_s6[26]","ar_s6[27]", "sigma_j_s6[27]","ar_s6[28]", "sigma_j_s6[28]",
            "ar_s6[29]", "sigma_j_s6[29]","ar_s6[30]", "sigma_j_s6[30]","ar_s6[31]", "sigma_j_s6[31]","ar_s6[32]", "sigma_j_s6[32]",
            "ar_s6[33]", "sigma_j_s6[33]","ar_s6[34]", "sigma_j_s6[34]","ar_s6[35]", "sigma_j_s6[35]","ar_s6[36]", "sigma_j_s6[36]",
            "ar_s6[37]", "sigma_j_s6[37]","ar_s6[38]", "sigma_j_s6[38]","ar_s6[39]", "sigma_j_s6[39]","ar_s6[40]", "sigma_j_s6[40]",
            "alpha_v1", "beta_v1", "sigma_v1",
            "ar_v1[1]","sigma_j_v1[1]","ar_v1[2]", "sigma_j_v1[2]","ar_v1[3]", "sigma_j_v1[3]","ar_v1[4]","sigma_j_v1[4]",
            "ar_v1[5]", "sigma_j_v1[5]","ar_v1[6]", "sigma_j_v1[6]","ar_v1[7]", "sigma_j_v1[7]","ar_v1[8]", "sigma_j_v1[8]",
            "ar_v1[9]", "sigma_j_v1[9]","ar_v1[10]", "sigma_j_v1[10]","ar_v1[11]", "sigma_j_v1[11]","ar_v1[12]", "sigma_j_v1[12]",
            "ar_v1[13]", "sigma_j_v1[13]","ar_v1[14]", "sigma_j_v1[14]", "ar_v1[15]", "sigma_j_v1[15]","ar_v1[16]","sigma_j_v1[16]",
            "ar_v1[17]", "sigma_j_v1[17]","ar_v1[18]", "sigma_j_v1[18]","ar_v1[19]", "sigma_j_v1[19]","ar_v1[20]", "sigma_j_v1[20]",
            "ar_v1[21]", "sigma_j_v1[21]","ar_v1[22]", "sigma_j_v1[22]","ar_v1[23]", "sigma_j_v1[23]","ar_v1[24]", "sigma_j_v1[24]",
            "ar_v1[25]", "sigma_j_v1[25]","ar_v1[26]", "sigma_j_v1[26]","ar_v1[27]", "sigma_j_v1[27]","ar_v1[28]", "sigma_j_v1[28]",
            "ar_v1[29]", "sigma_j_v1[29]","ar_v1[30]", "sigma_j_v1[30]","ar_v1[31]", "sigma_j_v1[31]","ar_v1[32]", "sigma_j_v1[32]",
            "ar_v1[33]", "sigma_j_v1[33]","ar_v1[34]", "sigma_j_v1[34]","ar_v1[35]", "sigma_j_v1[35]","ar_v1[36]", "sigma_j_v1[36]",
            "ar_v1[37]", "sigma_j_v1[37]","ar_v1[38]", "sigma_j_v1[38]","ar_v1[39]", "sigma_j_v1[39]","ar_v1[40]", "sigma_j_v1[40]",
            "alpha_v2", "beta_v2", "sigma_v2",
            "ar_v2[1]","sigma_j_v2[1]","ar_v2[2]", "sigma_j_v2[2]","ar_v2[3]", "sigma_j_v2[3]","ar_v2[4]","sigma_j_v2[4]",
            "ar_v2[5]", "sigma_j_v2[5]","ar_v2[6]", "sigma_j_v2[6]","ar_v2[7]", "sigma_j_v2[7]","ar_v2[8]", "sigma_j_v2[8]",
            "ar_v2[9]", "sigma_j_v2[9]","ar_v2[10]", "sigma_j_v2[10]","ar_v2[11]", "sigma_j_v2[11]","ar_v2[12]", "sigma_j_v2[12]",
            "ar_v2[13]", "sigma_j_v2[13]","ar_v2[14]", "sigma_j_v2[14]", "ar_v2[15]", "sigma_j_v2[15]","ar_v2[16]","sigma_j_v2[16]",
            "ar_v2[17]", "sigma_j_v2[17]","ar_v2[18]", "sigma_j_v2[18]","ar_v2[19]", "sigma_j_v2[19]","ar_v2[20]", "sigma_j_v2[20]",
            "ar_v2[21]", "sigma_j_v2[21]","ar_v2[22]", "sigma_j_v2[22]","ar_v2[23]", "sigma_j_v2[23]","ar_v2[24]", "sigma_j_v2[24]",
            "ar_v2[25]", "sigma_j_v2[25]","ar_v2[26]", "sigma_j_v2[26]","ar_v2[27]", "sigma_j_v2[27]","ar_v2[28]", "sigma_j_v2[28]",
            "ar_v2[29]", "sigma_j_v2[29]","ar_v2[30]", "sigma_j_v2[30]","ar_v2[31]", "sigma_j_v2[31]","ar_v2[32]", "sigma_j_v2[32]",
            "ar_v2[33]", "sigma_j_v2[33]","ar_v2[34]", "sigma_j_v2[34]","ar_v2[35]", "sigma_j_v2[35]","ar_v2[36]", "sigma_j_v2[36]",
            "ar_v2[37]", "sigma_j_v2[37]","ar_v2[38]", "sigma_j_v2[38]","ar_v2[39]", "sigma_j_v2[39]","ar_v2[40]", "sigma_j_v2[40]",
            "alpha_v3", "beta_v3", "sigma_v3",
            "ar_v3[1]","sigma_j_v3[1]","ar_v3[2]", "sigma_j_v3[2]","ar_v3[3]", "sigma_j_v3[3]","ar_v3[4]","sigma_j_v3[4]",
            "ar_v3[5]", "sigma_j_v3[5]","ar_v3[6]", "sigma_j_v3[6]","ar_v3[7]", "sigma_j_v3[7]","ar_v3[8]", "sigma_j_v3[8]",
            "ar_v3[9]", "sigma_j_v3[9]","ar_v3[10]", "sigma_j_v3[10]","ar_v3[11]", "sigma_j_v3[11]","ar_v3[12]", "sigma_j_v3[12]",
            "ar_v3[13]", "sigma_j_v3[13]","ar_v3[14]", "sigma_j_v3[14]", "ar_v3[15]", "sigma_j_v3[15]","ar_v3[16]","sigma_j_v3[16]",
            "ar_v3[17]", "sigma_j_v3[17]","ar_v3[18]", "sigma_j_v3[18]","ar_v3[19]", "sigma_j_v3[19]","ar_v3[20]", "sigma_j_v3[20]",
            "ar_v3[21]", "sigma_j_v3[21]","ar_v3[22]", "sigma_j_v3[22]","ar_v3[23]", "sigma_j_v3[23]","ar_v3[24]", "sigma_j_v3[24]",
            "ar_v3[25]", "sigma_j_v3[25]","ar_v3[26]", "sigma_j_v3[26]","ar_v3[27]", "sigma_j_v3[27]","ar_v3[28]", "sigma_j_v3[28]",
            "ar_v3[29]", "sigma_j_v3[29]","ar_v3[30]", "sigma_j_v3[30]","ar_v3[31]", "sigma_j_v3[31]","ar_v3[32]", "sigma_j_v3[32]",
            "ar_v3[33]", "sigma_j_v3[33]","ar_v3[34]", "sigma_j_v3[34]","ar_v3[35]", "sigma_j_v3[35]","ar_v3[36]", "sigma_j_v3[36]",
            "ar_v3[37]", "sigma_j_v3[37]","ar_v3[38]", "sigma_j_v3[38]","ar_v3[39]", "sigma_j_v3[39]","ar_v3[40]", "sigma_j_v3[40]",
            "alpha_v4", "beta_v4", "sigma_v4",
            "ar_v4[1]","sigma_j_v4[1]","ar_v4[2]", "sigma_j_v4[2]","ar_v4[3]", "sigma_j_v4[3]","ar_v4[4]","sigma_j_v4[4]",
            "ar_v4[5]", "sigma_j_v4[5]","ar_v4[6]", "sigma_j_v4[6]","ar_v4[7]", "sigma_j_v4[7]","ar_v4[8]", "sigma_j_v4[8]",
            "ar_v4[9]", "sigma_j_v4[9]","ar_v4[10]", "sigma_j_v4[10]","ar_v4[11]", "sigma_j_v4[11]","ar_v4[12]", "sigma_j_v4[12]",
            "ar_v4[13]", "sigma_j_v4[13]","ar_v4[14]", "sigma_j_v4[14]", "ar_v4[15]", "sigma_j_v4[15]","ar_v4[16]","sigma_j_v4[16]",
            "ar_v4[17]", "sigma_j_v4[17]","ar_v4[18]", "sigma_j_v4[18]","ar_v4[19]", "sigma_j_v4[19]","ar_v4[20]", "sigma_j_v4[20]",
            "ar_v4[21]", "sigma_j_v4[21]","ar_v4[22]", "sigma_j_v4[22]","ar_v4[23]", "sigma_j_v4[23]","ar_v4[24]", "sigma_j_v4[24]",
            "ar_v4[25]", "sigma_j_v4[25]","ar_v4[26]", "sigma_j_v4[26]","ar_v4[27]", "sigma_j_v4[27]","ar_v4[28]", "sigma_j_v4[28]",
            "ar_v4[29]", "sigma_j_v4[29]","ar_v4[30]", "sigma_j_v4[30]","ar_v4[31]", "sigma_j_v4[31]","ar_v4[32]", "sigma_j_v4[32]",
            "ar_v4[33]", "sigma_j_v4[33]","ar_v4[34]", "sigma_j_v4[34]","ar_v4[35]", "sigma_j_v4[35]","ar_v4[36]", "sigma_j_v4[36]",
            "ar_v4[37]", "sigma_j_v4[37]","ar_v4[38]", "sigma_j_v4[38]","ar_v4[39]", "sigma_j_v4[39]","ar_v4[40]", "sigma_j_v4[40]",
            "alpha_v5", "beta_v5", "sigma_v5",
            "ar_v5[1]","sigma_j_v5[1]","ar_v5[2]", "sigma_j_v5[2]","ar_v5[3]", "sigma_j_v5[3]","ar_v5[4]","sigma_j_v5[4]",
            "ar_v5[5]", "sigma_j_v5[5]","ar_v5[6]", "sigma_j_v5[6]","ar_v5[7]", "sigma_j_v5[7]","ar_v5[8]", "sigma_j_v5[8]",
            "ar_v5[9]", "sigma_j_v5[9]","ar_v5[10]", "sigma_j_v5[10]","ar_v5[11]", "sigma_j_v5[11]","ar_v5[12]", "sigma_j_v5[12]",
            "ar_v5[13]", "sigma_j_v5[13]","ar_v5[14]", "sigma_j_v5[14]", "ar_v5[15]", "sigma_j_v5[15]","ar_v5[16]","sigma_j_v5[16]",
            "ar_v5[17]", "sigma_j_v5[17]","ar_v5[18]", "sigma_j_v5[18]","ar_v5[19]", "sigma_j_v5[19]","ar_v5[20]", "sigma_j_v5[20]",
            "ar_v5[21]", "sigma_j_v5[21]","ar_v5[22]", "sigma_j_v5[22]","ar_v5[23]", "sigma_j_v5[23]","ar_v5[24]", "sigma_j_v5[24]",
            "ar_v5[25]", "sigma_j_v5[25]","ar_v5[26]", "sigma_j_v5[26]","ar_v5[27]", "sigma_j_v5[27]","ar_v5[28]", "sigma_j_v5[28]",
            "ar_v5[29]", "sigma_j_v5[29]","ar_v5[30]", "sigma_j_v5[30]","ar_v5[31]", "sigma_j_v5[31]","ar_v5[32]", "sigma_j_v5[32]",
            "ar_v5[33]", "sigma_j_v5[33]","ar_v5[34]", "sigma_j_v5[34]","ar_v5[35]", "sigma_j_v5[35]","ar_v5[36]", "sigma_j_v5[36]",
            "ar_v5[37]", "sigma_j_v5[37]","ar_v5[38]", "sigma_j_v5[38]","ar_v5[39]", "sigma_j_v5[39]","ar_v5[40]", "sigma_j_v5[40]",
            "alpha_y1", "beta_y1", "sigma_y1",
            "ar_y1[1]","sigma_j_y1[1]","ar_y1[2]", "sigma_j_y1[2]","ar_y1[3]", "sigma_j_y1[3]","ar_y1[4]","sigma_j_y1[4]",
            "ar_y1[5]", "sigma_j_y1[5]","ar_y1[6]", "sigma_j_y1[6]","ar_y1[7]", "sigma_j_y1[7]","ar_y1[8]", "sigma_j_y1[8]",
            "ar_y1[9]", "sigma_j_y1[9]","ar_y1[10]", "sigma_j_y1[10]","ar_y1[11]", "sigma_j_y1[11]","ar_y1[12]", "sigma_j_y1[12]",
            "ar_y1[13]", "sigma_j_y1[13]","ar_y1[14]", "sigma_j_y1[14]", "ar_y1[15]", "sigma_j_y1[15]","ar_y1[16]","sigma_j_y1[16]",
            "ar_y1[17]", "sigma_j_y1[17]","ar_y1[18]", "sigma_j_y1[18]","ar_y1[19]", "sigma_j_y1[19]","ar_y1[20]", "sigma_j_y1[20]",
            "ar_y1[21]", "sigma_j_y1[21]","ar_y1[22]", "sigma_j_y1[22]","ar_y1[23]", "sigma_j_y1[23]","ar_y1[24]", "sigma_j_y1[24]",
            "ar_y1[25]", "sigma_j_y1[25]","ar_y1[26]", "sigma_j_y1[26]","ar_y1[27]", "sigma_j_y1[27]","ar_y1[28]", "sigma_j_y1[28]",
            "ar_y1[29]", "sigma_j_y1[29]","ar_y1[30]", "sigma_j_y1[30]","ar_y1[31]", "sigma_j_y1[31]","ar_y1[32]", "sigma_j_y1[32]",
            "ar_y1[33]", "sigma_j_y1[33]","ar_y1[34]", "sigma_j_y1[34]","ar_y1[35]", "sigma_j_y1[35]","ar_y1[36]", "sigma_j_y1[36]",
            "ar_y1[37]", "sigma_j_y1[37]","ar_y1[38]", "sigma_j_y1[38]","ar_y1[39]", "sigma_j_y1[39]","ar_y1[40]", "sigma_j_y1[40]",
            "alpha_y2", "beta_y2", "sigma_y2",
            "ar_y2[1]","sigma_j_y2[1]","ar_y2[2]", "sigma_j_y2[2]","ar_y2[3]", "sigma_j_y2[3]","ar_y2[4]","sigma_j_y2[4]",
            "ar_y2[5]", "sigma_j_y2[5]","ar_y2[6]", "sigma_j_y2[6]","ar_y2[7]", "sigma_j_y2[7]","ar_y2[8]", "sigma_j_y2[8]",
            "ar_y2[9]", "sigma_j_y2[9]","ar_y2[10]", "sigma_j_y2[10]","ar_y2[11]", "sigma_j_y2[11]","ar_y2[12]", "sigma_j_y2[12]",
            "ar_y2[13]", "sigma_j_y2[13]","ar_y2[14]", "sigma_j_y2[14]", "ar_y2[15]", "sigma_j_y2[15]","ar_y2[16]","sigma_j_y2[16]",
            "ar_y2[17]", "sigma_j_y2[17]","ar_y2[18]", "sigma_j_y2[18]","ar_y2[19]", "sigma_j_y2[19]","ar_y2[20]", "sigma_j_y2[20]",
            "ar_y2[21]", "sigma_j_y2[21]","ar_y2[22]", "sigma_j_y2[22]","ar_y2[23]", "sigma_j_y2[23]","ar_y2[24]", "sigma_j_y2[24]",
            "ar_y2[25]", "sigma_j_y2[25]","ar_y2[26]", "sigma_j_y2[26]","ar_y2[27]", "sigma_j_y2[27]","ar_y2[28]", "sigma_j_y2[28]",
            "ar_y2[29]", "sigma_j_y2[29]","ar_y2[30]", "sigma_j_y2[30]","ar_y2[31]", "sigma_j_y2[31]","ar_y2[32]", "sigma_j_y2[32]",
            "ar_y2[33]", "sigma_j_y2[33]","ar_y2[34]", "sigma_j_y2[34]","ar_y2[35]", "sigma_j_y2[35]","ar_y2[36]", "sigma_j_y2[36]",
            "ar_y2[37]", "sigma_j_y2[37]","ar_y2[38]", "sigma_j_y2[38]","ar_y2[39]", "sigma_j_y2[39]","ar_y2[40]", "sigma_j_y2[40]",
            "alpha_y3", "beta_y3", "sigma_y3",
            "ar_y3[1]","sigma_j_y3[1]","ar_y3[2]", "sigma_j_y3[2]","ar_y3[3]", "sigma_j_y3[3]","ar_y3[4]","sigma_j_y3[4]",
            "ar_y3[5]", "sigma_j_y3[5]","ar_y3[6]", "sigma_j_y3[6]","ar_y3[7]", "sigma_j_y3[7]","ar_y3[8]", "sigma_j_y3[8]",
            "ar_y3[9]", "sigma_j_y3[9]","ar_y3[10]", "sigma_j_y3[10]","ar_y3[11]", "sigma_j_y3[11]","ar_y3[12]", "sigma_j_y3[12]",
            "ar_y3[13]", "sigma_j_y3[13]","ar_y3[14]", "sigma_j_y3[14]", "ar_y3[15]", "sigma_j_y3[15]","ar_y3[16]","sigma_j_y3[16]",
            "ar_y3[17]", "sigma_j_y3[17]","ar_y3[18]", "sigma_j_y3[18]","ar_y3[19]", "sigma_j_y3[19]","ar_y3[20]", "sigma_j_y3[20]",
            "ar_y3[21]", "sigma_j_y3[21]","ar_y3[22]", "sigma_j_y3[22]","ar_y3[23]", "sigma_j_y3[23]","ar_y3[24]", "sigma_j_y3[24]",
            "ar_y3[25]", "sigma_j_y3[25]","ar_y3[26]", "sigma_j_y3[26]","ar_y3[27]", "sigma_j_y3[27]","ar_y3[28]", "sigma_j_y3[28]",
            "ar_y3[29]", "sigma_j_y3[29]","ar_y3[30]", "sigma_j_y3[30]","ar_y3[31]", "sigma_j_y3[31]","ar_y3[32]", "sigma_j_y3[32]",
            "ar_y3[33]", "sigma_j_y3[33]","ar_y3[34]", "sigma_j_y3[34]","ar_y3[35]", "sigma_j_y3[35]","ar_y3[36]", "sigma_j_y3[36]",
            "ar_y3[37]", "sigma_j_y3[37]","ar_y3[38]", "sigma_j_y3[38]","ar_y3[39]", "sigma_j_y3[39]","ar_y3[40]", "sigma_j_y3[40]",
            "alpha_y1_ref", "beta1_y1_ref", "beta2_y1_ref","sigma_y1_ref",
            "alpha_y2_ref", "beta1_y2_ref", "beta2_y2_ref","sigma_y2_ref",
            "alpha_y3_ref", "beta1_y3_ref", "beta2_y3_ref","sigma_y3_ref")
fit_lm1 <- jags(data = jagsdata_s1,  parameters.to.save = params,model.file = lm1_jags,
                n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC = F)



lm1_mcmc <- as.mcmc(fit_lm1)
lm1_mcmc_combi <- as.mcmc(rbind(lm1_mcmc[[1]], lm1_mcmc[[2]], lm1_mcmc[[3]],lm1_mcmc[[4]],lm1_mcmc[[5]],lm1_mcmc[[6]],lm1_mcmc[[7]]))
lm1_mcmc_combi_rep <- do.call(rbind, rep(list(lm1_mcmc_combi), 1)) # replication

lm1_mcmc_gdp_pre <- do.call(cbind, rep(list("alpha_gdp_pre[1]"=lm1_mcmc_combi[,"alpha_gdp_pre[1]"],"alpha_gdp_pre[2]"=lm1_mcmc_combi[,"alpha_gdp_pre[2]"],
                                            "alpha_gdp_pre[3]"=lm1_mcmc_combi[,"alpha_gdp_pre[3]"],"alpha_gdp_pre[4]"=lm1_mcmc_combi[,"alpha_gdp_pre[4]"],
                                            "alpha_gdp_pre[5]"=lm1_mcmc_combi[,"alpha_gdp_pre[5]"],"alpha_gdp_pre[6]"=lm1_mcmc_combi[,"alpha_gdp_pre[6]"],
                                            "alpha_gdp_pre[7]"=lm1_mcmc_combi[,"alpha_gdp_pre[7]"], "alpha_gdp_pre[8]"=lm1_mcmc_combi[,"alpha_gdp_pre[8]"],
                                            "alpha_gdp_pre[9]"=lm1_mcmc_combi[,"alpha_gdp_pre[9]"],"alpha_gdp_pre[10]"=lm1_mcmc_combi[,"alpha_gdp_pre[10]"],
                                            "alpha_gdp_pre[11]"=lm1_mcmc_combi[,"alpha_gdp_pre[11]"],"alpha_gdp_pre[12]"=lm1_mcmc_combi[,"alpha_gdp_pre[12]"],
                                            "alpha_gdp_pre[13]"=lm1_mcmc_combi[,"alpha_gdp_pre[13]"],"alpha_gdp_pre[14]"=lm1_mcmc_combi[,"alpha_gdp_pre[14]"],
                                            "alpha_gdp_pre[15]"=lm1_mcmc_combi[,"alpha_gdp_pre[15]"],"alpha_gdp_pre[16]"=lm1_mcmc_combi[,"alpha_gdp_pre[16]"],
                                            "alpha_gdp_pre[17]"=lm1_mcmc_combi[,"alpha_gdp_pre[17]"], "alpha_gdp_pre[18]"=lm1_mcmc_combi[,"alpha_gdp_pre[18]"],
                                            "alpha_gdp_pre[19]"=lm1_mcmc_combi[,"alpha_gdp_pre[19]"],"alpha_gdp_pre[20]"=lm1_mcmc_combi[,"alpha_gdp_pre[20]"],
                                            "alpha_gdp_pre[21]"=lm1_mcmc_combi[,"alpha_gdp_pre[21]"],"alpha_gdp_pre[22]"=lm1_mcmc_combi[,"alpha_gdp_pre[22]"],
                                            "alpha_gdp_pre[23]"=lm1_mcmc_combi[,"alpha_gdp_pre[23]"],"alpha_gdp_pre[24]"=lm1_mcmc_combi[,"alpha_gdp_pre[24]"],
                                            "alpha_gdp_pre[25]"=lm1_mcmc_combi[,"alpha_gdp_pre[25]"],"alpha_gdp_pre[26]"=lm1_mcmc_combi[,"alpha_gdp_pre[26]"],
                                            "alpha_gdp_pre[27]"=lm1_mcmc_combi[,"alpha_gdp_pre[27]"], "alpha_gdp_pre[28]"=lm1_mcmc_combi[,"alpha_gdp_pre[28]"],
                                            "alpha_gdp_pre[29]"=lm1_mcmc_combi[,"alpha_gdp_pre[29]"],"alpha_gdp_pre[30]"=lm1_mcmc_combi[,"alpha_gdp_pre[30]"],
                                            "alpha_gdp_pre[31]"=lm1_mcmc_combi[,"alpha_gdp_pre[31]"],"alpha_gdp_pre[32]"=lm1_mcmc_combi[,"alpha_gdp_pre[32]"],
                                            "alpha_gdp_pre[33]"=lm1_mcmc_combi[,"alpha_gdp_pre[33]"],"alpha_gdp_pre[34]"=lm1_mcmc_combi[,"alpha_gdp_pre[34]"],
                                            "alpha_gdp_pre[35]"=lm1_mcmc_combi[,"alpha_gdp_pre[35]"],"alpha_gdp_pre[36]"=lm1_mcmc_combi[,"alpha_gdp_pre[36]"],
                                            "alpha_gdp_pre[37]"=lm1_mcmc_combi[,"alpha_gdp_pre[37]"], "alpha_gdp_pre[38]"=lm1_mcmc_combi[,"alpha_gdp_pre[38]"],
                                            "alpha_gdp_pre[39]"=lm1_mcmc_combi[,"alpha_gdp_pre[39]"],"alpha_gdp_pre[40]"=lm1_mcmc_combi[,"alpha_gdp_pre[40]"],
                                            "error1[1]"=lm1_mcmc_combi[,"error1[1]"],"error1[2]"=lm1_mcmc_combi[,"error1[2]"],"error1[3]"=lm1_mcmc_combi[,"error1[3]"],
                                            "error1[4]"=lm1_mcmc_combi[,"error1[4]"],"error1[5]"=lm1_mcmc_combi[,"error1[5]"],"error1[6]"=lm1_mcmc_combi[,"error1[6]"],
                                            "error1[7]"=lm1_mcmc_combi[,"error1[7]"],"error1[8]"=lm1_mcmc_combi[,"error1[8]"],"error1[9]"=lm1_mcmc_combi[,"error1[9]"],
                                            "error1[10]"=lm1_mcmc_combi[,"error1[10]"],"error1[11]"=lm1_mcmc_combi[,"error1[11]"],"error1[12]"=lm1_mcmc_combi[,"error1[12]"],
                                            "error1[13]"=lm1_mcmc_combi[,"error1[13]"],"error1[14]"=lm1_mcmc_combi[,"error1[14]"],"error1[15]"=lm1_mcmc_combi[,"error1[15]"],
                                            "error1[16]"=lm1_mcmc_combi[,"error1[16]"],"error1[17]"=lm1_mcmc_combi[,"error1[17]"], "error1[18]"=lm1_mcmc_combi[,"error1[18]"],
                                            "error1[19]"=lm1_mcmc_combi[,"error1[19]"],"error1[20]"=lm1_mcmc_combi[,"error1[20]"],"error1[21]"=lm1_mcmc_combi[,"error1[21]"],
                                            "error1[22]"=lm1_mcmc_combi[,"error1[22]"],"error1[23]"=lm1_mcmc_combi[,"error1[23]"],"error1[24]"=lm1_mcmc_combi[,"error1[24]"],
                                            "error1[25]"=lm1_mcmc_combi[,"error1[25]"],"error1[26]"=lm1_mcmc_combi[,"error1[26]"],"error1[27]"=lm1_mcmc_combi[,"error1[27]"],
                                            "error1[28]"=lm1_mcmc_combi[,"error1[28]"],"error1[29]"=lm1_mcmc_combi[,"error1[29]"],"error1[30]"=lm1_mcmc_combi[,"error1[30]"],
                                            "error1[31]"=lm1_mcmc_combi[,"error1[31]"],"error1[32]"=lm1_mcmc_combi[,"error1[32]"],"error1[33]"=lm1_mcmc_combi[,"error1[33]"],
                                            "error1[34]"=lm1_mcmc_combi[,"error1[34]"],"error1[35]"=lm1_mcmc_combi[,"error1[35]"],"error1[36]"=lm1_mcmc_combi[,"error1[36]"],
                                            "error1[37]"=lm1_mcmc_combi[,"error1[37]"],"error1[38]"=lm1_mcmc_combi[,"error1[38]"],"error1[39]"=lm1_mcmc_combi[,"error1[39]"],
                                            "error1[40]"=lm1_mcmc_combi[,"error1[40]"],
                                            "error2[1]"=lm1_mcmc_combi[,"error2[1]"],"error2[2]"=lm1_mcmc_combi[,"error2[2]"],"error2[3]"=lm1_mcmc_combi[,"error2[3]"],
                                            "error2[4]"=lm1_mcmc_combi[,"error2[4]"],"error2[5]"=lm1_mcmc_combi[,"error2[5]"],"error2[6]"=lm1_mcmc_combi[,"error2[6]"],
                                            "error2[7]"=lm1_mcmc_combi[,"error2[7]"],"error2[8]"=lm1_mcmc_combi[,"error2[8]"],"error2[9]"=lm1_mcmc_combi[,"error2[9]"],
                                            "error2[10]"=lm1_mcmc_combi[,"error2[10]"],"error2[11]"=lm1_mcmc_combi[,"error2[11]"],"error2[12]"=lm1_mcmc_combi[,"error2[12]"],
                                            "error2[13]"=lm1_mcmc_combi[,"error2[13]"],"error2[14]"=lm1_mcmc_combi[,"error2[14]"],"error2[15]"=lm1_mcmc_combi[,"error2[15]"],
                                            "error2[16]"=lm1_mcmc_combi[,"error2[16]"],"error2[17]"=lm1_mcmc_combi[,"error2[17]"], "error2[18]"=lm1_mcmc_combi[,"error2[18]"],
                                            "error2[19]"=lm1_mcmc_combi[,"error2[19]"],"error2[20]"=lm1_mcmc_combi[,"error2[20]"],"error2[21]"=lm1_mcmc_combi[,"error2[21]"],
                                            "error2[22]"=lm1_mcmc_combi[,"error2[22]"],"error2[23]"=lm1_mcmc_combi[,"error2[23]"],"error2[24]"=lm1_mcmc_combi[,"error2[24]"],
                                            "error2[25]"=lm1_mcmc_combi[,"error2[25]"],"error2[26]"=lm1_mcmc_combi[,"error2[26]"],"error2[27]"=lm1_mcmc_combi[,"error2[27]"],
                                            "error2[28]"=lm1_mcmc_combi[,"error2[28]"],"error2[29]"=lm1_mcmc_combi[,"error2[29]"],"error2[30]"=lm1_mcmc_combi[,"error2[30]"],
                                            "error2[31]"=lm1_mcmc_combi[,"error2[31]"],"error2[32]"=lm1_mcmc_combi[,"error2[32]"],"error2[33]"=lm1_mcmc_combi[,"error2[33]"],
                                            "error2[34]"=lm1_mcmc_combi[,"error2[34]"],"error2[35]"=lm1_mcmc_combi[,"error2[35]"],"error2[36]"=lm1_mcmc_combi[,"error2[36]"],
                                            "error2[37]"=lm1_mcmc_combi[,"error2[37]"], "error2[38]"=lm1_mcmc_combi[,"error2[38]"],"error2[39]"=lm1_mcmc_combi[,"error2[39]"],
                                            "error2[40]"=lm1_mcmc_combi[,"error2[40]"],
                                            "error3[1]"=lm1_mcmc_combi[,"error3[1]"],"error3[2]"=lm1_mcmc_combi[,"error3[2]"],"error3[3]"=lm1_mcmc_combi[,"error3[3]"],
                                            "error3[4]"=lm1_mcmc_combi[,"error3[4]"],"error3[5]"=lm1_mcmc_combi[,"error3[5]"],"error3[6]"=lm1_mcmc_combi[,"error3[6]"],
                                            "error3[7]"=lm1_mcmc_combi[,"error3[7]"], "error3[8]"=lm1_mcmc_combi[,"error3[8]"],"error3[9]"=lm1_mcmc_combi[,"error3[9]"],
                                            "error3[10]"=lm1_mcmc_combi[,"error3[10]"],"error3[11]"=lm1_mcmc_combi[,"error3[11]"],"error3[12]"=lm1_mcmc_combi[,"error3[12]"],
                                            "error3[13]"=lm1_mcmc_combi[,"error3[13]"],"error3[14]"=lm1_mcmc_combi[,"error3[14]"],"error3[15]"=lm1_mcmc_combi[,"error3[15]"],
                                            "error3[16]"=lm1_mcmc_combi[,"error3[16]"],"error3[17]"=lm1_mcmc_combi[,"error3[17]"], "error3[18]"=lm1_mcmc_combi[,"error3[18]"],
                                            "error3[19]"=lm1_mcmc_combi[,"error3[19]"], "error3[20]"=lm1_mcmc_combi[,"error3[20]"],"error3[21]"=lm1_mcmc_combi[,"error3[21]"],
                                            "error3[22]"=lm1_mcmc_combi[,"error3[22]"],"error3[23]"=lm1_mcmc_combi[,"error3[23]"],"error3[24]"=lm1_mcmc_combi[,"error3[24]"],
                                            "error3[25]"=lm1_mcmc_combi[,"error3[25]"],"error3[26]"=lm1_mcmc_combi[,"error3[26]"],"error3[27]"=lm1_mcmc_combi[,"error3[27]"],
                                            "error3[28]"=lm1_mcmc_combi[,"error3[28]"],"error3[29]"=lm1_mcmc_combi[,"error3[29]"],"error3[30]"=lm1_mcmc_combi[,"error3[30]"],
                                            "error3[31]"=lm1_mcmc_combi[,"error3[31]"],"error3[32]"=lm1_mcmc_combi[,"error3[32]"],"error3[33]"=lm1_mcmc_combi[,"error3[33]"],
                                            "error3[34]"=lm1_mcmc_combi[,"error3[34]"],"error3[35]"=lm1_mcmc_combi[,"error3[35]"],"error3[36]"=lm1_mcmc_combi[,"error3[36]"],
                                            "error3[37]"=lm1_mcmc_combi[,"error3[37]"], "error3[38]"=lm1_mcmc_combi[,"error3[38]"],"error3[39]"=lm1_mcmc_combi[,"error3[39]"],
                                            "error3[40]"=lm1_mcmc_combi[,"error3[40]"],
                                            "error4[1]"=lm1_mcmc_combi[,"error4[1]"],"error4[2]"=lm1_mcmc_combi[,"error4[2]"],"error4[3]"=lm1_mcmc_combi[,"error4[3]"],
                                            "error4[4]"=lm1_mcmc_combi[,"error4[4]"],"error4[5]"=lm1_mcmc_combi[,"error4[5]"],"error4[6]"=lm1_mcmc_combi[,"error4[6]"],
                                            "error4[7]"=lm1_mcmc_combi[,"error4[7]"], "error4[8]"=lm1_mcmc_combi[,"error4[8]"],"error4[9]"=lm1_mcmc_combi[,"error4[9]"],
                                            "error4[10]"=lm1_mcmc_combi[,"error4[10]"],"error4[11]"=lm1_mcmc_combi[,"error4[11]"],"error4[12]"=lm1_mcmc_combi[,"error4[12]"],
                                            "error4[13]"=lm1_mcmc_combi[,"error4[13]"],"error4[14]"=lm1_mcmc_combi[,"error4[14]"],"error4[15]"=lm1_mcmc_combi[,"error4[15]"],
                                            "error4[16]"=lm1_mcmc_combi[,"error4[16]"],"error4[17]"=lm1_mcmc_combi[,"error4[17]"], "error4[18]"=lm1_mcmc_combi[,"error4[18]"],
                                            "error4[19]"=lm1_mcmc_combi[,"error4[19]"], "error4[20]"=lm1_mcmc_combi[,"error4[20]"],"error4[21]"=lm1_mcmc_combi[,"error4[21]"],
                                            "error4[22]"=lm1_mcmc_combi[,"error4[22]"],"error4[23]"=lm1_mcmc_combi[,"error4[23]"],"error4[24]"=lm1_mcmc_combi[,"error4[24]"],
                                            "error4[25]"=lm1_mcmc_combi[,"error4[25]"],"error4[26]"=lm1_mcmc_combi[,"error4[26]"],"error4[27]"=lm1_mcmc_combi[,"error4[27]"],
                                            "error4[28]"=lm1_mcmc_combi[,"error4[28]"],"error4[29]"=lm1_mcmc_combi[,"error4[29]"],"error4[30]"=lm1_mcmc_combi[,"error4[30]"],
                                            "error4[31]"=lm1_mcmc_combi[,"error4[31]"],"error4[32]"=lm1_mcmc_combi[,"error4[32]"],"error4[33]"=lm1_mcmc_combi[,"error4[33]"],
                                            "error4[34]"=lm1_mcmc_combi[,"error4[34]"],"error4[35]"=lm1_mcmc_combi[,"error4[35]"],"error4[36]"=lm1_mcmc_combi[,"error4[36]"],
                                            "error4[37]"=lm1_mcmc_combi[,"error4[37]"], "error4[38]"=lm1_mcmc_combi[,"error4[38]"],"error4[39]"=lm1_mcmc_combi[,"error4[39]"],
                                            "error4[40]"=lm1_mcmc_combi[,"error4[40]"],
                                            "error5[1]"=lm1_mcmc_combi[,"error5[1]"],"error5[2]"=lm1_mcmc_combi[,"error5[2]"],"error5[3]"=lm1_mcmc_combi[,"error5[3]"],
                                            "error5[4]"=lm1_mcmc_combi[,"error5[4]"],"error5[5]"=lm1_mcmc_combi[,"error5[5]"],"error5[6]"=lm1_mcmc_combi[,"error5[6]"],
                                            "error5[7]"=lm1_mcmc_combi[,"error5[7]"], "error5[8]"=lm1_mcmc_combi[,"error5[8]"],"error5[9]"=lm1_mcmc_combi[,"error5[9]"],
                                            "error5[10]"=lm1_mcmc_combi[,"error5[10]"],"error5[11]"=lm1_mcmc_combi[,"error5[11]"],"error5[12]"=lm1_mcmc_combi[,"error5[12]"],
                                            "error5[13]"=lm1_mcmc_combi[,"error5[13]"],"error5[14]"=lm1_mcmc_combi[,"error5[14]"],"error5[15]"=lm1_mcmc_combi[,"error5[15]"],
                                            "error5[16]"=lm1_mcmc_combi[,"error5[16]"],"error5[17]"=lm1_mcmc_combi[,"error5[17]"], "error5[18]"=lm1_mcmc_combi[,"error5[18]"],
                                            "error5[19]"=lm1_mcmc_combi[,"error5[19]"], "error5[20]"=lm1_mcmc_combi[,"error5[20]"],"error5[21]"=lm1_mcmc_combi[,"error5[21]"],
                                            "error5[22]"=lm1_mcmc_combi[,"error5[22]"],"error5[23]"=lm1_mcmc_combi[,"error5[23]"],"error5[24]"=lm1_mcmc_combi[,"error5[24]"],
                                            "error5[25]"=lm1_mcmc_combi[,"error5[25]"],"error5[26]"=lm1_mcmc_combi[,"error5[26]"],"error5[27]"=lm1_mcmc_combi[,"error5[27]"],
                                            "error5[28]"=lm1_mcmc_combi[,"error5[28]"],"error5[29]"=lm1_mcmc_combi[,"error5[29]"],"error5[30]"=lm1_mcmc_combi[,"error5[30]"],
                                            "error5[31]"=lm1_mcmc_combi[,"error5[31]"],"error5[32]"=lm1_mcmc_combi[,"error5[32]"],"error5[33]"=lm1_mcmc_combi[,"error5[33]"],
                                            "error5[34]"=lm1_mcmc_combi[,"error5[34]"],"error5[35]"=lm1_mcmc_combi[,"error5[35]"],"error5[36]"=lm1_mcmc_combi[,"error5[36]"],
                                            "error5[37]"=lm1_mcmc_combi[,"error5[37]"], "error5[38]"=lm1_mcmc_combi[,"error5[38]"],"error5[39]"=lm1_mcmc_combi[,"error5[39]"],
                                            "error5[40]"=lm1_mcmc_combi[,"error5[40]"],
                                            "error6[1]"=lm1_mcmc_combi[,"error6[1]"],"error6[2]"=lm1_mcmc_combi[,"error6[2]"],"error6[3]"=lm1_mcmc_combi[,"error6[3]"],
                                            "error6[4]"=lm1_mcmc_combi[,"error6[4]"],"error6[5]"=lm1_mcmc_combi[,"error6[5]"],"error6[6]"=lm1_mcmc_combi[,"error6[6]"],
                                            "error6[7]"=lm1_mcmc_combi[,"error6[7]"], "error6[8]"=lm1_mcmc_combi[,"error6[8]"],"error6[9]"=lm1_mcmc_combi[,"error6[9]"],
                                            "error6[10]"=lm1_mcmc_combi[,"error6[10]"],"error6[11]"=lm1_mcmc_combi[,"error6[11]"],"error6[12]"=lm1_mcmc_combi[,"error6[12]"],
                                            "error6[13]"=lm1_mcmc_combi[,"error6[13]"],"error6[14]"=lm1_mcmc_combi[,"error6[14]"],"error6[15]"=lm1_mcmc_combi[,"error6[15]"],
                                            "error6[16]"=lm1_mcmc_combi[,"error6[16]"],"error6[17]"=lm1_mcmc_combi[,"error6[17]"], "error6[18]"=lm1_mcmc_combi[,"error6[18]"],
                                            "error6[19]"=lm1_mcmc_combi[,"error6[19]"], "error6[20]"=lm1_mcmc_combi[,"error6[20]"],"error6[21]"=lm1_mcmc_combi[,"error6[21]"],
                                            "error6[22]"=lm1_mcmc_combi[,"error6[22]"],"error6[23]"=lm1_mcmc_combi[,"error6[23]"],"error6[24]"=lm1_mcmc_combi[,"error6[24]"],
                                            "error6[25]"=lm1_mcmc_combi[,"error6[25]"],"error6[26]"=lm1_mcmc_combi[,"error6[26]"],"error6[27]"=lm1_mcmc_combi[,"error6[27]"],
                                            "error6[28]"=lm1_mcmc_combi[,"error6[28]"],"error6[29]"=lm1_mcmc_combi[,"error6[29]"],"error6[30]"=lm1_mcmc_combi[,"error6[30]"],
                                            "error6[31]"=lm1_mcmc_combi[,"error6[31]"],"error6[32]"=lm1_mcmc_combi[,"error6[32]"],"error6[33]"=lm1_mcmc_combi[,"error6[33]"],
                                            "error6[34]"=lm1_mcmc_combi[,"error6[34]"],"error6[35]"=lm1_mcmc_combi[,"error6[35]"],"error6[36]"=lm1_mcmc_combi[,"error6[36]"],
                                            "error6[37]"=lm1_mcmc_combi[,"error6[37]"], "error6[38]"=lm1_mcmc_combi[,"error6[38]"],"error6[39]"=lm1_mcmc_combi[,"error6[39]"],
                                            "error6[40]"=lm1_mcmc_combi[,"error6[40]"],
                                            "error7[1]"=lm1_mcmc_combi[,"error7[1]"],"error7[2]"=lm1_mcmc_combi[,"error7[2]"],"error7[3]"=lm1_mcmc_combi[,"error7[3]"],
                                            "error7[4]"=lm1_mcmc_combi[,"error7[4]"],"error7[5]"=lm1_mcmc_combi[,"error7[5]"],"error7[6]"=lm1_mcmc_combi[,"error7[6]"],
                                            "error7[7]"=lm1_mcmc_combi[,"error7[7]"], "error7[8]"=lm1_mcmc_combi[,"error7[8]"],"error7[9]"=lm1_mcmc_combi[,"error7[9]"],
                                            "error7[10]"=lm1_mcmc_combi[,"error7[10]"],"error7[11]"=lm1_mcmc_combi[,"error7[11]"],"error7[12]"=lm1_mcmc_combi[,"error7[12]"],
                                            "error7[13]"=lm1_mcmc_combi[,"error7[13]"],"error7[14]"=lm1_mcmc_combi[,"error7[14]"],"error7[15]"=lm1_mcmc_combi[,"error7[15]"],
                                            "error7[16]"=lm1_mcmc_combi[,"error7[16]"],"error7[17]"=lm1_mcmc_combi[,"error7[17]"], "error7[18]"=lm1_mcmc_combi[,"error7[18]"],
                                            "error7[19]"=lm1_mcmc_combi[,"error7[19]"], "error7[20]"=lm1_mcmc_combi[,"error7[20]"],"error7[21]"=lm1_mcmc_combi[,"error7[21]"],
                                            "error7[22]"=lm1_mcmc_combi[,"error7[22]"],"error7[23]"=lm1_mcmc_combi[,"error7[23]"],"error7[24]"=lm1_mcmc_combi[,"error7[24]"],
                                            "error7[25]"=lm1_mcmc_combi[,"error7[25]"],"error7[26]"=lm1_mcmc_combi[,"error7[26]"],"error7[27]"=lm1_mcmc_combi[,"error7[27]"],
                                            "error7[28]"=lm1_mcmc_combi[,"error7[28]"],"error7[29]"=lm1_mcmc_combi[,"error7[29]"],"error7[30]"=lm1_mcmc_combi[,"error7[30]"],
                                            "error7[31]"=lm1_mcmc_combi[,"error7[31]"],"error7[32]"=lm1_mcmc_combi[,"error7[32]"],"error7[33]"=lm1_mcmc_combi[,"error7[33]"],
                                            "error7[34]"=lm1_mcmc_combi[,"error7[34]"],"error7[35]"=lm1_mcmc_combi[,"error7[35]"],"error7[36]"=lm1_mcmc_combi[,"error7[36]"],
                                            "error7[37]"=lm1_mcmc_combi[,"error7[37]"], "error7[38]"=lm1_mcmc_combi[,"error7[38]"],"error7[39]"=lm1_mcmc_combi[,"error7[39]"],
                                            "error7[40]"=lm1_mcmc_combi[,"error7[40]"],
                                            "error8[1]"=lm1_mcmc_combi[,"error8[1]"],"error8[2]"=lm1_mcmc_combi[,"error8[2]"],"error8[3]"=lm1_mcmc_combi[,"error8[3]"],
                                            "error8[4]"=lm1_mcmc_combi[,"error8[4]"],"error8[5]"=lm1_mcmc_combi[,"error8[5]"],"error8[6]"=lm1_mcmc_combi[,"error8[6]"],
                                            "error8[7]"=lm1_mcmc_combi[,"error8[7]"], "error8[8]"=lm1_mcmc_combi[,"error8[8]"],"error8[9]"=lm1_mcmc_combi[,"error8[9]"],
                                            "error8[10]"=lm1_mcmc_combi[,"error8[10]"],"error8[11]"=lm1_mcmc_combi[,"error8[11]"],"error8[12]"=lm1_mcmc_combi[,"error8[12]"],
                                            "error8[13]"=lm1_mcmc_combi[,"error8[13]"],"error8[14]"=lm1_mcmc_combi[,"error8[14]"],"error8[15]"=lm1_mcmc_combi[,"error8[15]"],
                                            "error8[16]"=lm1_mcmc_combi[,"error8[16]"],"error8[17]"=lm1_mcmc_combi[,"error8[17]"], "error8[18]"=lm1_mcmc_combi[,"error8[18]"],
                                            "error8[19]"=lm1_mcmc_combi[,"error8[19]"], "error8[20]"=lm1_mcmc_combi[,"error8[20]"],"error8[21]"=lm1_mcmc_combi[,"error8[21]"],
                                            "error8[22]"=lm1_mcmc_combi[,"error8[22]"],"error8[23]"=lm1_mcmc_combi[,"error8[23]"],"error8[24]"=lm1_mcmc_combi[,"error8[24]"],
                                            "error8[25]"=lm1_mcmc_combi[,"error8[25]"],"error8[26]"=lm1_mcmc_combi[,"error8[26]"],"error8[27]"=lm1_mcmc_combi[,"error8[27]"],
                                            "error8[28]"=lm1_mcmc_combi[,"error8[28]"],"error8[29]"=lm1_mcmc_combi[,"error8[29]"],"error8[30]"=lm1_mcmc_combi[,"error8[30]"],
                                            "error8[31]"=lm1_mcmc_combi[,"error8[31]"],"error8[32]"=lm1_mcmc_combi[,"error8[32]"],"error8[33]"=lm1_mcmc_combi[,"error8[33]"],
                                            "error8[34]"=lm1_mcmc_combi[,"error8[34]"],"error8[35]"=lm1_mcmc_combi[,"error8[35]"],"error8[36]"=lm1_mcmc_combi[,"error8[36]"],
                                            "error8[37]"=lm1_mcmc_combi[,"error8[37]"], "error8[38]"=lm1_mcmc_combi[,"error8[38]"],"error8[39]"=lm1_mcmc_combi[,"error8[39]"],
                                            "error8[40]"=lm1_mcmc_combi[,"error8[40]"],
                                            "error9[1]"=lm1_mcmc_combi[,"error9[1]"],"error9[2]"=lm1_mcmc_combi[,"error9[2]"],"error9[3]"=lm1_mcmc_combi[,"error9[3]"],
                                            "error9[4]"=lm1_mcmc_combi[,"error9[4]"],"error9[5]"=lm1_mcmc_combi[,"error9[5]"],"error9[6]"=lm1_mcmc_combi[,"error9[6]"],
                                            "error9[7]"=lm1_mcmc_combi[,"error9[7]"], "error9[8]"=lm1_mcmc_combi[,"error9[8]"],"error9[9]"=lm1_mcmc_combi[,"error9[9]"],
                                            "error9[10]"=lm1_mcmc_combi[,"error9[10]"],"error9[11]"=lm1_mcmc_combi[,"error9[11]"],"error9[12]"=lm1_mcmc_combi[,"error9[12]"],
                                            "error9[13]"=lm1_mcmc_combi[,"error9[13]"],"error9[14]"=lm1_mcmc_combi[,"error9[14]"],"error9[15]"=lm1_mcmc_combi[,"error9[15]"],
                                            "error9[16]"=lm1_mcmc_combi[,"error9[16]"],"error9[17]"=lm1_mcmc_combi[,"error9[17]"], "error9[18]"=lm1_mcmc_combi[,"error9[18]"],
                                            "error9[19]"=lm1_mcmc_combi[,"error9[19]"], "error9[20]"=lm1_mcmc_combi[,"error9[20]"],"error9[21]"=lm1_mcmc_combi[,"error9[21]"],
                                            "error9[22]"=lm1_mcmc_combi[,"error9[22]"],"error9[23]"=lm1_mcmc_combi[,"error9[23]"],"error9[24]"=lm1_mcmc_combi[,"error9[24]"],
                                            "error9[25]"=lm1_mcmc_combi[,"error9[25]"],"error9[26]"=lm1_mcmc_combi[,"error9[26]"],"error9[27]"=lm1_mcmc_combi[,"error9[27]"],
                                            "error9[28]"=lm1_mcmc_combi[,"error9[28]"],"error9[29]"=lm1_mcmc_combi[,"error9[29]"],"error9[30]"=lm1_mcmc_combi[,"error9[30]"],
                                            "error9[31]"=lm1_mcmc_combi[,"error9[31]"],"error9[32]"=lm1_mcmc_combi[,"error9[32]"],"error9[33]"=lm1_mcmc_combi[,"error9[33]"],
                                            "error9[34]"=lm1_mcmc_combi[,"error9[34]"],"error9[35]"=lm1_mcmc_combi[,"error9[35]"],"error9[36]"=lm1_mcmc_combi[,"error9[36]"],
                                            "error9[37]"=lm1_mcmc_combi[,"error9[37]"], "error9[38]"=lm1_mcmc_combi[,"error9[38]"],"error9[39]"=lm1_mcmc_combi[,"error9[39]"],
                                            "error9[40]"=lm1_mcmc_combi[,"error9[40]"],
                                            "error10[1]"=lm1_mcmc_combi[,"error10[1]"],"error10[2]"=lm1_mcmc_combi[,"error10[2]"],"error10[3]"=lm1_mcmc_combi[,"error10[3]"],
                                            "error10[4]"=lm1_mcmc_combi[,"error10[4]"],"error10[5]"=lm1_mcmc_combi[,"error10[5]"],"error10[6]"=lm1_mcmc_combi[,"error10[6]"],
                                            "error10[7]"=lm1_mcmc_combi[,"error10[7]"], "error10[8]"=lm1_mcmc_combi[,"error10[8]"],"error10[9]"=lm1_mcmc_combi[,"error10[9]"],
                                            "error10[10]"=lm1_mcmc_combi[,"error10[10]"],"error10[11]"=lm1_mcmc_combi[,"error10[11]"],"error10[12]"=lm1_mcmc_combi[,"error10[12]"],
                                            "error10[13]"=lm1_mcmc_combi[,"error10[13]"],"error10[14]"=lm1_mcmc_combi[,"error10[14]"],"error10[15]"=lm1_mcmc_combi[,"error10[15]"],
                                            "error10[16]"=lm1_mcmc_combi[,"error10[16]"],"error10[17]"=lm1_mcmc_combi[,"error10[17]"],"error10[18]"=lm1_mcmc_combi[,"error10[18]"],
                                            "error10[19]"=lm1_mcmc_combi[,"error10[19]"],"error10[20]"=lm1_mcmc_combi[,"error10[20]"],"error10[21]"=lm1_mcmc_combi[,"error10[21]"],
                                            "error10[22]"=lm1_mcmc_combi[,"error10[22]"],"error10[23]"=lm1_mcmc_combi[,"error10[23]"],"error10[24]"=lm1_mcmc_combi[,"error10[24]"],
                                            "error10[25]"=lm1_mcmc_combi[,"error10[25]"],"error10[26]"=lm1_mcmc_combi[,"error10[26]"],"error10[27]"=lm1_mcmc_combi[,"error10[27]"],
                                            "error10[28]"=lm1_mcmc_combi[,"error10[28]"],"error10[29]"=lm1_mcmc_combi[,"error10[29]"],"error10[30]"=lm1_mcmc_combi[,"error10[30]"],
                                            "error10[31]"=lm1_mcmc_combi[,"error10[31]"],"error10[32]"=lm1_mcmc_combi[,"error10[32]"],"error10[33]"=lm1_mcmc_combi[,"error10[33]"],
                                            "error10[34]"=lm1_mcmc_combi[,"error10[34]"],"error10[35]"=lm1_mcmc_combi[,"error10[35]"],"error10[36]"=lm1_mcmc_combi[,"error10[36]"],
                                            "error10[37]"=lm1_mcmc_combi[,"error10[37]"],"error10[38]"=lm1_mcmc_combi[,"error10[38]"],"error10[39]"=lm1_mcmc_combi[,"error10[39]"],
                                            "error10[40]"=lm1_mcmc_combi[,"error10[40]"],
                                            "error11[1]"=lm1_mcmc_combi[,"error11[1]"],"error11[2]"=lm1_mcmc_combi[,"error11[2]"],"error11[3]"=lm1_mcmc_combi[,"error11[3]"],
                                            "error11[4]"=lm1_mcmc_combi[,"error11[4]"],"error11[5]"=lm1_mcmc_combi[,"error11[5]"],"error11[6]"=lm1_mcmc_combi[,"error11[6]"],
                                            "error11[7]"=lm1_mcmc_combi[,"error11[7]"], "error11[8]"=lm1_mcmc_combi[,"error11[8]"],"error11[9]"=lm1_mcmc_combi[,"error11[9]"],
                                            "error11[10]"=lm1_mcmc_combi[,"error11[10]"],"error11[11]"=lm1_mcmc_combi[,"error11[11]"],"error11[12]"=lm1_mcmc_combi[,"error11[12]"],
                                            "error11[13]"=lm1_mcmc_combi[,"error11[13]"],"error11[14]"=lm1_mcmc_combi[,"error11[14]"],"error11[15]"=lm1_mcmc_combi[,"error11[15]"],
                                            "error11[16]"=lm1_mcmc_combi[,"error11[16]"],"error11[17]"=lm1_mcmc_combi[,"error11[17]"],"error11[18]"=lm1_mcmc_combi[,"error11[18]"],
                                            "error11[19]"=lm1_mcmc_combi[,"error11[19]"],"error11[20]"=lm1_mcmc_combi[,"error11[20]"],"error11[21]"=lm1_mcmc_combi[,"error11[21]"],
                                            "error11[22]"=lm1_mcmc_combi[,"error11[22]"],"error11[23]"=lm1_mcmc_combi[,"error11[23]"],"error11[24]"=lm1_mcmc_combi[,"error11[24]"],
                                            "error11[25]"=lm1_mcmc_combi[,"error11[25]"],"error11[26]"=lm1_mcmc_combi[,"error11[26]"],"error11[27]"=lm1_mcmc_combi[,"error11[27]"],
                                            "error11[28]"=lm1_mcmc_combi[,"error11[28]"],"error11[29]"=lm1_mcmc_combi[,"error11[29]"],"error11[30]"=lm1_mcmc_combi[,"error11[30]"],
                                            "error11[31]"=lm1_mcmc_combi[,"error11[31]"],"error11[32]"=lm1_mcmc_combi[,"error11[32]"],"error11[33]"=lm1_mcmc_combi[,"error11[33]"],
                                            "error11[34]"=lm1_mcmc_combi[,"error11[34]"],"error11[35]"=lm1_mcmc_combi[,"error11[35]"],"error11[36]"=lm1_mcmc_combi[,"error11[36]"],
                                            "error11[37]"=lm1_mcmc_combi[,"error11[37]"],"error11[38]"=lm1_mcmc_combi[,"error11[38]"],"error11[39]"=lm1_mcmc_combi[,"error11[39]"],
                                            "error11[40]"=lm1_mcmc_combi[,"error11[40]"],
                                            "error12[1]"=lm1_mcmc_combi[,"error12[1]"],"error12[2]"=lm1_mcmc_combi[,"error12[2]"],"error12[3]"=lm1_mcmc_combi[,"error12[3]"],
                                            "error12[4]"=lm1_mcmc_combi[,"error12[4]"],"error12[5]"=lm1_mcmc_combi[,"error12[5]"],"error12[6]"=lm1_mcmc_combi[,"error12[6]"],
                                            "error12[7]"=lm1_mcmc_combi[,"error12[7]"], "error12[8]"=lm1_mcmc_combi[,"error12[8]"],"error12[9]"=lm1_mcmc_combi[,"error12[9]"],
                                            "error12[10]"=lm1_mcmc_combi[,"error12[10]"],"error12[11]"=lm1_mcmc_combi[,"error12[11]"],"error12[12]"=lm1_mcmc_combi[,"error12[12]"],
                                            "error12[13]"=lm1_mcmc_combi[,"error12[13]"],"error12[14]"=lm1_mcmc_combi[,"error12[14]"],"error12[15]"=lm1_mcmc_combi[,"error12[15]"],                                 "error12[16]"=lm1_mcmc_combi[,"error12[16]"],"error12[17]"=lm1_mcmc_combi[,"error12[17]"],"error12[18]"=lm1_mcmc_combi[,"error12[18]"],
                                            "error12[19]"=lm1_mcmc_combi[,"error12[19]"],"error12[20]"=lm1_mcmc_combi[,"error12[20]"],"error12[21]"=lm1_mcmc_combi[,"error12[21]"],
                                            "error12[22]"=lm1_mcmc_combi[,"error12[22]"],"error12[23]"=lm1_mcmc_combi[,"error12[23]"],"error12[24]"=lm1_mcmc_combi[,"error12[24]"],
                                            "error12[25]"=lm1_mcmc_combi[,"error12[25]"],"error12[26]"=lm1_mcmc_combi[,"error12[26]"],"error12[27]"=lm1_mcmc_combi[,"error12[27]"],
                                            "error12[28]"=lm1_mcmc_combi[,"error12[28]"],"error12[29]"=lm1_mcmc_combi[,"error12[29]"],"error12[30]"=lm1_mcmc_combi[,"error12[30]"],
                                            "error12[31]"=lm1_mcmc_combi[,"error12[31]"],"error12[32]"=lm1_mcmc_combi[,"error12[32]"],"error12[33]"=lm1_mcmc_combi[,"error12[33]"],
                                            "error12[34]"=lm1_mcmc_combi[,"error12[34]"],"error12[35]"=lm1_mcmc_combi[,"error12[35]"],"error12[36]"=lm1_mcmc_combi[,"error12[36]"],
                                            "error12[37]"=lm1_mcmc_combi[,"error12[37]"],"error12[38]"=lm1_mcmc_combi[,"error12[38]"],"error12[39]"=lm1_mcmc_combi[,"error12[39]"],
                                            "error12[40]"=lm1_mcmc_combi[,"error12[40]"],
                                            "error13[1]"=lm1_mcmc_combi[,"error13[1]"],"error13[2]"=lm1_mcmc_combi[,"error13[2]"],"error13[3]"=lm1_mcmc_combi[,"error13[3]"],
                                            "error13[4]"=lm1_mcmc_combi[,"error13[4]"],"error13[5]"=lm1_mcmc_combi[,"error13[5]"],"error13[6]"=lm1_mcmc_combi[,"error13[6]"],
                                            "error13[7]"=lm1_mcmc_combi[,"error13[7]"], "error13[8]"=lm1_mcmc_combi[,"error13[8]"],"error13[9]"=lm1_mcmc_combi[,"error13[9]"],
                                            "error13[10]"=lm1_mcmc_combi[,"error13[10]"],"error13[11]"=lm1_mcmc_combi[,"error13[11]"],"error13[12]"=lm1_mcmc_combi[,"error13[12]"],
                                            "error13[13]"=lm1_mcmc_combi[,"error13[13]"],"error13[14]"=lm1_mcmc_combi[,"error13[14]"],"error13[15]"=lm1_mcmc_combi[,"error13[15]"],
                                            "error13[16]"=lm1_mcmc_combi[,"error13[16]"],"error13[17]"=lm1_mcmc_combi[,"error13[17]"],"error13[18]"=lm1_mcmc_combi[,"error13[18]"],
                                            "error13[19]"=lm1_mcmc_combi[,"error13[19]"],"error13[20]"=lm1_mcmc_combi[,"error13[20]"],"error13[21]"=lm1_mcmc_combi[,"error13[21]"],
                                            "error13[22]"=lm1_mcmc_combi[,"error13[22]"],"error13[23]"=lm1_mcmc_combi[,"error13[23]"],"error13[24]"=lm1_mcmc_combi[,"error13[24]"],
                                            "error13[25]"=lm1_mcmc_combi[,"error13[25]"],"error13[26]"=lm1_mcmc_combi[,"error13[26]"],"error13[27]"=lm1_mcmc_combi[,"error13[27]"],
                                            "error13[28]"=lm1_mcmc_combi[,"error13[28]"],"error13[29]"=lm1_mcmc_combi[,"error13[29]"],"error13[30]"=lm1_mcmc_combi[,"error13[30]"],
                                            "error13[31]"=lm1_mcmc_combi[,"error13[31]"],"error13[32]"=lm1_mcmc_combi[,"error13[32]"],"error13[33]"=lm1_mcmc_combi[,"error13[33]"],
                                            "error13[34]"=lm1_mcmc_combi[,"error13[34]"],"error13[35]"=lm1_mcmc_combi[,"error13[35]"],"error13[36]"=lm1_mcmc_combi[,"error13[36]"],
                                            "error13[37]"=lm1_mcmc_combi[,"error13[37]"],"error13[38]"=lm1_mcmc_combi[,"error13[38]"],"error13[39]"=lm1_mcmc_combi[,"error13[39]"],
                                            "error13[40]"=lm1_mcmc_combi[,"error13[40]"],
                                            "error14[1]"=lm1_mcmc_combi[,"error14[1]"],"error14[2]"=lm1_mcmc_combi[,"error14[2]"],"error14[3]"=lm1_mcmc_combi[,"error14[3]"],
                                            "error14[4]"=lm1_mcmc_combi[,"error14[4]"],"error14[5]"=lm1_mcmc_combi[,"error14[5]"],"error14[6]"=lm1_mcmc_combi[,"error14[6]"],
                                            "error14[7]"=lm1_mcmc_combi[,"error14[7]"], "error14[8]"=lm1_mcmc_combi[,"error14[8]"],"error14[9]"=lm1_mcmc_combi[,"error14[9]"],
                                            "error14[10]"=lm1_mcmc_combi[,"error14[10]"],"error14[11]"=lm1_mcmc_combi[,"error14[11]"],"error14[12]"=lm1_mcmc_combi[,"error14[12]"],
                                            "error14[13]"=lm1_mcmc_combi[,"error14[13]"],"error14[14]"=lm1_mcmc_combi[,"error14[14]"],"error14[15]"=lm1_mcmc_combi[,"error14[15]"],
                                            "error14[16]"=lm1_mcmc_combi[,"error14[16]"],"error14[17]"=lm1_mcmc_combi[,"error14[17]"],"error14[18]"=lm1_mcmc_combi[,"error14[18]"],
                                            "error14[19]"=lm1_mcmc_combi[,"error14[19]"],"error14[20]"=lm1_mcmc_combi[,"error14[20]"],"error14[21]"=lm1_mcmc_combi[,"error14[21]"],
                                            "error14[22]"=lm1_mcmc_combi[,"error14[22]"],"error14[23]"=lm1_mcmc_combi[,"error14[23]"],"error14[24]"=lm1_mcmc_combi[,"error14[24]"],
                                            "error14[25]"=lm1_mcmc_combi[,"error14[25]"],"error14[26]"=lm1_mcmc_combi[,"error14[26]"],"error14[27]"=lm1_mcmc_combi[,"error14[27]"],
                                            "error14[28]"=lm1_mcmc_combi[,"error14[28]"],"error14[29]"=lm1_mcmc_combi[,"error14[29]"],"error14[30]"=lm1_mcmc_combi[,"error14[30]"],
                                            "error14[31]"=lm1_mcmc_combi[,"error14[31]"],"error14[32]"=lm1_mcmc_combi[,"error14[32]"],"error14[33]"=lm1_mcmc_combi[,"error14[33]"],
                                            "error14[34]"=lm1_mcmc_combi[,"error14[34]"],"error14[35]"=lm1_mcmc_combi[,"error14[35]"],"error14[36]"=lm1_mcmc_combi[,"error14[36]"],
                                            "error14[37]"=lm1_mcmc_combi[,"error14[37]"],"error14[38]"=lm1_mcmc_combi[,"error14[38]"],"error14[39]"=lm1_mcmc_combi[,"error14[39]"],
                                            "error14[40]"=lm1_mcmc_combi[,"error14[40]"],
                                            "error15[1]"=lm1_mcmc_combi[,"error15[1]"],"error15[2]"=lm1_mcmc_combi[,"error15[2]"],"error15[3]"=lm1_mcmc_combi[,"error15[3]"],
                                            "error15[4]"=lm1_mcmc_combi[,"error15[4]"],"error15[5]"=lm1_mcmc_combi[,"error15[5]"],"error15[6]"=lm1_mcmc_combi[,"error15[6]"],
                                            "error15[7]"=lm1_mcmc_combi[,"error15[7]"], "error15[8]"=lm1_mcmc_combi[,"error15[8]"],"error15[9]"=lm1_mcmc_combi[,"error15[9]"],
                                            "error15[10]"=lm1_mcmc_combi[,"error15[10]"],"error15[11]"=lm1_mcmc_combi[,"error15[11]"],"error15[12]"=lm1_mcmc_combi[,"error15[12]"],
                                            "error15[13]"=lm1_mcmc_combi[,"error15[13]"],"error15[14]"=lm1_mcmc_combi[,"error15[14]"],"error15[15]"=lm1_mcmc_combi[,"error15[15]"],
                                            "error15[16]"=lm1_mcmc_combi[,"error15[16]"],"error15[17]"=lm1_mcmc_combi[,"error15[17]"],"error15[18]"=lm1_mcmc_combi[,"error15[18]"],
                                            "error15[19]"=lm1_mcmc_combi[,"error15[19]"],"error15[20]"=lm1_mcmc_combi[,"error15[20]"],"error15[21]"=lm1_mcmc_combi[,"error15[21]"],
                                            "error15[22]"=lm1_mcmc_combi[,"error15[22]"],"error15[23]"=lm1_mcmc_combi[,"error15[23]"],"error15[24]"=lm1_mcmc_combi[,"error15[24]"],
                                            "error15[25]"=lm1_mcmc_combi[,"error15[25]"],"error15[26]"=lm1_mcmc_combi[,"error15[26]"],"error15[27]"=lm1_mcmc_combi[,"error15[27]"],
                                            "error15[28]"=lm1_mcmc_combi[,"error15[28]"],"error15[29]"=lm1_mcmc_combi[,"error15[29]"],"error15[30]"=lm1_mcmc_combi[,"error15[30]"],
                                            "error15[31]"=lm1_mcmc_combi[,"error15[31]"],"error15[32]"=lm1_mcmc_combi[,"error15[32]"],"error15[33]"=lm1_mcmc_combi[,"error15[33]"],
                                            "error15[34]"=lm1_mcmc_combi[,"error15[34]"],"error15[35]"=lm1_mcmc_combi[,"error15[35]"],"error15[36]"=lm1_mcmc_combi[,"error15[36]"],
                                            "error15[37]"=lm1_mcmc_combi[,"error15[37]"],"error15[38]"=lm1_mcmc_combi[,"error15[38]"],"error15[39]"=lm1_mcmc_combi[,"error15[39]"],
                                            "error15[40]"=lm1_mcmc_combi[,"error15[40]"],
                                            "error16[1]"=lm1_mcmc_combi[,"error16[1]"],"error16[2]"=lm1_mcmc_combi[,"error16[2]"],"error16[3]"=lm1_mcmc_combi[,"error16[3]"],
                                            "error16[4]"=lm1_mcmc_combi[,"error16[4]"],"error16[5]"=lm1_mcmc_combi[,"error16[5]"],"error16[6]"=lm1_mcmc_combi[,"error16[6]"],
                                            "error16[7]"=lm1_mcmc_combi[,"error16[7]"], "error16[8]"=lm1_mcmc_combi[,"error16[8]"],"error16[9]"=lm1_mcmc_combi[,"error16[9]"],
                                            "error16[10]"=lm1_mcmc_combi[,"error16[10]"],"error16[11]"=lm1_mcmc_combi[,"error16[11]"],"error16[12]"=lm1_mcmc_combi[,"error16[12]"],
                                            "error16[13]"=lm1_mcmc_combi[,"error16[13]"],"error16[14]"=lm1_mcmc_combi[,"error16[14]"],"error16[15]"=lm1_mcmc_combi[,"error16[15]"],
                                            "error16[16]"=lm1_mcmc_combi[,"error16[16]"],"error16[17]"=lm1_mcmc_combi[,"error16[17]"],"error16[18]"=lm1_mcmc_combi[,"error16[18]"],
                                            "error16[19]"=lm1_mcmc_combi[,"error16[19]"],"error16[20]"=lm1_mcmc_combi[,"error16[20]"],"error16[21]"=lm1_mcmc_combi[,"error16[21]"],
                                            "error16[22]"=lm1_mcmc_combi[,"error16[22]"],"error16[23]"=lm1_mcmc_combi[,"error16[23]"],"error16[24]"=lm1_mcmc_combi[,"error16[24]"],
                                            "error16[25]"=lm1_mcmc_combi[,"error16[25]"],"error16[26]"=lm1_mcmc_combi[,"error16[26]"],"error16[27]"=lm1_mcmc_combi[,"error16[27]"],
                                            "error16[28]"=lm1_mcmc_combi[,"error16[28]"],"error16[29]"=lm1_mcmc_combi[,"error16[29]"],"error16[30]"=lm1_mcmc_combi[,"error16[30]"],
                                            "error16[31]"=lm1_mcmc_combi[,"error16[31]"],"error16[32]"=lm1_mcmc_combi[,"error16[32]"],"error16[33]"=lm1_mcmc_combi[,"error16[33]"],
                                            "error16[34]"=lm1_mcmc_combi[,"error16[34]"],"error16[35]"=lm1_mcmc_combi[,"error16[35]"],"error16[36]"=lm1_mcmc_combi[,"error16[36]"],
                                            "error16[37]"=lm1_mcmc_combi[,"error16[37]"],"error16[38]"=lm1_mcmc_combi[,"error16[38]"],"error16[39]"=lm1_mcmc_combi[,"error16[39]"],
                                            "error16[40]"=lm1_mcmc_combi[,"error16[40]"],
                                            "error17[1]"=lm1_mcmc_combi[,"error17[1]"],"error17[2]"=lm1_mcmc_combi[,"error17[2]"],"error17[3]"=lm1_mcmc_combi[,"error17[3]"],
                                            "error17[4]"=lm1_mcmc_combi[,"error17[4]"],"error17[5]"=lm1_mcmc_combi[,"error17[5]"],"error17[6]"=lm1_mcmc_combi[,"error17[6]"],
                                            "error17[7]"=lm1_mcmc_combi[,"error17[7]"], "error17[8]"=lm1_mcmc_combi[,"error17[8]"],"error17[9]"=lm1_mcmc_combi[,"error17[9]"],
                                            "error17[10]"=lm1_mcmc_combi[,"error17[10]"],"error17[11]"=lm1_mcmc_combi[,"error17[11]"],"error17[12]"=lm1_mcmc_combi[,"error17[12]"],
                                            "error17[13]"=lm1_mcmc_combi[,"error17[13]"],"error17[14]"=lm1_mcmc_combi[,"error17[14]"],"error17[15]"=lm1_mcmc_combi[,"error17[15]"],
                                            "error17[16]"=lm1_mcmc_combi[,"error17[16]"],"error17[17]"=lm1_mcmc_combi[,"error17[17]"],"error17[18]"=lm1_mcmc_combi[,"error17[18]"],
                                            "error17[19]"=lm1_mcmc_combi[,"error17[19]"],"error17[20]"=lm1_mcmc_combi[,"error17[20]"],"error17[21]"=lm1_mcmc_combi[,"error17[21]"],
                                            "error17[22]"=lm1_mcmc_combi[,"error17[22]"],"error17[23]"=lm1_mcmc_combi[,"error17[23]"],"error17[24]"=lm1_mcmc_combi[,"error17[24]"],
                                            "error17[25]"=lm1_mcmc_combi[,"error17[25]"],"error17[26]"=lm1_mcmc_combi[,"error17[26]"],"error17[27]"=lm1_mcmc_combi[,"error17[27]"],
                                            "error17[28]"=lm1_mcmc_combi[,"error17[28]"],"error17[29]"=lm1_mcmc_combi[,"error17[29]"],"error17[30]"=lm1_mcmc_combi[,"error17[30]"],
                                            "error17[31]"=lm1_mcmc_combi[,"error17[31]"],"error17[32]"=lm1_mcmc_combi[,"error17[32]"],"error17[33]"=lm1_mcmc_combi[,"error17[33]"],
                                            "error17[34]"=lm1_mcmc_combi[,"error17[34]"],"error17[35]"=lm1_mcmc_combi[,"error17[35]"],"error17[36]"=lm1_mcmc_combi[,"error17[36]"],
                                            "error17[37]"=lm1_mcmc_combi[,"error17[37]"],"error17[38]"=lm1_mcmc_combi[,"error17[38]"],"error17[39]"=lm1_mcmc_combi[,"error17[39]"],
                                            "error17[40]"=lm1_mcmc_combi[,"error17[40]"],
                                            "error18[1]"=lm1_mcmc_combi[,"error18[1]"],"error18[2]"=lm1_mcmc_combi[,"error18[2]"],"error18[3]"=lm1_mcmc_combi[,"error18[3]"],
                                            "error18[4]"=lm1_mcmc_combi[,"error18[4]"],"error18[5]"=lm1_mcmc_combi[,"error18[5]"],"error18[6]"=lm1_mcmc_combi[,"error18[6]"],
                                            "error18[7]"=lm1_mcmc_combi[,"error18[7]"], "error18[8]"=lm1_mcmc_combi[,"error18[8]"],"error18[9]"=lm1_mcmc_combi[,"error18[9]"],
                                            "error18[10]"=lm1_mcmc_combi[,"error18[10]"],"error18[11]"=lm1_mcmc_combi[,"error18[11]"],"error18[12]"=lm1_mcmc_combi[,"error18[12]"],
                                            "error18[13]"=lm1_mcmc_combi[,"error18[13]"],"error18[14]"=lm1_mcmc_combi[,"error18[14]"],"error18[15]"=lm1_mcmc_combi[,"error18[15]"],
                                            "error18[16]"=lm1_mcmc_combi[,"error18[16]"],"error18[17]"=lm1_mcmc_combi[,"error18[17]"],"error18[18]"=lm1_mcmc_combi[,"error18[18]"],
                                            "error18[19]"=lm1_mcmc_combi[,"error18[19]"],"error18[20]"=lm1_mcmc_combi[,"error18[20]"],"error18[21]"=lm1_mcmc_combi[,"error18[21]"],
                                            "error18[22]"=lm1_mcmc_combi[,"error18[22]"],"error18[23]"=lm1_mcmc_combi[,"error18[23]"],"error18[24]"=lm1_mcmc_combi[,"error18[24]"],
                                            "error18[25]"=lm1_mcmc_combi[,"error18[25]"],"error18[26]"=lm1_mcmc_combi[,"error18[26]"],"error18[27]"=lm1_mcmc_combi[,"error18[27]"],
                                            "error18[28]"=lm1_mcmc_combi[,"error18[28]"],"error18[29]"=lm1_mcmc_combi[,"error18[29]"],"error18[30]"=lm1_mcmc_combi[,"error18[30]"],
                                            "error18[31]"=lm1_mcmc_combi[,"error18[31]"],"error18[32]"=lm1_mcmc_combi[,"error18[32]"],"error18[33]"=lm1_mcmc_combi[,"error18[33]"],
                                            "error18[34]"=lm1_mcmc_combi[,"error18[34]"],"error18[35]"=lm1_mcmc_combi[,"error18[35]"],"error18[36]"=lm1_mcmc_combi[,"error18[36]"],
                                            "error18[37]"=lm1_mcmc_combi[,"error18[37]"],"error18[38]"=lm1_mcmc_combi[,"error18[38]"],"error18[39]"=lm1_mcmc_combi[,"error18[39]"],
                                            "error18[40]"=lm1_mcmc_combi[,"error18[40]"],
                                            "error19[1]"=lm1_mcmc_combi[,"error19[1]"],"error19[2]"=lm1_mcmc_combi[,"error19[2]"],"error19[3]"=lm1_mcmc_combi[,"error19[3]"],
                                            "error19[4]"=lm1_mcmc_combi[,"error19[4]"],"error19[5]"=lm1_mcmc_combi[,"error19[5]"],"error19[6]"=lm1_mcmc_combi[,"error19[6]"],
                                            "error19[7]"=lm1_mcmc_combi[,"error19[7]"], "error19[8]"=lm1_mcmc_combi[,"error19[8]"],"error19[9]"=lm1_mcmc_combi[,"error19[9]"],
                                            "error19[10]"=lm1_mcmc_combi[,"error19[10]"],"error19[11]"=lm1_mcmc_combi[,"error19[11]"],"error19[12]"=lm1_mcmc_combi[,"error19[12]"],
                                            "error19[13]"=lm1_mcmc_combi[,"error19[13]"],"error19[14]"=lm1_mcmc_combi[,"error19[14]"],"error19[15]"=lm1_mcmc_combi[,"error19[15]"],
                                            "error19[16]"=lm1_mcmc_combi[,"error19[16]"],"error19[17]"=lm1_mcmc_combi[,"error19[17]"],"error19[18]"=lm1_mcmc_combi[,"error19[18]"],
                                            "error19[19]"=lm1_mcmc_combi[,"error19[19]"],"error19[20]"=lm1_mcmc_combi[,"error19[20]"],"error19[21]"=lm1_mcmc_combi[,"error19[21]"],
                                            "error19[22]"=lm1_mcmc_combi[,"error19[22]"],"error19[23]"=lm1_mcmc_combi[,"error19[23]"],"error19[24]"=lm1_mcmc_combi[,"error19[24]"],
                                            "error19[25]"=lm1_mcmc_combi[,"error19[25]"],"error19[26]"=lm1_mcmc_combi[,"error19[26]"],"error19[27]"=lm1_mcmc_combi[,"error19[27]"],
                                            "error19[28]"=lm1_mcmc_combi[,"error19[28]"],"error19[29]"=lm1_mcmc_combi[,"error19[29]"],"error19[30]"=lm1_mcmc_combi[,"error19[30]"],
                                            "error19[31]"=lm1_mcmc_combi[,"error19[31]"],"error19[32]"=lm1_mcmc_combi[,"error19[32]"],"error19[33]"=lm1_mcmc_combi[,"error19[33]"],
                                            "error19[34]"=lm1_mcmc_combi[,"error19[34]"],"error19[35]"=lm1_mcmc_combi[,"error19[35]"],"error19[36]"=lm1_mcmc_combi[,"error19[36]"],
                                            "error19[37]"=lm1_mcmc_combi[,"error19[37]"],"error19[38]"=lm1_mcmc_combi[,"error19[38]"],"error19[39]"=lm1_mcmc_combi[,"error19[39]"],
                                            "error19[40]"=lm1_mcmc_combi[,"error19[40]"],
                                            "error20[1]"=lm1_mcmc_combi[,"error20[1]"],"error20[2]"=lm1_mcmc_combi[,"error20[2]"],"error20[3]"=lm1_mcmc_combi[,"error20[3]"],
                                            "error20[4]"=lm1_mcmc_combi[,"error20[4]"],"error20[5]"=lm1_mcmc_combi[,"error20[5]"],"error20[6]"=lm1_mcmc_combi[,"error20[6]"],
                                            "error20[7]"=lm1_mcmc_combi[,"error20[7]"], "error20[8]"=lm1_mcmc_combi[,"error20[8]"],"error20[9]"=lm1_mcmc_combi[,"error20[9]"],
                                            "error20[10]"=lm1_mcmc_combi[,"error20[10]"],"error20[11]"=lm1_mcmc_combi[,"error20[11]"],"error20[12]"=lm1_mcmc_combi[,"error20[12]"],
                                            "error20[13]"=lm1_mcmc_combi[,"error20[13]"],"error20[14]"=lm1_mcmc_combi[,"error20[14]"],"error20[15]"=lm1_mcmc_combi[,"error20[15]"],  
                                            "error20[16]"=lm1_mcmc_combi[,"error20[16]"],"error20[17]"=lm1_mcmc_combi[,"error20[17]"],"error20[18]"=lm1_mcmc_combi[,"error20[18]"],
                                            "error20[19]"=lm1_mcmc_combi[,"error20[19]"],"error20[20]"=lm1_mcmc_combi[,"error20[20]"],"error20[21]"=lm1_mcmc_combi[,"error20[21]"],
                                            "error20[22]"=lm1_mcmc_combi[,"error20[22]"],"error20[23]"=lm1_mcmc_combi[,"error20[23]"],"error20[24]"=lm1_mcmc_combi[,"error20[24]"],
                                            "error20[25]"=lm1_mcmc_combi[,"error20[25]"],"error20[26]"=lm1_mcmc_combi[,"error20[26]"],"error20[27]"=lm1_mcmc_combi[,"error20[27]"],
                                            "error20[28]"=lm1_mcmc_combi[,"error20[28]"],"error20[29]"=lm1_mcmc_combi[,"error20[29]"],"error20[30]"=lm1_mcmc_combi[,"error20[30]"],
                                            "error20[31]"=lm1_mcmc_combi[,"error20[31]"],"error20[32]"=lm1_mcmc_combi[,"error20[32]"],"error20[33]"=lm1_mcmc_combi[,"error20[33]"],
                                            "error20[34]"=lm1_mcmc_combi[,"error20[34]"],"error20[35]"=lm1_mcmc_combi[,"error20[35]"],"error20[36]"=lm1_mcmc_combi[,"error20[36]"],
                                            "error20[37]"=lm1_mcmc_combi[,"error20[37]"],"error20[38]"=lm1_mcmc_combi[,"error20[38]"],"error20[39]"=lm1_mcmc_combi[,"error20[39]"],
                                            "error20[40]"=lm1_mcmc_combi[,"error20[40]"],
                                            "error21[1]"=lm1_mcmc_combi[,"error21[1]"],"error21[2]"=lm1_mcmc_combi[,"error21[2]"],"error21[3]"=lm1_mcmc_combi[,"error21[3]"],
                                            "error21[4]"=lm1_mcmc_combi[,"error21[4]"],"error21[5]"=lm1_mcmc_combi[,"error21[5]"],"error21[6]"=lm1_mcmc_combi[,"error21[6]"],
                                            "error21[7]"=lm1_mcmc_combi[,"error21[7]"], "error21[8]"=lm1_mcmc_combi[,"error21[8]"],"error21[9]"=lm1_mcmc_combi[,"error21[9]"],
                                            "error21[10]"=lm1_mcmc_combi[,"error21[10]"],"error21[11]"=lm1_mcmc_combi[,"error21[11]"],"error21[12]"=lm1_mcmc_combi[,"error21[12]"],
                                            "error21[13]"=lm1_mcmc_combi[,"error21[13]"],"error21[14]"=lm1_mcmc_combi[,"error21[14]"],"error21[15]"=lm1_mcmc_combi[,"error21[15]"],
                                            "error21[16]"=lm1_mcmc_combi[,"error21[16]"],"error21[17]"=lm1_mcmc_combi[,"error21[17]"],"error21[18]"=lm1_mcmc_combi[,"error21[18]"],
                                            "error21[19]"=lm1_mcmc_combi[,"error21[19]"],"error21[20]"=lm1_mcmc_combi[,"error21[20]"],"error21[21]"=lm1_mcmc_combi[,"error21[21]"],
                                            "error21[22]"=lm1_mcmc_combi[,"error21[22]"],"error21[23]"=lm1_mcmc_combi[,"error21[23]"],"error21[24]"=lm1_mcmc_combi[,"error21[24]"],
                                            "error21[25]"=lm1_mcmc_combi[,"error21[25]"],"error21[26]"=lm1_mcmc_combi[,"error21[26]"],"error21[27]"=lm1_mcmc_combi[,"error21[27]"],
                                            "error21[28]"=lm1_mcmc_combi[,"error21[28]"],"error21[29]"=lm1_mcmc_combi[,"error21[29]"],"error21[30]"=lm1_mcmc_combi[,"error21[30]"],
                                            "error21[31]"=lm1_mcmc_combi[,"error21[31]"],"error21[32]"=lm1_mcmc_combi[,"error21[32]"],"error21[33]"=lm1_mcmc_combi[,"error21[33]"],
                                            "error21[34]"=lm1_mcmc_combi[,"error21[34]"],"error21[35]"=lm1_mcmc_combi[,"error21[35]"],"error21[36]"=lm1_mcmc_combi[,"error21[36]"],
                                            "error21[37]"=lm1_mcmc_combi[,"error21[37]"],"error21[38]"=lm1_mcmc_combi[,"error21[38]"],"error21[39]"=lm1_mcmc_combi[,"error21[39]"],
                                            "error21[40]"=lm1_mcmc_combi[,"error21[40]"],
                                            "error22[1]"=lm1_mcmc_combi[,"error22[1]"],"error22[2]"=lm1_mcmc_combi[,"error22[2]"],"error22[3]"=lm1_mcmc_combi[,"error22[3]"],
                                            "error22[4]"=lm1_mcmc_combi[,"error22[4]"],"error22[5]"=lm1_mcmc_combi[,"error22[5]"],"error22[6]"=lm1_mcmc_combi[,"error22[6]"],
                                            "error22[7]"=lm1_mcmc_combi[,"error22[7]"], "error22[8]"=lm1_mcmc_combi[,"error22[8]"],"error22[9]"=lm1_mcmc_combi[,"error22[9]"],
                                            "error22[10]"=lm1_mcmc_combi[,"error22[10]"],"error22[11]"=lm1_mcmc_combi[,"error22[11]"],"error22[12]"=lm1_mcmc_combi[,"error22[12]"],
                                            "error22[13]"=lm1_mcmc_combi[,"error22[13]"],"error22[14]"=lm1_mcmc_combi[,"error22[14]"],"error22[15]"=lm1_mcmc_combi[,"error22[15]"],
                                            "error22[16]"=lm1_mcmc_combi[,"error22[16]"],"error22[17]"=lm1_mcmc_combi[,"error22[17]"],"error22[18]"=lm1_mcmc_combi[,"error22[18]"],
                                            "error22[19]"=lm1_mcmc_combi[,"error22[19]"],"error22[20]"=lm1_mcmc_combi[,"error22[20]"],"error22[21]"=lm1_mcmc_combi[,"error22[21]"],
                                            "error22[22]"=lm1_mcmc_combi[,"error22[22]"],"error22[23]"=lm1_mcmc_combi[,"error22[23]"],"error22[24]"=lm1_mcmc_combi[,"error22[24]"],
                                            "error22[25]"=lm1_mcmc_combi[,"error22[25]"],"error22[26]"=lm1_mcmc_combi[,"error22[26]"],"error22[27]"=lm1_mcmc_combi[,"error22[27]"],
                                            "error22[28]"=lm1_mcmc_combi[,"error22[28]"],"error22[29]"=lm1_mcmc_combi[,"error22[29]"],"error22[30]"=lm1_mcmc_combi[,"error22[30]"],
                                            "error22[31]"=lm1_mcmc_combi[,"error22[31]"],"error22[32]"=lm1_mcmc_combi[,"error22[32]"],"error22[33]"=lm1_mcmc_combi[,"error22[33]"],
                                            "error22[34]"=lm1_mcmc_combi[,"error22[34]"],"error22[35]"=lm1_mcmc_combi[,"error22[35]"],"error22[36]"=lm1_mcmc_combi[,"error22[36]"],
                                            "error22[37]"=lm1_mcmc_combi[,"error22[37]"],"error22[38]"=lm1_mcmc_combi[,"error22[38]"],"error22[39]"=lm1_mcmc_combi[,"error22[39]"],
                                            "error22[40]"=lm1_mcmc_combi[,"error22[40]"],
                                            "error23[1]"=lm1_mcmc_combi[,"error23[1]"],"error23[2]"=lm1_mcmc_combi[,"error23[2]"],"error23[3]"=lm1_mcmc_combi[,"error23[3]"],
                                            "error23[4]"=lm1_mcmc_combi[,"error23[4]"],"error23[5]"=lm1_mcmc_combi[,"error23[5]"],"error23[6]"=lm1_mcmc_combi[,"error23[6]"],
                                            "error23[7]"=lm1_mcmc_combi[,"error23[7]"], "error23[8]"=lm1_mcmc_combi[,"error23[8]"],"error23[9]"=lm1_mcmc_combi[,"error23[9]"],
                                            "error23[10]"=lm1_mcmc_combi[,"error23[10]"],"error23[11]"=lm1_mcmc_combi[,"error23[11]"],"error23[12]"=lm1_mcmc_combi[,"error23[12]"],
                                            "error23[13]"=lm1_mcmc_combi[,"error23[13]"],"error23[14]"=lm1_mcmc_combi[,"error23[14]"],"error23[15]"=lm1_mcmc_combi[,"error23[15]"],
                                            "error23[16]"=lm1_mcmc_combi[,"error23[16]"],"error23[17]"=lm1_mcmc_combi[,"error23[17]"],"error23[18]"=lm1_mcmc_combi[,"error23[18]"],
                                            "error23[19]"=lm1_mcmc_combi[,"error23[19]"],"error23[20]"=lm1_mcmc_combi[,"error23[20]"],"error23[21]"=lm1_mcmc_combi[,"error23[21]"],
                                            "error23[22]"=lm1_mcmc_combi[,"error23[22]"],"error23[23]"=lm1_mcmc_combi[,"error23[23]"],"error23[24]"=lm1_mcmc_combi[,"error23[24]"],
                                            "error23[25]"=lm1_mcmc_combi[,"error23[25]"],"error23[26]"=lm1_mcmc_combi[,"error23[26]"],"error23[27]"=lm1_mcmc_combi[,"error23[27]"],
                                            "error23[28]"=lm1_mcmc_combi[,"error23[28]"],"error23[29]"=lm1_mcmc_combi[,"error23[29]"],"error23[30]"=lm1_mcmc_combi[,"error23[30]"],
                                            "error23[31]"=lm1_mcmc_combi[,"error23[31]"],"error23[32]"=lm1_mcmc_combi[,"error23[32]"],"error23[33]"=lm1_mcmc_combi[,"error23[33]"],
                                            "error23[34]"=lm1_mcmc_combi[,"error23[34]"],"error23[35]"=lm1_mcmc_combi[,"error23[35]"],"error23[36]"=lm1_mcmc_combi[,"error23[36]"],
                                            "error23[37]"=lm1_mcmc_combi[,"error23[37]"],"error23[38]"=lm1_mcmc_combi[,"error23[38]"],"error23[39]"=lm1_mcmc_combi[,"error23[39]"],
                                            "error23[40]"=lm1_mcmc_combi[,"error23[40]"],
                                            "error24[1]"=lm1_mcmc_combi[,"error24[1]"],"error24[2]"=lm1_mcmc_combi[,"error24[2]"],"error24[3]"=lm1_mcmc_combi[,"error24[3]"],
                                            "error24[4]"=lm1_mcmc_combi[,"error24[4]"],"error24[5]"=lm1_mcmc_combi[,"error24[5]"],"error24[6]"=lm1_mcmc_combi[,"error24[6]"],
                                            "error24[7]"=lm1_mcmc_combi[,"error24[7]"], "error24[8]"=lm1_mcmc_combi[,"error24[8]"],"error24[9]"=lm1_mcmc_combi[,"error24[9]"],
                                            "error24[10]"=lm1_mcmc_combi[,"error24[10]"],"error24[11]"=lm1_mcmc_combi[,"error24[11]"],"error24[12]"=lm1_mcmc_combi[,"error24[12]"],
                                            "error24[13]"=lm1_mcmc_combi[,"error24[13]"],"error24[14]"=lm1_mcmc_combi[,"error24[14]"],"error24[15]"=lm1_mcmc_combi[,"error24[15]"],
                                            "error24[16]"=lm1_mcmc_combi[,"error24[16]"],"error24[17]"=lm1_mcmc_combi[,"error24[17]"],"error24[18]"=lm1_mcmc_combi[,"error24[18]"],
                                            "error24[19]"=lm1_mcmc_combi[,"error24[19]"],"error24[20]"=lm1_mcmc_combi[,"error24[20]"],"error24[21]"=lm1_mcmc_combi[,"error24[21]"],
                                            "error24[22]"=lm1_mcmc_combi[,"error24[22]"],"error24[23]"=lm1_mcmc_combi[,"error24[23]"],"error24[24]"=lm1_mcmc_combi[,"error24[24]"],
                                            "error24[25]"=lm1_mcmc_combi[,"error24[25]"],"error24[26]"=lm1_mcmc_combi[,"error24[26]"],"error24[27]"=lm1_mcmc_combi[,"error24[27]"],
                                            "error24[28]"=lm1_mcmc_combi[,"error24[28]"],"error24[29]"=lm1_mcmc_combi[,"error24[29]"],"error24[30]"=lm1_mcmc_combi[,"error24[30]"],
                                            "error24[31]"=lm1_mcmc_combi[,"error24[31]"],"error24[32]"=lm1_mcmc_combi[,"error24[32]"],"error24[33]"=lm1_mcmc_combi[,"error24[33]"],
                                            "error24[34]"=lm1_mcmc_combi[,"error24[34]"],"error24[35]"=lm1_mcmc_combi[,"error24[35]"],"error24[36]"=lm1_mcmc_combi[,"error24[36]"],
                                            "error24[37]"=lm1_mcmc_combi[,"error24[37]"],"error24[38]"=lm1_mcmc_combi[,"error24[38]"],"error24[39]"=lm1_mcmc_combi[,"error24[39]"],
                                            "error24[40]"=lm1_mcmc_combi[,"error24[40]"],
                                            "error25[1]"=lm1_mcmc_combi[,"error25[1]"],"error25[2]"=lm1_mcmc_combi[,"error25[2]"],"error25[3]"=lm1_mcmc_combi[,"error25[3]"],
                                            "error25[4]"=lm1_mcmc_combi[,"error25[4]"],"error25[5]"=lm1_mcmc_combi[,"error25[5]"],"error25[6]"=lm1_mcmc_combi[,"error25[6]"],
                                            "error25[7]"=lm1_mcmc_combi[,"error25[7]"], "error25[8]"=lm1_mcmc_combi[,"error25[8]"],"error25[9]"=lm1_mcmc_combi[,"error25[9]"],
                                            "error25[10]"=lm1_mcmc_combi[,"error25[10]"],"error25[11]"=lm1_mcmc_combi[,"error25[11]"],"error25[12]"=lm1_mcmc_combi[,"error25[12]"],
                                            "error25[13]"=lm1_mcmc_combi[,"error25[13]"],"error25[14]"=lm1_mcmc_combi[,"error25[14]"],"error25[15]"=lm1_mcmc_combi[,"error25[15]"],
                                            "error25[16]"=lm1_mcmc_combi[,"error25[16]"],"error25[17]"=lm1_mcmc_combi[,"error25[17]"],"error25[18]"=lm1_mcmc_combi[,"error25[18]"],
                                            "error25[19]"=lm1_mcmc_combi[,"error25[19]"],"error25[20]"=lm1_mcmc_combi[,"error25[20]"],"error25[21]"=lm1_mcmc_combi[,"error25[21]"],
                                            "error25[22]"=lm1_mcmc_combi[,"error25[22]"],"error25[23]"=lm1_mcmc_combi[,"error25[23]"],"error25[24]"=lm1_mcmc_combi[,"error25[24]"],
                                            "error25[25]"=lm1_mcmc_combi[,"error25[25]"],"error25[26]"=lm1_mcmc_combi[,"error25[26]"],"error25[27]"=lm1_mcmc_combi[,"error25[27]"],
                                            "error25[28]"=lm1_mcmc_combi[,"error25[28]"],"error25[29]"=lm1_mcmc_combi[,"error25[29]"],"error25[30]"=lm1_mcmc_combi[,"error25[30]"],
                                            "error25[31]"=lm1_mcmc_combi[,"error25[31]"],"error25[32]"=lm1_mcmc_combi[,"error25[32]"],"error25[33]"=lm1_mcmc_combi[,"error25[33]"],
                                            "error25[34]"=lm1_mcmc_combi[,"error25[34]"],"error25[35]"=lm1_mcmc_combi[,"error25[35]"],"error25[36]"=lm1_mcmc_combi[,"error25[36]"],
                                            "error25[37]"=lm1_mcmc_combi[,"error25[37]"],"error25[38]"=lm1_mcmc_combi[,"error25[38]"],"error25[39]"=lm1_mcmc_combi[,"error25[39]"],
                                            "error25[40]"=lm1_mcmc_combi[,"error25[40]"],
                                            "error26[1]"=lm1_mcmc_combi[,"error26[1]"],"error26[2]"=lm1_mcmc_combi[,"error26[2]"],"error26[3]"=lm1_mcmc_combi[,"error26[3]"],
                                            "error26[4]"=lm1_mcmc_combi[,"error26[4]"],"error26[5]"=lm1_mcmc_combi[,"error26[5]"],"error26[6]"=lm1_mcmc_combi[,"error26[6]"],
                                            "error26[7]"=lm1_mcmc_combi[,"error26[7]"], "error26[8]"=lm1_mcmc_combi[,"error26[8]"],"error26[9]"=lm1_mcmc_combi[,"error26[9]"],
                                            "error26[10]"=lm1_mcmc_combi[,"error26[10]"],"error26[11]"=lm1_mcmc_combi[,"error26[11]"],"error26[12]"=lm1_mcmc_combi[,"error26[12]"],
                                            "error26[13]"=lm1_mcmc_combi[,"error26[13]"],"error26[14]"=lm1_mcmc_combi[,"error26[14]"],"error26[15]"=lm1_mcmc_combi[,"error26[15]"],
                                            "error26[16]"=lm1_mcmc_combi[,"error26[16]"],"error26[17]"=lm1_mcmc_combi[,"error26[17]"],"error26[18]"=lm1_mcmc_combi[,"error26[18]"],
                                            "error26[19]"=lm1_mcmc_combi[,"error26[19]"],"error26[20]"=lm1_mcmc_combi[,"error26[20]"],"error26[21]"=lm1_mcmc_combi[,"error26[21]"],
                                            "error26[22]"=lm1_mcmc_combi[,"error26[22]"],"error26[23]"=lm1_mcmc_combi[,"error26[23]"],"error26[24]"=lm1_mcmc_combi[,"error26[24]"],
                                            "error26[25]"=lm1_mcmc_combi[,"error26[25]"],"error26[26]"=lm1_mcmc_combi[,"error26[26]"],"error26[27]"=lm1_mcmc_combi[,"error26[27]"],
                                            "error26[28]"=lm1_mcmc_combi[,"error26[28]"],"error26[29]"=lm1_mcmc_combi[,"error26[29]"],"error26[30]"=lm1_mcmc_combi[,"error26[30]"],
                                            "error26[31]"=lm1_mcmc_combi[,"error26[31]"],"error26[32]"=lm1_mcmc_combi[,"error26[32]"],"error26[33]"=lm1_mcmc_combi[,"error26[33]"],
                                            "error26[34]"=lm1_mcmc_combi[,"error26[34]"],"error26[35]"=lm1_mcmc_combi[,"error26[35]"],"error26[36]"=lm1_mcmc_combi[,"error26[36]"],
                                            "error26[37]"=lm1_mcmc_combi[,"error26[37]"],"error26[38]"=lm1_mcmc_combi[,"error26[38]"],"error26[39]"=lm1_mcmc_combi[,"error26[39]"],
                                            "error26[40]"=lm1_mcmc_combi[,"error26[40]"],
                                            "error27[1]"=lm1_mcmc_combi[,"error27[1]"],"error27[2]"=lm1_mcmc_combi[,"error27[2]"],"error27[3]"=lm1_mcmc_combi[,"error27[3]"],
                                            "error27[4]"=lm1_mcmc_combi[,"error27[4]"],"error27[5]"=lm1_mcmc_combi[,"error27[5]"],"error27[6]"=lm1_mcmc_combi[,"error27[6]"],
                                            "error27[7]"=lm1_mcmc_combi[,"error27[7]"], "error27[8]"=lm1_mcmc_combi[,"error27[8]"],"error27[9]"=lm1_mcmc_combi[,"error27[9]"],
                                            "error27[10]"=lm1_mcmc_combi[,"error27[10]"],"error27[11]"=lm1_mcmc_combi[,"error27[11]"],"error27[12]"=lm1_mcmc_combi[,"error27[12]"],
                                            "error27[13]"=lm1_mcmc_combi[,"error27[13]"],"error27[14]"=lm1_mcmc_combi[,"error27[14]"],"error27[15]"=lm1_mcmc_combi[,"error27[15]"],
                                            "error27[16]"=lm1_mcmc_combi[,"error27[16]"],"error27[17]"=lm1_mcmc_combi[,"error27[17]"],"error27[18]"=lm1_mcmc_combi[,"error27[18]"],
                                            "error27[19]"=lm1_mcmc_combi[,"error27[19]"],"error27[20]"=lm1_mcmc_combi[,"error27[20]"],"error27[21]"=lm1_mcmc_combi[,"error27[21]"],
                                            "error27[22]"=lm1_mcmc_combi[,"error27[22]"],"error27[23]"=lm1_mcmc_combi[,"error27[23]"],"error27[24]"=lm1_mcmc_combi[,"error27[24]"],
                                            "error27[25]"=lm1_mcmc_combi[,"error27[25]"],"error27[26]"=lm1_mcmc_combi[,"error27[26]"],"error27[27]"=lm1_mcmc_combi[,"error27[27]"],
                                            "error27[28]"=lm1_mcmc_combi[,"error27[28]"],"error27[29]"=lm1_mcmc_combi[,"error27[29]"],"error27[30]"=lm1_mcmc_combi[,"error27[30]"],
                                            "error27[31]"=lm1_mcmc_combi[,"error27[31]"],"error27[32]"=lm1_mcmc_combi[,"error27[32]"],"error27[33]"=lm1_mcmc_combi[,"error27[33]"],
                                            "error27[34]"=lm1_mcmc_combi[,"error27[34]"],"error27[35]"=lm1_mcmc_combi[,"error27[35]"],"error27[36]"=lm1_mcmc_combi[,"error27[36]"],
                                            "error27[37]"=lm1_mcmc_combi[,"error27[37]"],"error27[38]"=lm1_mcmc_combi[,"error27[38]"],"error27[39]"=lm1_mcmc_combi[,"error27[39]"],
                                            "error27[40]"=lm1_mcmc_combi[,"error27[40]"],
                                            "error28[1]"=lm1_mcmc_combi[,"error28[1]"],"error28[2]"=lm1_mcmc_combi[,"error28[2]"],"error28[3]"=lm1_mcmc_combi[,"error28[3]"],
                                            "error28[4]"=lm1_mcmc_combi[,"error28[4]"],"error28[5]"=lm1_mcmc_combi[,"error28[5]"],"error28[6]"=lm1_mcmc_combi[,"error28[6]"],
                                            "error28[7]"=lm1_mcmc_combi[,"error28[7]"], "error28[8]"=lm1_mcmc_combi[,"error28[8]"],"error28[9]"=lm1_mcmc_combi[,"error28[9]"],
                                            "error28[10]"=lm1_mcmc_combi[,"error28[10]"],"error28[11]"=lm1_mcmc_combi[,"error28[11]"],"error28[12]"=lm1_mcmc_combi[,"error28[12]"],
                                            "error28[13]"=lm1_mcmc_combi[,"error28[13]"],"error28[14]"=lm1_mcmc_combi[,"error28[14]"],"error28[15]"=lm1_mcmc_combi[,"error28[15]"],
                                            "error28[16]"=lm1_mcmc_combi[,"error28[16]"],"error28[17]"=lm1_mcmc_combi[,"error28[17]"],"error28[18]"=lm1_mcmc_combi[,"error28[18]"],
                                            "error28[19]"=lm1_mcmc_combi[,"error28[19]"],"error28[20]"=lm1_mcmc_combi[,"error28[20]"],"error28[21]"=lm1_mcmc_combi[,"error28[21]"],
                                            "error28[22]"=lm1_mcmc_combi[,"error28[22]"],"error28[23]"=lm1_mcmc_combi[,"error28[23]"],"error28[24]"=lm1_mcmc_combi[,"error28[24]"],
                                            "error28[25]"=lm1_mcmc_combi[,"error28[25]"],"error28[26]"=lm1_mcmc_combi[,"error28[26]"],"error28[27]"=lm1_mcmc_combi[,"error28[27]"],
                                            "error28[28]"=lm1_mcmc_combi[,"error28[28]"],"error28[29]"=lm1_mcmc_combi[,"error28[29]"],"error28[30]"=lm1_mcmc_combi[,"error28[30]"],
                                            "error28[31]"=lm1_mcmc_combi[,"error28[31]"],"error28[32]"=lm1_mcmc_combi[,"error28[32]"],"error28[33]"=lm1_mcmc_combi[,"error28[33]"],
                                            "error28[34]"=lm1_mcmc_combi[,"error28[34]"],"error28[35]"=lm1_mcmc_combi[,"error28[35]"],"error28[36]"=lm1_mcmc_combi[,"error28[36]"],
                                            "error28[37]"=lm1_mcmc_combi[,"error28[37]"],"error28[38]"=lm1_mcmc_combi[,"error28[38]"],"error28[39]"=lm1_mcmc_combi[,"error28[39]"],
                                            "error28[40]"=lm1_mcmc_combi[,"error28[40]"],
                                            "error29[1]"=lm1_mcmc_combi[,"error29[1]"],"error29[2]"=lm1_mcmc_combi[,"error29[2]"],"error29[3]"=lm1_mcmc_combi[,"error29[3]"],
                                            "error29[4]"=lm1_mcmc_combi[,"error29[4]"],"error29[5]"=lm1_mcmc_combi[,"error29[5]"],"error29[6]"=lm1_mcmc_combi[,"error29[6]"],
                                            "error29[7]"=lm1_mcmc_combi[,"error29[7]"], "error29[8]"=lm1_mcmc_combi[,"error29[8]"],"error29[9]"=lm1_mcmc_combi[,"error29[9]"],
                                            "error29[10]"=lm1_mcmc_combi[,"error29[10]"],"error29[11]"=lm1_mcmc_combi[,"error29[11]"],"error29[12]"=lm1_mcmc_combi[,"error29[12]"],
                                            "error29[13]"=lm1_mcmc_combi[,"error29[13]"],"error29[14]"=lm1_mcmc_combi[,"error29[14]"],"error29[15]"=lm1_mcmc_combi[,"error29[15]"],
                                            "error29[16]"=lm1_mcmc_combi[,"error29[16]"],"error29[17]"=lm1_mcmc_combi[,"error29[17]"],"error29[18]"=lm1_mcmc_combi[,"error29[18]"],
                                            "error29[19]"=lm1_mcmc_combi[,"error29[19]"],"error29[20]"=lm1_mcmc_combi[,"error29[20]"],"error29[21]"=lm1_mcmc_combi[,"error29[21]"],
                                            "error29[22]"=lm1_mcmc_combi[,"error29[22]"],"error29[23]"=lm1_mcmc_combi[,"error29[23]"],"error29[24]"=lm1_mcmc_combi[,"error29[24]"],
                                            "error29[25]"=lm1_mcmc_combi[,"error29[25]"],"error29[26]"=lm1_mcmc_combi[,"error29[26]"],"error29[27]"=lm1_mcmc_combi[,"error29[27]"],
                                            "error29[28]"=lm1_mcmc_combi[,"error29[28]"],"error29[29]"=lm1_mcmc_combi[,"error29[29]"],"error29[30]"=lm1_mcmc_combi[,"error29[30]"],
                                            "error29[31]"=lm1_mcmc_combi[,"error29[31]"],"error29[32]"=lm1_mcmc_combi[,"error29[32]"],"error29[33]"=lm1_mcmc_combi[,"error29[33]"],
                                            "error29[34]"=lm1_mcmc_combi[,"error29[34]"],"error29[35]"=lm1_mcmc_combi[,"error29[35]"],"error29[36]"=lm1_mcmc_combi[,"error29[36]"],
                                            "error29[37]"=lm1_mcmc_combi[,"error29[37]"],"error29[38]"=lm1_mcmc_combi[,"error29[38]"],"error29[39]"=lm1_mcmc_combi[,"error29[39]"],
                                            "error29[40]"=lm1_mcmc_combi[,"error29[40]"],
                                            "error30[1]"=lm1_mcmc_combi[,"error30[1]"],"error30[2]"=lm1_mcmc_combi[,"error30[2]"],"error30[3]"=lm1_mcmc_combi[,"error30[3]"],
                                            "error30[4]"=lm1_mcmc_combi[,"error30[4]"],"error30[5]"=lm1_mcmc_combi[,"error30[5]"],"error30[6]"=lm1_mcmc_combi[,"error30[6]"],
                                            "error30[7]"=lm1_mcmc_combi[,"error30[7]"], "error30[8]"=lm1_mcmc_combi[,"error30[8]"],"error30[9]"=lm1_mcmc_combi[,"error30[9]"],
                                            "error30[10]"=lm1_mcmc_combi[,"error30[10]"],"error30[11]"=lm1_mcmc_combi[,"error30[11]"],"error30[12]"=lm1_mcmc_combi[,"error30[12]"],
                                            "error30[13]"=lm1_mcmc_combi[,"error30[13]"],"error30[14]"=lm1_mcmc_combi[,"error30[14]"],"error30[15]"=lm1_mcmc_combi[,"error30[15]"],
                                            "error30[16]"=lm1_mcmc_combi[,"error30[16]"],"error30[17]"=lm1_mcmc_combi[,"error30[17]"],"error30[18]"=lm1_mcmc_combi[,"error30[18]"],
                                            "error30[19]"=lm1_mcmc_combi[,"error30[19]"],"error30[20]"=lm1_mcmc_combi[,"error30[20]"],"error30[21]"=lm1_mcmc_combi[,"error30[21]"],
                                            "error30[22]"=lm1_mcmc_combi[,"error30[22]"],"error30[23]"=lm1_mcmc_combi[,"error30[23]"],"error30[24]"=lm1_mcmc_combi[,"error30[24]"],
                                            "error30[25]"=lm1_mcmc_combi[,"error30[25]"],"error30[26]"=lm1_mcmc_combi[,"error30[26]"],"error30[27]"=lm1_mcmc_combi[,"error30[27]"],
                                            "error30[28]"=lm1_mcmc_combi[,"error30[28]"],"error30[29]"=lm1_mcmc_combi[,"error30[29]"],"error30[30]"=lm1_mcmc_combi[,"error30[30]"],
                                            "error30[31]"=lm1_mcmc_combi[,"error30[31]"],"error30[32]"=lm1_mcmc_combi[,"error30[32]"],"error30[33]"=lm1_mcmc_combi[,"error30[33]"],
                                            "error30[34]"=lm1_mcmc_combi[,"error30[34]"],"error30[35]"=lm1_mcmc_combi[,"error30[35]"],"error30[36]"=lm1_mcmc_combi[,"error30[36]"],
                                            "error30[37]"=lm1_mcmc_combi[,"error30[37]"],"error30[38]"=lm1_mcmc_combi[,"error30[38]"],"error30[39]"=lm1_mcmc_combi[,"error30[39]"],
                                            "error30[40]"=lm1_mcmc_combi[,"error30[40]"],
                                            "error31[1]"=lm1_mcmc_combi[,"error31[1]"],"error31[2]"=lm1_mcmc_combi[,"error31[2]"],"error31[3]"=lm1_mcmc_combi[,"error31[3]"],
                                            "error31[4]"=lm1_mcmc_combi[,"error31[4]"],"error31[5]"=lm1_mcmc_combi[,"error31[5]"],"error31[6]"=lm1_mcmc_combi[,"error31[6]"],
                                            "error31[7]"=lm1_mcmc_combi[,"error31[7]"], "error31[8]"=lm1_mcmc_combi[,"error31[8]"],"error31[9]"=lm1_mcmc_combi[,"error31[9]"],
                                            "error31[10]"=lm1_mcmc_combi[,"error31[10]"],"error31[11]"=lm1_mcmc_combi[,"error31[11]"],"error31[12]"=lm1_mcmc_combi[,"error31[12]"],
                                            "error31[13]"=lm1_mcmc_combi[,"error31[13]"],"error31[14]"=lm1_mcmc_combi[,"error31[14]"],"error31[15]"=lm1_mcmc_combi[,"error31[15]"],
                                            "error31[16]"=lm1_mcmc_combi[,"error31[16]"],"error31[17]"=lm1_mcmc_combi[,"error31[17]"],"error31[18]"=lm1_mcmc_combi[,"error31[18]"],
                                            "error31[19]"=lm1_mcmc_combi[,"error31[19]"],"error31[20]"=lm1_mcmc_combi[,"error31[20]"],"error31[21]"=lm1_mcmc_combi[,"error31[21]"],
                                            "error31[22]"=lm1_mcmc_combi[,"error31[22]"],"error31[23]"=lm1_mcmc_combi[,"error31[23]"],"error31[24]"=lm1_mcmc_combi[,"error31[24]"],
                                            "error31[25]"=lm1_mcmc_combi[,"error31[25]"],"error31[26]"=lm1_mcmc_combi[,"error31[26]"],"error31[27]"=lm1_mcmc_combi[,"error31[27]"],
                                            "error31[28]"=lm1_mcmc_combi[,"error31[28]"],"error31[29]"=lm1_mcmc_combi[,"error31[29]"],"error31[30]"=lm1_mcmc_combi[,"error31[30]"],
                                            "error31[31]"=lm1_mcmc_combi[,"error31[31]"],"error31[32]"=lm1_mcmc_combi[,"error31[32]"],"error31[33]"=lm1_mcmc_combi[,"error31[33]"],
                                            "error31[34]"=lm1_mcmc_combi[,"error31[34]"],"error31[35]"=lm1_mcmc_combi[,"error31[35]"],"error31[36]"=lm1_mcmc_combi[,"error31[36]"],
                                            "error31[37]"=lm1_mcmc_combi[,"error31[37]"],"error31[38]"=lm1_mcmc_combi[,"error31[38]"],"error31[39]"=lm1_mcmc_combi[,"error31[39]"],
                                            "error31[40]"=lm1_mcmc_combi[,"error31[40]"],
                                            "error32[1]"=lm1_mcmc_combi[,"error32[1]"],"error32[2]"=lm1_mcmc_combi[,"error32[2]"],"error32[3]"=lm1_mcmc_combi[,"error32[3]"],
                                            "error32[4]"=lm1_mcmc_combi[,"error32[4]"],"error32[5]"=lm1_mcmc_combi[,"error32[5]"],"error32[6]"=lm1_mcmc_combi[,"error32[6]"],
                                            "error32[7]"=lm1_mcmc_combi[,"error32[7]"], "error32[8]"=lm1_mcmc_combi[,"error32[8]"],"error32[9]"=lm1_mcmc_combi[,"error32[9]"],
                                            "error32[10]"=lm1_mcmc_combi[,"error32[10]"],"error32[11]"=lm1_mcmc_combi[,"error32[11]"],"error32[12]"=lm1_mcmc_combi[,"error32[12]"],
                                            "error32[13]"=lm1_mcmc_combi[,"error32[13]"],"error32[14]"=lm1_mcmc_combi[,"error32[14]"],"error32[15]"=lm1_mcmc_combi[,"error32[15]"],
                                            "error32[16]"=lm1_mcmc_combi[,"error32[16]"],"error32[17]"=lm1_mcmc_combi[,"error32[17]"],"error32[18]"=lm1_mcmc_combi[,"error32[18]"],
                                            "error32[19]"=lm1_mcmc_combi[,"error32[19]"],"error32[20]"=lm1_mcmc_combi[,"error32[20]"],"error32[21]"=lm1_mcmc_combi[,"error32[21]"],
                                            "error32[22]"=lm1_mcmc_combi[,"error32[22]"],"error32[23]"=lm1_mcmc_combi[,"error32[23]"],"error32[24]"=lm1_mcmc_combi[,"error32[24]"],
                                            "error32[25]"=lm1_mcmc_combi[,"error32[25]"],"error32[26]"=lm1_mcmc_combi[,"error32[26]"],"error32[27]"=lm1_mcmc_combi[,"error32[27]"],
                                            "error32[28]"=lm1_mcmc_combi[,"error32[28]"],"error32[29]"=lm1_mcmc_combi[,"error32[29]"],"error32[30]"=lm1_mcmc_combi[,"error32[30]"],
                                            "error32[31]"=lm1_mcmc_combi[,"error32[31]"],"error32[32]"=lm1_mcmc_combi[,"error32[32]"],"error32[33]"=lm1_mcmc_combi[,"error32[33]"],
                                            "error32[34]"=lm1_mcmc_combi[,"error32[34]"],"error32[35]"=lm1_mcmc_combi[,"error32[35]"],"error32[36]"=lm1_mcmc_combi[,"error32[36]"],
                                            "error32[37]"=lm1_mcmc_combi[,"error32[37]"],"error32[38]"=lm1_mcmc_combi[,"error32[38]"],"error32[39]"=lm1_mcmc_combi[,"error32[39]"],
                                            "error32[40]"=lm1_mcmc_combi[,"error32[40]"]), 1)) 
lm1_mcmc_fh <- do.call(cbind, rep(list("ar_fh1[1]"=lm1_mcmc_combi[,"ar_fh1[1]"], "ar_fh1[2]"=lm1_mcmc_combi[,"ar_fh1[2]"],"ar_fh1[3]"=lm1_mcmc_combi[,"ar_fh1[3]"],
                                       "ar_fh1[4]"=lm1_mcmc_combi[,"ar_fh1[4]"], "ar_fh1[5]"=lm1_mcmc_combi[,"ar_fh1[5]"],"ar_fh1[6]"=lm1_mcmc_combi[,"ar_fh1[6]"],
                                       "ar_fh1[7]"=lm1_mcmc_combi[,"ar_fh1[7]"], "ar_fh1[8]"=lm1_mcmc_combi[,"ar_fh1[8]"],"ar_fh1[9]"=lm1_mcmc_combi[,"ar_fh1[9]"], 
                                       "ar_fh1[10]"=lm1_mcmc_combi[,"ar_fh1[10]"],"ar_fh1[11]"=lm1_mcmc_combi[,"ar_fh1[11]"],"ar_fh1[12]"=lm1_mcmc_combi[,"ar_fh1[12]"],
                                       "ar_fh1[13]"=lm1_mcmc_combi[,"ar_fh1[13]"],"ar_fh1[14]"=lm1_mcmc_combi[,"ar_fh1[14]"],"ar_fh1[15]"=lm1_mcmc_combi[,"ar_fh1[15]"], 
                                       "ar_fh1[16]"=lm1_mcmc_combi[,"ar_fh1[16]"],"ar_fh1[17]"=lm1_mcmc_combi[,"ar_fh1[17]"],"ar_fh1[18]"=lm1_mcmc_combi[,"ar_fh1[18]"],
                                       "ar_fh1[19]"=lm1_mcmc_combi[,"ar_fh1[19]"],"ar_fh1[20]"=lm1_mcmc_combi[,"ar_fh1[20]"],"ar_fh1[21]"=lm1_mcmc_combi[,"ar_fh1[21]"],
                                       "ar_fh1[22]"=lm1_mcmc_combi[,"ar_fh1[22]"],"ar_fh1[23]"=lm1_mcmc_combi[,"ar_fh1[23]"],"ar_fh1[24]"=lm1_mcmc_combi[,"ar_fh1[24]"],
                                       "ar_fh1[25]"=lm1_mcmc_combi[,"ar_fh1[25]"],"ar_fh1[26]"=lm1_mcmc_combi[,"ar_fh1[26]"],"ar_fh1[27]"=lm1_mcmc_combi[,"ar_fh1[27]"], 
                                       "ar_fh1[28]"=lm1_mcmc_combi[,"ar_fh1[28]"],"ar_fh1[29]"=lm1_mcmc_combi[,"ar_fh1[29]"],"ar_fh1[30]"=lm1_mcmc_combi[,"ar_fh1[30]"],
                                       "ar_fh1[31]"=lm1_mcmc_combi[,"ar_fh1[31]"],"ar_fh1[32]"=lm1_mcmc_combi[,"ar_fh1[32]"],"ar_fh1[33]"=lm1_mcmc_combi[,"ar_fh1[33]"],
                                       "ar_fh1[34]"=lm1_mcmc_combi[,"ar_fh1[34]"],"ar_fh1[35]"=lm1_mcmc_combi[,"ar_fh1[35]"],"ar_fh1[36]"=lm1_mcmc_combi[,"ar_fh1[36]"],
                                       "ar_fh1[37]"=lm1_mcmc_combi[,"ar_fh1[37]"],"ar_fh1[38]"=lm1_mcmc_combi[,"ar_fh1[38]"],"ar_fh1[39]"=lm1_mcmc_combi[,"ar_fh1[39]"], 
                                       "ar_fh1[40]"=lm1_mcmc_combi[,"ar_fh1[40]"],
                                       "ar_fh2[1]"=lm1_mcmc_combi[,"ar_fh1[1]"], "ar_fh2[2]"=lm1_mcmc_combi[,"ar_fh2[2]"],"ar_fh2[3]"=lm1_mcmc_combi[,"ar_fh2[3]"],"ar_fh2[4]"=lm1_mcmc_combi[,"ar_fh2[4]"],
                                       "ar_fh2[5]"=lm1_mcmc_combi[,"ar_fh1[5]"], "ar_fh2[6]"=lm1_mcmc_combi[,"ar_fh2[6]"],"ar_fh2[7]"=lm1_mcmc_combi[,"ar_fh2[7]"], "ar_fh2[8]"=lm1_mcmc_combi[,"ar_fh2[8]"],
                                       "ar_fh2[9]"=lm1_mcmc_combi[,"ar_fh1[9]"], "ar_fh2[10]"=lm1_mcmc_combi[,"ar_fh2[10]"],"ar_fh2[11]"=lm1_mcmc_combi[,"ar_fh2[11]"],"ar_fh2[12]"=lm1_mcmc_combi[,"ar_fh2[12]"],
                                       "ar_fh2[13]"=lm1_mcmc_combi[,"ar_fh1[13]"],"ar_fh2[14]"=lm1_mcmc_combi[,"ar_fh2[14]"],"ar_fh2[15]"=lm1_mcmc_combi[,"ar_fh2[15]"],
                                       "ar_fh2[16]"=lm1_mcmc_combi[,"ar_fh2[16]"],"ar_fh2[17]"=lm1_mcmc_combi[,"ar_fh2[17]"],"ar_fh2[18]"=lm1_mcmc_combi[,"ar_fh2[18]"],
                                       "ar_fh2[19]"=lm1_mcmc_combi[,"ar_fh2[19]"],"ar_fh2[20]"=lm1_mcmc_combi[,"ar_fh2[20]"],"ar_fh2[21]"=lm1_mcmc_combi[,"ar_fh2[21]"],
                                       "ar_fh2[22]"=lm1_mcmc_combi[,"ar_fh2[22]"],"ar_fh2[23]"=lm1_mcmc_combi[,"ar_fh2[23]"],"ar_fh2[24]"=lm1_mcmc_combi[,"ar_fh2[24]"],
                                       "ar_fh2[25]"=lm1_mcmc_combi[,"ar_fh2[25]"],"ar_fh2[26]"=lm1_mcmc_combi[,"ar_fh2[26]"],"ar_fh2[27]"=lm1_mcmc_combi[,"ar_fh2[27]"], 
                                       "ar_fh2[28]"=lm1_mcmc_combi[,"ar_fh2[28]"],"ar_fh2[29]"=lm1_mcmc_combi[,"ar_fh2[29]"],"ar_fh2[30]"=lm1_mcmc_combi[,"ar_fh2[30]"],
                                       "ar_fh2[31]"=lm1_mcmc_combi[,"ar_fh2[31]"],"ar_fh2[32]"=lm1_mcmc_combi[,"ar_fh2[32]"],"ar_fh2[33]"=lm1_mcmc_combi[,"ar_fh2[33]"],
                                       "ar_fh2[34]"=lm1_mcmc_combi[,"ar_fh2[34]"],"ar_fh2[35]"=lm1_mcmc_combi[,"ar_fh2[35]"],"ar_fh2[36]"=lm1_mcmc_combi[,"ar_fh2[36]"],
                                       "ar_fh2[37]"=lm1_mcmc_combi[,"ar_fh2[37]"],"ar_fh2[38]"=lm1_mcmc_combi[,"ar_fh2[38]"],"ar_fh2[39]"=lm1_mcmc_combi[,"ar_fh2[39]"], 
                                       "ar_fh2[40]"=lm1_mcmc_combi[,"ar_fh2[40]"],
                                       "ar_fh3[1]"=lm1_mcmc_combi[,"ar_fh3[1]"], "ar_fh3[2]"=lm1_mcmc_combi[,"ar_fh3[2]"],"ar_fh3[3]"=lm1_mcmc_combi[,"ar_fh3[3]"],"ar_fh3[4]"=lm1_mcmc_combi[,"ar_fh3[4]"],
                                       "ar_fh3[5]"=lm1_mcmc_combi[,"ar_fh3[5]"], "ar_fh3[6]"=lm1_mcmc_combi[,"ar_fh3[6]"],"ar_fh3[7]"=lm1_mcmc_combi[,"ar_fh3[7]"], "ar_fh3[8]"=lm1_mcmc_combi[,"ar_fh3[8]"],
                                       "ar_fh3[9]"=lm1_mcmc_combi[,"ar_fh3[9]"], "ar_fh3[10]"=lm1_mcmc_combi[,"ar_fh3[10]"],"ar_fh3[11]"=lm1_mcmc_combi[,"ar_fh3[11]"],"ar_fh3[12]"=lm1_mcmc_combi[,"ar_fh3[12]"],
                                       "ar_fh3[13]"=lm1_mcmc_combi[,"ar_fh3[13]"],"ar_fh3[14]"=lm1_mcmc_combi[,"ar_fh3[14]"],"ar_fh3[15]"=lm1_mcmc_combi[,"ar_fh3[15]"],
                                       "ar_fh3[16]"=lm1_mcmc_combi[,"ar_fh3[16]"],"ar_fh3[17]"=lm1_mcmc_combi[,"ar_fh3[17]"],"ar_fh3[18]"=lm1_mcmc_combi[,"ar_fh3[18]"],
                                       "ar_fh3[19]"=lm1_mcmc_combi[,"ar_fh3[19]"],"ar_fh3[20]"=lm1_mcmc_combi[,"ar_fh3[20]"],"ar_fh3[21]"=lm1_mcmc_combi[,"ar_fh3[21]"],
                                       "ar_fh3[22]"=lm1_mcmc_combi[,"ar_fh3[22]"],"ar_fh3[23]"=lm1_mcmc_combi[,"ar_fh3[23]"],"ar_fh3[24]"=lm1_mcmc_combi[,"ar_fh3[24]"],
                                       "ar_fh3[25]"=lm1_mcmc_combi[,"ar_fh3[25]"],"ar_fh3[26]"=lm1_mcmc_combi[,"ar_fh3[26]"],"ar_fh3[27]"=lm1_mcmc_combi[,"ar_fh3[27]"], 
                                       "ar_fh3[28]"=lm1_mcmc_combi[,"ar_fh3[28]"],"ar_fh3[29]"=lm1_mcmc_combi[,"ar_fh3[29]"],"ar_fh3[30]"=lm1_mcmc_combi[,"ar_fh3[30]"],
                                       "ar_fh3[31]"=lm1_mcmc_combi[,"ar_fh3[31]"],"ar_fh3[32]"=lm1_mcmc_combi[,"ar_fh3[32]"],"ar_fh3[33]"=lm1_mcmc_combi[,"ar_fh3[33]"],
                                       "ar_fh3[34]"=lm1_mcmc_combi[,"ar_fh3[34]"],"ar_fh3[35]"=lm1_mcmc_combi[,"ar_fh3[35]"],"ar_fh3[36]"=lm1_mcmc_combi[,"ar_fh3[36]"],
                                       "ar_fh3[37]"=lm1_mcmc_combi[,"ar_fh3[37]"],"ar_fh3[38]"=lm1_mcmc_combi[,"ar_fh3[38]"],"ar_fh3[39]"=lm1_mcmc_combi[,"ar_fh3[39]"], 
                                       "ar_fh3[40]"=lm1_mcmc_combi[,"ar_fh3[40]"],
                                       "ar_fh4[1]"=lm1_mcmc_combi[,"ar_fh4[1]"], "ar_fh4[2]"=lm1_mcmc_combi[,"ar_fh4[2]"],"ar_fh4[3]"=lm1_mcmc_combi[,"ar_fh4[3]"],"ar_fh4[4]"=lm1_mcmc_combi[,"ar_fh4[4]"],
                                       "ar_fh4[5]"=lm1_mcmc_combi[,"ar_fh4[5]"], "ar_fh4[6]"=lm1_mcmc_combi[,"ar_fh4[6]"],"ar_fh4[7]"=lm1_mcmc_combi[,"ar_fh4[7]"], "ar_fh4[8]"=lm1_mcmc_combi[,"ar_fh4[8]"],
                                       "ar_fh4[9]"=lm1_mcmc_combi[,"ar_fh4[9]"], "ar_fh4[10]"=lm1_mcmc_combi[,"ar_fh4[10]"],"ar_fh4[11]"=lm1_mcmc_combi[,"ar_fh4[11]"],"ar_fh4[12]"=lm1_mcmc_combi[,"ar_fh4[12]"],
                                       "ar_fh4[13]"=lm1_mcmc_combi[,"ar_fh4[13]"],"ar_fh4[14]"=lm1_mcmc_combi[,"ar_fh4[14]"],"ar_fh4[15]"=lm1_mcmc_combi[,"ar_fh4[15]"],
                                       "ar_fh4[16]"=lm1_mcmc_combi[,"ar_fh4[16]"],"ar_fh4[17]"=lm1_mcmc_combi[,"ar_fh4[17]"],"ar_fh4[18]"=lm1_mcmc_combi[,"ar_fh4[18]"],
                                       "ar_fh4[19]"=lm1_mcmc_combi[,"ar_fh4[19]"],"ar_fh4[20]"=lm1_mcmc_combi[,"ar_fh4[20]"],"ar_fh4[21]"=lm1_mcmc_combi[,"ar_fh4[21]"],
                                       "ar_fh4[22]"=lm1_mcmc_combi[,"ar_fh4[22]"],"ar_fh4[23]"=lm1_mcmc_combi[,"ar_fh4[23]"],"ar_fh4[24]"=lm1_mcmc_combi[,"ar_fh4[24]"],
                                       "ar_fh4[25]"=lm1_mcmc_combi[,"ar_fh4[25]"],"ar_fh4[26]"=lm1_mcmc_combi[,"ar_fh4[26]"],"ar_fh4[27]"=lm1_mcmc_combi[,"ar_fh4[27]"], 
                                       "ar_fh4[28]"=lm1_mcmc_combi[,"ar_fh4[28]"],"ar_fh4[29]"=lm1_mcmc_combi[,"ar_fh4[29]"],"ar_fh4[30]"=lm1_mcmc_combi[,"ar_fh4[30]"],
                                       "ar_fh4[31]"=lm1_mcmc_combi[,"ar_fh4[31]"],"ar_fh4[32]"=lm1_mcmc_combi[,"ar_fh4[32]"],"ar_fh4[33]"=lm1_mcmc_combi[,"ar_fh4[33]"],
                                       "ar_fh4[34]"=lm1_mcmc_combi[,"ar_fh4[34]"],"ar_fh4[35]"=lm1_mcmc_combi[,"ar_fh4[35]"],"ar_fh4[36]"=lm1_mcmc_combi[,"ar_fh4[36]"],
                                       "ar_fh4[37]"=lm1_mcmc_combi[,"ar_fh4[37]"],"ar_fh4[38]"=lm1_mcmc_combi[,"ar_fh4[38]"],"ar_fh4[39]"=lm1_mcmc_combi[,"ar_fh4[39]"], 
                                       "ar_fh4[40]"=lm1_mcmc_combi[,"ar_fh4[40]"],
                                       "ar_fh5[1]"=lm1_mcmc_combi[,"ar_fh5[1]"], "ar_fh5[2]"=lm1_mcmc_combi[,"ar_fh5[2]"],"ar_fh5[3]"=lm1_mcmc_combi[,"ar_fh5[3]"],"ar_fh5[4]"=lm1_mcmc_combi[,"ar_fh5[4]"],
                                       "ar_fh5[5]"=lm1_mcmc_combi[,"ar_fh5[5]"], "ar_fh5[6]"=lm1_mcmc_combi[,"ar_fh5[6]"],"ar_fh5[7]"=lm1_mcmc_combi[,"ar_fh5[7]"], "ar_fh5[8]"=lm1_mcmc_combi[,"ar_fh5[8]"],
                                       "ar_fh5[9]"=lm1_mcmc_combi[,"ar_fh5[9]"], "ar_fh5[10]"=lm1_mcmc_combi[,"ar_fh5[10]"],"ar_fh5[11]"=lm1_mcmc_combi[,"ar_fh5[11]"],"ar_fh5[12]"=lm1_mcmc_combi[,"ar_fh5[12]"],
                                       "ar_fh5[13]"=lm1_mcmc_combi[,"ar_fh5[13]"],"ar_fh5[14]"=lm1_mcmc_combi[,"ar_fh5[14]"],"ar_fh5[15]"=lm1_mcmc_combi[,"ar_fh5[15]"],
                                       "ar_fh5[16]"=lm1_mcmc_combi[,"ar_fh5[16]"],"ar_fh5[17]"=lm1_mcmc_combi[,"ar_fh5[17]"],"ar_fh5[18]"=lm1_mcmc_combi[,"ar_fh5[18]"],
                                       "ar_fh5[19]"=lm1_mcmc_combi[,"ar_fh5[19]"],"ar_fh5[20]"=lm1_mcmc_combi[,"ar_fh5[20]"],"ar_fh5[21]"=lm1_mcmc_combi[,"ar_fh5[21]"],
                                       "ar_fh5[22]"=lm1_mcmc_combi[,"ar_fh5[22]"],"ar_fh5[23]"=lm1_mcmc_combi[,"ar_fh5[23]"],"ar_fh5[24]"=lm1_mcmc_combi[,"ar_fh5[24]"],
                                       "ar_fh5[25]"=lm1_mcmc_combi[,"ar_fh5[25]"],"ar_fh5[26]"=lm1_mcmc_combi[,"ar_fh5[26]"],"ar_fh5[27]"=lm1_mcmc_combi[,"ar_fh5[27]"], 
                                       "ar_fh5[28]"=lm1_mcmc_combi[,"ar_fh5[28]"],"ar_fh5[29]"=lm1_mcmc_combi[,"ar_fh5[29]"],"ar_fh5[30]"=lm1_mcmc_combi[,"ar_fh5[30]"],
                                       "ar_fh5[31]"=lm1_mcmc_combi[,"ar_fh5[31]"],"ar_fh5[32]"=lm1_mcmc_combi[,"ar_fh5[32]"],"ar_fh5[33]"=lm1_mcmc_combi[,"ar_fh5[33]"],
                                       "ar_fh5[34]"=lm1_mcmc_combi[,"ar_fh5[34]"],"ar_fh5[35]"=lm1_mcmc_combi[,"ar_fh5[35]"],"ar_fh5[36]"=lm1_mcmc_combi[,"ar_fh5[36]"],
                                       "ar_fh5[37]"=lm1_mcmc_combi[,"ar_fh5[37]"],"ar_fh5[38]"=lm1_mcmc_combi[,"ar_fh5[38]"],"ar_fh5[39]"=lm1_mcmc_combi[,"ar_fh5[39]"], 
                                       "ar_fh5[40]"=lm1_mcmc_combi[,"ar_fh5[40]"],
                                       "ar_fh6[1]"=lm1_mcmc_combi[,"ar_fh6[1]"], "ar_fh6[2]"=lm1_mcmc_combi[,"ar_fh6[2]"],"ar_fh6[3]"=lm1_mcmc_combi[,"ar_fh6[3]"],"ar_fh6[4]"=lm1_mcmc_combi[,"ar_fh6[4]"],
                                       "ar_fh6[5]"=lm1_mcmc_combi[,"ar_fh6[5]"], "ar_fh6[6]"=lm1_mcmc_combi[,"ar_fh6[6]"],"ar_fh6[7]"=lm1_mcmc_combi[,"ar_fh6[7]"], "ar_fh6[8]"=lm1_mcmc_combi[,"ar_fh6[8]"],
                                       "ar_fh6[9]"=lm1_mcmc_combi[,"ar_fh6[9]"], "ar_fh6[10]"=lm1_mcmc_combi[,"ar_fh6[10]"],"ar_fh6[11]"=lm1_mcmc_combi[,"ar_fh6[11]"],"ar_fh6[12]"=lm1_mcmc_combi[,"ar_fh6[12]"],
                                       "ar_fh6[13]"=lm1_mcmc_combi[,"ar_fh6[13]"],"ar_fh6[14]"=lm1_mcmc_combi[,"ar_fh6[14]"],"ar_fh6[15]"=lm1_mcmc_combi[,"ar_fh6[15]"],    
                                       "ar_fh6[16]"=lm1_mcmc_combi[,"ar_fh6[16]"],"ar_fh6[17]"=lm1_mcmc_combi[,"ar_fh6[17]"],"ar_fh6[18]"=lm1_mcmc_combi[,"ar_fh6[18]"],
                                       "ar_fh6[19]"=lm1_mcmc_combi[,"ar_fh6[19]"],"ar_fh6[20]"=lm1_mcmc_combi[,"ar_fh6[20]"],"ar_fh6[21]"=lm1_mcmc_combi[,"ar_fh6[21]"],
                                       "ar_fh6[22]"=lm1_mcmc_combi[,"ar_fh6[22]"],"ar_fh6[23]"=lm1_mcmc_combi[,"ar_fh6[23]"],"ar_fh6[24]"=lm1_mcmc_combi[,"ar_fh6[24]"],
                                       "ar_fh6[25]"=lm1_mcmc_combi[,"ar_fh6[25]"],"ar_fh6[26]"=lm1_mcmc_combi[,"ar_fh6[26]"],"ar_fh6[27]"=lm1_mcmc_combi[,"ar_fh6[27]"], 
                                       "ar_fh6[28]"=lm1_mcmc_combi[,"ar_fh6[28]"],"ar_fh6[29]"=lm1_mcmc_combi[,"ar_fh6[29]"],"ar_fh6[30]"=lm1_mcmc_combi[,"ar_fh6[30]"],
                                       "ar_fh6[31]"=lm1_mcmc_combi[,"ar_fh6[31]"],"ar_fh6[32]"=lm1_mcmc_combi[,"ar_fh6[32]"],"ar_fh6[33]"=lm1_mcmc_combi[,"ar_fh6[33]"],
                                       "ar_fh6[34]"=lm1_mcmc_combi[,"ar_fh6[34]"],"ar_fh6[35]"=lm1_mcmc_combi[,"ar_fh6[35]"],"ar_fh6[36]"=lm1_mcmc_combi[,"ar_fh6[36]"],
                                       "ar_fh6[37]"=lm1_mcmc_combi[,"ar_fh6[37]"],"ar_fh6[38]"=lm1_mcmc_combi[,"ar_fh6[38]"],"ar_fh6[39]"=lm1_mcmc_combi[,"ar_fh6[39]"], 
                                       "ar_fh6[40]"=lm1_mcmc_combi[,"ar_fh6[40]"],                                  
                                       "sigma_j_fh1[1]"=lm1_mcmc_combi[,"sigma_j_fh1[1]"],"sigma_j_fh1[2]"=lm1_mcmc_combi[,"sigma_j_fh1[2]"],"sigma_j_fh1[3]"=lm1_mcmc_combi[,"sigma_j_fh1[3]"],"sigma_j_fh1[4]"=lm1_mcmc_combi[,"sigma_j_fh1[4]"],
                                       "sigma_j_fh1[5]"=lm1_mcmc_combi[,"sigma_j_fh1[5]"],"sigma_j_fh1[6]"=lm1_mcmc_combi[,"sigma_j_fh1[6]"],"sigma_j_fh1[7]"=lm1_mcmc_combi[,"sigma_j_fh1[7]"],"sigma_j_fh1[8]"=lm1_mcmc_combi[,"sigma_j_fh1[8]"],
                                       "sigma_j_fh1[9]"=lm1_mcmc_combi[,"sigma_j_fh1[9]"],"sigma_j_fh1[10]"=lm1_mcmc_combi[,"sigma_j_fh1[10]"],"sigma_j_fh1[11]"=lm1_mcmc_combi[,"sigma_j_fh1[11]"],"sigma_j_fh1[12]"=lm1_mcmc_combi[,"sigma_j_fh1[12]"],
                                       "sigma_j_fh1[13]"=lm1_mcmc_combi[,"sigma_j_fh1[13]"],"sigma_j_fh1[14]"=lm1_mcmc_combi[,"sigma_j_fh1[14]"],"sigma_j_fh1[15]"=lm1_mcmc_combi[,"sigma_j_fh1[15]"],
                                       "sigma_j_fh1[16]"=lm1_mcmc_combi[,"sigma_j_fh1[16]"],"sigma_j_fh1[17]"=lm1_mcmc_combi[,"sigma_j_fh1[17]"],"sigma_j_fh1[18]"=lm1_mcmc_combi[,"sigma_j_fh1[18]"],
                                       "sigma_j_fh1[19]"=lm1_mcmc_combi[,"sigma_j_fh1[19]"],"sigma_j_fh1[20]"=lm1_mcmc_combi[,"sigma_j_fh1[20]"], "sigma_j_fh1[21]"=lm1_mcmc_combi[,"sigma_j_fh1[21]"],
                                       "sigma_j_fh1[22]"=lm1_mcmc_combi[,"sigma_j_fh1[22]"],"sigma_j_fh1[23]"=lm1_mcmc_combi[,"sigma_j_fh1[23]"],"sigma_j_fh1[24]"=lm1_mcmc_combi[,"sigma_j_fh1[24]"],
                                       "sigma_j_fh1[25]"=lm1_mcmc_combi[,"sigma_j_fh1[25]"],"sigma_j_fh1[26]"=lm1_mcmc_combi[,"sigma_j_fh1[26]"],"sigma_j_fh1[27]"=lm1_mcmc_combi[,"sigma_j_fh1[27]"],
                                       "sigma_j_fh1[28]"=lm1_mcmc_combi[,"sigma_j_fh1[28]"],"sigma_j_fh1[29]"=lm1_mcmc_combi[,"sigma_j_fh1[29]"],"sigma_j_fh1[30]"=lm1_mcmc_combi[,"sigma_j_fh1[30]"],
                                       "sigma_j_fh1[31]"=lm1_mcmc_combi[,"sigma_j_fh1[31]"],"sigma_j_fh1[32]"=lm1_mcmc_combi[,"sigma_j_fh1[32]"],"sigma_j_fh1[33]"=lm1_mcmc_combi[,"sigma_j_fh1[33]"],
                                       "sigma_j_fh1[34]"=lm1_mcmc_combi[,"sigma_j_fh1[34]"],"sigma_j_fh1[35]"=lm1_mcmc_combi[,"sigma_j_fh1[35]"],"sigma_j_fh1[36]"=lm1_mcmc_combi[,"sigma_j_fh1[36]"],
                                       "sigma_j_fh1[37]"=lm1_mcmc_combi[,"sigma_j_fh1[37]"],"sigma_j_fh1[38]"=lm1_mcmc_combi[,"sigma_j_fh1[38]"],"sigma_j_fh1[39]"=lm1_mcmc_combi[,"sigma_j_fh1[39]"],
                                       "sigma_j_fh1[40]"=lm1_mcmc_combi[,"sigma_j_fh1[40]"],
                                       "sigma_j_fh2[1]"=lm1_mcmc_combi[,"sigma_j_fh2[1]"],"sigma_j_fh2[2]"=lm1_mcmc_combi[,"sigma_j_fh2[2]"],"sigma_j_fh2[3]"=lm1_mcmc_combi[,"sigma_j_fh2[3]"],"sigma_j_fh2[4]"=lm1_mcmc_combi[,"sigma_j_fh2[4]"],
                                       "sigma_j_fh2[5]"=lm1_mcmc_combi[,"sigma_j_fh2[5]"],"sigma_j_fh2[6]"=lm1_mcmc_combi[,"sigma_j_fh2[6]"],"sigma_j_fh2[7]"=lm1_mcmc_combi[,"sigma_j_fh2[7]"],"sigma_j_fh2[8]"=lm1_mcmc_combi[,"sigma_j_fh2[8]"],
                                       "sigma_j_fh2[9]"=lm1_mcmc_combi[,"sigma_j_fh2[9]"],"sigma_j_fh2[10]"=lm1_mcmc_combi[,"sigma_j_fh2[10]"],"sigma_j_fh2[11]"=lm1_mcmc_combi[,"sigma_j_fh2[11]"],"sigma_j_fh2[12]"=lm1_mcmc_combi[,"sigma_j_fh2[12]"],
                                       "sigma_j_fh2[13]"=lm1_mcmc_combi[,"sigma_j_fh2[13]"],"sigma_j_fh2[14]"=lm1_mcmc_combi[,"sigma_j_fh2[14]"],"sigma_j_fh2[15]"=lm1_mcmc_combi[,"sigma_j_fh2[15]"],
                                       "sigma_j_fh2[16]"=lm1_mcmc_combi[,"sigma_j_fh2[16]"],"sigma_j_fh2[17]"=lm1_mcmc_combi[,"sigma_j_fh2[17]"],"sigma_j_fh2[18]"=lm1_mcmc_combi[,"sigma_j_fh2[18]"],
                                       "sigma_j_fh2[19]"=lm1_mcmc_combi[,"sigma_j_fh2[19]"],"sigma_j_fh2[20]"=lm1_mcmc_combi[,"sigma_j_fh2[20]"],"sigma_j_fh2[21]"=lm1_mcmc_combi[,"sigma_j_fh2[21]"],
                                       "sigma_j_fh2[22]"=lm1_mcmc_combi[,"sigma_j_fh2[22]"],"sigma_j_fh2[23]"=lm1_mcmc_combi[,"sigma_j_fh2[23]"],"sigma_j_fh2[24]"=lm1_mcmc_combi[,"sigma_j_fh2[24]"],
                                       "sigma_j_fh2[25]"=lm1_mcmc_combi[,"sigma_j_fh2[25]"],"sigma_j_fh2[26]"=lm1_mcmc_combi[,"sigma_j_fh2[26]"],"sigma_j_fh2[27]"=lm1_mcmc_combi[,"sigma_j_fh2[27]"],
                                       "sigma_j_fh2[28]"=lm1_mcmc_combi[,"sigma_j_fh2[28]"],"sigma_j_fh2[29]"=lm1_mcmc_combi[,"sigma_j_fh2[29]"],"sigma_j_fh2[30]"=lm1_mcmc_combi[,"sigma_j_fh2[30]"],
                                       "sigma_j_fh2[31]"=lm1_mcmc_combi[,"sigma_j_fh2[31]"],"sigma_j_fh2[32]"=lm1_mcmc_combi[,"sigma_j_fh2[32]"],"sigma_j_fh2[33]"=lm1_mcmc_combi[,"sigma_j_fh2[33]"],
                                       "sigma_j_fh2[34]"=lm1_mcmc_combi[,"sigma_j_fh2[34]"],"sigma_j_fh2[35]"=lm1_mcmc_combi[,"sigma_j_fh2[35]"],"sigma_j_fh2[36]"=lm1_mcmc_combi[,"sigma_j_fh2[36]"],
                                       "sigma_j_fh2[37]"=lm1_mcmc_combi[,"sigma_j_fh2[37]"],"sigma_j_fh2[38]"=lm1_mcmc_combi[,"sigma_j_fh2[38]"],"sigma_j_fh2[39]"=lm1_mcmc_combi[,"sigma_j_fh2[39]"],
                                       "sigma_j_fh2[40]"=lm1_mcmc_combi[,"sigma_j_fh2[40]"],
                                       "sigma_j_fh3[1]"=lm1_mcmc_combi[,"sigma_j_fh3[1]"],"sigma_j_fh3[2]"=lm1_mcmc_combi[,"sigma_j_fh3[2]"],"sigma_j_fh3[3]"=lm1_mcmc_combi[,"sigma_j_fh3[3]"],"sigma_j_fh3[4]"=lm1_mcmc_combi[,"sigma_j_fh3[4]"],
                                       "sigma_j_fh3[5]"=lm1_mcmc_combi[,"sigma_j_fh3[5]"],"sigma_j_fh3[6]"=lm1_mcmc_combi[,"sigma_j_fh3[6]"],"sigma_j_fh3[7]"=lm1_mcmc_combi[,"sigma_j_fh3[7]"],"sigma_j_fh3[8]"=lm1_mcmc_combi[,"sigma_j_fh3[8]"],
                                       "sigma_j_fh3[9]"=lm1_mcmc_combi[,"sigma_j_fh3[9]"],"sigma_j_fh3[10]"=lm1_mcmc_combi[,"sigma_j_fh3[10]"],"sigma_j_fh3[11]"=lm1_mcmc_combi[,"sigma_j_fh3[11]"],"sigma_j_fh3[12]"=lm1_mcmc_combi[,"sigma_j_fh3[12]"],
                                       "sigma_j_fh3[13]"=lm1_mcmc_combi[,"sigma_j_fh3[13]"],"sigma_j_fh3[14]"=lm1_mcmc_combi[,"sigma_j_fh3[14]"],"sigma_j_fh3[15]"=lm1_mcmc_combi[,"sigma_j_fh3[15]"],
                                       "sigma_j_fh3[16]"=lm1_mcmc_combi[,"sigma_j_fh3[16]"],"sigma_j_fh3[17]"=lm1_mcmc_combi[,"sigma_j_fh3[17]"],"sigma_j_fh3[18]"=lm1_mcmc_combi[,"sigma_j_fh3[18]"],
                                       "sigma_j_fh3[19]"=lm1_mcmc_combi[,"sigma_j_fh3[19]"],"sigma_j_fh3[20]"=lm1_mcmc_combi[,"sigma_j_fh3[20]"],"sigma_j_fh3[21]"=lm1_mcmc_combi[,"sigma_j_fh3[21]"],
                                       "sigma_j_fh3[22]"=lm1_mcmc_combi[,"sigma_j_fh3[22]"],"sigma_j_fh3[23]"=lm1_mcmc_combi[,"sigma_j_fh3[23]"],"sigma_j_fh3[24]"=lm1_mcmc_combi[,"sigma_j_fh3[24]"],
                                       "sigma_j_fh3[25]"=lm1_mcmc_combi[,"sigma_j_fh3[25]"],"sigma_j_fh3[26]"=lm1_mcmc_combi[,"sigma_j_fh3[26]"],"sigma_j_fh3[27]"=lm1_mcmc_combi[,"sigma_j_fh3[27]"],
                                       "sigma_j_fh3[28]"=lm1_mcmc_combi[,"sigma_j_fh3[28]"],"sigma_j_fh3[29]"=lm1_mcmc_combi[,"sigma_j_fh3[29]"],"sigma_j_fh3[30]"=lm1_mcmc_combi[,"sigma_j_fh3[30]"],
                                       "sigma_j_fh3[31]"=lm1_mcmc_combi[,"sigma_j_fh3[31]"],"sigma_j_fh3[32]"=lm1_mcmc_combi[,"sigma_j_fh3[32]"],"sigma_j_fh3[33]"=lm1_mcmc_combi[,"sigma_j_fh3[33]"],
                                       "sigma_j_fh3[34]"=lm1_mcmc_combi[,"sigma_j_fh3[34]"],"sigma_j_fh3[35]"=lm1_mcmc_combi[,"sigma_j_fh3[35]"],"sigma_j_fh3[36]"=lm1_mcmc_combi[,"sigma_j_fh3[36]"],
                                       "sigma_j_fh3[37]"=lm1_mcmc_combi[,"sigma_j_fh3[37]"],"sigma_j_fh3[38]"=lm1_mcmc_combi[,"sigma_j_fh3[38]"],"sigma_j_fh3[39]"=lm1_mcmc_combi[,"sigma_j_fh3[39]"],
                                       "sigma_j_fh3[40]"=lm1_mcmc_combi[,"sigma_j_fh3[40]"],
                                       "sigma_j_fh4[1]"=lm1_mcmc_combi[,"sigma_j_fh4[1]"],"sigma_j_fh4[2]"=lm1_mcmc_combi[,"sigma_j_fh4[2]"],"sigma_j_fh4[3]"=lm1_mcmc_combi[,"sigma_j_fh4[3]"],"sigma_j_fh4[4]"=lm1_mcmc_combi[,"sigma_j_fh4[4]"],
                                       "sigma_j_fh4[5]"=lm1_mcmc_combi[,"sigma_j_fh4[5]"],"sigma_j_fh4[6]"=lm1_mcmc_combi[,"sigma_j_fh4[6]"],"sigma_j_fh4[7]"=lm1_mcmc_combi[,"sigma_j_fh4[7]"],"sigma_j_fh4[8]"=lm1_mcmc_combi[,"sigma_j_fh4[8]"],
                                       "sigma_j_fh4[9]"=lm1_mcmc_combi[,"sigma_j_fh4[9]"],"sigma_j_fh4[10]"=lm1_mcmc_combi[,"sigma_j_fh4[10]"],"sigma_j_fh4[11]"=lm1_mcmc_combi[,"sigma_j_fh4[11]"],"sigma_j_fh4[12]"=lm1_mcmc_combi[,"sigma_j_fh4[12]"],
                                       "sigma_j_fh4[13]"=lm1_mcmc_combi[,"sigma_j_fh4[13]"],"sigma_j_fh4[14]"=lm1_mcmc_combi[,"sigma_j_fh4[14]"],"sigma_j_fh4[15]"=lm1_mcmc_combi[,"sigma_j_fh4[15]"],
                                       "sigma_j_fh4[16]"=lm1_mcmc_combi[,"sigma_j_fh4[16]"],"sigma_j_fh4[17]"=lm1_mcmc_combi[,"sigma_j_fh4[17]"],"sigma_j_fh4[18]"=lm1_mcmc_combi[,"sigma_j_fh4[18]"],
                                       "sigma_j_fh4[19]"=lm1_mcmc_combi[,"sigma_j_fh4[19]"],"sigma_j_fh4[20]"=lm1_mcmc_combi[,"sigma_j_fh4[20]"],"sigma_j_fh4[21]"=lm1_mcmc_combi[,"sigma_j_fh4[21]"],
                                       "sigma_j_fh4[22]"=lm1_mcmc_combi[,"sigma_j_fh4[22]"],"sigma_j_fh4[23]"=lm1_mcmc_combi[,"sigma_j_fh4[23]"],"sigma_j_fh4[24]"=lm1_mcmc_combi[,"sigma_j_fh4[24]"],
                                       "sigma_j_fh4[25]"=lm1_mcmc_combi[,"sigma_j_fh4[25]"],"sigma_j_fh4[26]"=lm1_mcmc_combi[,"sigma_j_fh4[26]"],"sigma_j_fh4[27]"=lm1_mcmc_combi[,"sigma_j_fh4[27]"],
                                       "sigma_j_fh4[28]"=lm1_mcmc_combi[,"sigma_j_fh4[28]"],"sigma_j_fh4[29]"=lm1_mcmc_combi[,"sigma_j_fh4[29]"],"sigma_j_fh4[30]"=lm1_mcmc_combi[,"sigma_j_fh4[30]"],
                                       "sigma_j_fh4[31]"=lm1_mcmc_combi[,"sigma_j_fh4[31]"],"sigma_j_fh4[32]"=lm1_mcmc_combi[,"sigma_j_fh4[32]"],"sigma_j_fh4[33]"=lm1_mcmc_combi[,"sigma_j_fh4[33]"],
                                       "sigma_j_fh4[34]"=lm1_mcmc_combi[,"sigma_j_fh4[34]"],"sigma_j_fh4[35]"=lm1_mcmc_combi[,"sigma_j_fh4[35]"],"sigma_j_fh4[36]"=lm1_mcmc_combi[,"sigma_j_fh4[36]"],
                                       "sigma_j_fh4[37]"=lm1_mcmc_combi[,"sigma_j_fh4[37]"],"sigma_j_fh4[38]"=lm1_mcmc_combi[,"sigma_j_fh4[38]"],"sigma_j_fh4[39]"=lm1_mcmc_combi[,"sigma_j_fh4[39]"],
                                       "sigma_j_fh4[40]"=lm1_mcmc_combi[,"sigma_j_fh4[40]"],
                                       "sigma_j_fh5[1]"=lm1_mcmc_combi[,"sigma_j_fh5[1]"],"sigma_j_fh5[2]"=lm1_mcmc_combi[,"sigma_j_fh5[2]"],"sigma_j_fh5[3]"=lm1_mcmc_combi[,"sigma_j_fh5[3]"],"sigma_j_fh5[4]"=lm1_mcmc_combi[,"sigma_j_fh5[4]"],
                                       "sigma_j_fh5[5]"=lm1_mcmc_combi[,"sigma_j_fh5[5]"],"sigma_j_fh5[6]"=lm1_mcmc_combi[,"sigma_j_fh5[6]"],"sigma_j_fh5[7]"=lm1_mcmc_combi[,"sigma_j_fh5[7]"],"sigma_j_fh5[8]"=lm1_mcmc_combi[,"sigma_j_fh5[8]"],
                                       "sigma_j_fh5[9]"=lm1_mcmc_combi[,"sigma_j_fh5[9]"],"sigma_j_fh5[10]"=lm1_mcmc_combi[,"sigma_j_fh5[10]"],"sigma_j_fh5[11]"=lm1_mcmc_combi[,"sigma_j_fh5[11]"],"sigma_j_fh5[12]"=lm1_mcmc_combi[,"sigma_j_fh5[12]"],
                                       "sigma_j_fh5[13]"=lm1_mcmc_combi[,"sigma_j_fh5[13]"],"sigma_j_fh5[14]"=lm1_mcmc_combi[,"sigma_j_fh5[14]"],"sigma_j_fh5[15]"=lm1_mcmc_combi[,"sigma_j_fh5[15]"],
                                       "sigma_j_fh5[16]"=lm1_mcmc_combi[,"sigma_j_fh5[16]"],"sigma_j_fh5[17]"=lm1_mcmc_combi[,"sigma_j_fh5[17]"],"sigma_j_fh5[18]"=lm1_mcmc_combi[,"sigma_j_fh5[18]"],
                                       "sigma_j_fh5[19]"=lm1_mcmc_combi[,"sigma_j_fh5[19]"],"sigma_j_fh5[20]"=lm1_mcmc_combi[,"sigma_j_fh5[20]"],"sigma_j_fh5[21]"=lm1_mcmc_combi[,"sigma_j_fh5[21]"],
                                       "sigma_j_fh5[22]"=lm1_mcmc_combi[,"sigma_j_fh5[22]"],"sigma_j_fh5[23]"=lm1_mcmc_combi[,"sigma_j_fh5[23]"],"sigma_j_fh5[24]"=lm1_mcmc_combi[,"sigma_j_fh5[24]"],
                                       "sigma_j_fh5[25]"=lm1_mcmc_combi[,"sigma_j_fh5[25]"],"sigma_j_fh5[26]"=lm1_mcmc_combi[,"sigma_j_fh5[26]"],"sigma_j_fh5[27]"=lm1_mcmc_combi[,"sigma_j_fh5[27]"],
                                       "sigma_j_fh5[28]"=lm1_mcmc_combi[,"sigma_j_fh5[28]"],"sigma_j_fh5[29]"=lm1_mcmc_combi[,"sigma_j_fh5[29]"],"sigma_j_fh5[30]"=lm1_mcmc_combi[,"sigma_j_fh5[30]"],
                                       "sigma_j_fh5[31]"=lm1_mcmc_combi[,"sigma_j_fh5[31]"],"sigma_j_fh5[32]"=lm1_mcmc_combi[,"sigma_j_fh5[32]"],"sigma_j_fh5[33]"=lm1_mcmc_combi[,"sigma_j_fh5[33]"],
                                       "sigma_j_fh5[34]"=lm1_mcmc_combi[,"sigma_j_fh5[34]"],"sigma_j_fh5[35]"=lm1_mcmc_combi[,"sigma_j_fh5[35]"],"sigma_j_fh5[36]"=lm1_mcmc_combi[,"sigma_j_fh5[36]"],
                                       "sigma_j_fh5[37]"=lm1_mcmc_combi[,"sigma_j_fh5[37]"],"sigma_j_fh5[38]"=lm1_mcmc_combi[,"sigma_j_fh5[38]"],"sigma_j_fh5[39]"=lm1_mcmc_combi[,"sigma_j_fh5[39]"],
                                       "sigma_j_fh5[40]"=lm1_mcmc_combi[,"sigma_j_fh5[40]"],
                                       "sigma_j_fh6[1]"=lm1_mcmc_combi[,"sigma_j_fh6[1]"],"sigma_j_fh6[2]"=lm1_mcmc_combi[,"sigma_j_fh6[2]"],"sigma_j_fh6[3]"=lm1_mcmc_combi[,"sigma_j_fh6[3]"],"sigma_j_fh6[4]"=lm1_mcmc_combi[,"sigma_j_fh6[4]"],
                                       "sigma_j_fh6[5]"=lm1_mcmc_combi[,"sigma_j_fh6[5]"],"sigma_j_fh6[6]"=lm1_mcmc_combi[,"sigma_j_fh6[6]"],"sigma_j_fh6[7]"=lm1_mcmc_combi[,"sigma_j_fh6[7]"],"sigma_j_fh6[8]"=lm1_mcmc_combi[,"sigma_j_fh6[8]"],
                                       "sigma_j_fh6[9]"=lm1_mcmc_combi[,"sigma_j_fh6[9]"],"sigma_j_fh6[10]"=lm1_mcmc_combi[,"sigma_j_fh6[10]"],"sigma_j_fh6[11]"=lm1_mcmc_combi[,"sigma_j_fh6[11]"],"sigma_j_fh6[12]"=lm1_mcmc_combi[,"sigma_j_fh6[12]"],
                                       "sigma_j_fh6[13]"=lm1_mcmc_combi[,"sigma_j_fh6[13]"],"sigma_j_fh6[14]"=lm1_mcmc_combi[,"sigma_j_fh6[14]"],"sigma_j_fh6[15]"=lm1_mcmc_combi[,"sigma_j_fh6[15]"],
                                       "sigma_j_fh6[16]"=lm1_mcmc_combi[,"sigma_j_fh6[16]"],"sigma_j_fh6[17]"=lm1_mcmc_combi[,"sigma_j_fh6[17]"],"sigma_j_fh6[18]"=lm1_mcmc_combi[,"sigma_j_fh6[18]"],
                                       "sigma_j_fh6[19]"=lm1_mcmc_combi[,"sigma_j_fh6[19]"],"sigma_j_fh6[20]"=lm1_mcmc_combi[,"sigma_j_fh6[20]"],"sigma_j_fh6[21]"=lm1_mcmc_combi[,"sigma_j_fh6[21]"],
                                       "sigma_j_fh6[22]"=lm1_mcmc_combi[,"sigma_j_fh6[22]"],"sigma_j_fh6[23]"=lm1_mcmc_combi[,"sigma_j_fh6[23]"],"sigma_j_fh6[24]"=lm1_mcmc_combi[,"sigma_j_fh6[24]"],
                                       "sigma_j_fh6[25]"=lm1_mcmc_combi[,"sigma_j_fh6[25]"],"sigma_j_fh6[26]"=lm1_mcmc_combi[,"sigma_j_fh6[26]"],"sigma_j_fh6[27]"=lm1_mcmc_combi[,"sigma_j_fh6[27]"],
                                       "sigma_j_fh6[28]"=lm1_mcmc_combi[,"sigma_j_fh6[28]"],"sigma_j_fh6[29]"=lm1_mcmc_combi[,"sigma_j_fh6[29]"],"sigma_j_fh6[30]"=lm1_mcmc_combi[,"sigma_j_fh6[30]"],
                                       "sigma_j_fh6[31]"=lm1_mcmc_combi[,"sigma_j_fh6[31]"],"sigma_j_fh6[32]"=lm1_mcmc_combi[,"sigma_j_fh6[32]"],"sigma_j_fh6[33]"=lm1_mcmc_combi[,"sigma_j_fh6[33]"],
                                       "sigma_j_fh6[34]"=lm1_mcmc_combi[,"sigma_j_fh6[34]"],"sigma_j_fh6[35]"=lm1_mcmc_combi[,"sigma_j_fh6[35]"],"sigma_j_fh6[36]"=lm1_mcmc_combi[,"sigma_j_fh6[36]"],
                                       "sigma_j_fh6[37]"=lm1_mcmc_combi[,"sigma_j_fh6[37]"],"sigma_j_fh6[38]"=lm1_mcmc_combi[,"sigma_j_fh6[38]"],"sigma_j_fh6[39]"=lm1_mcmc_combi[,"sigma_j_fh6[39]"],
                                       "sigma_j_fh6[40]"=lm1_mcmc_combi[,"sigma_j_fh6[40]"]), 1)) 

lm1_mcmc_fv <- do.call(cbind, rep(list("ar_fv1[1]"=lm1_mcmc_combi[,"ar_fv1[1]"], "ar_fv1[2]"=lm1_mcmc_combi[,"ar_fv1[2]"],"ar_fv1[3]"=lm1_mcmc_combi[,"ar_fv1[3]"],
                                       "ar_fv1[4]"=lm1_mcmc_combi[,"ar_fv1[4]"],"ar_fv1[5]"=lm1_mcmc_combi[,"ar_fv1[5]"], "ar_fv1[6]"=lm1_mcmc_combi[,"ar_fv1[6]"],
                                       "ar_fv1[7]"=lm1_mcmc_combi[,"ar_fv1[7]"], "ar_fv1[8]"=lm1_mcmc_combi[,"ar_fv1[8]"],"ar_fv1[9]"=lm1_mcmc_combi[,"ar_fv1[9]"], 
                                       "ar_fv1[10]"=lm1_mcmc_combi[,"ar_fv1[10]"],"ar_fv1[11]"=lm1_mcmc_combi[,"ar_fv1[11]"],"ar_fv1[12]"=lm1_mcmc_combi[,"ar_fv1[12]"],
                                       "ar_fv1[13]"=lm1_mcmc_combi[,"ar_fv1[13]"],"ar_fv1[14]"=lm1_mcmc_combi[,"ar_fv1[14]"],"ar_fv1[15]"=lm1_mcmc_combi[,"ar_fv1[15]"], 
                                       "ar_fv1[16]"=lm1_mcmc_combi[,"ar_fv1[16]"],"ar_fv1[17]"=lm1_mcmc_combi[,"ar_fv1[17]"], "ar_fv1[18]"=lm1_mcmc_combi[,"ar_fv1[18]"],
                                       "ar_fv1[19]"=lm1_mcmc_combi[,"ar_fv1[19]"],"ar_fv1[20]"=lm1_mcmc_combi[,"ar_fv1[20]"],"ar_fv1[21]"=lm1_mcmc_combi[,"ar_fv1[21]"],
                                       "ar_fv1[22]"=lm1_mcmc_combi[,"ar_fv1[22]"],"ar_fv1[23]"=lm1_mcmc_combi[,"ar_fv1[23]"],"ar_fv1[24]"=lm1_mcmc_combi[,"ar_fv1[24]"],
                                       "ar_fv1[25]"=lm1_mcmc_combi[,"ar_fv1[25]"],"ar_fv1[26]"=lm1_mcmc_combi[,"ar_fv1[26]"],"ar_fv1[27]"=lm1_mcmc_combi[,"ar_fv1[27]"], 
                                       "ar_fv1[28]"=lm1_mcmc_combi[,"ar_fv1[28]"],"ar_fv1[29]"=lm1_mcmc_combi[,"ar_fv1[29]"], "ar_fv1[30]"=lm1_mcmc_combi[,"ar_fv1[30]"],
                                       "ar_fv1[31]"=lm1_mcmc_combi[,"ar_fv1[31]"], "ar_fv1[32]"=lm1_mcmc_combi[,"ar_fv1[32]"],"ar_fv1[33]"=lm1_mcmc_combi[,"ar_fv1[33]"],
                                       "ar_fv1[34]"=lm1_mcmc_combi[,"ar_fv1[34]"],"ar_fv1[35]"=lm1_mcmc_combi[,"ar_fv1[35]"],"ar_fv1[36]"=lm1_mcmc_combi[,"ar_fv1[36]"],
                                       "ar_fv1[37]"=lm1_mcmc_combi[,"ar_fv1[37]"], "ar_fv1[38]"=lm1_mcmc_combi[,"ar_fv1[38]"],"ar_fv1[39]"=lm1_mcmc_combi[,"ar_fv1[39]"], 
                                       "ar_fv1[40]"=lm1_mcmc_combi[,"ar_fv1[40]"],
                                       "ar_fv2[1]"=lm1_mcmc_combi[,"ar_fv1[1]"], "ar_fv2[2]"=lm1_mcmc_combi[,"ar_fv2[2]"],"ar_fv2[3]"=lm1_mcmc_combi[,"ar_fv2[3]"],"ar_fv2[4]"=lm1_mcmc_combi[,"ar_fv2[4]"],
                                       "ar_fv2[5]"=lm1_mcmc_combi[,"ar_fv1[5]"], "ar_fv2[6]"=lm1_mcmc_combi[,"ar_fv2[6]"],"ar_fv2[7]"=lm1_mcmc_combi[,"ar_fv2[7]"], "ar_fv2[8]"=lm1_mcmc_combi[,"ar_fv2[8]"],
                                       "ar_fv2[9]"=lm1_mcmc_combi[,"ar_fv1[9]"], "ar_fv2[10]"=lm1_mcmc_combi[,"ar_fv2[10]"],"ar_fv2[11]"=lm1_mcmc_combi[,"ar_fv2[11]"],"ar_fv2[12]"=lm1_mcmc_combi[,"ar_fv2[12]"],
                                       "ar_fv2[13]"=lm1_mcmc_combi[,"ar_fv1[13]"],"ar_fv2[14]"=lm1_mcmc_combi[,"ar_fv2[14]"],"ar_fv2[15]"=lm1_mcmc_combi[,"ar_fv2[15]"],
                                       "ar_fv2[16]"=lm1_mcmc_combi[,"ar_fv2[16]"],"ar_fv2[17]"=lm1_mcmc_combi[,"ar_fv2[17]"],"ar_fv2[18]"=lm1_mcmc_combi[,"ar_fv2[18]"],
                                       "ar_fv2[19]"=lm1_mcmc_combi[,"ar_fv2[19]"],"ar_fv2[20]"=lm1_mcmc_combi[,"ar_fv2[20]"],"ar_fv2[21]"=lm1_mcmc_combi[,"ar_fv2[21]"],
                                       "ar_fv2[22]"=lm1_mcmc_combi[,"ar_fv2[22]"],"ar_fv2[23]"=lm1_mcmc_combi[,"ar_fv2[23]"],"ar_fv2[24]"=lm1_mcmc_combi[,"ar_fv2[24]"],
                                       "ar_fv2[25]"=lm1_mcmc_combi[,"ar_fv2[25]"],"ar_fv2[26]"=lm1_mcmc_combi[,"ar_fv2[26]"],"ar_fv2[27]"=lm1_mcmc_combi[,"ar_fv2[27]"], 
                                       "ar_fv2[28]"=lm1_mcmc_combi[,"ar_fv2[28]"],"ar_fv2[29]"=lm1_mcmc_combi[,"ar_fv2[29]"],"ar_fv2[30]"=lm1_mcmc_combi[,"ar_fv2[30]"],
                                       "ar_fv2[31]"=lm1_mcmc_combi[,"ar_fv2[31]"],"ar_fv2[32]"=lm1_mcmc_combi[,"ar_fv2[32]"],"ar_fv2[33]"=lm1_mcmc_combi[,"ar_fv2[33]"],
                                       "ar_fv2[34]"=lm1_mcmc_combi[,"ar_fv2[34]"],"ar_fv2[35]"=lm1_mcmc_combi[,"ar_fv2[35]"],"ar_fv2[36]"=lm1_mcmc_combi[,"ar_fv2[36]"],
                                       "ar_fv2[37]"=lm1_mcmc_combi[,"ar_fv2[37]"],"ar_fv2[38]"=lm1_mcmc_combi[,"ar_fv2[38]"],"ar_fv2[39]"=lm1_mcmc_combi[,"ar_fv2[39]"], 
                                       "ar_fv2[40]"=lm1_mcmc_combi[,"ar_fv2[40]"],
                                       "ar_fv3[1]"=lm1_mcmc_combi[,"ar_fv3[1]"], "ar_fv3[2]"=lm1_mcmc_combi[,"ar_fv3[2]"],"ar_fv3[3]"=lm1_mcmc_combi[,"ar_fv3[3]"],"ar_fv3[4]"=lm1_mcmc_combi[,"ar_fv3[4]"],
                                       "ar_fv3[5]"=lm1_mcmc_combi[,"ar_fv3[5]"], "ar_fv3[6]"=lm1_mcmc_combi[,"ar_fv3[6]"],"ar_fv3[7]"=lm1_mcmc_combi[,"ar_fv3[7]"], "ar_fv3[8]"=lm1_mcmc_combi[,"ar_fv3[8]"],
                                       "ar_fv3[9]"=lm1_mcmc_combi[,"ar_fv3[9]"], "ar_fv3[10]"=lm1_mcmc_combi[,"ar_fv3[10]"],"ar_fv3[11]"=lm1_mcmc_combi[,"ar_fv3[11]"],"ar_fv3[12]"=lm1_mcmc_combi[,"ar_fv3[12]"],
                                       "ar_fv3[13]"=lm1_mcmc_combi[,"ar_fv3[13]"],"ar_fv3[14]"=lm1_mcmc_combi[,"ar_fv3[14]"],"ar_fv3[15]"=lm1_mcmc_combi[,"ar_fv3[15]"],
                                       "ar_fv3[16]"=lm1_mcmc_combi[,"ar_fv3[16]"],"ar_fv3[17]"=lm1_mcmc_combi[,"ar_fv3[17]"],"ar_fv3[18]"=lm1_mcmc_combi[,"ar_fv3[18]"],
                                       "ar_fv3[19]"=lm1_mcmc_combi[,"ar_fv3[19]"],"ar_fv3[20]"=lm1_mcmc_combi[,"ar_fv3[20]"],"ar_fv3[21]"=lm1_mcmc_combi[,"ar_fv3[21]"],
                                       "ar_fv3[22]"=lm1_mcmc_combi[,"ar_fv3[22]"],"ar_fv3[23]"=lm1_mcmc_combi[,"ar_fv3[23]"],"ar_fv3[24]"=lm1_mcmc_combi[,"ar_fv3[24]"],
                                       "ar_fv3[25]"=lm1_mcmc_combi[,"ar_fv3[25]"],"ar_fv3[26]"=lm1_mcmc_combi[,"ar_fv3[26]"],"ar_fv3[27]"=lm1_mcmc_combi[,"ar_fv3[27]"], 
                                       "ar_fv3[28]"=lm1_mcmc_combi[,"ar_fv3[28]"],"ar_fv3[29]"=lm1_mcmc_combi[,"ar_fv3[29]"],"ar_fv3[30]"=lm1_mcmc_combi[,"ar_fv3[30]"],
                                       "ar_fv3[31]"=lm1_mcmc_combi[,"ar_fv3[31]"],"ar_fv3[32]"=lm1_mcmc_combi[,"ar_fv3[32]"],"ar_fv3[33]"=lm1_mcmc_combi[,"ar_fv3[33]"],
                                       "ar_fv3[34]"=lm1_mcmc_combi[,"ar_fv3[34]"],"ar_fv3[35]"=lm1_mcmc_combi[,"ar_fv3[35]"],"ar_fv3[36]"=lm1_mcmc_combi[,"ar_fv3[36]"],
                                       "ar_fv3[37]"=lm1_mcmc_combi[,"ar_fv3[37]"],"ar_fv3[38]"=lm1_mcmc_combi[,"ar_fv3[38]"],"ar_fv3[39]"=lm1_mcmc_combi[,"ar_fv3[39]"], 
                                       "ar_fv3[40]"=lm1_mcmc_combi[,"ar_fv3[40]"],
                                       "ar_fv4[1]"=lm1_mcmc_combi[,"ar_fv4[1]"], "ar_fv4[2]"=lm1_mcmc_combi[,"ar_fv4[2]"],"ar_fv4[3]"=lm1_mcmc_combi[,"ar_fv4[3]"],"ar_fv4[4]"=lm1_mcmc_combi[,"ar_fv4[4]"],
                                       "ar_fv4[5]"=lm1_mcmc_combi[,"ar_fv4[5]"], "ar_fv4[6]"=lm1_mcmc_combi[,"ar_fv4[6]"],"ar_fv4[7]"=lm1_mcmc_combi[,"ar_fv4[7]"], "ar_fv4[8]"=lm1_mcmc_combi[,"ar_fv4[8]"],
                                       "ar_fv4[9]"=lm1_mcmc_combi[,"ar_fv4[9]"], "ar_fv4[10]"=lm1_mcmc_combi[,"ar_fv4[10]"],"ar_fv4[11]"=lm1_mcmc_combi[,"ar_fv4[11]"],"ar_fv4[12]"=lm1_mcmc_combi[,"ar_fv4[12]"],
                                       "ar_fv4[13]"=lm1_mcmc_combi[,"ar_fv4[13]"],"ar_fv4[14]"=lm1_mcmc_combi[,"ar_fv4[14]"],"ar_fv4[15]"=lm1_mcmc_combi[,"ar_fv4[15]"],
                                       "ar_fv4[16]"=lm1_mcmc_combi[,"ar_fv4[16]"],"ar_fv4[17]"=lm1_mcmc_combi[,"ar_fv4[17]"],"ar_fv4[18]"=lm1_mcmc_combi[,"ar_fv4[18]"],
                                       "ar_fv4[19]"=lm1_mcmc_combi[,"ar_fv4[19]"],"ar_fv4[20]"=lm1_mcmc_combi[,"ar_fv4[20]"],"ar_fv4[21]"=lm1_mcmc_combi[,"ar_fv4[21]"],
                                       "ar_fv4[22]"=lm1_mcmc_combi[,"ar_fv4[22]"],"ar_fv4[23]"=lm1_mcmc_combi[,"ar_fv4[23]"],"ar_fv4[24]"=lm1_mcmc_combi[,"ar_fv4[24]"],
                                       "ar_fv4[25]"=lm1_mcmc_combi[,"ar_fv4[25]"],"ar_fv4[26]"=lm1_mcmc_combi[,"ar_fv4[26]"],"ar_fv4[27]"=lm1_mcmc_combi[,"ar_fv4[27]"], 
                                       "ar_fv4[28]"=lm1_mcmc_combi[,"ar_fv4[28]"],"ar_fv4[29]"=lm1_mcmc_combi[,"ar_fv4[29]"],"ar_fv4[30]"=lm1_mcmc_combi[,"ar_fv4[30]"],
                                       "ar_fv4[31]"=lm1_mcmc_combi[,"ar_fv4[31]"],"ar_fv4[32]"=lm1_mcmc_combi[,"ar_fv4[32]"],"ar_fv4[33]"=lm1_mcmc_combi[,"ar_fv4[33]"],
                                       "ar_fv4[34]"=lm1_mcmc_combi[,"ar_fv4[34]"],"ar_fv4[35]"=lm1_mcmc_combi[,"ar_fv4[35]"],"ar_fv4[36]"=lm1_mcmc_combi[,"ar_fv4[36]"],
                                       "ar_fv4[37]"=lm1_mcmc_combi[,"ar_fv4[37]"],"ar_fv4[38]"=lm1_mcmc_combi[,"ar_fv4[38]"],"ar_fv4[39]"=lm1_mcmc_combi[,"ar_fv4[39]"], 
                                       "ar_fv4[40]"=lm1_mcmc_combi[,"ar_fv4[40]"],
                                       "ar_fv5[1]"=lm1_mcmc_combi[,"ar_fv5[1]"], "ar_fv5[2]"=lm1_mcmc_combi[,"ar_fv5[2]"],"ar_fv5[3]"=lm1_mcmc_combi[,"ar_fv5[3]"],"ar_fv5[4]"=lm1_mcmc_combi[,"ar_fv5[4]"],
                                       "ar_fv5[5]"=lm1_mcmc_combi[,"ar_fv5[5]"], "ar_fv5[6]"=lm1_mcmc_combi[,"ar_fv5[6]"],"ar_fv5[7]"=lm1_mcmc_combi[,"ar_fv5[7]"], "ar_fv5[8]"=lm1_mcmc_combi[,"ar_fv5[8]"],
                                       "ar_fv5[9]"=lm1_mcmc_combi[,"ar_fv5[9]"], "ar_fv5[10]"=lm1_mcmc_combi[,"ar_fv5[10]"],"ar_fv5[11]"=lm1_mcmc_combi[,"ar_fv5[11]"],"ar_fv5[12]"=lm1_mcmc_combi[,"ar_fv5[12]"],
                                       "ar_fv5[13]"=lm1_mcmc_combi[,"ar_fv5[13]"],"ar_fv5[14]"=lm1_mcmc_combi[,"ar_fv5[14]"],"ar_fv5[15]"=lm1_mcmc_combi[,"ar_fv5[15]"],
                                       "ar_fv5[16]"=lm1_mcmc_combi[,"ar_fv5[16]"],"ar_fv5[17]"=lm1_mcmc_combi[,"ar_fv5[17]"],"ar_fv5[18]"=lm1_mcmc_combi[,"ar_fv5[18]"],
                                       "ar_fv5[19]"=lm1_mcmc_combi[,"ar_fv5[19]"],"ar_fv5[20]"=lm1_mcmc_combi[,"ar_fv5[20]"],"ar_fv5[21]"=lm1_mcmc_combi[,"ar_fv5[21]"],
                                       "ar_fv5[22]"=lm1_mcmc_combi[,"ar_fv5[22]"],"ar_fv5[23]"=lm1_mcmc_combi[,"ar_fv5[23]"],"ar_fv5[24]"=lm1_mcmc_combi[,"ar_fv5[24]"],
                                       "ar_fv5[25]"=lm1_mcmc_combi[,"ar_fv5[25]"],"ar_fv5[26]"=lm1_mcmc_combi[,"ar_fv5[26]"],"ar_fv5[27]"=lm1_mcmc_combi[,"ar_fv5[27]"], 
                                       "ar_fv5[28]"=lm1_mcmc_combi[,"ar_fv5[28]"],"ar_fv5[29]"=lm1_mcmc_combi[,"ar_fv5[29]"],"ar_fv5[30]"=lm1_mcmc_combi[,"ar_fv5[30]"],
                                       "ar_fv5[31]"=lm1_mcmc_combi[,"ar_fv5[31]"],"ar_fv5[32]"=lm1_mcmc_combi[,"ar_fv5[32]"],"ar_fv5[33]"=lm1_mcmc_combi[,"ar_fv5[33]"],
                                       "ar_fv5[34]"=lm1_mcmc_combi[,"ar_fv5[34]"],"ar_fv5[35]"=lm1_mcmc_combi[,"ar_fv5[35]"],"ar_fv5[36]"=lm1_mcmc_combi[,"ar_fv5[36]"],
                                       "ar_fv5[37]"=lm1_mcmc_combi[,"ar_fv5[37]"],"ar_fv5[38]"=lm1_mcmc_combi[,"ar_fv5[38]"],"ar_fv5[39]"=lm1_mcmc_combi[,"ar_fv5[39]"], 
                                       "ar_fv5[40]"=lm1_mcmc_combi[,"ar_fv5[40]"],
                                       "sigma_j_fv1[1]"=lm1_mcmc_combi[,"sigma_j_fv1[1]"],"sigma_j_fv1[2]"=lm1_mcmc_combi[,"sigma_j_fv1[2]"],"sigma_j_fv1[3]"=lm1_mcmc_combi[,"sigma_j_fv1[3]"],"sigma_j_fv1[4]"=lm1_mcmc_combi[,"sigma_j_fv1[4]"],
                                       "sigma_j_fv1[5]"=lm1_mcmc_combi[,"sigma_j_fv1[5]"],"sigma_j_fv1[6]"=lm1_mcmc_combi[,"sigma_j_fv1[6]"],"sigma_j_fv1[7]"=lm1_mcmc_combi[,"sigma_j_fv1[7]"],"sigma_j_fv1[8]"=lm1_mcmc_combi[,"sigma_j_fv1[8]"],
                                       "sigma_j_fv1[9]"=lm1_mcmc_combi[,"sigma_j_fv1[9]"],"sigma_j_fv1[10]"=lm1_mcmc_combi[,"sigma_j_fv1[10]"],"sigma_j_fv1[11]"=lm1_mcmc_combi[,"sigma_j_fv1[11]"],"sigma_j_fv1[12]"=lm1_mcmc_combi[,"sigma_j_fv1[12]"],
                                       "sigma_j_fv1[13]"=lm1_mcmc_combi[,"sigma_j_fv1[13]"],"sigma_j_fv1[14]"=lm1_mcmc_combi[,"sigma_j_fv1[14]"],"sigma_j_fv1[15]"=lm1_mcmc_combi[,"sigma_j_fv1[15]"],
                                       "sigma_j_fv1[16]"=lm1_mcmc_combi[,"sigma_j_fv1[16]"],"sigma_j_fv1[17]"=lm1_mcmc_combi[,"sigma_j_fv1[17]"],"sigma_j_fv1[18]"=lm1_mcmc_combi[,"sigma_j_fv1[18]"],
                                       "sigma_j_fv1[19]"=lm1_mcmc_combi[,"sigma_j_fv1[19]"],"sigma_j_fv1[20]"=lm1_mcmc_combi[,"sigma_j_fv1[20]"], "sigma_j_fv1[21]"=lm1_mcmc_combi[,"sigma_j_fv1[21]"],
                                       "sigma_j_fv1[22]"=lm1_mcmc_combi[,"sigma_j_fv1[22]"],"sigma_j_fv1[23]"=lm1_mcmc_combi[,"sigma_j_fv1[23]"],"sigma_j_fv1[24]"=lm1_mcmc_combi[,"sigma_j_fv1[24]"],
                                       "sigma_j_fv1[25]"=lm1_mcmc_combi[,"sigma_j_fv1[25]"],"sigma_j_fv1[26]"=lm1_mcmc_combi[,"sigma_j_fv1[26]"],"sigma_j_fv1[27]"=lm1_mcmc_combi[,"sigma_j_fv1[27]"],
                                       "sigma_j_fv1[28]"=lm1_mcmc_combi[,"sigma_j_fv1[28]"],"sigma_j_fv1[29]"=lm1_mcmc_combi[,"sigma_j_fv1[29]"],"sigma_j_fv1[30]"=lm1_mcmc_combi[,"sigma_j_fv1[30]"],
                                       "sigma_j_fv1[31]"=lm1_mcmc_combi[,"sigma_j_fv1[31]"],"sigma_j_fv1[32]"=lm1_mcmc_combi[,"sigma_j_fv1[32]"],"sigma_j_fv1[33]"=lm1_mcmc_combi[,"sigma_j_fv1[33]"],
                                       "sigma_j_fv1[34]"=lm1_mcmc_combi[,"sigma_j_fv1[34]"],"sigma_j_fv1[35]"=lm1_mcmc_combi[,"sigma_j_fv1[35]"],"sigma_j_fv1[36]"=lm1_mcmc_combi[,"sigma_j_fv1[36]"],
                                       "sigma_j_fv1[37]"=lm1_mcmc_combi[,"sigma_j_fv1[37]"],"sigma_j_fv1[38]"=lm1_mcmc_combi[,"sigma_j_fv1[38]"],"sigma_j_fv1[39]"=lm1_mcmc_combi[,"sigma_j_fv1[39]"],
                                       "sigma_j_fv1[40]"=lm1_mcmc_combi[,"sigma_j_fv1[40]"],
                                       "sigma_j_fv2[1]"=lm1_mcmc_combi[,"sigma_j_fv2[1]"],"sigma_j_fv2[2]"=lm1_mcmc_combi[,"sigma_j_fv2[2]"],"sigma_j_fv2[3]"=lm1_mcmc_combi[,"sigma_j_fv2[3]"],"sigma_j_fv2[4]"=lm1_mcmc_combi[,"sigma_j_fv2[4]"],
                                       "sigma_j_fv2[5]"=lm1_mcmc_combi[,"sigma_j_fv2[5]"],"sigma_j_fv2[6]"=lm1_mcmc_combi[,"sigma_j_fv2[6]"],"sigma_j_fv2[7]"=lm1_mcmc_combi[,"sigma_j_fv2[7]"],"sigma_j_fv2[8]"=lm1_mcmc_combi[,"sigma_j_fv2[8]"],
                                       "sigma_j_fv2[9]"=lm1_mcmc_combi[,"sigma_j_fv2[9]"],"sigma_j_fv2[10]"=lm1_mcmc_combi[,"sigma_j_fv2[10]"],"sigma_j_fv2[11]"=lm1_mcmc_combi[,"sigma_j_fv2[11]"],"sigma_j_fv2[12]"=lm1_mcmc_combi[,"sigma_j_fv2[12]"],
                                       "sigma_j_fv2[13]"=lm1_mcmc_combi[,"sigma_j_fv2[13]"],"sigma_j_fv2[14]"=lm1_mcmc_combi[,"sigma_j_fv2[14]"],"sigma_j_fv2[15]"=lm1_mcmc_combi[,"sigma_j_fv2[15]"],
                                       "sigma_j_fv2[16]"=lm1_mcmc_combi[,"sigma_j_fv2[16]"],"sigma_j_fv2[17]"=lm1_mcmc_combi[,"sigma_j_fv2[17]"],"sigma_j_fv2[18]"=lm1_mcmc_combi[,"sigma_j_fv2[18]"],
                                       "sigma_j_fv2[19]"=lm1_mcmc_combi[,"sigma_j_fv2[19]"],"sigma_j_fv2[20]"=lm1_mcmc_combi[,"sigma_j_fv2[20]"],"sigma_j_fv2[21]"=lm1_mcmc_combi[,"sigma_j_fv2[21]"],
                                       "sigma_j_fv2[22]"=lm1_mcmc_combi[,"sigma_j_fv2[22]"],"sigma_j_fv2[23]"=lm1_mcmc_combi[,"sigma_j_fv2[23]"],"sigma_j_fv2[24]"=lm1_mcmc_combi[,"sigma_j_fv2[24]"],
                                       "sigma_j_fv2[25]"=lm1_mcmc_combi[,"sigma_j_fv2[25]"],"sigma_j_fv2[26]"=lm1_mcmc_combi[,"sigma_j_fv2[26]"],"sigma_j_fv2[27]"=lm1_mcmc_combi[,"sigma_j_fv2[27]"],
                                       "sigma_j_fv2[28]"=lm1_mcmc_combi[,"sigma_j_fv2[28]"],"sigma_j_fv2[29]"=lm1_mcmc_combi[,"sigma_j_fv2[29]"],"sigma_j_fv2[30]"=lm1_mcmc_combi[,"sigma_j_fv2[30]"],
                                       "sigma_j_fv2[31]"=lm1_mcmc_combi[,"sigma_j_fv2[31]"],"sigma_j_fv2[32]"=lm1_mcmc_combi[,"sigma_j_fv2[32]"],"sigma_j_fv2[33]"=lm1_mcmc_combi[,"sigma_j_fv2[33]"],
                                       "sigma_j_fv2[34]"=lm1_mcmc_combi[,"sigma_j_fv2[34]"],"sigma_j_fv2[35]"=lm1_mcmc_combi[,"sigma_j_fv2[35]"],"sigma_j_fv2[36]"=lm1_mcmc_combi[,"sigma_j_fv2[36]"],
                                       "sigma_j_fv2[37]"=lm1_mcmc_combi[,"sigma_j_fv2[37]"],"sigma_j_fv2[38]"=lm1_mcmc_combi[,"sigma_j_fv2[38]"],"sigma_j_fv2[39]"=lm1_mcmc_combi[,"sigma_j_fv2[39]"],
                                       "sigma_j_fv2[40]"=lm1_mcmc_combi[,"sigma_j_fv2[40]"],
                                       "sigma_j_fv3[1]"=lm1_mcmc_combi[,"sigma_j_fv3[1]"],"sigma_j_fv3[2]"=lm1_mcmc_combi[,"sigma_j_fv3[2]"],"sigma_j_fv3[3]"=lm1_mcmc_combi[,"sigma_j_fv3[3]"],"sigma_j_fv3[4]"=lm1_mcmc_combi[,"sigma_j_fv3[4]"],
                                       "sigma_j_fv3[5]"=lm1_mcmc_combi[,"sigma_j_fv3[5]"],"sigma_j_fv3[6]"=lm1_mcmc_combi[,"sigma_j_fv3[6]"],"sigma_j_fv3[7]"=lm1_mcmc_combi[,"sigma_j_fv3[7]"],"sigma_j_fv3[8]"=lm1_mcmc_combi[,"sigma_j_fv3[8]"],
                                       "sigma_j_fv3[9]"=lm1_mcmc_combi[,"sigma_j_fv3[9]"],"sigma_j_fv3[10]"=lm1_mcmc_combi[,"sigma_j_fv3[10]"],"sigma_j_fv3[11]"=lm1_mcmc_combi[,"sigma_j_fv3[11]"],"sigma_j_fv3[12]"=lm1_mcmc_combi[,"sigma_j_fv3[12]"],
                                       "sigma_j_fv3[13]"=lm1_mcmc_combi[,"sigma_j_fv3[13]"],"sigma_j_fv3[14]"=lm1_mcmc_combi[,"sigma_j_fv3[14]"],"sigma_j_fv3[15]"=lm1_mcmc_combi[,"sigma_j_fv3[15]"],
                                       "sigma_j_fv3[16]"=lm1_mcmc_combi[,"sigma_j_fv3[16]"],"sigma_j_fv3[17]"=lm1_mcmc_combi[,"sigma_j_fv3[17]"],"sigma_j_fv3[18]"=lm1_mcmc_combi[,"sigma_j_fv3[18]"],
                                       "sigma_j_fv3[19]"=lm1_mcmc_combi[,"sigma_j_fv3[19]"],"sigma_j_fv3[20]"=lm1_mcmc_combi[,"sigma_j_fv3[20]"],"sigma_j_fv3[21]"=lm1_mcmc_combi[,"sigma_j_fv3[21]"],
                                       "sigma_j_fv3[22]"=lm1_mcmc_combi[,"sigma_j_fv3[22]"],"sigma_j_fv3[23]"=lm1_mcmc_combi[,"sigma_j_fv3[23]"],"sigma_j_fv3[24]"=lm1_mcmc_combi[,"sigma_j_fv3[24]"],
                                       "sigma_j_fv3[25]"=lm1_mcmc_combi[,"sigma_j_fv3[25]"],"sigma_j_fv3[26]"=lm1_mcmc_combi[,"sigma_j_fv3[26]"],"sigma_j_fv3[27]"=lm1_mcmc_combi[,"sigma_j_fv3[27]"],
                                       "sigma_j_fv3[28]"=lm1_mcmc_combi[,"sigma_j_fv3[28]"],"sigma_j_fv3[29]"=lm1_mcmc_combi[,"sigma_j_fv3[29]"],"sigma_j_fv3[30]"=lm1_mcmc_combi[,"sigma_j_fv3[30]"],
                                       "sigma_j_fv3[31]"=lm1_mcmc_combi[,"sigma_j_fv3[31]"],"sigma_j_fv3[32]"=lm1_mcmc_combi[,"sigma_j_fv3[32]"],"sigma_j_fv3[33]"=lm1_mcmc_combi[,"sigma_j_fv3[33]"],
                                       "sigma_j_fv3[34]"=lm1_mcmc_combi[,"sigma_j_fv3[34]"],"sigma_j_fv3[35]"=lm1_mcmc_combi[,"sigma_j_fv3[35]"],"sigma_j_fv3[36]"=lm1_mcmc_combi[,"sigma_j_fv3[36]"],
                                       "sigma_j_fv3[37]"=lm1_mcmc_combi[,"sigma_j_fv3[37]"],"sigma_j_fv3[38]"=lm1_mcmc_combi[,"sigma_j_fv3[38]"],"sigma_j_fv3[39]"=lm1_mcmc_combi[,"sigma_j_fv3[39]"],
                                       "sigma_j_fv3[40]"=lm1_mcmc_combi[,"sigma_j_fv3[40]"],
                                       "sigma_j_fv4[1]"=lm1_mcmc_combi[,"sigma_j_fv4[1]"],"sigma_j_fv4[2]"=lm1_mcmc_combi[,"sigma_j_fv4[2]"],"sigma_j_fv4[3]"=lm1_mcmc_combi[,"sigma_j_fv4[3]"],"sigma_j_fv4[4]"=lm1_mcmc_combi[,"sigma_j_fv4[4]"],
                                       "sigma_j_fv4[5]"=lm1_mcmc_combi[,"sigma_j_fv4[5]"],"sigma_j_fv4[6]"=lm1_mcmc_combi[,"sigma_j_fv4[6]"],"sigma_j_fv4[7]"=lm1_mcmc_combi[,"sigma_j_fv4[7]"],"sigma_j_fv4[8]"=lm1_mcmc_combi[,"sigma_j_fv4[8]"],
                                       "sigma_j_fv4[9]"=lm1_mcmc_combi[,"sigma_j_fv4[9]"],"sigma_j_fv4[10]"=lm1_mcmc_combi[,"sigma_j_fv4[10]"],"sigma_j_fv4[11]"=lm1_mcmc_combi[,"sigma_j_fv4[11]"],"sigma_j_fv4[12]"=lm1_mcmc_combi[,"sigma_j_fv4[12]"],
                                       "sigma_j_fv4[13]"=lm1_mcmc_combi[,"sigma_j_fv4[13]"],"sigma_j_fv4[14]"=lm1_mcmc_combi[,"sigma_j_fv4[14]"],"sigma_j_fv4[15]"=lm1_mcmc_combi[,"sigma_j_fv4[15]"],
                                       "sigma_j_fv4[16]"=lm1_mcmc_combi[,"sigma_j_fv4[16]"],"sigma_j_fv4[17]"=lm1_mcmc_combi[,"sigma_j_fv4[17]"],"sigma_j_fv4[18]"=lm1_mcmc_combi[,"sigma_j_fv4[18]"],
                                       "sigma_j_fv4[19]"=lm1_mcmc_combi[,"sigma_j_fv4[19]"],"sigma_j_fv4[20]"=lm1_mcmc_combi[,"sigma_j_fv4[20]"],"sigma_j_fv4[21]"=lm1_mcmc_combi[,"sigma_j_fv4[21]"],
                                       "sigma_j_fv4[22]"=lm1_mcmc_combi[,"sigma_j_fv4[22]"],"sigma_j_fv4[23]"=lm1_mcmc_combi[,"sigma_j_fv4[23]"],"sigma_j_fv4[24]"=lm1_mcmc_combi[,"sigma_j_fv4[24]"],
                                       "sigma_j_fv4[25]"=lm1_mcmc_combi[,"sigma_j_fv4[25]"],"sigma_j_fv4[26]"=lm1_mcmc_combi[,"sigma_j_fv4[26]"],"sigma_j_fv4[27]"=lm1_mcmc_combi[,"sigma_j_fv4[27]"],
                                       "sigma_j_fv4[28]"=lm1_mcmc_combi[,"sigma_j_fv4[28]"],"sigma_j_fv4[29]"=lm1_mcmc_combi[,"sigma_j_fv4[29]"],"sigma_j_fv4[30]"=lm1_mcmc_combi[,"sigma_j_fv4[30]"],
                                       "sigma_j_fv4[31]"=lm1_mcmc_combi[,"sigma_j_fv4[31]"],"sigma_j_fv4[32]"=lm1_mcmc_combi[,"sigma_j_fv4[32]"],"sigma_j_fv4[33]"=lm1_mcmc_combi[,"sigma_j_fv4[33]"],
                                       "sigma_j_fv4[34]"=lm1_mcmc_combi[,"sigma_j_fv4[34]"],"sigma_j_fv4[35]"=lm1_mcmc_combi[,"sigma_j_fv4[35]"],"sigma_j_fv4[36]"=lm1_mcmc_combi[,"sigma_j_fv4[36]"],
                                       "sigma_j_fv4[37]"=lm1_mcmc_combi[,"sigma_j_fv4[37]"],"sigma_j_fv4[38]"=lm1_mcmc_combi[,"sigma_j_fv4[38]"],"sigma_j_fv4[39]"=lm1_mcmc_combi[,"sigma_j_fv4[39]"],
                                       "sigma_j_fv4[40]"=lm1_mcmc_combi[,"sigma_j_fv4[40]"],
                                       "sigma_j_fv5[1]"=lm1_mcmc_combi[,"sigma_j_fv5[1]"],"sigma_j_fv5[2]"=lm1_mcmc_combi[,"sigma_j_fv5[2]"],"sigma_j_fv5[3]"=lm1_mcmc_combi[,"sigma_j_fv5[3]"],"sigma_j_fv5[4]"=lm1_mcmc_combi[,"sigma_j_fv5[4]"],
                                       "sigma_j_fv5[5]"=lm1_mcmc_combi[,"sigma_j_fv5[5]"],"sigma_j_fv5[6]"=lm1_mcmc_combi[,"sigma_j_fv5[6]"],"sigma_j_fv5[7]"=lm1_mcmc_combi[,"sigma_j_fv5[7]"],"sigma_j_fv5[8]"=lm1_mcmc_combi[,"sigma_j_fv5[8]"],
                                       "sigma_j_fv5[9]"=lm1_mcmc_combi[,"sigma_j_fv5[9]"],"sigma_j_fv5[10]"=lm1_mcmc_combi[,"sigma_j_fv5[10]"],"sigma_j_fv5[11]"=lm1_mcmc_combi[,"sigma_j_fv5[11]"],"sigma_j_fv5[12]"=lm1_mcmc_combi[,"sigma_j_fv5[12]"],
                                       "sigma_j_fv5[13]"=lm1_mcmc_combi[,"sigma_j_fv5[13]"],"sigma_j_fv5[14]"=lm1_mcmc_combi[,"sigma_j_fv5[14]"],"sigma_j_fv5[15]"=lm1_mcmc_combi[,"sigma_j_fv5[15]"],
                                       "sigma_j_fv5[16]"=lm1_mcmc_combi[,"sigma_j_fv5[16]"],"sigma_j_fv5[17]"=lm1_mcmc_combi[,"sigma_j_fv5[17]"],"sigma_j_fv5[18]"=lm1_mcmc_combi[,"sigma_j_fv5[18]"],
                                       "sigma_j_fv5[19]"=lm1_mcmc_combi[,"sigma_j_fv5[19]"],"sigma_j_fv5[20]"=lm1_mcmc_combi[,"sigma_j_fv5[20]"],"sigma_j_fv5[21]"=lm1_mcmc_combi[,"sigma_j_fv5[21]"],
                                       "sigma_j_fv5[22]"=lm1_mcmc_combi[,"sigma_j_fv5[22]"],"sigma_j_fv5[23]"=lm1_mcmc_combi[,"sigma_j_fv5[23]"],"sigma_j_fv5[24]"=lm1_mcmc_combi[,"sigma_j_fv5[24]"],
                                       "sigma_j_fv5[25]"=lm1_mcmc_combi[,"sigma_j_fv5[25]"],"sigma_j_fv5[26]"=lm1_mcmc_combi[,"sigma_j_fv5[26]"],"sigma_j_fv5[27]"=lm1_mcmc_combi[,"sigma_j_fv5[27]"],
                                       "sigma_j_fv5[28]"=lm1_mcmc_combi[,"sigma_j_fv5[28]"],"sigma_j_fv5[29]"=lm1_mcmc_combi[,"sigma_j_fv5[29]"],"sigma_j_fv5[30]"=lm1_mcmc_combi[,"sigma_j_fv5[30]"],
                                       "sigma_j_fv5[31]"=lm1_mcmc_combi[,"sigma_j_fv5[31]"],"sigma_j_fv5[32]"=lm1_mcmc_combi[,"sigma_j_fv5[32]"],"sigma_j_fv5[33]"=lm1_mcmc_combi[,"sigma_j_fv5[33]"],
                                       "sigma_j_fv5[34]"=lm1_mcmc_combi[,"sigma_j_fv5[34]"],"sigma_j_fv5[35]"=lm1_mcmc_combi[,"sigma_j_fv5[35]"],"sigma_j_fv5[36]"=lm1_mcmc_combi[,"sigma_j_fv5[36]"],
                                       "sigma_j_fv5[37]"=lm1_mcmc_combi[,"sigma_j_fv5[37]"],"sigma_j_fv5[38]"=lm1_mcmc_combi[,"sigma_j_fv5[38]"],"sigma_j_fv5[39]"=lm1_mcmc_combi[,"sigma_j_fv5[39]"],
                                       "sigma_j_fv5[40]"=lm1_mcmc_combi[,"sigma_j_fv5[40]"]), 1)) 
lm1_mcmc_fg <- do.call(cbind, rep(list("ar_fg1[1]"=lm1_mcmc_combi[,"ar_fg1[1]"], "ar_fg1[2]"=lm1_mcmc_combi[,"ar_fg1[2]"],"ar_fg1[3]"=lm1_mcmc_combi[,"ar_fg1[3]"],
                                       "ar_fg1[4]"=lm1_mcmc_combi[,"ar_fg1[4]"],"ar_fg1[5]"=lm1_mcmc_combi[,"ar_fg1[5]"], "ar_fg1[6]"=lm1_mcmc_combi[,"ar_fg1[6]"],
                                       "ar_fg1[7]"=lm1_mcmc_combi[,"ar_fg1[7]"], "ar_fg1[8]"=lm1_mcmc_combi[,"ar_fg1[8]"],"ar_fg1[9]"=lm1_mcmc_combi[,"ar_fg1[9]"], 
                                       "ar_fg1[10]"=lm1_mcmc_combi[,"ar_fg1[10]"],"ar_fg1[11]"=lm1_mcmc_combi[,"ar_fg1[11]"],"ar_fg1[12]"=lm1_mcmc_combi[,"ar_fg1[12]"],
                                       "ar_fg1[13]"=lm1_mcmc_combi[,"ar_fg1[13]"],"ar_fg1[14]"=lm1_mcmc_combi[,"ar_fg1[14]"],"ar_fg1[15]"=lm1_mcmc_combi[,"ar_fg1[15]"], 
                                       "ar_fg1[16]"=lm1_mcmc_combi[,"ar_fg1[16]"],"ar_fg1[17]"=lm1_mcmc_combi[,"ar_fg1[17]"], "ar_fg1[18]"=lm1_mcmc_combi[,"ar_fg1[18]"],
                                       "ar_fg1[19]"=lm1_mcmc_combi[,"ar_fg1[19]"], "ar_fg1[20]"=lm1_mcmc_combi[,"ar_fg1[20]"],"ar_fg1[21]"=lm1_mcmc_combi[,"ar_fg1[21]"],
                                       "ar_fg1[22]"=lm1_mcmc_combi[,"ar_fg1[22]"],"ar_fg1[23]"=lm1_mcmc_combi[,"ar_fg1[23]"],"ar_fg1[24]"=lm1_mcmc_combi[,"ar_fg1[24]"],
                                       "ar_fg1[25]"=lm1_mcmc_combi[,"ar_fg1[25]"],"ar_fg1[26]"=lm1_mcmc_combi[,"ar_fg1[26]"],"ar_fg1[27]"=lm1_mcmc_combi[,"ar_fg1[27]"], 
                                       "ar_fg1[28]"=lm1_mcmc_combi[,"ar_fg1[28]"],"ar_fg1[29]"=lm1_mcmc_combi[,"ar_fg1[29]"], "ar_fg1[30]"=lm1_mcmc_combi[,"ar_fg1[30]"],
                                       "ar_fg1[31]"=lm1_mcmc_combi[,"ar_fg1[31]"], "ar_fg1[32]"=lm1_mcmc_combi[,"ar_fg1[32]"],"ar_fg1[33]"=lm1_mcmc_combi[,"ar_fg1[33]"],
                                       "ar_fg1[34]"=lm1_mcmc_combi[,"ar_fg1[34]"],"ar_fg1[35]"=lm1_mcmc_combi[,"ar_fg1[35]"],"ar_fg1[36]"=lm1_mcmc_combi[,"ar_fg1[36]"],
                                       "ar_fg1[37]"=lm1_mcmc_combi[,"ar_fg1[37]"], "ar_fg1[38]"=lm1_mcmc_combi[,"ar_fg1[38]"],"ar_fg1[39]"=lm1_mcmc_combi[,"ar_fg1[39]"], 
                                       "ar_fg1[40]"=lm1_mcmc_combi[,"ar_fg1[40]"],
                                       "sigma_j_fg1[1]"=lm1_mcmc_combi[,"sigma_j_fg1[1]"],"sigma_j_fg1[2]"=lm1_mcmc_combi[,"sigma_j_fg1[2]"],"sigma_j_fg1[3]"=lm1_mcmc_combi[,"sigma_j_fg1[3]"],"sigma_j_fg1[4]"=lm1_mcmc_combi[,"sigma_j_fg1[4]"],
                                       "sigma_j_fg1[5]"=lm1_mcmc_combi[,"sigma_j_fg1[5]"],"sigma_j_fg1[6]"=lm1_mcmc_combi[,"sigma_j_fg1[6]"],"sigma_j_fg1[7]"=lm1_mcmc_combi[,"sigma_j_fg1[7]"],"sigma_j_fg1[8]"=lm1_mcmc_combi[,"sigma_j_fg1[8]"],
                                       "sigma_j_fg1[9]"=lm1_mcmc_combi[,"sigma_j_fg1[9]"],"sigma_j_fg1[10]"=lm1_mcmc_combi[,"sigma_j_fg1[10]"],"sigma_j_fg1[11]"=lm1_mcmc_combi[,"sigma_j_fg1[11]"],"sigma_j_fg1[12]"=lm1_mcmc_combi[,"sigma_j_fg1[12]"],
                                       "sigma_j_fg1[13]"=lm1_mcmc_combi[,"sigma_j_fg1[13]"],"sigma_j_fg1[14]"=lm1_mcmc_combi[,"sigma_j_fg1[14]"],"sigma_j_fg1[15]"=lm1_mcmc_combi[,"sigma_j_fg1[15]"],
                                       "sigma_j_fg1[16]"=lm1_mcmc_combi[,"sigma_j_fg1[16]"],"sigma_j_fg1[17]"=lm1_mcmc_combi[,"sigma_j_fg1[17]"],"sigma_j_fg1[18]"=lm1_mcmc_combi[,"sigma_j_fg1[18]"],
                                       "sigma_j_fg1[19]"=lm1_mcmc_combi[,"sigma_j_fg1[19]"],"sigma_j_fg1[20]"=lm1_mcmc_combi[,"sigma_j_fg1[20]"], "sigma_j_fg1[21]"=lm1_mcmc_combi[,"sigma_j_fg1[21]"],
                                       "sigma_j_fg1[22]"=lm1_mcmc_combi[,"sigma_j_fg1[22]"],"sigma_j_fg1[23]"=lm1_mcmc_combi[,"sigma_j_fg1[23]"],"sigma_j_fg1[24]"=lm1_mcmc_combi[,"sigma_j_fg1[24]"],
                                       "sigma_j_fg1[25]"=lm1_mcmc_combi[,"sigma_j_fg1[25]"],"sigma_j_fg1[26]"=lm1_mcmc_combi[,"sigma_j_fg1[26]"],"sigma_j_fg1[27]"=lm1_mcmc_combi[,"sigma_j_fg1[27]"],
                                       "sigma_j_fg1[28]"=lm1_mcmc_combi[,"sigma_j_fg1[28]"],"sigma_j_fg1[29]"=lm1_mcmc_combi[,"sigma_j_fg1[29]"],"sigma_j_fg1[30]"=lm1_mcmc_combi[,"sigma_j_fg1[30]"],
                                       "sigma_j_fg1[31]"=lm1_mcmc_combi[,"sigma_j_fg1[31]"],"sigma_j_fg1[32]"=lm1_mcmc_combi[,"sigma_j_fg1[32]"],"sigma_j_fg1[33]"=lm1_mcmc_combi[,"sigma_j_fg1[33]"],
                                       "sigma_j_fg1[34]"=lm1_mcmc_combi[,"sigma_j_fg1[34]"],"sigma_j_fg1[35]"=lm1_mcmc_combi[,"sigma_j_fg1[35]"],"sigma_j_fg1[36]"=lm1_mcmc_combi[,"sigma_j_fg1[36]"],
                                       "sigma_j_fg1[37]"=lm1_mcmc_combi[,"sigma_j_fg1[37]"],"sigma_j_fg1[38]"=lm1_mcmc_combi[,"sigma_j_fg1[38]"],"sigma_j_fg1[39]"=lm1_mcmc_combi[,"sigma_j_fg1[39]"],
                                       "sigma_j_fg1[40]"=lm1_mcmc_combi[,"sigma_j_fg1[40]"]), 1)) 

lm1_mcmc_s <- do.call(cbind, rep(list("ar_s1[1]"=lm1_mcmc_combi[,"ar_s1[1]"], "ar_s1[2]"=lm1_mcmc_combi[,"ar_s1[2]"],"ar_s1[3]"=lm1_mcmc_combi[,"ar_s1[3]"],"ar_s1[4]"=lm1_mcmc_combi[,"ar_s1[4]"],
                                      "ar_s1[5]"=lm1_mcmc_combi[,"ar_s1[5]"], "ar_s1[6]"=lm1_mcmc_combi[,"ar_s1[6]"],"ar_s1[7]"=lm1_mcmc_combi[,"ar_s1[7]"], "ar_s1[8]"=lm1_mcmc_combi[,"ar_s1[8]"],
                                      "ar_s1[9]"=lm1_mcmc_combi[,"ar_s1[9]"], "ar_s1[10]"=lm1_mcmc_combi[,"ar_s1[10]"],"ar_s1[11]"=lm1_mcmc_combi[,"ar_s1[11]"],"ar_s1[12]"=lm1_mcmc_combi[,"ar_s1[12]"],
                                      "ar_s1[13]"=lm1_mcmc_combi[,"ar_s1[13]"],"ar_s1[14]"=lm1_mcmc_combi[,"ar_s1[14]"],"ar_s1[15]"=lm1_mcmc_combi[,"ar_s1[15]"],
                                      "ar_s1[16]"=lm1_mcmc_combi[,"ar_s1[16]"],"ar_s1[17]"=lm1_mcmc_combi[,"ar_s1[17]"],"ar_s1[18]"=lm1_mcmc_combi[,"ar_s1[18]"],
                                      "ar_s1[19]"=lm1_mcmc_combi[,"ar_s1[19]"],"ar_s1[20]"=lm1_mcmc_combi[,"ar_s1[20]"],"ar_s1[21]"=lm1_mcmc_combi[,"ar_s1[21]"],
                                      "ar_s1[22]"=lm1_mcmc_combi[,"ar_s1[22]"],"ar_s1[23]"=lm1_mcmc_combi[,"ar_s1[23]"],"ar_s1[24]"=lm1_mcmc_combi[,"ar_s1[24]"],
                                      "ar_s1[25]"=lm1_mcmc_combi[,"ar_s1[25]"],"ar_s1[26]"=lm1_mcmc_combi[,"ar_s1[26]"],"ar_s1[27]"=lm1_mcmc_combi[,"ar_s1[27]"], 
                                      "ar_s1[28]"=lm1_mcmc_combi[,"ar_s1[28]"],"ar_s1[29]"=lm1_mcmc_combi[,"ar_s1[29]"],"ar_s1[30]"=lm1_mcmc_combi[,"ar_s1[30]"],
                                      "ar_s1[31]"=lm1_mcmc_combi[,"ar_s1[31]"],"ar_s1[32]"=lm1_mcmc_combi[,"ar_s1[32]"],"ar_s1[33]"=lm1_mcmc_combi[,"ar_s1[33]"],
                                      "ar_s1[34]"=lm1_mcmc_combi[,"ar_s1[34]"],"ar_s1[35]"=lm1_mcmc_combi[,"ar_s1[35]"],"ar_s1[36]"=lm1_mcmc_combi[,"ar_s1[36]"],
                                      "ar_s1[37]"=lm1_mcmc_combi[,"ar_s1[37]"],"ar_s1[38]"=lm1_mcmc_combi[,"ar_s1[38]"],"ar_s1[39]"=lm1_mcmc_combi[,"ar_s1[39]"], 
                                      "ar_s1[40]"=lm1_mcmc_combi[,"ar_s1[40]"],
                                      "ar_s2[1]"=lm1_mcmc_combi[,"ar_s2[1]"], "ar_s2[2]"=lm1_mcmc_combi[,"ar_s2[2]"],"ar_s2[3]"=lm1_mcmc_combi[,"ar_s2[3]"],"ar_s2[4]"=lm1_mcmc_combi[,"ar_s2[4]"],
                                      "ar_s2[5]"=lm1_mcmc_combi[,"ar_s2[5]"], "ar_s2[6]"=lm1_mcmc_combi[,"ar_s2[6]"],"ar_s2[7]"=lm1_mcmc_combi[,"ar_s2[7]"], "ar_s2[8]"=lm1_mcmc_combi[,"ar_s2[8]"],
                                      "ar_s2[9]"=lm1_mcmc_combi[,"ar_s2[9]"], "ar_s2[10]"=lm1_mcmc_combi[,"ar_s2[10]"],"ar_s2[11]"=lm1_mcmc_combi[,"ar_s2[11]"],"ar_s2[12]"=lm1_mcmc_combi[,"ar_s2[12]"],
                                      "ar_s2[13]"=lm1_mcmc_combi[,"ar_s2[13]"],"ar_s2[14]"=lm1_mcmc_combi[,"ar_s2[14]"],"ar_s2[15]"=lm1_mcmc_combi[,"ar_s2[15]"],
                                      "ar_s2[16]"=lm1_mcmc_combi[,"ar_s2[16]"],"ar_s2[17]"=lm1_mcmc_combi[,"ar_s2[17]"],"ar_s2[18]"=lm1_mcmc_combi[,"ar_s2[18]"],
                                      "ar_s2[19]"=lm1_mcmc_combi[,"ar_s2[19]"],"ar_s2[20]"=lm1_mcmc_combi[,"ar_s2[20]"],"ar_s2[21]"=lm1_mcmc_combi[,"ar_s2[21]"],
                                      "ar_s2[22]"=lm1_mcmc_combi[,"ar_s2[22]"],"ar_s2[23]"=lm1_mcmc_combi[,"ar_s2[23]"],"ar_s2[24]"=lm1_mcmc_combi[,"ar_s2[24]"],
                                      "ar_s2[25]"=lm1_mcmc_combi[,"ar_s2[25]"],"ar_s2[26]"=lm1_mcmc_combi[,"ar_s2[26]"],"ar_s2[27]"=lm1_mcmc_combi[,"ar_s2[27]"], 
                                      "ar_s2[28]"=lm1_mcmc_combi[,"ar_s2[28]"],"ar_s2[29]"=lm1_mcmc_combi[,"ar_s2[29]"],"ar_s2[30]"=lm1_mcmc_combi[,"ar_s2[30]"],
                                      "ar_s2[31]"=lm1_mcmc_combi[,"ar_s2[31]"],"ar_s2[32]"=lm1_mcmc_combi[,"ar_s2[32]"],"ar_s2[33]"=lm1_mcmc_combi[,"ar_s2[33]"],
                                      "ar_s2[34]"=lm1_mcmc_combi[,"ar_s2[34]"],"ar_s2[35]"=lm1_mcmc_combi[,"ar_s2[35]"],"ar_s2[36]"=lm1_mcmc_combi[,"ar_s2[36]"],
                                      "ar_s2[37]"=lm1_mcmc_combi[,"ar_s2[37]"],"ar_s2[38]"=lm1_mcmc_combi[,"ar_s2[38]"],"ar_s2[39]"=lm1_mcmc_combi[,"ar_s2[39]"], 
                                      "ar_s2[40]"=lm1_mcmc_combi[,"ar_s2[40]"],
                                      "ar_s3[1]"=lm1_mcmc_combi[,"ar_s3[1]"], "ar_s3[2]"=lm1_mcmc_combi[,"ar_s3[2]"],"ar_s3[3]"=lm1_mcmc_combi[,"ar_s3[3]"],"ar_s3[4]"=lm1_mcmc_combi[,"ar_s3[4]"],
                                      "ar_s3[5]"=lm1_mcmc_combi[,"ar_s3[5]"], "ar_s3[6]"=lm1_mcmc_combi[,"ar_s3[6]"],"ar_s3[7]"=lm1_mcmc_combi[,"ar_s3[7]"], "ar_s3[8]"=lm1_mcmc_combi[,"ar_s3[8]"],
                                      "ar_s3[9]"=lm1_mcmc_combi[,"ar_s3[9]"], "ar_s3[10]"=lm1_mcmc_combi[,"ar_s3[10]"],"ar_s3[11]"=lm1_mcmc_combi[,"ar_s3[11]"],"ar_s3[12]"=lm1_mcmc_combi[,"ar_s3[12]"],
                                      "ar_s3[13]"=lm1_mcmc_combi[,"ar_s3[13]"],"ar_s3[14]"=lm1_mcmc_combi[,"ar_s3[14]"],"ar_s3[15]"=lm1_mcmc_combi[,"ar_s3[15]"],
                                      "ar_s3[16]"=lm1_mcmc_combi[,"ar_s3[16]"],"ar_s3[17]"=lm1_mcmc_combi[,"ar_s3[17]"],"ar_s3[18]"=lm1_mcmc_combi[,"ar_s3[18]"],
                                      "ar_s3[19]"=lm1_mcmc_combi[,"ar_s3[19]"],"ar_s3[20]"=lm1_mcmc_combi[,"ar_s3[20]"],"ar_s3[21]"=lm1_mcmc_combi[,"ar_s3[21]"],
                                      "ar_s3[22]"=lm1_mcmc_combi[,"ar_s3[22]"],"ar_s3[23]"=lm1_mcmc_combi[,"ar_s3[23]"],"ar_s3[24]"=lm1_mcmc_combi[,"ar_s3[24]"],
                                      "ar_s3[25]"=lm1_mcmc_combi[,"ar_s3[25]"],"ar_s3[26]"=lm1_mcmc_combi[,"ar_s3[26]"],"ar_s3[27]"=lm1_mcmc_combi[,"ar_s3[27]"], 
                                      "ar_s3[28]"=lm1_mcmc_combi[,"ar_s3[28]"],"ar_s3[29]"=lm1_mcmc_combi[,"ar_s3[29]"],"ar_s3[30]"=lm1_mcmc_combi[,"ar_s3[30]"],
                                      "ar_s3[31]"=lm1_mcmc_combi[,"ar_s3[31]"],"ar_s3[32]"=lm1_mcmc_combi[,"ar_s3[32]"],"ar_s3[33]"=lm1_mcmc_combi[,"ar_s3[33]"],
                                      "ar_s3[34]"=lm1_mcmc_combi[,"ar_s3[34]"],"ar_s3[35]"=lm1_mcmc_combi[,"ar_s3[35]"],"ar_s3[36]"=lm1_mcmc_combi[,"ar_s3[36]"],
                                      "ar_s3[37]"=lm1_mcmc_combi[,"ar_s3[37]"],"ar_s3[38]"=lm1_mcmc_combi[,"ar_s3[38]"],"ar_s3[39]"=lm1_mcmc_combi[,"ar_s3[39]"], 
                                      "ar_s3[40]"=lm1_mcmc_combi[,"ar_s3[40]"],
                                      "ar_s4[1]"=lm1_mcmc_combi[,"ar_s4[1]"], "ar_s4[2]"=lm1_mcmc_combi[,"ar_s4[2]"],"ar_s4[3]"=lm1_mcmc_combi[,"ar_s4[3]"],"ar_s4[4]"=lm1_mcmc_combi[,"ar_s4[4]"],
                                      "ar_s4[5]"=lm1_mcmc_combi[,"ar_s4[5]"], "ar_s4[6]"=lm1_mcmc_combi[,"ar_s4[6]"],"ar_s4[7]"=lm1_mcmc_combi[,"ar_s4[7]"], "ar_s4[8]"=lm1_mcmc_combi[,"ar_s4[8]"],
                                      "ar_s4[9]"=lm1_mcmc_combi[,"ar_s4[9]"], "ar_s4[10]"=lm1_mcmc_combi[,"ar_s4[10]"],"ar_s4[11]"=lm1_mcmc_combi[,"ar_s4[11]"],"ar_s4[12]"=lm1_mcmc_combi[,"ar_s4[12]"],
                                      "ar_s4[13]"=lm1_mcmc_combi[,"ar_s4[13]"],"ar_s4[14]"=lm1_mcmc_combi[,"ar_s4[14]"],"ar_s4[15]"=lm1_mcmc_combi[,"ar_s4[15]"],
                                      "ar_s4[16]"=lm1_mcmc_combi[,"ar_s4[16]"],"ar_s4[17]"=lm1_mcmc_combi[,"ar_s4[17]"],"ar_s4[18]"=lm1_mcmc_combi[,"ar_s4[18]"],
                                      "ar_s4[19]"=lm1_mcmc_combi[,"ar_s4[19]"],"ar_s4[20]"=lm1_mcmc_combi[,"ar_s4[20]"],"ar_s4[21]"=lm1_mcmc_combi[,"ar_s4[21]"],
                                      "ar_s4[22]"=lm1_mcmc_combi[,"ar_s4[22]"],"ar_s4[23]"=lm1_mcmc_combi[,"ar_s4[23]"],"ar_s4[24]"=lm1_mcmc_combi[,"ar_s4[24]"],
                                      "ar_s4[25]"=lm1_mcmc_combi[,"ar_s4[25]"],"ar_s4[26]"=lm1_mcmc_combi[,"ar_s4[26]"],"ar_s4[27]"=lm1_mcmc_combi[,"ar_s4[27]"], 
                                      "ar_s4[28]"=lm1_mcmc_combi[,"ar_s4[28]"],"ar_s4[29]"=lm1_mcmc_combi[,"ar_s4[29]"],"ar_s4[30]"=lm1_mcmc_combi[,"ar_s4[30]"],
                                      "ar_s4[31]"=lm1_mcmc_combi[,"ar_s4[31]"],"ar_s4[32]"=lm1_mcmc_combi[,"ar_s4[32]"],"ar_s4[33]"=lm1_mcmc_combi[,"ar_s4[33]"],
                                      "ar_s4[34]"=lm1_mcmc_combi[,"ar_s4[34]"],"ar_s4[35]"=lm1_mcmc_combi[,"ar_s4[35]"],"ar_s4[36]"=lm1_mcmc_combi[,"ar_s4[36]"],
                                      "ar_s4[37]"=lm1_mcmc_combi[,"ar_s4[37]"],"ar_s4[38]"=lm1_mcmc_combi[,"ar_s4[38]"],"ar_s4[39]"=lm1_mcmc_combi[,"ar_s4[39]"], 
                                      "ar_s4[40]"=lm1_mcmc_combi[,"ar_s4[40]"],
                                      "ar_s5[1]"=lm1_mcmc_combi[,"ar_s5[1]"], "ar_s5[2]"=lm1_mcmc_combi[,"ar_s5[2]"],"ar_s5[3]"=lm1_mcmc_combi[,"ar_s5[3]"],"ar_s5[4]"=lm1_mcmc_combi[,"ar_s5[4]"],
                                      "ar_s5[5]"=lm1_mcmc_combi[,"ar_s5[5]"], "ar_s5[6]"=lm1_mcmc_combi[,"ar_s5[6]"],"ar_s5[7]"=lm1_mcmc_combi[,"ar_s5[7]"], "ar_s5[8]"=lm1_mcmc_combi[,"ar_s5[8]"],
                                      "ar_s5[9]"=lm1_mcmc_combi[,"ar_s5[9]"], "ar_s5[10]"=lm1_mcmc_combi[,"ar_s5[10]"],"ar_s5[11]"=lm1_mcmc_combi[,"ar_s5[11]"],"ar_s5[12]"=lm1_mcmc_combi[,"ar_s5[12]"],
                                      "ar_s5[13]"=lm1_mcmc_combi[,"ar_s5[13]"],"ar_s5[14]"=lm1_mcmc_combi[,"ar_s5[14]"],"ar_s5[15]"=lm1_mcmc_combi[,"ar_s5[15]"],
                                      "ar_s5[16]"=lm1_mcmc_combi[,"ar_s5[16]"],"ar_s5[17]"=lm1_mcmc_combi[,"ar_s5[17]"],"ar_s5[18]"=lm1_mcmc_combi[,"ar_s5[18]"],
                                      "ar_s5[19]"=lm1_mcmc_combi[,"ar_s5[19]"],"ar_s5[20]"=lm1_mcmc_combi[,"ar_s5[20]"],"ar_s5[21]"=lm1_mcmc_combi[,"ar_s5[21]"],
                                      "ar_s5[22]"=lm1_mcmc_combi[,"ar_s5[22]"],"ar_s5[23]"=lm1_mcmc_combi[,"ar_s5[23]"],"ar_s5[24]"=lm1_mcmc_combi[,"ar_s5[24]"],
                                      "ar_s5[25]"=lm1_mcmc_combi[,"ar_s5[25]"],"ar_s5[26]"=lm1_mcmc_combi[,"ar_s5[26]"],"ar_s5[27]"=lm1_mcmc_combi[,"ar_s5[27]"], 
                                      "ar_s5[28]"=lm1_mcmc_combi[,"ar_s5[28]"],"ar_s5[29]"=lm1_mcmc_combi[,"ar_s5[29]"],"ar_s5[30]"=lm1_mcmc_combi[,"ar_s5[30]"],
                                      "ar_s5[31]"=lm1_mcmc_combi[,"ar_s5[31]"],"ar_s5[32]"=lm1_mcmc_combi[,"ar_s5[32]"],"ar_s5[33]"=lm1_mcmc_combi[,"ar_s5[33]"],
                                      "ar_s5[34]"=lm1_mcmc_combi[,"ar_s5[34]"],"ar_s5[35]"=lm1_mcmc_combi[,"ar_s5[35]"],"ar_s5[36]"=lm1_mcmc_combi[,"ar_s5[36]"],
                                      "ar_s5[37]"=lm1_mcmc_combi[,"ar_s5[37]"],"ar_s5[38]"=lm1_mcmc_combi[,"ar_s5[38]"],"ar_s5[39]"=lm1_mcmc_combi[,"ar_s5[39]"], 
                                      "ar_s5[40]"=lm1_mcmc_combi[,"ar_s5[40]"],
                                      "ar_s6[1]"=lm1_mcmc_combi[,"ar_s6[1]"], "ar_s6[2]"=lm1_mcmc_combi[,"ar_s6[2]"],"ar_s6[3]"=lm1_mcmc_combi[,"ar_s6[3]"],"ar_s6[4]"=lm1_mcmc_combi[,"ar_s6[4]"],
                                      "ar_s6[5]"=lm1_mcmc_combi[,"ar_s6[5]"], "ar_s6[6]"=lm1_mcmc_combi[,"ar_s6[6]"],"ar_s6[7]"=lm1_mcmc_combi[,"ar_s6[7]"], "ar_s6[8]"=lm1_mcmc_combi[,"ar_s6[8]"],
                                      "ar_s6[9]"=lm1_mcmc_combi[,"ar_s6[9]"], "ar_s6[10]"=lm1_mcmc_combi[,"ar_s6[10]"],"ar_s6[11]"=lm1_mcmc_combi[,"ar_s6[11]"],"ar_s6[12]"=lm1_mcmc_combi[,"ar_s6[12]"],
                                      "ar_s6[13]"=lm1_mcmc_combi[,"ar_s6[13]"],"ar_s6[14]"=lm1_mcmc_combi[,"ar_s6[14]"],"ar_s6[15]"=lm1_mcmc_combi[,"ar_s6[15]"],                 "ar_s6[16]"=lm1_mcmc_combi[,"ar_s6[16]"],"ar_s6[17]"=lm1_mcmc_combi[,"ar_s6[17]"],"ar_s6[18]"=lm1_mcmc_combi[,"ar_s6[18]"],
                                      "ar_s6[19]"=lm1_mcmc_combi[,"ar_s6[19]"],"ar_s6[20]"=lm1_mcmc_combi[,"ar_s6[20]"],"ar_s6[21]"=lm1_mcmc_combi[,"ar_s6[21]"],
                                      "ar_s6[22]"=lm1_mcmc_combi[,"ar_s6[22]"],"ar_s6[23]"=lm1_mcmc_combi[,"ar_s6[23]"],"ar_s6[24]"=lm1_mcmc_combi[,"ar_s6[24]"],
                                      "ar_s6[25]"=lm1_mcmc_combi[,"ar_s6[25]"],"ar_s6[26]"=lm1_mcmc_combi[,"ar_s6[26]"],"ar_s6[27]"=lm1_mcmc_combi[,"ar_s6[27]"], 
                                      "ar_s6[28]"=lm1_mcmc_combi[,"ar_s6[28]"],"ar_s6[29]"=lm1_mcmc_combi[,"ar_s6[29]"],"ar_s6[30]"=lm1_mcmc_combi[,"ar_s6[30]"],
                                      "ar_s6[31]"=lm1_mcmc_combi[,"ar_s6[31]"],"ar_s6[32]"=lm1_mcmc_combi[,"ar_s6[32]"],"ar_s6[33]"=lm1_mcmc_combi[,"ar_s6[33]"],
                                      "ar_s6[34]"=lm1_mcmc_combi[,"ar_s6[34]"],"ar_s6[35]"=lm1_mcmc_combi[,"ar_s6[35]"],"ar_s6[36]"=lm1_mcmc_combi[,"ar_s6[36]"],
                                      "ar_s6[37]"=lm1_mcmc_combi[,"ar_s6[37]"],"ar_s6[38]"=lm1_mcmc_combi[,"ar_s6[38]"],"ar_s6[39]"=lm1_mcmc_combi[,"ar_s6[39]"], 
                                      "ar_s6[40]"=lm1_mcmc_combi[,"ar_s6[40]"],                           
                                      "sigma_j_s1[1]"=lm1_mcmc_combi[,"sigma_j_s1[1]"],"sigma_j_s1[2]"=lm1_mcmc_combi[,"sigma_j_s1[2]"],"sigma_j_s1[3]"=lm1_mcmc_combi[,"sigma_j_s1[3]"],"sigma_j_s1[4]"=lm1_mcmc_combi[,"sigma_j_s1[4]"],
                                      "sigma_j_s1[5]"=lm1_mcmc_combi[,"sigma_j_s1[5]"],"sigma_j_s1[6]"=lm1_mcmc_combi[,"sigma_j_s1[6]"],"sigma_j_s1[7]"=lm1_mcmc_combi[,"sigma_j_s1[7]"],"sigma_j_s1[8]"=lm1_mcmc_combi[,"sigma_j_s1[8]"],
                                      "sigma_j_s1[9]"=lm1_mcmc_combi[,"sigma_j_s1[9]"],"sigma_j_s1[10]"=lm1_mcmc_combi[,"sigma_j_s1[10]"],"sigma_j_s1[11]"=lm1_mcmc_combi[,"sigma_j_s1[11]"],"sigma_j_s1[12]"=lm1_mcmc_combi[,"sigma_j_s1[12]"],
                                      "sigma_j_s1[13]"=lm1_mcmc_combi[,"sigma_j_s1[13]"],"sigma_j_s1[14]"=lm1_mcmc_combi[,"sigma_j_s1[14]"],"sigma_j_s1[15]"=lm1_mcmc_combi[,"sigma_j_s1[15]"],
                                      "sigma_j_s1[16]"=lm1_mcmc_combi[,"sigma_j_s1[16]"],"sigma_j_s1[17]"=lm1_mcmc_combi[,"sigma_j_s1[17]"],"sigma_j_s1[18]"=lm1_mcmc_combi[,"sigma_j_s1[18]"],
                                      "sigma_j_s1[19]"=lm1_mcmc_combi[,"sigma_j_s1[19]"],"sigma_j_s1[20]"=lm1_mcmc_combi[,"sigma_j_s1[20]"], "sigma_j_s1[21]"=lm1_mcmc_combi[,"sigma_j_s1[21]"],
                                      "sigma_j_s1[22]"=lm1_mcmc_combi[,"sigma_j_s1[22]"],"sigma_j_s1[23]"=lm1_mcmc_combi[,"sigma_j_s1[23]"],"sigma_j_s1[24]"=lm1_mcmc_combi[,"sigma_j_s1[24]"],
                                      "sigma_j_s1[25]"=lm1_mcmc_combi[,"sigma_j_s1[25]"],"sigma_j_s1[26]"=lm1_mcmc_combi[,"sigma_j_s1[26]"],"sigma_j_s1[27]"=lm1_mcmc_combi[,"sigma_j_s1[27]"],
                                      "sigma_j_s1[28]"=lm1_mcmc_combi[,"sigma_j_s1[28]"],"sigma_j_s1[29]"=lm1_mcmc_combi[,"sigma_j_s1[29]"],"sigma_j_s1[30]"=lm1_mcmc_combi[,"sigma_j_s1[30]"],
                                      "sigma_j_s1[31]"=lm1_mcmc_combi[,"sigma_j_s1[31]"],"sigma_j_s1[32]"=lm1_mcmc_combi[,"sigma_j_s1[32]"],"sigma_j_s1[33]"=lm1_mcmc_combi[,"sigma_j_s1[33]"],
                                      "sigma_j_s1[34]"=lm1_mcmc_combi[,"sigma_j_s1[34]"],"sigma_j_s1[35]"=lm1_mcmc_combi[,"sigma_j_s1[35]"],"sigma_j_s1[36]"=lm1_mcmc_combi[,"sigma_j_s1[36]"],
                                      "sigma_j_s1[37]"=lm1_mcmc_combi[,"sigma_j_s1[37]"],"sigma_j_s1[38]"=lm1_mcmc_combi[,"sigma_j_s1[38]"],"sigma_j_s1[39]"=lm1_mcmc_combi[,"sigma_j_s1[39]"],
                                      "sigma_j_s1[40]"=lm1_mcmc_combi[,"sigma_j_s1[40]"],
                                      "sigma_j_s2[1]"=lm1_mcmc_combi[,"sigma_j_s2[1]"],"sigma_j_s2[2]"=lm1_mcmc_combi[,"sigma_j_s2[2]"],"sigma_j_s2[3]"=lm1_mcmc_combi[,"sigma_j_s2[3]"],"sigma_j_s2[4]"=lm1_mcmc_combi[,"sigma_j_s2[4]"],
                                      "sigma_j_s2[5]"=lm1_mcmc_combi[,"sigma_j_s2[5]"],"sigma_j_s2[6]"=lm1_mcmc_combi[,"sigma_j_s2[6]"],"sigma_j_s2[7]"=lm1_mcmc_combi[,"sigma_j_s2[7]"],"sigma_j_s2[8]"=lm1_mcmc_combi[,"sigma_j_s2[8]"],
                                      "sigma_j_s2[9]"=lm1_mcmc_combi[,"sigma_j_s2[9]"],"sigma_j_s2[10]"=lm1_mcmc_combi[,"sigma_j_s2[10]"],"sigma_j_s2[11]"=lm1_mcmc_combi[,"sigma_j_s2[11]"],"sigma_j_s2[12]"=lm1_mcmc_combi[,"sigma_j_s2[12]"],
                                      "sigma_j_s2[13]"=lm1_mcmc_combi[,"sigma_j_s2[13]"],"sigma_j_s2[14]"=lm1_mcmc_combi[,"sigma_j_s2[14]"],"sigma_j_s2[15]"=lm1_mcmc_combi[,"sigma_j_s2[15]"],
                                      "sigma_j_s2[16]"=lm1_mcmc_combi[,"sigma_j_s2[16]"],"sigma_j_s2[17]"=lm1_mcmc_combi[,"sigma_j_s2[17]"],"sigma_j_s2[18]"=lm1_mcmc_combi[,"sigma_j_s2[18]"],
                                      "sigma_j_s2[19]"=lm1_mcmc_combi[,"sigma_j_s2[19]"],"sigma_j_s2[20]"=lm1_mcmc_combi[,"sigma_j_s2[20]"], "sigma_j_s2[21]"=lm1_mcmc_combi[,"sigma_j_s2[21]"],
                                      "sigma_j_s2[22]"=lm1_mcmc_combi[,"sigma_j_s2[22]"],"sigma_j_s2[23]"=lm1_mcmc_combi[,"sigma_j_s2[23]"],"sigma_j_s2[24]"=lm1_mcmc_combi[,"sigma_j_s2[24]"],
                                      "sigma_j_s2[25]"=lm1_mcmc_combi[,"sigma_j_s2[25]"],"sigma_j_s2[26]"=lm1_mcmc_combi[,"sigma_j_s2[26]"],"sigma_j_s2[27]"=lm1_mcmc_combi[,"sigma_j_s2[27]"],
                                      "sigma_j_s2[28]"=lm1_mcmc_combi[,"sigma_j_s2[28]"],"sigma_j_s2[29]"=lm1_mcmc_combi[,"sigma_j_s2[29]"],"sigma_j_s2[30]"=lm1_mcmc_combi[,"sigma_j_s2[30]"],
                                      "sigma_j_s2[31]"=lm1_mcmc_combi[,"sigma_j_s2[31]"],"sigma_j_s2[32]"=lm1_mcmc_combi[,"sigma_j_s2[32]"],"sigma_j_s2[33]"=lm1_mcmc_combi[,"sigma_j_s2[33]"],
                                      "sigma_j_s2[34]"=lm1_mcmc_combi[,"sigma_j_s2[34]"],"sigma_j_s2[35]"=lm1_mcmc_combi[,"sigma_j_s2[35]"],"sigma_j_s2[36]"=lm1_mcmc_combi[,"sigma_j_s2[36]"],
                                      "sigma_j_s2[37]"=lm1_mcmc_combi[,"sigma_j_s2[37]"],"sigma_j_s2[38]"=lm1_mcmc_combi[,"sigma_j_s2[38]"],"sigma_j_s2[39]"=lm1_mcmc_combi[,"sigma_j_s2[39]"],
                                      "sigma_j_s2[40]"=lm1_mcmc_combi[,"sigma_j_s2[40]"],
                                      "sigma_j_s3[1]"=lm1_mcmc_combi[,"sigma_j_s3[1]"],"sigma_j_s3[2]"=lm1_mcmc_combi[,"sigma_j_s3[2]"],"sigma_j_s3[3]"=lm1_mcmc_combi[,"sigma_j_s3[3]"],"sigma_j_s3[4]"=lm1_mcmc_combi[,"sigma_j_s3[4]"],
                                      "sigma_j_s3[5]"=lm1_mcmc_combi[,"sigma_j_s3[5]"],"sigma_j_s3[6]"=lm1_mcmc_combi[,"sigma_j_s3[6]"],"sigma_j_s3[7]"=lm1_mcmc_combi[,"sigma_j_s3[7]"],"sigma_j_s3[8]"=lm1_mcmc_combi[,"sigma_j_s3[8]"],
                                      "sigma_j_s3[9]"=lm1_mcmc_combi[,"sigma_j_s3[9]"],"sigma_j_s3[10]"=lm1_mcmc_combi[,"sigma_j_s3[10]"],"sigma_j_s3[11]"=lm1_mcmc_combi[,"sigma_j_s3[11]"],"sigma_j_s3[12]"=lm1_mcmc_combi[,"sigma_j_s3[12]"],
                                      "sigma_j_s3[13]"=lm1_mcmc_combi[,"sigma_j_s3[13]"],"sigma_j_s3[14]"=lm1_mcmc_combi[,"sigma_j_s3[14]"],"sigma_j_s3[15]"=lm1_mcmc_combi[,"sigma_j_s3[15]"],
                                      "sigma_j_s3[16]"=lm1_mcmc_combi[,"sigma_j_s3[16]"],"sigma_j_s3[17]"=lm1_mcmc_combi[,"sigma_j_s3[17]"],"sigma_j_s3[18]"=lm1_mcmc_combi[,"sigma_j_s3[18]"],
                                      "sigma_j_s3[19]"=lm1_mcmc_combi[,"sigma_j_s3[19]"],"sigma_j_s3[20]"=lm1_mcmc_combi[,"sigma_j_s3[20]"], "sigma_j_s3[21]"=lm1_mcmc_combi[,"sigma_j_s3[21]"],
                                      "sigma_j_s3[22]"=lm1_mcmc_combi[,"sigma_j_s3[22]"],"sigma_j_s3[23]"=lm1_mcmc_combi[,"sigma_j_s3[23]"],"sigma_j_s3[24]"=lm1_mcmc_combi[,"sigma_j_s3[24]"],
                                      "sigma_j_s3[25]"=lm1_mcmc_combi[,"sigma_j_s3[25]"],"sigma_j_s3[26]"=lm1_mcmc_combi[,"sigma_j_s3[26]"],"sigma_j_s3[27]"=lm1_mcmc_combi[,"sigma_j_s3[27]"],
                                      "sigma_j_s3[28]"=lm1_mcmc_combi[,"sigma_j_s3[28]"],"sigma_j_s3[29]"=lm1_mcmc_combi[,"sigma_j_s3[29]"],"sigma_j_s3[30]"=lm1_mcmc_combi[,"sigma_j_s3[30]"],
                                      "sigma_j_s3[31]"=lm1_mcmc_combi[,"sigma_j_s3[31]"],"sigma_j_s3[32]"=lm1_mcmc_combi[,"sigma_j_s3[32]"],"sigma_j_s3[33]"=lm1_mcmc_combi[,"sigma_j_s3[33]"],
                                      "sigma_j_s3[34]"=lm1_mcmc_combi[,"sigma_j_s3[34]"],"sigma_j_s3[35]"=lm1_mcmc_combi[,"sigma_j_s3[35]"],"sigma_j_s3[36]"=lm1_mcmc_combi[,"sigma_j_s3[36]"],
                                      "sigma_j_s3[37]"=lm1_mcmc_combi[,"sigma_j_s3[37]"],"sigma_j_s3[38]"=lm1_mcmc_combi[,"sigma_j_s3[38]"],"sigma_j_s3[39]"=lm1_mcmc_combi[,"sigma_j_s3[39]"],
                                      "sigma_j_s3[40]"=lm1_mcmc_combi[,"sigma_j_s3[40]"],
                                      "sigma_j_s4[1]"=lm1_mcmc_combi[,"sigma_j_s4[1]"],"sigma_j_s4[2]"=lm1_mcmc_combi[,"sigma_j_s4[2]"],"sigma_j_s4[3]"=lm1_mcmc_combi[,"sigma_j_s4[3]"],"sigma_j_s4[4]"=lm1_mcmc_combi[,"sigma_j_s4[4]"],
                                      "sigma_j_s4[5]"=lm1_mcmc_combi[,"sigma_j_s4[5]"],"sigma_j_s4[6]"=lm1_mcmc_combi[,"sigma_j_s4[6]"],"sigma_j_s4[7]"=lm1_mcmc_combi[,"sigma_j_s4[7]"],"sigma_j_s4[8]"=lm1_mcmc_combi[,"sigma_j_s4[8]"],
                                      "sigma_j_s4[9]"=lm1_mcmc_combi[,"sigma_j_s4[9]"],"sigma_j_s4[10]"=lm1_mcmc_combi[,"sigma_j_s4[10]"],"sigma_j_s4[11]"=lm1_mcmc_combi[,"sigma_j_s4[11]"],"sigma_j_s4[12]"=lm1_mcmc_combi[,"sigma_j_s4[12]"],
                                      "sigma_j_s4[13]"=lm1_mcmc_combi[,"sigma_j_s4[13]"],"sigma_j_s4[14]"=lm1_mcmc_combi[,"sigma_j_s4[14]"],"sigma_j_s4[15]"=lm1_mcmc_combi[,"sigma_j_s4[15]"],
                                      "sigma_j_s4[16]"=lm1_mcmc_combi[,"sigma_j_s4[16]"],"sigma_j_s4[17]"=lm1_mcmc_combi[,"sigma_j_s4[17]"],"sigma_j_s4[18]"=lm1_mcmc_combi[,"sigma_j_s4[18]"],
                                      "sigma_j_s4[19]"=lm1_mcmc_combi[,"sigma_j_s4[19]"],"sigma_j_s4[20]"=lm1_mcmc_combi[,"sigma_j_s4[20]"], "sigma_j_s4[21]"=lm1_mcmc_combi[,"sigma_j_s4[21]"],
                                      "sigma_j_s4[22]"=lm1_mcmc_combi[,"sigma_j_s4[22]"],"sigma_j_s4[23]"=lm1_mcmc_combi[,"sigma_j_s4[23]"],"sigma_j_s4[24]"=lm1_mcmc_combi[,"sigma_j_s4[24]"],
                                      "sigma_j_s4[25]"=lm1_mcmc_combi[,"sigma_j_s4[25]"],"sigma_j_s4[26]"=lm1_mcmc_combi[,"sigma_j_s4[26]"],"sigma_j_s4[27]"=lm1_mcmc_combi[,"sigma_j_s4[27]"],
                                      "sigma_j_s4[28]"=lm1_mcmc_combi[,"sigma_j_s4[28]"],"sigma_j_s4[29]"=lm1_mcmc_combi[,"sigma_j_s4[29]"],"sigma_j_s4[30]"=lm1_mcmc_combi[,"sigma_j_s4[30]"],
                                      "sigma_j_s4[31]"=lm1_mcmc_combi[,"sigma_j_s4[31]"],"sigma_j_s4[32]"=lm1_mcmc_combi[,"sigma_j_s4[32]"],"sigma_j_s4[33]"=lm1_mcmc_combi[,"sigma_j_s4[33]"],
                                      "sigma_j_s4[34]"=lm1_mcmc_combi[,"sigma_j_s4[34]"],"sigma_j_s4[35]"=lm1_mcmc_combi[,"sigma_j_s4[35]"],"sigma_j_s4[36]"=lm1_mcmc_combi[,"sigma_j_s4[36]"],
                                      "sigma_j_s4[37]"=lm1_mcmc_combi[,"sigma_j_s4[37]"],"sigma_j_s4[38]"=lm1_mcmc_combi[,"sigma_j_s4[38]"],"sigma_j_s4[39]"=lm1_mcmc_combi[,"sigma_j_s4[39]"],
                                      "sigma_j_s4[40]"=lm1_mcmc_combi[,"sigma_j_s4[40]"],
                                      "sigma_j_s5[1]"=lm1_mcmc_combi[,"sigma_j_s5[1]"],"sigma_j_s5[2]"=lm1_mcmc_combi[,"sigma_j_s5[2]"],"sigma_j_s5[3]"=lm1_mcmc_combi[,"sigma_j_s5[3]"],"sigma_j_s5[4]"=lm1_mcmc_combi[,"sigma_j_s5[4]"],
                                      "sigma_j_s5[5]"=lm1_mcmc_combi[,"sigma_j_s5[5]"],"sigma_j_s5[6]"=lm1_mcmc_combi[,"sigma_j_s5[6]"],"sigma_j_s5[7]"=lm1_mcmc_combi[,"sigma_j_s5[7]"],"sigma_j_s5[8]"=lm1_mcmc_combi[,"sigma_j_s5[8]"],
                                      "sigma_j_s5[9]"=lm1_mcmc_combi[,"sigma_j_s5[9]"],"sigma_j_s5[10]"=lm1_mcmc_combi[,"sigma_j_s5[10]"],"sigma_j_s5[11]"=lm1_mcmc_combi[,"sigma_j_s5[11]"],"sigma_j_s5[12]"=lm1_mcmc_combi[,"sigma_j_s5[12]"],
                                      "sigma_j_s5[13]"=lm1_mcmc_combi[,"sigma_j_s5[13]"],"sigma_j_s5[14]"=lm1_mcmc_combi[,"sigma_j_s5[14]"],"sigma_j_s5[15]"=lm1_mcmc_combi[,"sigma_j_s5[15]"],
                                      "sigma_j_s5[16]"=lm1_mcmc_combi[,"sigma_j_s5[16]"],"sigma_j_s5[17]"=lm1_mcmc_combi[,"sigma_j_s5[17]"],"sigma_j_s5[18]"=lm1_mcmc_combi[,"sigma_j_s5[18]"],
                                      "sigma_j_s5[19]"=lm1_mcmc_combi[,"sigma_j_s5[19]"],"sigma_j_s5[20]"=lm1_mcmc_combi[,"sigma_j_s5[20]"], "sigma_j_s5[21]"=lm1_mcmc_combi[,"sigma_j_s5[21]"],
                                      "sigma_j_s5[22]"=lm1_mcmc_combi[,"sigma_j_s5[22]"],"sigma_j_s5[23]"=lm1_mcmc_combi[,"sigma_j_s5[23]"],"sigma_j_s5[24]"=lm1_mcmc_combi[,"sigma_j_s5[24]"],
                                      "sigma_j_s5[25]"=lm1_mcmc_combi[,"sigma_j_s5[25]"],"sigma_j_s5[26]"=lm1_mcmc_combi[,"sigma_j_s5[26]"],"sigma_j_s5[27]"=lm1_mcmc_combi[,"sigma_j_s5[27]"],
                                      "sigma_j_s5[28]"=lm1_mcmc_combi[,"sigma_j_s5[28]"],"sigma_j_s5[29]"=lm1_mcmc_combi[,"sigma_j_s5[29]"],"sigma_j_s5[30]"=lm1_mcmc_combi[,"sigma_j_s5[30]"],
                                      "sigma_j_s5[31]"=lm1_mcmc_combi[,"sigma_j_s5[31]"],"sigma_j_s5[32]"=lm1_mcmc_combi[,"sigma_j_s5[32]"],"sigma_j_s5[33]"=lm1_mcmc_combi[,"sigma_j_s5[33]"],
                                      "sigma_j_s5[34]"=lm1_mcmc_combi[,"sigma_j_s5[34]"],"sigma_j_s5[35]"=lm1_mcmc_combi[,"sigma_j_s5[35]"],"sigma_j_s5[36]"=lm1_mcmc_combi[,"sigma_j_s5[36]"],
                                      "sigma_j_s5[37]"=lm1_mcmc_combi[,"sigma_j_s5[37]"],"sigma_j_s5[38]"=lm1_mcmc_combi[,"sigma_j_s5[38]"],"sigma_j_s5[39]"=lm1_mcmc_combi[,"sigma_j_s5[39]"],
                                      "sigma_j_s5[40]"=lm1_mcmc_combi[,"sigma_j_s5[40]"],
                                      "sigma_j_s6[1]"=lm1_mcmc_combi[,"sigma_j_s6[1]"],"sigma_j_s6[2]"=lm1_mcmc_combi[,"sigma_j_s6[2]"],"sigma_j_s6[3]"=lm1_mcmc_combi[,"sigma_j_s6[3]"],"sigma_j_s6[4]"=lm1_mcmc_combi[,"sigma_j_s6[4]"],
                                      "sigma_j_s6[5]"=lm1_mcmc_combi[,"sigma_j_s6[5]"],"sigma_j_s6[6]"=lm1_mcmc_combi[,"sigma_j_s6[6]"],"sigma_j_s6[7]"=lm1_mcmc_combi[,"sigma_j_s6[7]"],"sigma_j_s6[8]"=lm1_mcmc_combi[,"sigma_j_s6[8]"],
                                      "sigma_j_s6[9]"=lm1_mcmc_combi[,"sigma_j_s6[9]"],"sigma_j_s6[10]"=lm1_mcmc_combi[,"sigma_j_s6[10]"],"sigma_j_s6[11]"=lm1_mcmc_combi[,"sigma_j_s6[11]"],"sigma_j_s6[12]"=lm1_mcmc_combi[,"sigma_j_s6[12]"],
                                      "sigma_j_s6[13]"=lm1_mcmc_combi[,"sigma_j_s6[13]"],"sigma_j_s6[14]"=lm1_mcmc_combi[,"sigma_j_s6[14]"],"sigma_j_s6[15]"=lm1_mcmc_combi[,"sigma_j_s6[15]"],
                                      "sigma_j_s6[16]"=lm1_mcmc_combi[,"sigma_j_s6[16]"],"sigma_j_s6[17]"=lm1_mcmc_combi[,"sigma_j_s6[17]"],"sigma_j_s6[18]"=lm1_mcmc_combi[,"sigma_j_s6[18]"],
                                      "sigma_j_s6[19]"=lm1_mcmc_combi[,"sigma_j_s6[19]"],"sigma_j_s6[20]"=lm1_mcmc_combi[,"sigma_j_s6[20]"], "sigma_j_s6[21]"=lm1_mcmc_combi[,"sigma_j_s6[21]"],
                                      "sigma_j_s6[22]"=lm1_mcmc_combi[,"sigma_j_s6[22]"],"sigma_j_s6[23]"=lm1_mcmc_combi[,"sigma_j_s6[23]"],"sigma_j_s6[24]"=lm1_mcmc_combi[,"sigma_j_s6[24]"],
                                      "sigma_j_s6[25]"=lm1_mcmc_combi[,"sigma_j_s6[25]"],"sigma_j_s6[26]"=lm1_mcmc_combi[,"sigma_j_s6[26]"],"sigma_j_s6[27]"=lm1_mcmc_combi[,"sigma_j_s6[27]"],
                                      "sigma_j_s6[28]"=lm1_mcmc_combi[,"sigma_j_s6[28]"],"sigma_j_s6[29]"=lm1_mcmc_combi[,"sigma_j_s6[29]"],"sigma_j_s6[30]"=lm1_mcmc_combi[,"sigma_j_s6[30]"],
                                      "sigma_j_s6[31]"=lm1_mcmc_combi[,"sigma_j_s6[31]"],"sigma_j_s6[32]"=lm1_mcmc_combi[,"sigma_j_s6[32]"],"sigma_j_s6[33]"=lm1_mcmc_combi[,"sigma_j_s6[33]"],
                                      "sigma_j_s6[34]"=lm1_mcmc_combi[,"sigma_j_s6[34]"],"sigma_j_s6[35]"=lm1_mcmc_combi[,"sigma_j_s6[35]"],"sigma_j_s6[36]"=lm1_mcmc_combi[,"sigma_j_s6[36]"],
                                      "sigma_j_s6[37]"=lm1_mcmc_combi[,"sigma_j_s6[37]"],"sigma_j_s6[38]"=lm1_mcmc_combi[,"sigma_j_s6[38]"],"sigma_j_s6[39]"=lm1_mcmc_combi[,"sigma_j_s6[39]"],
                                      "sigma_j_s6[40]"=lm1_mcmc_combi[,"sigma_j_s6[40]"]), 1)) 
lm1_mcmc_v <- do.call(cbind, rep(list("ar_v1[1]"=lm1_mcmc_combi[,"ar_v1[1]"], "ar_v1[2]"=lm1_mcmc_combi[,"ar_v1[2]"],"ar_v1[3]"=lm1_mcmc_combi[,"ar_v1[3]"],"ar_v1[4]"=lm1_mcmc_combi[,"ar_v1[4]"],
                                      "ar_v1[5]"=lm1_mcmc_combi[,"ar_v1[5]"], "ar_v1[6]"=lm1_mcmc_combi[,"ar_v1[6]"],"ar_v1[7]"=lm1_mcmc_combi[,"ar_v1[7]"], "ar_v1[8]"=lm1_mcmc_combi[,"ar_v1[8]"],
                                      "ar_v1[9]"=lm1_mcmc_combi[,"ar_v1[9]"], "ar_v1[10]"=lm1_mcmc_combi[,"ar_v1[10]"],"ar_v1[11]"=lm1_mcmc_combi[,"ar_v1[11]"],"ar_v1[12]"=lm1_mcmc_combi[,"ar_v1[12]"],
                                      "ar_v1[13]"=lm1_mcmc_combi[,"ar_v1[13]"],"ar_v1[14]"=lm1_mcmc_combi[,"ar_v1[14]"],"ar_v1[15]"=lm1_mcmc_combi[,"ar_v1[15]"],
                                      "ar_v1[16]"=lm1_mcmc_combi[,"ar_v1[16]"],"ar_v1[17]"=lm1_mcmc_combi[,"ar_v1[17]"],"ar_v1[18]"=lm1_mcmc_combi[,"ar_v1[18]"],
                                      "ar_v1[19]"=lm1_mcmc_combi[,"ar_v1[19]"],"ar_v1[20]"=lm1_mcmc_combi[,"ar_v1[20]"],"ar_v1[21]"=lm1_mcmc_combi[,"ar_v1[21]"],
                                      "ar_v1[22]"=lm1_mcmc_combi[,"ar_v1[22]"],"ar_v1[23]"=lm1_mcmc_combi[,"ar_v1[23]"],"ar_v1[24]"=lm1_mcmc_combi[,"ar_v1[24]"],
                                      "ar_v1[25]"=lm1_mcmc_combi[,"ar_v1[25]"],"ar_v1[26]"=lm1_mcmc_combi[,"ar_v1[26]"],"ar_v1[27]"=lm1_mcmc_combi[,"ar_v1[27]"], 
                                      "ar_v1[28]"=lm1_mcmc_combi[,"ar_v1[28]"],"ar_v1[29]"=lm1_mcmc_combi[,"ar_v1[29]"],"ar_v1[30]"=lm1_mcmc_combi[,"ar_v1[30]"],
                                      "ar_v1[31]"=lm1_mcmc_combi[,"ar_v1[31]"],"ar_v1[32]"=lm1_mcmc_combi[,"ar_v1[32]"],"ar_v1[33]"=lm1_mcmc_combi[,"ar_v1[33]"],
                                      "ar_v1[34]"=lm1_mcmc_combi[,"ar_v1[34]"],"ar_v1[35]"=lm1_mcmc_combi[,"ar_v1[35]"],"ar_v1[36]"=lm1_mcmc_combi[,"ar_v1[36]"],
                                      "ar_v1[37]"=lm1_mcmc_combi[,"ar_v1[37]"],"ar_v1[38]"=lm1_mcmc_combi[,"ar_v1[38]"],"ar_v1[39]"=lm1_mcmc_combi[,"ar_v1[39]"], 
                                      "ar_v1[40]"=lm1_mcmc_combi[,"ar_v1[40]"],
                                      "ar_v2[1]"=lm1_mcmc_combi[,"ar_v2[1]"], "ar_v2[2]"=lm1_mcmc_combi[,"ar_v2[2]"],"ar_v2[3]"=lm1_mcmc_combi[,"ar_v2[3]"],"ar_v2[4]"=lm1_mcmc_combi[,"ar_v2[4]"],
                                      "ar_v2[5]"=lm1_mcmc_combi[,"ar_v2[5]"], "ar_v2[6]"=lm1_mcmc_combi[,"ar_v2[6]"],"ar_v2[7]"=lm1_mcmc_combi[,"ar_v2[7]"], "ar_v2[8]"=lm1_mcmc_combi[,"ar_v2[8]"],
                                      "ar_v2[9]"=lm1_mcmc_combi[,"ar_v2[9]"], "ar_v2[10]"=lm1_mcmc_combi[,"ar_v2[10]"],"ar_v2[11]"=lm1_mcmc_combi[,"ar_v2[11]"],"ar_v2[12]"=lm1_mcmc_combi[,"ar_v2[12]"],
                                      "ar_v2[13]"=lm1_mcmc_combi[,"ar_v2[13]"],"ar_v2[14]"=lm1_mcmc_combi[,"ar_v2[14]"],"ar_v2[15]"=lm1_mcmc_combi[,"ar_v2[15]"],
                                      "ar_v2[16]"=lm1_mcmc_combi[,"ar_v2[16]"],"ar_v2[17]"=lm1_mcmc_combi[,"ar_v2[17]"],"ar_v2[18]"=lm1_mcmc_combi[,"ar_v2[18]"],
                                      "ar_v2[19]"=lm1_mcmc_combi[,"ar_v2[19]"],"ar_v2[20]"=lm1_mcmc_combi[,"ar_v2[20]"],"ar_v2[21]"=lm1_mcmc_combi[,"ar_v2[21]"],
                                      "ar_v2[22]"=lm1_mcmc_combi[,"ar_v2[22]"],"ar_v2[23]"=lm1_mcmc_combi[,"ar_v2[23]"],"ar_v2[24]"=lm1_mcmc_combi[,"ar_v2[24]"],
                                      "ar_v2[25]"=lm1_mcmc_combi[,"ar_v2[25]"],"ar_v2[26]"=lm1_mcmc_combi[,"ar_v2[26]"],"ar_v2[27]"=lm1_mcmc_combi[,"ar_v2[27]"], 
                                      "ar_v2[28]"=lm1_mcmc_combi[,"ar_v2[28]"],"ar_v2[29]"=lm1_mcmc_combi[,"ar_v2[29]"],"ar_v2[30]"=lm1_mcmc_combi[,"ar_v2[30]"],
                                      "ar_v2[31]"=lm1_mcmc_combi[,"ar_v2[31]"],"ar_v2[32]"=lm1_mcmc_combi[,"ar_v2[32]"],"ar_v2[33]"=lm1_mcmc_combi[,"ar_v2[33]"],
                                      "ar_v2[34]"=lm1_mcmc_combi[,"ar_v2[34]"],"ar_v2[35]"=lm1_mcmc_combi[,"ar_v2[35]"],"ar_v2[36]"=lm1_mcmc_combi[,"ar_v2[36]"],
                                      "ar_v2[37]"=lm1_mcmc_combi[,"ar_v2[37]"],"ar_v2[38]"=lm1_mcmc_combi[,"ar_v2[38]"],"ar_v2[39]"=lm1_mcmc_combi[,"ar_v2[39]"], 
                                      "ar_v2[40]"=lm1_mcmc_combi[,"ar_v2[40]"],
                                      "ar_v3[1]"=lm1_mcmc_combi[,"ar_v3[1]"], "ar_v3[2]"=lm1_mcmc_combi[,"ar_v3[2]"],"ar_v3[3]"=lm1_mcmc_combi[,"ar_v3[3]"],"ar_v3[4]"=lm1_mcmc_combi[,"ar_v3[4]"],
                                      "ar_v3[5]"=lm1_mcmc_combi[,"ar_v3[5]"], "ar_v3[6]"=lm1_mcmc_combi[,"ar_v3[6]"],"ar_v3[7]"=lm1_mcmc_combi[,"ar_v3[7]"], "ar_v3[8]"=lm1_mcmc_combi[,"ar_v3[8]"],
                                      "ar_v3[9]"=lm1_mcmc_combi[,"ar_v3[9]"], "ar_v3[10]"=lm1_mcmc_combi[,"ar_v3[10]"],"ar_v3[11]"=lm1_mcmc_combi[,"ar_v3[11]"],"ar_v3[12]"=lm1_mcmc_combi[,"ar_v3[12]"],
                                      "ar_v3[13]"=lm1_mcmc_combi[,"ar_v3[13]"],"ar_v3[14]"=lm1_mcmc_combi[,"ar_v3[14]"],"ar_v3[15]"=lm1_mcmc_combi[,"ar_v3[15]"],
                                      "ar_v3[16]"=lm1_mcmc_combi[,"ar_v3[16]"],"ar_v3[17]"=lm1_mcmc_combi[,"ar_v3[17]"],"ar_v3[18]"=lm1_mcmc_combi[,"ar_v3[18]"],
                                      "ar_v3[19]"=lm1_mcmc_combi[,"ar_v3[19]"],"ar_v3[20]"=lm1_mcmc_combi[,"ar_v3[20]"],"ar_v3[21]"=lm1_mcmc_combi[,"ar_v3[21]"],
                                      "ar_v3[22]"=lm1_mcmc_combi[,"ar_v3[22]"],"ar_v3[23]"=lm1_mcmc_combi[,"ar_v3[23]"],"ar_v3[24]"=lm1_mcmc_combi[,"ar_v3[24]"],
                                      "ar_v3[25]"=lm1_mcmc_combi[,"ar_v3[25]"],"ar_v3[26]"=lm1_mcmc_combi[,"ar_v3[26]"],"ar_v3[27]"=lm1_mcmc_combi[,"ar_v3[27]"], 
                                      "ar_v3[28]"=lm1_mcmc_combi[,"ar_v3[28]"],"ar_v3[29]"=lm1_mcmc_combi[,"ar_v3[29]"],"ar_v3[30]"=lm1_mcmc_combi[,"ar_v3[30]"],
                                      "ar_v3[31]"=lm1_mcmc_combi[,"ar_v3[31]"],"ar_v3[32]"=lm1_mcmc_combi[,"ar_v3[32]"],"ar_v3[33]"=lm1_mcmc_combi[,"ar_v3[33]"],
                                      "ar_v3[34]"=lm1_mcmc_combi[,"ar_v3[34]"],"ar_v3[35]"=lm1_mcmc_combi[,"ar_v3[35]"],"ar_v3[36]"=lm1_mcmc_combi[,"ar_v3[36]"],
                                      "ar_v3[37]"=lm1_mcmc_combi[,"ar_v3[37]"],"ar_v3[38]"=lm1_mcmc_combi[,"ar_v3[38]"],"ar_v3[39]"=lm1_mcmc_combi[,"ar_v3[39]"], 
                                      "ar_v3[40]"=lm1_mcmc_combi[,"ar_v3[40]"],
                                      "ar_v4[1]"=lm1_mcmc_combi[,"ar_v4[1]"], "ar_v4[2]"=lm1_mcmc_combi[,"ar_v4[2]"],"ar_v4[3]"=lm1_mcmc_combi[,"ar_v4[3]"],"ar_v4[4]"=lm1_mcmc_combi[,"ar_v4[4]"],
                                      "ar_v4[5]"=lm1_mcmc_combi[,"ar_v4[5]"], "ar_v4[6]"=lm1_mcmc_combi[,"ar_v4[6]"],"ar_v4[7]"=lm1_mcmc_combi[,"ar_v4[7]"], "ar_v4[8]"=lm1_mcmc_combi[,"ar_v4[8]"],
                                      "ar_v4[9]"=lm1_mcmc_combi[,"ar_v4[9]"], "ar_v4[10]"=lm1_mcmc_combi[,"ar_v4[10]"],"ar_v4[11]"=lm1_mcmc_combi[,"ar_v4[11]"],"ar_v4[12]"=lm1_mcmc_combi[,"ar_v4[12]"],
                                      "ar_v4[13]"=lm1_mcmc_combi[,"ar_v4[13]"],"ar_v4[14]"=lm1_mcmc_combi[,"ar_v4[14]"],"ar_v4[15]"=lm1_mcmc_combi[,"ar_v4[15]"],
                                      "ar_v4[16]"=lm1_mcmc_combi[,"ar_v4[16]"],"ar_v4[17]"=lm1_mcmc_combi[,"ar_v4[17]"],"ar_v4[18]"=lm1_mcmc_combi[,"ar_v4[18]"],
                                      "ar_v4[19]"=lm1_mcmc_combi[,"ar_v4[19]"],"ar_v4[20]"=lm1_mcmc_combi[,"ar_v4[20]"],"ar_v4[21]"=lm1_mcmc_combi[,"ar_v4[21]"],
                                      "ar_v4[22]"=lm1_mcmc_combi[,"ar_v4[22]"],"ar_v4[23]"=lm1_mcmc_combi[,"ar_v4[23]"],"ar_v4[24]"=lm1_mcmc_combi[,"ar_v4[24]"],
                                      "ar_v4[25]"=lm1_mcmc_combi[,"ar_v4[25]"],"ar_v4[26]"=lm1_mcmc_combi[,"ar_v4[26]"],"ar_v4[27]"=lm1_mcmc_combi[,"ar_v4[27]"], 
                                      "ar_v4[28]"=lm1_mcmc_combi[,"ar_v4[28]"],"ar_v4[29]"=lm1_mcmc_combi[,"ar_v4[29]"],"ar_v4[30]"=lm1_mcmc_combi[,"ar_v4[30]"],
                                      "ar_v4[31]"=lm1_mcmc_combi[,"ar_v4[31]"],"ar_v4[32]"=lm1_mcmc_combi[,"ar_v4[32]"],"ar_v4[33]"=lm1_mcmc_combi[,"ar_v4[33]"],
                                      "ar_v4[34]"=lm1_mcmc_combi[,"ar_v4[34]"],"ar_v4[35]"=lm1_mcmc_combi[,"ar_v4[35]"],"ar_v4[36]"=lm1_mcmc_combi[,"ar_v4[36]"],
                                      "ar_v4[37]"=lm1_mcmc_combi[,"ar_v4[37]"],"ar_v4[38]"=lm1_mcmc_combi[,"ar_v4[38]"],"ar_v4[39]"=lm1_mcmc_combi[,"ar_v4[39]"], 
                                      "ar_v4[40]"=lm1_mcmc_combi[,"ar_v4[40]"],
                                      "ar_v5[1]"=lm1_mcmc_combi[,"ar_v5[1]"], "ar_v5[2]"=lm1_mcmc_combi[,"ar_v5[2]"],"ar_v5[3]"=lm1_mcmc_combi[,"ar_v5[3]"],"ar_v5[4]"=lm1_mcmc_combi[,"ar_v5[4]"],
                                      "ar_v5[5]"=lm1_mcmc_combi[,"ar_v5[5]"], "ar_v5[6]"=lm1_mcmc_combi[,"ar_v5[6]"],"ar_v5[7]"=lm1_mcmc_combi[,"ar_v5[7]"], "ar_v5[8]"=lm1_mcmc_combi[,"ar_v5[8]"],
                                      "ar_v5[9]"=lm1_mcmc_combi[,"ar_v5[9]"], "ar_v5[10]"=lm1_mcmc_combi[,"ar_v5[10]"],"ar_v5[11]"=lm1_mcmc_combi[,"ar_v5[11]"],"ar_v5[12]"=lm1_mcmc_combi[,"ar_v5[12]"],
                                      "ar_v5[13]"=lm1_mcmc_combi[,"ar_v5[13]"],"ar_v5[14]"=lm1_mcmc_combi[,"ar_v5[14]"],"ar_v5[15]"=lm1_mcmc_combi[,"ar_v5[15]"],
                                      "ar_v5[16]"=lm1_mcmc_combi[,"ar_v5[16]"],"ar_v5[17]"=lm1_mcmc_combi[,"ar_v5[17]"],"ar_v5[18]"=lm1_mcmc_combi[,"ar_v5[18]"],
                                      "ar_v5[19]"=lm1_mcmc_combi[,"ar_v5[19]"],"ar_v5[20]"=lm1_mcmc_combi[,"ar_v5[20]"],"ar_v5[21]"=lm1_mcmc_combi[,"ar_v5[21]"],
                                      "ar_v5[22]"=lm1_mcmc_combi[,"ar_v5[22]"],"ar_v5[23]"=lm1_mcmc_combi[,"ar_v5[23]"],"ar_v5[24]"=lm1_mcmc_combi[,"ar_v5[24]"],
                                      "ar_v5[25]"=lm1_mcmc_combi[,"ar_v5[25]"],"ar_v5[26]"=lm1_mcmc_combi[,"ar_v5[26]"],"ar_v5[27]"=lm1_mcmc_combi[,"ar_v5[27]"], 
                                      "ar_v5[28]"=lm1_mcmc_combi[,"ar_v5[28]"],"ar_v5[29]"=lm1_mcmc_combi[,"ar_v5[29]"],"ar_v5[30]"=lm1_mcmc_combi[,"ar_v5[30]"],
                                      "ar_v5[31]"=lm1_mcmc_combi[,"ar_v5[31]"],"ar_v5[32]"=lm1_mcmc_combi[,"ar_v5[32]"],"ar_v5[33]"=lm1_mcmc_combi[,"ar_v5[33]"],
                                      "ar_v5[34]"=lm1_mcmc_combi[,"ar_v5[34]"],"ar_v5[35]"=lm1_mcmc_combi[,"ar_v5[35]"],"ar_v5[36]"=lm1_mcmc_combi[,"ar_v5[36]"],
                                      "ar_v5[37]"=lm1_mcmc_combi[,"ar_v5[37]"],"ar_v5[38]"=lm1_mcmc_combi[,"ar_v5[38]"],"ar_v5[39]"=lm1_mcmc_combi[,"ar_v5[39]"], 
                                      "ar_v5[40]"=lm1_mcmc_combi[,"ar_v5[40]"],
                                      "sigma_j_v1[1]"=lm1_mcmc_combi[,"sigma_j_v1[1]"],"sigma_j_v1[2]"=lm1_mcmc_combi[,"sigma_j_v1[2]"],"sigma_j_v1[3]"=lm1_mcmc_combi[,"sigma_j_v1[3]"],"sigma_j_v1[4]"=lm1_mcmc_combi[,"sigma_j_v1[4]"],
                                      "sigma_j_v1[5]"=lm1_mcmc_combi[,"sigma_j_v1[5]"],"sigma_j_v1[6]"=lm1_mcmc_combi[,"sigma_j_v1[6]"],"sigma_j_v1[7]"=lm1_mcmc_combi[,"sigma_j_v1[7]"],"sigma_j_v1[8]"=lm1_mcmc_combi[,"sigma_j_v1[8]"],
                                      "sigma_j_v1[9]"=lm1_mcmc_combi[,"sigma_j_v1[9]"],"sigma_j_v1[10]"=lm1_mcmc_combi[,"sigma_j_v1[10]"],"sigma_j_v1[11]"=lm1_mcmc_combi[,"sigma_j_v1[11]"],"sigma_j_v1[12]"=lm1_mcmc_combi[,"sigma_j_v1[12]"],
                                      "sigma_j_v1[13]"=lm1_mcmc_combi[,"sigma_j_v1[13]"],"sigma_j_v1[14]"=lm1_mcmc_combi[,"sigma_j_v1[14]"],"sigma_j_v1[15]"=lm1_mcmc_combi[,"sigma_j_v1[15]"],
                                      "sigma_j_v1[16]"=lm1_mcmc_combi[,"sigma_j_v1[16]"],"sigma_j_v1[17]"=lm1_mcmc_combi[,"sigma_j_v1[17]"],"sigma_j_v1[18]"=lm1_mcmc_combi[,"sigma_j_v1[18]"],
                                      "sigma_j_v1[19]"=lm1_mcmc_combi[,"sigma_j_v1[19]"],"sigma_j_v1[20]"=lm1_mcmc_combi[,"sigma_j_v1[20]"],"sigma_j_v1[21]"=lm1_mcmc_combi[,"sigma_j_v1[21]"],
                                      "sigma_j_v1[22]"=lm1_mcmc_combi[,"sigma_j_v1[22]"],"sigma_j_v1[23]"=lm1_mcmc_combi[,"sigma_j_v1[23]"],"sigma_j_v1[24]"=lm1_mcmc_combi[,"sigma_j_v1[24]"],
                                      "sigma_j_v1[25]"=lm1_mcmc_combi[,"sigma_j_v1[25]"],"sigma_j_v1[26]"=lm1_mcmc_combi[,"sigma_j_v1[26]"],"sigma_j_v1[27]"=lm1_mcmc_combi[,"sigma_j_v1[27]"],
                                      "sigma_j_v1[28]"=lm1_mcmc_combi[,"sigma_j_v1[28]"],"sigma_j_v1[29]"=lm1_mcmc_combi[,"sigma_j_v1[29]"],"sigma_j_v1[30]"=lm1_mcmc_combi[,"sigma_j_v1[30]"],
                                      "sigma_j_v1[31]"=lm1_mcmc_combi[,"sigma_j_v1[31]"],"sigma_j_v1[32]"=lm1_mcmc_combi[,"sigma_j_v1[32]"],"sigma_j_v1[33]"=lm1_mcmc_combi[,"sigma_j_v1[33]"],
                                      "sigma_j_v1[34]"=lm1_mcmc_combi[,"sigma_j_v1[34]"],"sigma_j_v1[35]"=lm1_mcmc_combi[,"sigma_j_v1[35]"],"sigma_j_v1[36]"=lm1_mcmc_combi[,"sigma_j_v1[36]"],
                                      "sigma_j_v1[37]"=lm1_mcmc_combi[,"sigma_j_v1[37]"],"sigma_j_v1[38]"=lm1_mcmc_combi[,"sigma_j_v1[38]"],"sigma_j_v1[39]"=lm1_mcmc_combi[,"sigma_j_v1[39]"],
                                      "sigma_j_v1[40]"=lm1_mcmc_combi[,"sigma_j_v1[40]"],
                                      "sigma_j_v2[1]"=lm1_mcmc_combi[,"sigma_j_v2[1]"],"sigma_j_v2[2]"=lm1_mcmc_combi[,"sigma_j_v2[2]"],"sigma_j_v2[3]"=lm1_mcmc_combi[,"sigma_j_v2[3]"],"sigma_j_v2[4]"=lm1_mcmc_combi[,"sigma_j_v2[4]"],
                                      "sigma_j_v2[5]"=lm1_mcmc_combi[,"sigma_j_v2[5]"],"sigma_j_v2[6]"=lm1_mcmc_combi[,"sigma_j_v2[6]"],"sigma_j_v2[7]"=lm1_mcmc_combi[,"sigma_j_v2[7]"],"sigma_j_v2[8]"=lm1_mcmc_combi[,"sigma_j_v2[8]"],
                                      "sigma_j_v2[9]"=lm1_mcmc_combi[,"sigma_j_v2[9]"],"sigma_j_v2[10]"=lm1_mcmc_combi[,"sigma_j_v2[10]"],"sigma_j_v2[11]"=lm1_mcmc_combi[,"sigma_j_v2[11]"],"sigma_j_v2[12]"=lm1_mcmc_combi[,"sigma_j_v2[12]"],
                                      "sigma_j_v2[13]"=lm1_mcmc_combi[,"sigma_j_v2[13]"],"sigma_j_v2[14]"=lm1_mcmc_combi[,"sigma_j_v2[14]"],"sigma_j_v2[15]"=lm1_mcmc_combi[,"sigma_j_v2[15]"],
                                      "sigma_j_v2[16]"=lm1_mcmc_combi[,"sigma_j_v2[16]"],"sigma_j_v2[17]"=lm1_mcmc_combi[,"sigma_j_v2[17]"],"sigma_j_v2[18]"=lm1_mcmc_combi[,"sigma_j_v2[18]"],
                                      "sigma_j_v2[19]"=lm1_mcmc_combi[,"sigma_j_v2[19]"],"sigma_j_v2[20]"=lm1_mcmc_combi[,"sigma_j_v2[20]"],"sigma_j_v2[21]"=lm1_mcmc_combi[,"sigma_j_v2[21]"],
                                      "sigma_j_v2[22]"=lm1_mcmc_combi[,"sigma_j_v2[22]"],"sigma_j_v2[23]"=lm1_mcmc_combi[,"sigma_j_v2[23]"],"sigma_j_v2[24]"=lm1_mcmc_combi[,"sigma_j_v2[24]"],
                                      "sigma_j_v2[25]"=lm1_mcmc_combi[,"sigma_j_v2[25]"],"sigma_j_v2[26]"=lm1_mcmc_combi[,"sigma_j_v2[26]"],"sigma_j_v2[27]"=lm1_mcmc_combi[,"sigma_j_v2[27]"],
                                      "sigma_j_v2[28]"=lm1_mcmc_combi[,"sigma_j_v2[28]"],"sigma_j_v2[29]"=lm1_mcmc_combi[,"sigma_j_v2[29]"],"sigma_j_v2[30]"=lm1_mcmc_combi[,"sigma_j_v2[30]"],
                                      "sigma_j_v2[31]"=lm1_mcmc_combi[,"sigma_j_v2[31]"],"sigma_j_v2[32]"=lm1_mcmc_combi[,"sigma_j_v2[32]"],"sigma_j_v2[33]"=lm1_mcmc_combi[,"sigma_j_v2[33]"],
                                      "sigma_j_v2[34]"=lm1_mcmc_combi[,"sigma_j_v2[34]"],"sigma_j_v2[35]"=lm1_mcmc_combi[,"sigma_j_v2[35]"],"sigma_j_v2[36]"=lm1_mcmc_combi[,"sigma_j_v2[36]"],
                                      "sigma_j_v2[37]"=lm1_mcmc_combi[,"sigma_j_v2[37]"],"sigma_j_v2[38]"=lm1_mcmc_combi[,"sigma_j_v2[38]"],"sigma_j_v2[39]"=lm1_mcmc_combi[,"sigma_j_v2[39]"],
                                      "sigma_j_v2[40]"=lm1_mcmc_combi[,"sigma_j_v2[40]"],
                                      "sigma_j_v3[1]"=lm1_mcmc_combi[,"sigma_j_v3[1]"],"sigma_j_v3[2]"=lm1_mcmc_combi[,"sigma_j_v3[2]"],"sigma_j_v3[3]"=lm1_mcmc_combi[,"sigma_j_v3[3]"],"sigma_j_v3[4]"=lm1_mcmc_combi[,"sigma_j_v3[4]"],
                                      "sigma_j_v3[5]"=lm1_mcmc_combi[,"sigma_j_v3[5]"],"sigma_j_v3[6]"=lm1_mcmc_combi[,"sigma_j_v3[6]"],"sigma_j_v3[7]"=lm1_mcmc_combi[,"sigma_j_v3[7]"],"sigma_j_v3[8]"=lm1_mcmc_combi[,"sigma_j_v3[8]"],
                                      "sigma_j_v3[9]"=lm1_mcmc_combi[,"sigma_j_v3[9]"],"sigma_j_v3[10]"=lm1_mcmc_combi[,"sigma_j_v3[10]"],"sigma_j_v3[11]"=lm1_mcmc_combi[,"sigma_j_v3[11]"],"sigma_j_v3[12]"=lm1_mcmc_combi[,"sigma_j_v3[12]"],
                                      "sigma_j_v3[13]"=lm1_mcmc_combi[,"sigma_j_v3[13]"],"sigma_j_v3[14]"=lm1_mcmc_combi[,"sigma_j_v3[14]"],"sigma_j_v3[15]"=lm1_mcmc_combi[,"sigma_j_v3[15]"],
                                      "sigma_j_v3[16]"=lm1_mcmc_combi[,"sigma_j_v3[16]"],"sigma_j_v3[17]"=lm1_mcmc_combi[,"sigma_j_v3[17]"],"sigma_j_v3[18]"=lm1_mcmc_combi[,"sigma_j_v3[18]"],
                                      "sigma_j_v3[19]"=lm1_mcmc_combi[,"sigma_j_v3[19]"],"sigma_j_v3[20]"=lm1_mcmc_combi[,"sigma_j_v3[20]"],"sigma_j_v3[21]"=lm1_mcmc_combi[,"sigma_j_v3[21]"],
                                      "sigma_j_v3[22]"=lm1_mcmc_combi[,"sigma_j_v3[22]"],"sigma_j_v3[23]"=lm1_mcmc_combi[,"sigma_j_v3[23]"],"sigma_j_v3[24]"=lm1_mcmc_combi[,"sigma_j_v3[24]"],
                                      "sigma_j_v3[25]"=lm1_mcmc_combi[,"sigma_j_v3[25]"],"sigma_j_v3[26]"=lm1_mcmc_combi[,"sigma_j_v3[26]"],"sigma_j_v3[27]"=lm1_mcmc_combi[,"sigma_j_v3[27]"],
                                      "sigma_j_v3[28]"=lm1_mcmc_combi[,"sigma_j_v3[28]"],"sigma_j_v3[29]"=lm1_mcmc_combi[,"sigma_j_v3[29]"],"sigma_j_v3[30]"=lm1_mcmc_combi[,"sigma_j_v3[30]"],
                                      "sigma_j_v3[31]"=lm1_mcmc_combi[,"sigma_j_v3[31]"],"sigma_j_v3[32]"=lm1_mcmc_combi[,"sigma_j_v3[32]"],"sigma_j_v3[33]"=lm1_mcmc_combi[,"sigma_j_v3[33]"],
                                      "sigma_j_v3[34]"=lm1_mcmc_combi[,"sigma_j_v3[34]"],"sigma_j_v3[35]"=lm1_mcmc_combi[,"sigma_j_v3[35]"],"sigma_j_v3[36]"=lm1_mcmc_combi[,"sigma_j_v3[36]"],
                                      "sigma_j_v3[37]"=lm1_mcmc_combi[,"sigma_j_v3[37]"],"sigma_j_v3[38]"=lm1_mcmc_combi[,"sigma_j_v3[38]"],"sigma_j_v3[39]"=lm1_mcmc_combi[,"sigma_j_v3[39]"],
                                      "sigma_j_v3[40]"=lm1_mcmc_combi[,"sigma_j_v3[40]"],
                                      "sigma_j_v4[1]"=lm1_mcmc_combi[,"sigma_j_v4[1]"],"sigma_j_v4[2]"=lm1_mcmc_combi[,"sigma_j_v4[2]"],"sigma_j_v4[3]"=lm1_mcmc_combi[,"sigma_j_v4[3]"],"sigma_j_v4[4]"=lm1_mcmc_combi[,"sigma_j_v4[4]"],
                                      "sigma_j_v4[5]"=lm1_mcmc_combi[,"sigma_j_v4[5]"],"sigma_j_v4[6]"=lm1_mcmc_combi[,"sigma_j_v4[6]"],"sigma_j_v4[7]"=lm1_mcmc_combi[,"sigma_j_v4[7]"],"sigma_j_v4[8]"=lm1_mcmc_combi[,"sigma_j_v4[8]"],
                                      "sigma_j_v4[9]"=lm1_mcmc_combi[,"sigma_j_v4[9]"],"sigma_j_v4[10]"=lm1_mcmc_combi[,"sigma_j_v4[10]"],"sigma_j_v4[11]"=lm1_mcmc_combi[,"sigma_j_v4[11]"],"sigma_j_v4[12]"=lm1_mcmc_combi[,"sigma_j_v4[12]"],
                                      "sigma_j_v4[13]"=lm1_mcmc_combi[,"sigma_j_v4[13]"],"sigma_j_v4[14]"=lm1_mcmc_combi[,"sigma_j_v4[14]"],"sigma_j_v4[15]"=lm1_mcmc_combi[,"sigma_j_v4[15]"],
                                      "sigma_j_v4[16]"=lm1_mcmc_combi[,"sigma_j_v4[16]"],"sigma_j_v4[17]"=lm1_mcmc_combi[,"sigma_j_v4[17]"],"sigma_j_v4[18]"=lm1_mcmc_combi[,"sigma_j_v4[18]"],
                                      "sigma_j_v4[19]"=lm1_mcmc_combi[,"sigma_j_v4[19]"],"sigma_j_v4[20]"=lm1_mcmc_combi[,"sigma_j_v4[20]"],"sigma_j_v4[21]"=lm1_mcmc_combi[,"sigma_j_v4[21]"],
                                      "sigma_j_v4[22]"=lm1_mcmc_combi[,"sigma_j_v4[22]"],"sigma_j_v4[23]"=lm1_mcmc_combi[,"sigma_j_v4[23]"],"sigma_j_v4[24]"=lm1_mcmc_combi[,"sigma_j_v4[24]"],
                                      "sigma_j_v4[25]"=lm1_mcmc_combi[,"sigma_j_v4[25]"],"sigma_j_v4[26]"=lm1_mcmc_combi[,"sigma_j_v4[26]"],"sigma_j_v4[27]"=lm1_mcmc_combi[,"sigma_j_v4[27]"],
                                      "sigma_j_v4[28]"=lm1_mcmc_combi[,"sigma_j_v4[28]"],"sigma_j_v4[29]"=lm1_mcmc_combi[,"sigma_j_v4[29]"],"sigma_j_v4[30]"=lm1_mcmc_combi[,"sigma_j_v4[30]"],
                                      "sigma_j_v4[31]"=lm1_mcmc_combi[,"sigma_j_v4[31]"],"sigma_j_v4[32]"=lm1_mcmc_combi[,"sigma_j_v4[32]"],"sigma_j_v4[33]"=lm1_mcmc_combi[,"sigma_j_v4[33]"],
                                      "sigma_j_v4[34]"=lm1_mcmc_combi[,"sigma_j_v4[34]"],"sigma_j_v4[35]"=lm1_mcmc_combi[,"sigma_j_v4[35]"],"sigma_j_v4[36]"=lm1_mcmc_combi[,"sigma_j_v4[36]"],
                                      "sigma_j_v4[37]"=lm1_mcmc_combi[,"sigma_j_v4[37]"],"sigma_j_v4[38]"=lm1_mcmc_combi[,"sigma_j_v4[38]"],"sigma_j_v4[39]"=lm1_mcmc_combi[,"sigma_j_v4[39]"],
                                      "sigma_j_v4[40]"=lm1_mcmc_combi[,"sigma_j_v4[40]"],
                                      "sigma_j_v5[1]"=lm1_mcmc_combi[,"sigma_j_v5[1]"],"sigma_j_v5[2]"=lm1_mcmc_combi[,"sigma_j_v5[2]"],"sigma_j_v5[3]"=lm1_mcmc_combi[,"sigma_j_v5[3]"],"sigma_j_v5[4]"=lm1_mcmc_combi[,"sigma_j_v5[4]"],
                                      "sigma_j_v5[5]"=lm1_mcmc_combi[,"sigma_j_v5[5]"],"sigma_j_v5[6]"=lm1_mcmc_combi[,"sigma_j_v5[6]"],"sigma_j_v5[7]"=lm1_mcmc_combi[,"sigma_j_v5[7]"],"sigma_j_v5[8]"=lm1_mcmc_combi[,"sigma_j_v5[8]"],
                                      "sigma_j_v5[9]"=lm1_mcmc_combi[,"sigma_j_v5[9]"],"sigma_j_v5[10]"=lm1_mcmc_combi[,"sigma_j_v5[10]"],"sigma_j_v5[11]"=lm1_mcmc_combi[,"sigma_j_v5[11]"],"sigma_j_v5[12]"=lm1_mcmc_combi[,"sigma_j_v5[12]"],
                                      "sigma_j_v5[13]"=lm1_mcmc_combi[,"sigma_j_v5[13]"],"sigma_j_v5[14]"=lm1_mcmc_combi[,"sigma_j_v5[14]"],"sigma_j_v5[15]"=lm1_mcmc_combi[,"sigma_j_v5[15]"],
                                      "sigma_j_v5[16]"=lm1_mcmc_combi[,"sigma_j_v5[16]"],"sigma_j_v5[17]"=lm1_mcmc_combi[,"sigma_j_v5[17]"],"sigma_j_v5[18]"=lm1_mcmc_combi[,"sigma_j_v5[18]"],
                                      "sigma_j_v5[19]"=lm1_mcmc_combi[,"sigma_j_v5[19]"],"sigma_j_v5[20]"=lm1_mcmc_combi[,"sigma_j_v5[20]"],"sigma_j_v5[21]"=lm1_mcmc_combi[,"sigma_j_v5[21]"],
                                      "sigma_j_v5[22]"=lm1_mcmc_combi[,"sigma_j_v5[22]"],"sigma_j_v5[23]"=lm1_mcmc_combi[,"sigma_j_v5[23]"],"sigma_j_v5[24]"=lm1_mcmc_combi[,"sigma_j_v5[24]"],
                                      "sigma_j_v5[25]"=lm1_mcmc_combi[,"sigma_j_v5[25]"],"sigma_j_v5[26]"=lm1_mcmc_combi[,"sigma_j_v5[26]"],"sigma_j_v5[27]"=lm1_mcmc_combi[,"sigma_j_v5[27]"],
                                      "sigma_j_v5[28]"=lm1_mcmc_combi[,"sigma_j_v5[28]"],"sigma_j_v5[29]"=lm1_mcmc_combi[,"sigma_j_v5[29]"],"sigma_j_v5[30]"=lm1_mcmc_combi[,"sigma_j_v5[30]"],
                                      "sigma_j_v5[31]"=lm1_mcmc_combi[,"sigma_j_v5[31]"],"sigma_j_v5[32]"=lm1_mcmc_combi[,"sigma_j_v5[32]"],"sigma_j_v5[33]"=lm1_mcmc_combi[,"sigma_j_v5[33]"],
                                      "sigma_j_v5[34]"=lm1_mcmc_combi[,"sigma_j_v5[34]"],"sigma_j_v5[35]"=lm1_mcmc_combi[,"sigma_j_v5[35]"],"sigma_j_v5[36]"=lm1_mcmc_combi[,"sigma_j_v5[36]"],
                                      "sigma_j_v5[37]"=lm1_mcmc_combi[,"sigma_j_v5[37]"],"sigma_j_v5[38]"=lm1_mcmc_combi[,"sigma_j_v5[38]"],"sigma_j_v5[39]"=lm1_mcmc_combi[,"sigma_j_v5[39]"],
                                      "sigma_j_v5[40]"=lm1_mcmc_combi[,"sigma_j_v5[40]"]), 1))                                               
lm1_mcmc_y <- do.call(cbind, rep(list("ar_y1[1]"=lm1_mcmc_combi[,"ar_y1[1]"], "ar_y1[2]"=lm1_mcmc_combi[,"ar_y1[2]"],"ar_y1[3]"=lm1_mcmc_combi[,"ar_y1[3]"],"ar_y1[4]"=lm1_mcmc_combi[,"ar_y1[4]"], "ar_y1[5]"=lm1_mcmc_combi[,"ar_y1[5]"], "ar_y1[6]"=lm1_mcmc_combi[,"ar_y1[6]"],"ar_y1[7]"=lm1_mcmc_combi[,"ar_y1[7]"], "ar_y1[8]"=lm1_mcmc_combi[,"ar_y1[8]"], "ar_y1[9]"=lm1_mcmc_combi[,"ar_y1[9]"], "ar_y1[10]"=lm1_mcmc_combi[,"ar_y1[10]"],"ar_y1[11]"=lm1_mcmc_combi[,"ar_y1[11]"],"ar_y1[12]"=lm1_mcmc_combi[,"ar_y1[12]"],"ar_y1[13]"=lm1_mcmc_combi[,"ar_y1[13]"],"ar_y1[14]"=lm1_mcmc_combi[,"ar_y1[14]"],"ar_y1[15]"=lm1_mcmc_combi[,"ar_y1[15]"],
                                      "ar_y1[16]"=lm1_mcmc_combi[,"ar_y1[16]"],"ar_y1[17]"=lm1_mcmc_combi[,"ar_y1[17]"],"ar_y1[18]"=lm1_mcmc_combi[,"ar_y1[18]"],
                                      "ar_y1[19]"=lm1_mcmc_combi[,"ar_y1[19]"],"ar_y1[20]"=lm1_mcmc_combi[,"ar_y1[20]"],"ar_y1[21]"=lm1_mcmc_combi[,"ar_y1[21]"], 
                                      "ar_y1[22]"=lm1_mcmc_combi[,"ar_y1[22]"],"ar_y1[23]"=lm1_mcmc_combi[,"ar_y1[23]"],"ar_y1[24]"=lm1_mcmc_combi[,"ar_y1[24]"],
                                      "ar_y1[25]"=lm1_mcmc_combi[,"ar_y1[25]"],"ar_y1[26]"=lm1_mcmc_combi[,"ar_y1[26]"],"ar_y1[27]"=lm1_mcmc_combi[,"ar_y1[27]"],
                                      "ar_y1[28]"=lm1_mcmc_combi[,"ar_y1[28]"],"ar_y1[29]"=lm1_mcmc_combi[,"ar_y1[29]"],"ar_y1[30]"=lm1_mcmc_combi[,"ar_y1[30]"],
                                      "ar_y1[31]"=lm1_mcmc_combi[,"ar_y1[31]"],"ar_y1[32]"=lm1_mcmc_combi[,"ar_y1[32]"],"ar_y1[33]"=lm1_mcmc_combi[,"ar_y1[33]"],
                                      "ar_y1[34]"=lm1_mcmc_combi[,"ar_y1[34]"],"ar_y1[35]"=lm1_mcmc_combi[,"ar_y1[35]"], "ar_y1[36]"=lm1_mcmc_combi[,"ar_y1[36]"],
                                      "ar_y1[37]"=lm1_mcmc_combi[,"ar_y1[37]"],"ar_y1[38]"=lm1_mcmc_combi[,"ar_y1[38]"],"ar_y1[39]"=lm1_mcmc_combi[,"ar_y1[39]"],
                                      "ar_y1[40]"=lm1_mcmc_combi[,"ar_y1[40]"],
                                      "ar_y2[1]"=lm1_mcmc_combi[,"ar_y2[1]"], "ar_y2[2]"=lm1_mcmc_combi[,"ar_y2[2]"],"ar_y2[3]"=lm1_mcmc_combi[,"ar_y2[3]"],"ar_y2[4]"=lm1_mcmc_combi[,"ar_y2[4]"],
                                      "ar_y2[5]"=lm1_mcmc_combi[,"ar_y2[5]"], "ar_y2[6]"=lm1_mcmc_combi[,"ar_y2[6]"],"ar_y2[7]"=lm1_mcmc_combi[,"ar_y2[7]"], "ar_y2[8]"=lm1_mcmc_combi[,"ar_y2[8]"],
                                      "ar_y2[9]"=lm1_mcmc_combi[,"ar_y2[9]"], "ar_y2[10]"=lm1_mcmc_combi[,"ar_y2[10]"],"ar_y2[11]"=lm1_mcmc_combi[,"ar_y2[11]"],"ar_y2[12]"=lm1_mcmc_combi[,"ar_y2[12]"],
                                      "ar_y2[13]"=lm1_mcmc_combi[,"ar_y2[13]"],"ar_y2[14]"=lm1_mcmc_combi[,"ar_y2[14]"],"ar_y2[15]"=lm1_mcmc_combi[,"ar_y2[15]"],
                                      "ar_y2[16]"=lm1_mcmc_combi[,"ar_y2[16]"],"ar_y2[17]"=lm1_mcmc_combi[,"ar_y2[17]"],"ar_y2[18]"=lm1_mcmc_combi[,"ar_y2[18]"],
                                      "ar_y2[19]"=lm1_mcmc_combi[,"ar_y2[19]"],"ar_y2[20]"=lm1_mcmc_combi[,"ar_y2[20]"],"ar_y2[21]"=lm1_mcmc_combi[,"ar_y2[21]"], 
                                      "ar_y2[22]"=lm1_mcmc_combi[,"ar_y2[22]"],"ar_y2[23]"=lm1_mcmc_combi[,"ar_y2[23]"],"ar_y2[24]"=lm1_mcmc_combi[,"ar_y2[24]"],
                                      "ar_y2[25]"=lm1_mcmc_combi[,"ar_y2[25]"],"ar_y2[26]"=lm1_mcmc_combi[,"ar_y2[26]"],"ar_y2[27]"=lm1_mcmc_combi[,"ar_y2[27]"],
                                      "ar_y2[28]"=lm1_mcmc_combi[,"ar_y2[28]"],"ar_y2[29]"=lm1_mcmc_combi[,"ar_y2[29]"],"ar_y2[30]"=lm1_mcmc_combi[,"ar_y2[30]"],
                                      "ar_y2[31]"=lm1_mcmc_combi[,"ar_y2[31]"],"ar_y2[32]"=lm1_mcmc_combi[,"ar_y2[32]"],"ar_y2[33]"=lm1_mcmc_combi[,"ar_y2[33]"],
                                      "ar_y2[34]"=lm1_mcmc_combi[,"ar_y2[34]"],"ar_y2[35]"=lm1_mcmc_combi[,"ar_y2[35]"], "ar_y2[36]"=lm1_mcmc_combi[,"ar_y2[36]"],
                                      "ar_y2[37]"=lm1_mcmc_combi[,"ar_y2[37]"],"ar_y2[38]"=lm1_mcmc_combi[,"ar_y2[38]"],"ar_y2[39]"=lm1_mcmc_combi[,"ar_y2[39]"],
                                      "ar_y2[40]"=lm1_mcmc_combi[,"ar_y2[40]"],
                                      "ar_y3[1]"=lm1_mcmc_combi[,"ar_y3[1]"], "ar_y3[2]"=lm1_mcmc_combi[,"ar_y3[2]"],"ar_y3[3]"=lm1_mcmc_combi[,"ar_y3[3]"],"ar_y3[4]"=lm1_mcmc_combi[,"ar_y3[4]"],
                                      "ar_y3[5]"=lm1_mcmc_combi[,"ar_y3[5]"], "ar_y3[6]"=lm1_mcmc_combi[,"ar_y3[6]"],"ar_y3[7]"=lm1_mcmc_combi[,"ar_y3[7]"], "ar_y3[8]"=lm1_mcmc_combi[,"ar_y3[8]"],
                                      "ar_y3[9]"=lm1_mcmc_combi[,"ar_y3[9]"], "ar_y3[10]"=lm1_mcmc_combi[,"ar_y3[10]"],"ar_y3[11]"=lm1_mcmc_combi[,"ar_y3[11]"],"ar_y3[12]"=lm1_mcmc_combi[,"ar_y3[12]"],
                                      "ar_y3[13]"=lm1_mcmc_combi[,"ar_y3[13]"],"ar_y3[14]"=lm1_mcmc_combi[,"ar_y3[14]"],"ar_y3[15]"=lm1_mcmc_combi[,"ar_y3[15]"],
                                      "ar_y3[16]"=lm1_mcmc_combi[,"ar_y3[16]"],"ar_y3[17]"=lm1_mcmc_combi[,"ar_y3[17]"],"ar_y3[18]"=lm1_mcmc_combi[,"ar_y3[18]"],
                                      "ar_y3[19]"=lm1_mcmc_combi[,"ar_y3[19]"],"ar_y3[20]"=lm1_mcmc_combi[,"ar_y3[20]"],"ar_y3[21]"=lm1_mcmc_combi[,"ar_y3[21]"], 
                                      "ar_y3[22]"=lm1_mcmc_combi[,"ar_y3[22]"],"ar_y3[23]"=lm1_mcmc_combi[,"ar_y3[23]"],"ar_y3[24]"=lm1_mcmc_combi[,"ar_y3[24]"],
                                      "ar_y3[25]"=lm1_mcmc_combi[,"ar_y3[25]"],"ar_y3[26]"=lm1_mcmc_combi[,"ar_y3[26]"],"ar_y3[27]"=lm1_mcmc_combi[,"ar_y3[27]"],
                                      "ar_y3[28]"=lm1_mcmc_combi[,"ar_y3[28]"],"ar_y3[29]"=lm1_mcmc_combi[,"ar_y3[29]"],"ar_y3[30]"=lm1_mcmc_combi[,"ar_y3[30]"],
                                      "ar_y3[31]"=lm1_mcmc_combi[,"ar_y3[31]"],"ar_y3[32]"=lm1_mcmc_combi[,"ar_y3[32]"],"ar_y3[33]"=lm1_mcmc_combi[,"ar_y3[33]"],
                                      "ar_y3[34]"=lm1_mcmc_combi[,"ar_y3[34]"],"ar_y3[35]"=lm1_mcmc_combi[,"ar_y3[35]"], "ar_y3[36]"=lm1_mcmc_combi[,"ar_y3[36]"],
                                      "ar_y3[37]"=lm1_mcmc_combi[,"ar_y3[37]"],"ar_y3[38]"=lm1_mcmc_combi[,"ar_y3[38]"],"ar_y3[39]"=lm1_mcmc_combi[,"ar_y3[39]"],
                                      "ar_y3[40]"=lm1_mcmc_combi[,"ar_y3[40]"],
                                      "sigma_j_y1[1]"=lm1_mcmc_combi[,"sigma_j_y1[1]"],"sigma_j_y1[2]"=lm1_mcmc_combi[,"sigma_j_y1[2]"],"sigma_j_y1[3]"=lm1_mcmc_combi[,"sigma_j_y1[3]"],"sigma_j_y1[4]"=lm1_mcmc_combi[,"sigma_j_y1[4]"],
                                      "sigma_j_y1[5]"=lm1_mcmc_combi[,"sigma_j_y1[5]"],"sigma_j_y1[6]"=lm1_mcmc_combi[,"sigma_j_y1[6]"],"sigma_j_y1[7]"=lm1_mcmc_combi[,"sigma_j_y1[7]"],"sigma_j_y1[8]"=lm1_mcmc_combi[,"sigma_j_y1[8]"],
                                      "sigma_j_y1[9]"=lm1_mcmc_combi[,"sigma_j_y1[9]"],"sigma_j_y1[10]"=lm1_mcmc_combi[,"sigma_j_y1[10]"],"sigma_j_y1[11]"=lm1_mcmc_combi[,"sigma_j_y1[11]"],"sigma_j_y1[12]"=lm1_mcmc_combi[,"sigma_j_y1[12]"],
                                      "sigma_j_y1[13]"=lm1_mcmc_combi[,"sigma_j_y1[13]"],"sigma_j_y1[14]"=lm1_mcmc_combi[,"sigma_j_y1[14]"],"sigma_j_y1[15]"=lm1_mcmc_combi[,"sigma_j_y1[15]"],
                                      "sigma_j_y1[16]"=lm1_mcmc_combi[,"sigma_j_y1[16]"],"sigma_j_y1[17]"=lm1_mcmc_combi[,"sigma_j_y1[17]"],"sigma_j_y1[18]"=lm1_mcmc_combi[,"sigma_j_y1[18]"],
                                      "sigma_j_y1[19]"=lm1_mcmc_combi[,"sigma_j_y1[19]"],"sigma_j_y1[20]"=lm1_mcmc_combi[,"sigma_j_y1[20]"],"sigma_j_y1[21]"=lm1_mcmc_combi[,"sigma_j_y1[21]"],
                                      "sigma_j_y1[22]"=lm1_mcmc_combi[,"sigma_j_y1[22]"],"sigma_j_y1[23]"=lm1_mcmc_combi[,"sigma_j_y1[23]"],"sigma_j_y1[24]"=lm1_mcmc_combi[,"sigma_j_y1[24]"],
                                      "sigma_j_y1[25]"=lm1_mcmc_combi[,"sigma_j_y1[25]"],"sigma_j_y1[26]"=lm1_mcmc_combi[,"sigma_j_y1[26]"],"sigma_j_y1[27]"=lm1_mcmc_combi[,"sigma_j_y1[27]"],
                                      "sigma_j_y1[28]"=lm1_mcmc_combi[,"sigma_j_y1[28]"],"sigma_j_y1[29]"=lm1_mcmc_combi[,"sigma_j_y1[29]"],"sigma_j_y1[30]"=lm1_mcmc_combi[,"sigma_j_y1[30]"], 
                                      "sigma_j_y1[31]"=lm1_mcmc_combi[,"sigma_j_y1[31]"],"sigma_j_y1[32]"=lm1_mcmc_combi[,"sigma_j_y1[32]"],"sigma_j_y1[33]"=lm1_mcmc_combi[,"sigma_j_y1[33]"],
                                      "sigma_j_y1[34]"=lm1_mcmc_combi[,"sigma_j_y1[34]"],"sigma_j_y1[35]"=lm1_mcmc_combi[,"sigma_j_y1[35]"],"sigma_j_y1[36]"=lm1_mcmc_combi[,"sigma_j_y1[36]"],
                                      "sigma_j_y1[37]"=lm1_mcmc_combi[,"sigma_j_y1[37]"],"sigma_j_y1[38]"=lm1_mcmc_combi[,"sigma_j_y1[38]"],"sigma_j_y1[39]"=lm1_mcmc_combi[,"sigma_j_y1[39]"],
                                      "sigma_j_y1[40]"=lm1_mcmc_combi[,"sigma_j_y1[40]"],
                                      "sigma_j_y2[1]"=lm1_mcmc_combi[,"sigma_j_y2[1]"],"sigma_j_y2[2]"=lm1_mcmc_combi[,"sigma_j_y2[2]"],"sigma_j_y2[3]"=lm1_mcmc_combi[,"sigma_j_y2[3]"],"sigma_j_y2[4]"=lm1_mcmc_combi[,"sigma_j_y2[4]"],
                                      "sigma_j_y2[5]"=lm1_mcmc_combi[,"sigma_j_y2[5]"],"sigma_j_y2[6]"=lm1_mcmc_combi[,"sigma_j_y2[6]"],"sigma_j_y2[7]"=lm1_mcmc_combi[,"sigma_j_y2[7]"],"sigma_j_y2[8]"=lm1_mcmc_combi[,"sigma_j_y2[8]"],
                                      "sigma_j_y2[9]"=lm1_mcmc_combi[,"sigma_j_y2[9]"],"sigma_j_y2[10]"=lm1_mcmc_combi[,"sigma_j_y2[10]"],"sigma_j_y2[11]"=lm1_mcmc_combi[,"sigma_j_y2[11]"],"sigma_j_y2[12]"=lm1_mcmc_combi[,"sigma_j_y2[12]"],
                                      "sigma_j_y2[13]"=lm1_mcmc_combi[,"sigma_j_y2[13]"],"sigma_j_y2[14]"=lm1_mcmc_combi[,"sigma_j_y2[14]"],"sigma_j_y2[15]"=lm1_mcmc_combi[,"sigma_j_y2[15]"],
                                      "sigma_j_y2[16]"=lm1_mcmc_combi[,"sigma_j_y2[16]"],"sigma_j_y2[17]"=lm1_mcmc_combi[,"sigma_j_y2[17]"],"sigma_j_y2[18]"=lm1_mcmc_combi[,"sigma_j_y2[18]"],
                                      "sigma_j_y2[19]"=lm1_mcmc_combi[,"sigma_j_y2[19]"],"sigma_j_y2[20]"=lm1_mcmc_combi[,"sigma_j_y2[20]"],"sigma_j_y2[21]"=lm1_mcmc_combi[,"sigma_j_y2[21]"],
                                      "sigma_j_y2[22]"=lm1_mcmc_combi[,"sigma_j_y2[22]"],"sigma_j_y2[23]"=lm1_mcmc_combi[,"sigma_j_y2[23]"],"sigma_j_y2[24]"=lm1_mcmc_combi[,"sigma_j_y2[24]"],
                                      "sigma_j_y2[25]"=lm1_mcmc_combi[,"sigma_j_y2[25]"],"sigma_j_y2[26]"=lm1_mcmc_combi[,"sigma_j_y2[26]"],"sigma_j_y2[27]"=lm1_mcmc_combi[,"sigma_j_y2[27]"],
                                      "sigma_j_y2[28]"=lm1_mcmc_combi[,"sigma_j_y2[28]"],"sigma_j_y2[29]"=lm1_mcmc_combi[,"sigma_j_y2[29]"],"sigma_j_y2[30]"=lm1_mcmc_combi[,"sigma_j_y2[30]"], 
                                      "sigma_j_y2[31]"=lm1_mcmc_combi[,"sigma_j_y2[31]"],"sigma_j_y2[32]"=lm1_mcmc_combi[,"sigma_j_y2[32]"],"sigma_j_y2[33]"=lm1_mcmc_combi[,"sigma_j_y2[33]"],
                                      "sigma_j_y2[34]"=lm1_mcmc_combi[,"sigma_j_y2[34]"],"sigma_j_y2[35]"=lm1_mcmc_combi[,"sigma_j_y2[35]"],"sigma_j_y2[36]"=lm1_mcmc_combi[,"sigma_j_y2[36]"],
                                      "sigma_j_y2[37]"=lm1_mcmc_combi[,"sigma_j_y2[37]"],"sigma_j_y2[38]"=lm1_mcmc_combi[,"sigma_j_y2[38]"],"sigma_j_y2[39]"=lm1_mcmc_combi[,"sigma_j_y2[39]"],
                                      "sigma_j_y2[40]"=lm1_mcmc_combi[,"sigma_j_y2[40]"],
                                      "sigma_j_y3[1]"=lm1_mcmc_combi[,"sigma_j_y3[1]"],"sigma_j_y3[2]"=lm1_mcmc_combi[,"sigma_j_y3[2]"],"sigma_j_y3[3]"=lm1_mcmc_combi[,"sigma_j_y3[3]"],"sigma_j_y3[4]"=lm1_mcmc_combi[,"sigma_j_y3[4]"],
                                      "sigma_j_y3[5]"=lm1_mcmc_combi[,"sigma_j_y3[5]"],"sigma_j_y3[6]"=lm1_mcmc_combi[,"sigma_j_y3[6]"],"sigma_j_y3[7]"=lm1_mcmc_combi[,"sigma_j_y3[7]"],"sigma_j_y3[8]"=lm1_mcmc_combi[,"sigma_j_y3[8]"],
                                      "sigma_j_y3[9]"=lm1_mcmc_combi[,"sigma_j_y3[9]"],"sigma_j_y3[10]"=lm1_mcmc_combi[,"sigma_j_y3[10]"],"sigma_j_y3[11]"=lm1_mcmc_combi[,"sigma_j_y3[11]"],"sigma_j_y3[12]"=lm1_mcmc_combi[,"sigma_j_y3[12]"],
                                      "sigma_j_y3[13]"=lm1_mcmc_combi[,"sigma_j_y3[13]"],"sigma_j_y3[14]"=lm1_mcmc_combi[,"sigma_j_y3[14]"],"sigma_j_y3[15]"=lm1_mcmc_combi[,"sigma_j_y3[15]"],
                                      "sigma_j_y3[16]"=lm1_mcmc_combi[,"sigma_j_y3[16]"],"sigma_j_y3[17]"=lm1_mcmc_combi[,"sigma_j_y3[17]"],"sigma_j_y3[18]"=lm1_mcmc_combi[,"sigma_j_y3[18]"],
                                      "sigma_j_y3[19]"=lm1_mcmc_combi[,"sigma_j_y3[19]"],"sigma_j_y3[20]"=lm1_mcmc_combi[,"sigma_j_y3[20]"],"sigma_j_y3[21]"=lm1_mcmc_combi[,"sigma_j_y3[21]"],
                                      "sigma_j_y3[22]"=lm1_mcmc_combi[,"sigma_j_y3[22]"],"sigma_j_y3[23]"=lm1_mcmc_combi[,"sigma_j_y3[23]"],"sigma_j_y3[24]"=lm1_mcmc_combi[,"sigma_j_y3[24]"],
                                      "sigma_j_y3[25]"=lm1_mcmc_combi[,"sigma_j_y3[25]"],"sigma_j_y3[26]"=lm1_mcmc_combi[,"sigma_j_y3[26]"],"sigma_j_y3[27]"=lm1_mcmc_combi[,"sigma_j_y3[27]"],
                                      "sigma_j_y3[28]"=lm1_mcmc_combi[,"sigma_j_y3[28]"],"sigma_j_y3[29]"=lm1_mcmc_combi[,"sigma_j_y3[29]"],"sigma_j_y3[30]"=lm1_mcmc_combi[,"sigma_j_y3[30]"], 
                                      "sigma_j_y3[31]"=lm1_mcmc_combi[,"sigma_j_y3[31]"],"sigma_j_y3[32]"=lm1_mcmc_combi[,"sigma_j_y3[32]"],"sigma_j_y3[33]"=lm1_mcmc_combi[,"sigma_j_y3[33]"],
                                      "sigma_j_y3[34]"=lm1_mcmc_combi[,"sigma_j_y3[34]"],"sigma_j_y3[35]"=lm1_mcmc_combi[,"sigma_j_y3[35]"],"sigma_j_y3[36]"=lm1_mcmc_combi[,"sigma_j_y3[36]"],
                                      "sigma_j_y3[37]"=lm1_mcmc_combi[,"sigma_j_y3[37]"],"sigma_j_y3[38]"=lm1_mcmc_combi[,"sigma_j_y3[38]"],"sigma_j_y3[39]"=lm1_mcmc_combi[,"sigma_j_y3[39]"],
                                      "sigma_j_y3[40]"=lm1_mcmc_combi[,"sigma_j_y3[40]"]), 1))  
#projection of target level for emission intensity
j=33
pred_gdp <- matrix(NA, nrow = nrow(lm1_mcmc_gdp_pre),ncol=n_2050)
for (i in 1:nrow(lm1_mcmc_gdp_pre)){
  for (t in 1:24){pred_gdp[i,t] <- exp(x[t,j])}
  for (t in 25:25){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*1]}
    else {pred_gdp[i,t]<-0}}
  for (t in 26:26){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*2]}
    else {pred_gdp[i,t]<-0}}
  for (t in 27:27){
    if (pred_gdp[i,t-1] >- 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*3]}
    else {pred_gdp[i,t]<-0}}
  for (t in 28:28){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*4]}
    else {pred_gdp[i,t]<-0}}
  for (t in 29:29){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*5]}
    else {pred_gdp[i,t]<-0}}
  for (t in 30:30){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*6]}
    else {pred_gdp[i,t]<-0}}
  for (t in 31:31){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*7]}
    else {pred_gdp[i,t]<-0}}
  for (t in 32:32){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*8]}
    else {pred_gdp[i,t]<-0}}
  for (t in 33:33){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*9]}
    else {pred_gdp[i,t]<-0}}
  for (t in 34:34){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*10]}
    else {pred_gdp[i,t]<-0}}
  for (t in 35:35){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*11]}
    else {pred_gdp[i,t]<-0}}
  for (t in 36:36){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*12]}
    else {pred_gdp[i,t]<-0}}
  for (t in 37:37){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*13]}
    else {pred_gdp[i,t]<-0}}
  for (t in 38:38){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*14]}
    else {pred_gdp[i,t]<-0}}
  for (t in 39:39){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*15]}
    else {pred_gdp[i,t]<-0}}
  for (t in 40:40){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*16]}
    else {pred_gdp[i,t]<-0}}
  for (t in 41:41){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*17]}
    else {pred_gdp[i,t]<-0}}
  for (t in 42:42){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*18]}
    else {pred_gdp[i,t]<-0}}
  for (t in 43:43){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*19]}
    else {pred_gdp[i,t]<-0}}
  for (t in 44:44){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*20]}
    else {pred_gdp[i,t]<-0}}
  for (t in 45:45){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*21]}
    else {pred_gdp[i,t]<-0}}
  for (t in 46:46){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*22]}
    else {pred_gdp[i,t]<-0}}
  for (t in 47:47){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*23]}
    else {pred_gdp[i,t]<-0}}
  for (t in 48:48){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*24]}
    else {pred_gdp[i,t]<-0}}
  for (t in 49:49){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*25]}
    else {pred_gdp[i,t]<-0}}
  for (t in 50:50){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*26]}
    else {pred_gdp[i,t]<-0}}
  for (t in 51:51){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*27]}
    else {pred_gdp[i,t]<-0}}
  for (t in 52:52){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*28]}
    else {pred_gdp[i,t]<-0}}
  for (t in 53:53){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*29]}
    else {pred_gdp[i,t]<-0}}
  for (t in 54:54){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*30]}
    else {pred_gdp[i,t]<-0}}
  for (t in 55:55){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*31]}
    else {pred_gdp[i,t]<-0}}
  for (t in 56:56){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*32]}
    else {pred_gdp[i,t]<-0}}
}
for (i in 1:nrow(lm1_mcmc_gdp_pre)){
  for (t in 1:n_2050){
    if (pred_gdp[i,t] > 0){pred_gdp[i,t] <- pred_gdp[i,t]}
    else {pred_gdp[i,t]<-10^(-50)}
  }
}

pred_eihr1 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)
pred_eihr2 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)
pred_eihr3 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)
pred_eihr4 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)
pred_eihr5 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)
pred_eihr6 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)

pred_eivr1 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)
pred_eivr2 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)
pred_eivr3 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)
pred_eivr4 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)
pred_eivr5 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)

pred_eigr1 <- matrix(NA, nrow = nrow(lm1_mcmc_ref), ncol=n_2050)

for (i in 1:nrow(lm1_mcmc_ref)){
  for (t in 1:n_2050){
    pred_eihr1[i,t] <- (lm1_mcmc_ref[i,1]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,7])+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,13])
    pred_eihr3[i,t] <- (lm1_mcmc_ref[i,3]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,9]) +rnorm(1,mean=0,sd=lm1_mcmc_ref[i,15])
    pred_eihr6[i,t] <- (lm1_mcmc_ref[i,6]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,12]) +rnorm(1,mean=0,sd=lm1_mcmc_ref[i,18])
    pred_eivr2[i,t] <- (lm1_mcmc_ref[i,20]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,25]) +rnorm(1,mean=0,sd=lm1_mcmc_ref[i,30])
    pred_eivr3[i,t] <- (lm1_mcmc_ref[i,21]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,26]) +rnorm(1,mean=0,sd=lm1_mcmc_ref[i,31])
    pred_eivr4[i,t] <- (lm1_mcmc_ref[i,22]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,27]) +rnorm(1,mean=0,sd=lm1_mcmc_ref[i,32])
  }
}

j=20
pred_gdp <- matrix(NA, nrow = nrow(lm1_mcmc_gdp_pre),ncol=n_2050)
for (i in 1:nrow(lm1_mcmc_gdp_pre)){
  for (t in 1:24){pred_gdp[i,t] <- exp(x[t,j])}
  for (t in 25:25){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*1]}
    else {pred_gdp[i,t]<-0}}
  for (t in 26:26){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*2]}
    else {pred_gdp[i,t]<-0}}
  for (t in 27:27){
    if (pred_gdp[i,t-1] >- 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*3]}
    else {pred_gdp[i,t]<-0}}
  for (t in 28:28){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*4]}
    else {pred_gdp[i,t]<-0}}
  for (t in 29:29){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*5]}
    else {pred_gdp[i,t]<-0}}
  for (t in 30:30){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*6]}
    else {pred_gdp[i,t]<-0}}
  for (t in 31:31){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*7]}
    else {pred_gdp[i,t]<-0}}
  for (t in 32:32){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*8]}
    else {pred_gdp[i,t]<-0}}
  for (t in 33:33){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*9]}
    else {pred_gdp[i,t]<-0}}
  for (t in 34:34){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*10]}
    else {pred_gdp[i,t]<-0}}
  for (t in 35:35){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*11]}
    else {pred_gdp[i,t]<-0}}
  for (t in 36:36){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*12]}
    else {pred_gdp[i,t]<-0}}
  for (t in 37:37){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*13]}
    else {pred_gdp[i,t]<-0}}
  for (t in 38:38){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*14]}
    else {pred_gdp[i,t]<-0}}
  for (t in 39:39){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*15]}
    else {pred_gdp[i,t]<-0}}
  for (t in 40:40){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*16]}
    else {pred_gdp[i,t]<-0}}
  for (t in 41:41){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*17]}
    else {pred_gdp[i,t]<-0}}
  for (t in 42:42){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*18]}
    else {pred_gdp[i,t]<-0}}
  for (t in 43:43){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*19]}
    else {pred_gdp[i,t]<-0}}
  for (t in 44:44){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*20]}
    else {pred_gdp[i,t]<-0}}
  for (t in 45:45){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*21]}
    else {pred_gdp[i,t]<-0}}
  for (t in 46:46){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*22]}
    else {pred_gdp[i,t]<-0}}
  for (t in 47:47){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*23]}
    else {pred_gdp[i,t]<-0}}
  for (t in 48:48){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*24]}
    else {pred_gdp[i,t]<-0}}
  for (t in 49:49){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*25]}
    else {pred_gdp[i,t]<-0}}
  for (t in 50:50){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*26]}
    else {pred_gdp[i,t]<-0}}
  for (t in 51:51){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*27]}
    else {pred_gdp[i,t]<-0}}
  for (t in 52:52){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*28]}
    else {pred_gdp[i,t]<-0}}
  for (t in 53:53){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*29]}
    else {pred_gdp[i,t]<-0}}
  for (t in 54:54){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*30]}
    else {pred_gdp[i,t]<-0}}
  for (t in 55:55){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*31]}
    else {pred_gdp[i,t]<-0}}
  for (t in 56:56){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*32]}
    else {pred_gdp[i,t]<-0}}
}
for (i in 1:nrow(lm1_mcmc_gdp_pre)){
  for (t in 1:n_2050){
    if (pred_gdp[i,t] > 0){pred_gdp[i,t] <- pred_gdp[i,t]}
    else {pred_gdp[i,t]<-10^(-50)}
  }
}
for (i in 1:nrow(lm1_mcmc_ref)){
  for (t in 1:n_2050){
    pred_eihr4[i,t] <- (lm1_mcmc_ref[i,4]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,10]) +rnorm(1,mean=0,sd=lm1_mcmc_ref[i,16])
    pred_eivr5[i,t] <- (lm1_mcmc_ref[i,23]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,28])+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,33])
    pred_eigr1[i,t] <- (lm1_mcmc_ref[i,34]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,35]) +rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36])
  }
}

j=34
pred_gdp <- matrix(NA, nrow = nrow(lm1_mcmc_gdp_pre),ncol=n_2050)
for (i in 1:nrow(lm1_mcmc_gdp_pre)){
  for (t in 1:24){pred_gdp[i,t] <- exp(x[t,j])}
  for (t in 25:25){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*1]}
    else {pred_gdp[i,t]<-0}}
  for (t in 26:26){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*2]}
    else {pred_gdp[i,t]<-0}}
  for (t in 27:27){
    if (pred_gdp[i,t-1] >- 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*3]}
    else {pred_gdp[i,t]<-0}}
  for (t in 28:28){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*4]}
    else {pred_gdp[i,t]<-0}}
  for (t in 29:29){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*5]}
    else {pred_gdp[i,t]<-0}}
  for (t in 30:30){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*6]}
    else {pred_gdp[i,t]<-0}}
  for (t in 31:31){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*7]}
    else {pred_gdp[i,t]<-0}}
  for (t in 32:32){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*8]}
    else {pred_gdp[i,t]<-0}}
  for (t in 33:33){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*9]}
    else {pred_gdp[i,t]<-0}}
  for (t in 34:34){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*10]}
    else {pred_gdp[i,t]<-0}}
  for (t in 35:35){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*11]}
    else {pred_gdp[i,t]<-0}}
  for (t in 36:36){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*12]}
    else {pred_gdp[i,t]<-0}}
  for (t in 37:37){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*13]}
    else {pred_gdp[i,t]<-0}}
  for (t in 38:38){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*14]}
    else {pred_gdp[i,t]<-0}}
  for (t in 39:39){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*15]}
    else {pred_gdp[i,t]<-0}}
  for (t in 40:40){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*16]}
    else {pred_gdp[i,t]<-0}}
  for (t in 41:41){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*17]}
    else {pred_gdp[i,t]<-0}}
  for (t in 42:42){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*18]}
    else {pred_gdp[i,t]<-0}}
  for (t in 43:43){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*19]}
    else {pred_gdp[i,t]<-0}}
  for (t in 44:44){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*20]}
    else {pred_gdp[i,t]<-0}}
  for (t in 45:45){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*21]}
    else {pred_gdp[i,t]<-0}}
  for (t in 46:46){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*22]}
    else {pred_gdp[i,t]<-0}}
  for (t in 47:47){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*23]}
    else {pred_gdp[i,t]<-0}}
  for (t in 48:48){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*24]}
    else {pred_gdp[i,t]<-0}}
  for (t in 49:49){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*25]}
    else {pred_gdp[i,t]<-0}}
  for (t in 50:50){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*26]}
    else {pred_gdp[i,t]<-0}}
  for (t in 51:51){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*27]}
    else {pred_gdp[i,t]<-0}}
  for (t in 52:52){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*28]}
    else {pred_gdp[i,t]<-0}}
  for (t in 53:53){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*29]}
    else {pred_gdp[i,t]<-0}}
  for (t in 54:54){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*30]}
    else {pred_gdp[i,t]<-0}}
  for (t in 55:55){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*31]}
    else {pred_gdp[i,t]<-0}}
  for (t in 56:56){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*32]}
    else {pred_gdp[i,t]<-0}}
}
for (i in 1:nrow(lm1_mcmc_gdp_pre)){
  for (t in 1:n_2050){
    if (pred_gdp[i,t] > 0){pred_gdp[i,t] <- pred_gdp[i,t]}
    else {pred_gdp[i,t]<-10^(-50)}
  }
}
for (i in 1:nrow(lm1_mcmc_ref)){
  for (t in 1:n_2050){
    pred_eihr2[i,t] <- (lm1_mcmc_ref[i,2]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,8])+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,14])
    pred_eihr5[i,t] <- (lm1_mcmc_ref[i,5]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,11]) +rnorm(1,mean=0,sd=lm1_mcmc_ref[i,17])
  }
}

j=35
pred_gdp <- matrix(NA, nrow = nrow(lm1_mcmc_gdp_pre),ncol=n_2050)
for (i in 1:nrow(lm1_mcmc_gdp_pre)){
  for (t in 1:24){pred_gdp[i,t] <- exp(x[t,j])}
  for (t in 25:25){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*1]}
    else {pred_gdp[i,t]<-0}}
  for (t in 26:26){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*2]}
    else {pred_gdp[i,t]<-0}}
  for (t in 27:27){
    if (pred_gdp[i,t-1] >- 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*3]}
    else {pred_gdp[i,t]<-0}}
  for (t in 28:28){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*4]}
    else {pred_gdp[i,t]<-0}}
  for (t in 29:29){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*5]}
    else {pred_gdp[i,t]<-0}}
  for (t in 30:30){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*6]}
    else {pred_gdp[i,t]<-0}}
  for (t in 31:31){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*7]}
    else {pred_gdp[i,t]<-0}}
  for (t in 32:32){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*8]}
    else {pred_gdp[i,t]<-0}}
  for (t in 33:33){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*9]}
    else {pred_gdp[i,t]<-0}}
  for (t in 34:34){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*10]}
    else {pred_gdp[i,t]<-0}}
  for (t in 35:35){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*11]}
    else {pred_gdp[i,t]<-0}}
  for (t in 36:36){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*12]}
    else {pred_gdp[i,t]<-0}}
  for (t in 37:37){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*13]}
    else {pred_gdp[i,t]<-0}}
  for (t in 38:38){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*14]}
    else {pred_gdp[i,t]<-0}}
  for (t in 39:39){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*15]}
    else {pred_gdp[i,t]<-0}}
  for (t in 40:40){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*16]}
    else {pred_gdp[i,t]<-0}}
  for (t in 41:41){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*17]}
    else {pred_gdp[i,t]<-0}}
  for (t in 42:42){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*18]}
    else {pred_gdp[i,t]<-0}}
  for (t in 43:43){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*19]}
    else {pred_gdp[i,t]<-0}}
  for (t in 44:44){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*20]}
    else {pred_gdp[i,t]<-0}}
  for (t in 45:45){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*21]}
    else {pred_gdp[i,t]<-0}}
  for (t in 46:46){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*22]}
    else {pred_gdp[i,t]<-0}}
  for (t in 47:47){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*23]}
    else {pred_gdp[i,t]<-0}}
  for (t in 48:48){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*24]}
    else {pred_gdp[i,t]<-0}}
  for (t in 49:49){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*25]}
    else {pred_gdp[i,t]<-0}}
  for (t in 50:50){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*26]}
    else {pred_gdp[i,t]<-0}}
  for (t in 51:51){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*27]}
    else {pred_gdp[i,t]<-0}}
  for (t in 52:52){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*28]}
    else {pred_gdp[i,t]<-0}}
  for (t in 53:53){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*29]}
    else {pred_gdp[i,t]<-0}}
  for (t in 54:54){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*30]}
    else {pred_gdp[i,t]<-0}}
  for (t in 55:55){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*31]}
    else {pred_gdp[i,t]<-0}}
  for (t in 56:56){
    if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*32]}
    else {pred_gdp[i,t]<-0}}
}
for (i in 1:nrow(lm1_mcmc_gdp_pre)){
  for (t in 1:n_2050){
    if (pred_gdp[i,t] > 0){pred_gdp[i,t] <- pred_gdp[i,t]}
    else {pred_gdp[i,t]<-10^(-50)}
  }
}

for (i in 1:nrow(lm1_mcmc_ref)){
  for (t in 1:n_2050){
    pred_eivr1[i,t] <- (lm1_mcmc_ref[i,19]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,24])+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,29])
  }
}
#projection of emissions for 24 EU countries
#historical emissions
his_e= as.data.frame(matrix(NA, nrow =13,ncol=40))
his_e[1,]<-fh[21,]*so[21,]*yo[21,]*ro[21,]*gdp[21,]*p[21,]/1000
his_e[2,]<-fh[42,]*so[42,]*yo[21,]*ro[21,]*gdp[21,]*p[21,]/1000
his_e[3,]<-fh[63,]*so[63,]*yo[21,]*ro[21,]*gdp[21,]*p[21,]/1000
his_e[4,]<-fh[84,]*so[84,]*yo[21,]*ro[21,]*gdp[21,]*p[21,]/1000
his_e[5,]<-fh[105,]*so[105,]*yo[21,]*ro[21,]*gdp[21,]*p[21,]/1000
his_e[6,]<-fh[126,]*so[126,]*yo[21,]*ro[21,]*gdp[21,]*p[21,]/1000

his_e[7,]<-fv[21,]*vo[21,]*yo[50,]*ro[21,]*gdp[21,]*p[21,]/1000
his_e[8,]<-fv[42,]*vo[42,]*yo[50,]*ro[21,]*gdp[21,]*p[21,]/1000
his_e[9,]<-fv[63,]*vo[63,]*yo[50,]*ro[21,]*gdp[21,]*p[21,]/1000
his_e[10,]<-fv[84,]*vo[84,]*yo[50,]*ro[21,]*gdp[21,]*p[21,]/1000
his_e[11,]<-fv[105,]*vo[105,]*yo[50,]*ro[21,]*gdp[21,]*p[21,]/1000

his_e[12,]<-fg[21,]*go[21,]*yo[79,]*ro[21,]*gdp[21,]*p[21,]/1000

his_e[13,]<-his_e[1,]+his_e[2,]+his_e[3,]+his_e[4,]+his_e[5,]+his_e[6,]+his_e[7,]+his_e[8,]+his_e[9,]+his_e[10,]+his_e[11,]+his_e[12,]
for (j in 1:40){
  for (i in 1:nrow(his_e)){
    if (his_e[i,j] >- 0){his_e[i,j]=his_e[i,j]}
    else {his_e[i,j]<-10^(-50)}
  }
}


his_o_40= as.data.frame(matrix(NA, nrow =7,ncol=40))
his_o_40[1,]<-(so[21,]*yo[21,])*ro[21,]*gdp[21,]*p[21,]
his_o_40[2,]<-(so[42,]*yo[21,])*ro[21,]*gdp[21,]*p[21,]
his_o_40[3,]<-(so[63,]*yo[21,]+vo[21,]*yo[50,])*ro[21,]*gdp[21,]*p[21,]
his_o_40[4,]<-(so[84,]*yo[21,]+vo[42,]*yo[50,])*ro[21,]*gdp[21,]*p[21,]
his_o_40[5,]<-(so[105,]*yo[21,]+vo[63,]*yo[50,])*ro[21,]*gdp[21,]*p[21,]
his_o_40[6,]<-(so[126,]*yo[21,]+vo[84,]*yo[50,]+go[21,]*yo[79,])*ro[21,]*gdp[21,]*p[21,]
his_o_40[7,]<-(vo[105,]*yo[50,])*ro[21,]*gdp[21,]*p[21,]

his_e_40= as.data.frame(matrix(NA, nrow =7,ncol=40))
his_e_40[1,]<-(his_e[1,])
his_e_40[2,]<-(his_e[2,])
his_e_40[3,]<-(his_e[3,]+his_e[7,])
his_e_40[4,]<-(his_e[4,]+his_e[8,])
his_e_40[5,]<-(his_e[5,]+his_e[9,])
his_e_40[6,]<-(his_e[6,]+his_e[10,]+his_e[12,])
his_e_40[7,]<-(his_e[11,])

his_o= as.data.frame(matrix(NA, nrow =7,ncol=18))
his_o[,1]<-his_o_40[,1]+his_o_40[,2]+his_o_40[,3]+his_o_40[,4]+his_o_40[,5]+his_o_40[,6]+his_o_40[,7]+his_o_40[,8]+his_o_40[,9]+his_o_40[,10]+his_o_40[,11]+his_o_40[,12]+his_o_40[,13]+his_o_40[,14]+his_o_40[,15]+his_o_40[,16]+his_o_40[,17]+his_o_40[,18]+his_o_40[,19]+his_o_40[,20]+his_o_40[,21]+his_o_40[,22]+his_o_40[,23]+his_o_40[,24]
for (j in 2:17){
  his_o[,j]=his_o_40[,j+23]
}
his_o[,18]<-his_o_40[,1]+his_o_40[,2]+his_o_40[,3]+his_o_40[,4]+his_o_40[,5]+his_o_40[,6]+his_o_40[,7]+his_o_40[,8]+his_o_40[,9]+his_o_40[,10]+his_o_40[,11]+his_o_40[,12]+his_o_40[,13]+his_o_40[,14]+his_o_40[,15]+his_o_40[,16]+his_o_40[,17]+his_o_40[,18]+his_o_40[,19]+his_o_40[,20]+his_o_40[,21]+his_o_40[,22]+his_o_40[,23]+his_o_40[,24]+his_o_40[,25]+his_o_40[,26]+his_o_40[,27]+his_o_40[,28]+his_o_40[,29]+his_o_40[,30]+his_o_40[,31]+his_o_40[,32]+his_o_40[,33]+his_o_40[,34]+his_o_40[,35]+his_o_40[,36]+his_o_40[,37]+his_o_40[,38]+his_o_40[,39]+his_o_40[,40]

his_e_17= as.data.frame(matrix(NA, nrow =7,ncol=18))
his_e_17[,1]<-his_e_40[,1]+his_e_40[,2]+his_e_40[,3]+his_e_40[,4]+his_e_40[,5]+his_e_40[,6]+his_e_40[,7]+his_e_40[,8]+his_e_40[,9]+his_e_40[,10]+his_e_40[,11]+his_e_40[,12]+his_e_40[,13]+his_e_40[,14]+his_e_40[,15]+his_e_40[,16]+his_e_40[,17]+his_e_40[,18]+his_e_40[,19]+his_e_40[,20]+his_e_40[,21]+his_e_40[,22]+his_e_40[,23]+his_e_40[,24]
for (j in 2:17){
  his_e_17[,j]=his_e_40[,j+23]
}  
his_e_17[,18]<-his_e_40[,1]+his_e_40[,2]+his_e_40[,3]+his_e_40[,4]+his_e_40[,5]+his_e_40[,6]+his_e_40[,7]+his_e_40[,8]+his_e_40[,9]+his_e_40[,10]+his_e_40[,11]+his_e_40[,12]+his_e_40[,13]+his_e_40[,14]+his_e_40[,15]+his_e_40[,16]+his_e_40[,17]+his_e_40[,18]+his_e_40[,19]+his_e_40[,20]+his_e_40[,21]+his_e_40[,22]+his_e_40[,23]+his_e_40[,24]+his_e_40[,25]+his_e_40[,26]+his_e_40[,27]+his_e_40[,28]+his_e_40[,29]+his_e_40[,30]+his_e_40[,31]+his_e_40[,32]+his_e_40[,33]+his_e_40[,34]+his_e_40[,35]+his_e_40[,36]+his_e_40[,37]+his_e_40[,38]+his_e_40[,39]+his_e_40[,40]


his_ei= as.data.frame(matrix(NA, nrow =7,ncol=18))
for (i in 1:7){
  his_ei[i,]<-his_e_17[i,]/his_o[i,]*1000
}

#projecting emissions
pred_e_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_e1_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_e2_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_e3_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_e4_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_e5_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_e6_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_e7_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_e_s_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_e_v_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_e_g_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)

pred_ei1_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei2_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei3_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei4_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei5_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei6_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei7_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei_euro<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)

fac_f1_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_f2_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_f3_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_f4_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_f5_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_f6_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_f7_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_s1_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_s2_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_s3_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_s4_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_s5_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_s6_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_s7_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_f_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_s_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_y_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_r_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_gdp_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)
fac_p_reg<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=18)

for (j in 1:24){
  pred_p <- matrix(NA, nrow = nrow(lm1_mcmc_p),ncol=n_2050)
  for (i in 1:nrow(lm1_mcmc_p)){
    for (t in 1:26){pred_p[i,t] <- exp(x1[t,j])}
    for (t in 27:56){
      if (pred_p[i,t-1] > 0){pred_p[i,t]<-exp(rnorm(1,mean=0,sd=sd(lm1_mcmc_p[,j+80]))+log(pred_p[i,t-1])+d1[t-1,j])+rnorm(1,mean=0,sd=lm1_mcmc_p[i,j+40])}
      else {pred_p[i,t]<-0}}
  }
  for (i in 1:nrow(lm1_mcmc_p)){
    for (t in 1:n_2050){
      if (pred_p[i,t] > 0){pred_p[i,t] <- pred_p[i,t]}
      else {pred_p[i,t]<-10^(-50)}
    }
  }
  
  pred_gdp <- matrix(NA, nrow = nrow(lm1_mcmc_gdp_pre),ncol=n_2050)
  for (i in 1:nrow(lm1_mcmc_gdp_pre)){
    for (t in 1:24){pred_gdp[i,t] <- exp(x[t,j])}
    for (t in 25:25){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*1]}
      else {pred_gdp[i,t]<-0}}
    for (t in 26:26){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*2]}
      else {pred_gdp[i,t]<-0}}
    for (t in 27:27){
      if (pred_gdp[i,t-1] >- 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*3]}
      else {pred_gdp[i,t]<-0}}
    for (t in 28:28){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*4]}
      else {pred_gdp[i,t]<-0}}
    for (t in 29:29){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*5]}
      else {pred_gdp[i,t]<-0}}
    for (t in 30:30){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*6]}
      else {pred_gdp[i,t]<-0}}
    for (t in 31:31){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*7]}
      else {pred_gdp[i,t]<-0}}
    for (t in 32:32){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*8]}
      else {pred_gdp[i,t]<-0}}
    for (t in 33:33){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*9]}
      else {pred_gdp[i,t]<-0}}
    for (t in 34:34){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*10]}
      else {pred_gdp[i,t]<-0}}
    for (t in 35:35){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*11]}
      else {pred_gdp[i,t]<-0}}
    for (t in 36:36){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*12]}
      else {pred_gdp[i,t]<-0}}
    for (t in 37:37){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*13]}
      else {pred_gdp[i,t]<-0}}
    for (t in 38:38){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*14]}
      else {pred_gdp[i,t]<-0}}
    for (t in 39:39){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*15]}
      else {pred_gdp[i,t]<-0}}
    for (t in 40:40){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*16]}
      else {pred_gdp[i,t]<-0}}
    for (t in 41:41){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*17]}
      else {pred_gdp[i,t]<-0}}
    for (t in 42:42){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*18]}
      else {pred_gdp[i,t]<-0}}
    for (t in 43:43){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*19]}
      else {pred_gdp[i,t]<-0}}
    for (t in 44:44){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*20]}
      else {pred_gdp[i,t]<-0}}
    for (t in 45:45){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*21]}
      else {pred_gdp[i,t]<-0}}
    for (t in 46:46){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*22]}
      else {pred_gdp[i,t]<-0}}
    for (t in 47:47){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*23]}
      else {pred_gdp[i,t]<-0}}
    for (t in 48:48){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*24]}
      else {pred_gdp[i,t]<-0}}
    for (t in 49:49){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*25]}
      else {pred_gdp[i,t]<-0}}
    for (t in 50:50){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*26]}
      else {pred_gdp[i,t]<-0}}
    for (t in 51:51){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*27]}
      else {pred_gdp[i,t]<-0}}
    for (t in 52:52){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*28]}
      else {pred_gdp[i,t]<-0}}
    for (t in 53:53){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*29]}
      else {pred_gdp[i,t]<-0}}
    for (t in 54:54){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*30]}
      else {pred_gdp[i,t]<-0}}
    for (t in 55:55){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*31]}
      else {pred_gdp[i,t]<-0}}
    for (t in 56:56){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*32]}
      else {pred_gdp[i,t]<-0}}
  }
  for (i in 1:nrow(lm1_mcmc_gdp_pre)){
    for (t in 1:n_2050){
      if (pred_gdp[i,t] > 0){pred_gdp[i,t] <- pred_gdp[i,t]}
      else {pred_gdp[i,t]<-10^(-50)}
    }
  }
  
  pred_r <- matrix(NA, nrow = nrow(lm1_mcmc_p),ncol=n_2050)
  for (i in 1:nrow(lm1_mcmc_r)){
    for (t in 1:n_2015){pred_r[i,t] <- x2[t,j]}
    for (t in 22:56){
      if (pred_r[i,t-1] > 0){pred_r[i,t]<-rnorm(1,mean=0,sd=lm1_mcmc_r[i,j+80])+pred_r[i,t-1]+d2[t-1,j]+rnorm(1,mean=0,sd=lm1_mcmc_r[i,j+40])}
      else {pred_r[i,t]<-0}}
  }
  for (i in 1:nrow(lm1_mcmc_r)){
    for (t in 1:n_2050){
      if (pred_r[i,t] > 0){pred_r[i,t] <- pred_r[i,t]}
      else {pred_r[i,t]<-10^(-50)}
    }
  }
  
  pred_eih1 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih2 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih3 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih4 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih5 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih6 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv1 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv2 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv3 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv4 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv5 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eig1 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 2:21){
      pred_eih1[i,t-1] <- eih[t,j]
      pred_eih2[i,t-1] <- eih[t+21,j]
      pred_eih3[i,t-1] <- eih[t+21*2,j] 
      pred_eih4[i,t-1] <- eih[t+21*3,j]       
      pred_eih5[i,t-1] <- eih[t+21*4,j] 
      pred_eih6[i,t-1] <- eih[t+21*5,j] 
      
      pred_eiv1[i,t-1] <- eiv[t,j]
      pred_eiv2[i,t-1] <- eiv[t+21,j]
      pred_eiv3[i,t-1] <- eiv[t+21*2,j] 
      pred_eiv4[i,t-1] <- eiv[t+21*3,j]       
      pred_eiv5[i,t-1] <- eiv[t+21*4,j] 
      pred_eig1[i,t-1] <- eig[t,j]
    }
    for (t in 22:n_2050){
      pred_eih1[i,t-1] <-(pred_eihr1[i,t]+(pred_eih1[i,t-2]-pred_eihr1[i,t-1])*lm1_mcmc_fh[i,j+40*0])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240])
      pred_eih2[i,t-1] <-(pred_eihr2[i,t]+(pred_eih2[i,t-2]-pred_eihr2[i,t-1])*lm1_mcmc_fh[i,j+40*1])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*1])
      pred_eih3[i,t-1] <-(pred_eihr3[i,t]+(pred_eih3[i,t-2]-pred_eihr3[i,t-1])*lm1_mcmc_fh[i,j+40*2])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*2])
      pred_eih4[i,t-1] <-(pred_eihr4[i,t]+(pred_eih4[i,t-2]-pred_eihr4[i,t-1])*lm1_mcmc_fh[i,j+40*3])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*3])
      pred_eih5[i,t-1] <-(pred_eihr5[i,t]+(pred_eih5[i,t-2]-pred_eihr5[i,t-1])*lm1_mcmc_fh[i,j+40*4])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*4])
      pred_eih6[i,t-1] <-(pred_eihr6[i,t]+(pred_eih6[i,t-2]-pred_eihr6[i,t-1])*lm1_mcmc_fh[i,j+40*5])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*5])
      pred_eiv1[i,t-1] <-(pred_eivr1[i,t]+(pred_eiv1[i,t-2]-pred_eivr1[i,t-1])*lm1_mcmc_fv[i,j+40*0])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200])
      pred_eiv2[i,t-1] <-(pred_eivr2[i,t]+(pred_eiv2[i,t-2]-pred_eivr2[i,t-1])*lm1_mcmc_fv[i,j+40*1])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*1])
      pred_eiv3[i,t-1] <-(pred_eivr3[i,t]+(pred_eiv3[i,t-2]-pred_eivr3[i,t-1])*lm1_mcmc_fv[i,j+40*2])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*2])
      pred_eiv4[i,t-1] <-(pred_eivr4[i,t]+(pred_eiv4[i,t-2]-pred_eivr4[i,t-1])*lm1_mcmc_fv[i,j+40*3])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*3])
      pred_eiv5[i,t-1] <-(pred_eivr5[i,t]+(pred_eiv5[i,t-2]-pred_eivr5[i,t-1])*lm1_mcmc_fv[i,j+40*4])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*4])
      pred_eig1[i,t-1] <-(pred_eigr1[i,t]+(pred_eig1[i,t-2]-pred_eigr1[i,t-1])*lm1_mcmc_fg[i,j+40*0])+rnorm(1,mean=0,sd=lm1_mcmc_fg[i,j+40])
    }
    for (t in 2:n_2050){
      pred_eih1[i,t-1]=exp(pred_eih1[i,t-1])
      pred_eih2[i,t-1]=exp(pred_eih2[i,t-1])
      pred_eih3[i,t-1]=exp(pred_eih3[i,t-1])
      pred_eih4[i,t-1]=exp(pred_eih4[i,t-1])
      pred_eih5[i,t-1]=exp(pred_eih5[i,t-1])
      pred_eih6[i,t-1]=exp(pred_eih6[i,t-1])
      pred_eiv1[i,t-1]=exp(pred_eiv1[i,t-1])
      pred_eiv2[i,t-1]=exp(pred_eiv2[i,t-1])
      pred_eiv3[i,t-1]=exp(pred_eiv3[i,t-1])
      pred_eiv4[i,t-1]=exp(pred_eiv4[i,t-1])
      pred_eiv5[i,t-1]=exp(pred_eiv5[i,t-1])
      pred_eig1[i,t-1]=exp(pred_eig1[i,t-1])
    }
    for (t in 2:n_2050){
      if (pred_eih1[i,t-1] >- 0){
        pred_eih1[i,t-1]<-pred_eih1[i,t-1]}
      else {pred_eih1[i,t-1]<-10^(-50)}
      if (pred_eih2[i,t-1] >- 0){
        pred_eih2[i,t-1]<-pred_eih2[i,t-1]}
      else {pred_eih2[i,t-1]<-10^(-50)}
      if (pred_eih3[i,t-1] >- 0){
        pred_eih3[i,t-1]<-pred_eih3[i,t-1]}
      else {pred_eih3[i,t-1]<-10^(-50)}
      if (pred_eih4[i,t-1] >- 0){
        pred_eih4[i,t-1]<-pred_eih4[i,t-1]}
      else {pred_eih4[i,t-1]<-10^(-50)}
      if (pred_eih5[i,t-1] >- 0){
        pred_eih5[i,t-1]=pred_eih5[i,t-1]}
      else {pred_eih5[i,t-1]<-10^(-50)}
      if (pred_eih6[i,t-1] >- 0){
        pred_eih6[i,t-1]<-pred_eih6[i,t-1]}
      else {pred_eih6[i,t-1]<-10^(-50)}
      if (pred_eiv1[i,t-1] >- 0){
        pred_eiv1[i,t-1]<-pred_eiv1[i,t-1]}
      else {pred_eiv1[i,t-1]<-10^(-50)}
      if (pred_eiv2[i,t-1] >- 0){
        pred_eiv2[i,t-1]<-pred_eiv2[i,t-1]}
      else {pred_eiv2[i,t-1]<-10^(-50)}
      if (pred_eiv3[i,t-1] >- 0){
        pred_eiv3[i,t-1]<-pred_eiv3[i,t-1]}
      else {pred_eiv3[i,t-1]<-10^(-50)}
      if (pred_eiv4[i,t-1] >- 0){
        pred_eiv4[i,t-1]<-pred_eiv4[i,t-1]}
      else {pred_eiv4[i,t-1]<-10^(-50)}
      if (pred_eiv5[i,t-1] >- 0){
        pred_eiv5[i,t-1]=pred_eiv5[i,t-1]}
      else {pred_eiv5[i,t-1]<-10^(-50)}
      if (pred_eig1[i,t-1] >- 0){
        pred_eig1[i,t-1]<-pred_eig1[i,t-1]}
      else {pred_eig1[i,t-1]<-10^(-50)}
    }
  }
  
  pred_sr1 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s1 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr2 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s2 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr3 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s3 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr4 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s4 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr5 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s5 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr6 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s6 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 1:n_2050){
      pred_sr1[i,t] <- lm1_mcmc_ref[i,36+1]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+7]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+13])
      pred_sr2[i,t] <- lm1_mcmc_ref[i,36+2]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+8]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+14])
      pred_sr3[i,t] <- lm1_mcmc_ref[i,36+3]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+9]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+15])
      pred_sr4[i,t] <- lm1_mcmc_ref[i,36+4]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+10]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+16])
      pred_sr5[i,t] <- lm1_mcmc_ref[i,36+5]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+11]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+17])
      pred_sr6[i,t] <- lm1_mcmc_ref[i,36+6]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+12]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+18])
    }
    for (t in 2:21){
      pred_s1[i,t-1] <- s[t,j]
      pred_s2[i,t-1] <- s[t+21*1,j]
      pred_s3[i,t-1] <- s[t+21*2,j]
      pred_s4[i,t-1] <- s[t+21*3,j]
      pred_s5[i,t-1] <- s[t+21*4,j]
      pred_s6[i,t-1] <- s[t+21*5,j]
    }
    for (t in 22:n_2050){
      pred_s1[i,t-1] <-  pred_sr1[i,t]+(pred_s1[i,t-2]-pred_sr1[i,t-1])*lm1_mcmc_s[i,j]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240])
      pred_s2[i,t-1] <-  pred_sr2[i,t]+(pred_s2[i,t-2]-pred_sr2[i,t-1])*lm1_mcmc_s[i,j+40*1]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*1])
      pred_s3[i,t-1] <-  pred_sr3[i,t]+(pred_s3[i,t-2]-pred_sr3[i,t-1])*lm1_mcmc_s[i,j+40*2]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*2])
      pred_s4[i,t-1] <-  pred_sr4[i,t]+(pred_s4[i,t-2]-pred_sr4[i,t-1])*lm1_mcmc_s[i,j+40*3]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*3])
      pred_s5[i,t-1] <-  pred_sr5[i,t]+(pred_s5[i,t-2]-pred_sr5[i,t-1])*lm1_mcmc_s[i,j+40*4]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*4])
      pred_s6[i,t-1] <-  pred_sr6[i,t]+(pred_s6[i,t-2]-pred_sr6[i,t-1])*lm1_mcmc_s[i,j+40*5]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*5])
    }
    for (t in 2:n_2050){
      pred_s1[i,t-1]=exp(pred_s1[i,t-1])
      pred_s2[i,t-1]=exp(pred_s2[i,t-1])
      pred_s3[i,t-1]=exp(pred_s3[i,t-1])
      pred_s4[i,t-1]=exp(pred_s4[i,t-1])
      pred_s5[i,t-1]=exp(pred_s5[i,t-1])
      pred_s6[i,t-1]=exp(pred_s6[i,t-1])
    }
  }
  pred_norm_s1<-pred_s1/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s2<-pred_s2/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s3<-pred_s3/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s4<-pred_s4/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s5<-pred_s5/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s6<-pred_s6/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s1[is.na(pred_norm_s1)] <- 0
  pred_norm_s2[is.na(pred_norm_s2)] <- 0
  pred_norm_s3[is.na(pred_norm_s3)] <- 0
  pred_norm_s4[is.na(pred_norm_s4)] <- 0
  pred_norm_s5[is.na(pred_norm_s5)] <- 0
  pred_norm_s6[is.na(pred_norm_s6)] <- 0
  for (i in 1:nrow(lm1_mcmc_ref)){
    for (t in 2:n_2050){
      if (pred_norm_s1[i,t-1] > 0){
        pred_norm_s1[i,t-1]<-pred_norm_s1[i,t-1]}
      else {pred_norm_s1[i,t-1]<-10^(-50)}
      if (pred_norm_s2[i,t-1] > 0){
        pred_norm_s2[i,t-1]<-pred_norm_s2[i,t-1]}
      else {pred_norm_s2[i,t-1]<-10^(-50)}
      if (pred_norm_s3[i,t-1] >- 0){
        pred_norm_s3[i,t-1]=pred_norm_s3[i,t-1]}
      else {pred_norm_s3[i,t-1]<-10^(-50)}
      if (pred_norm_s4[i,t-1] > 0){
        pred_norm_s4[i,t-1]<-pred_norm_s4[i,t-1]}
      else {pred_norm_s4[i,t-1]<-10^(-50)}
      if (pred_norm_s5[i,t-1] > 0){
        pred_norm_s5[i,t-1]<-pred_norm_s5[i,t-1]}
      else {pred_norm_s5[i,t-1]<-10^(-50)}
      if (pred_norm_s6[i,t-1] > 0){
        pred_norm_s6[i,t-1]<-pred_norm_s6[i,t-1]}
      else {pred_norm_s6[i,t-1]<-10^(-50)}
    }
  } 
  
  pred_vr1 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v1 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr2 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v2 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr3 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v3 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr4 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v4 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr5 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v5 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 1:n_2050){
      pred_vr1[i,t] <- lm1_mcmc_ref[i,54+1]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+6]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+11])
      pred_vr2[i,t] <- lm1_mcmc_ref[i,54+2]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+7]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+12])
      pred_vr3[i,t] <- lm1_mcmc_ref[i,54+3]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+8]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+13])
      pred_vr4[i,t] <- lm1_mcmc_ref[i,54+4]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+9]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+14])
      pred_vr5[i,t] <- lm1_mcmc_ref[i,54+5]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+10]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+15])
    }
    for (t in 2:21){
      pred_v1[i,t-1] <- v[t,j]
      pred_v2[i,t-1] <- v[t+21*1,j]
      pred_v3[i,t-1] <- v[t+21*2,j]
      pred_v4[i,t-1] <- v[t+21*3,j]
      pred_v5[i,t-1] <- v[t+21*4,j]
    }
    for (t in 22:n_2050){
      pred_v1[i,t-1] <-  pred_vr1[i,t]+(pred_v1[i,t-2]-pred_vr1[i,t-1])*lm1_mcmc_v[i,j]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200])
      pred_v2[i,t-1] <-  pred_vr2[i,t]+(pred_v2[i,t-2]-pred_vr2[i,t-1])*lm1_mcmc_v[i,j+40*1]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*1])
      pred_v3[i,t-1] <-  pred_vr3[i,t]+(pred_v3[i,t-2]-pred_vr3[i,t-1])*lm1_mcmc_v[i,j+40*2]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*2])
      pred_v4[i,t-1] <-  pred_vr4[i,t]+(pred_v4[i,t-2]-pred_vr4[i,t-1])*lm1_mcmc_v[i,j+40*3]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*3])
      pred_v5[i,t-1] <-  pred_vr5[i,t]+(pred_v5[i,t-2]-pred_vr5[i,t-1])*lm1_mcmc_v[i,j+40*4]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*4])
    }
    for (t in 2:n_2050){
      pred_v1[i,t-1]=exp(pred_v1[i,t-1])
      pred_v2[i,t-1]=exp(pred_v2[i,t-1])
      pred_v3[i,t-1]=exp(pred_v3[i,t-1])
      pred_v4[i,t-1]=exp(pred_v4[i,t-1])
      pred_v5[i,t-1]=exp(pred_v5[i,t-1])
    }
  }
  pred_norm_v1<-pred_v1/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v2<-pred_v2/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v3<-pred_v3/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v4<-pred_v4/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v5<-pred_v5/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v1[is.na(pred_norm_v1)] <- 0
  pred_norm_v2[is.na(pred_norm_v2)] <- 0
  pred_norm_v3[is.na(pred_norm_v3)] <- 0
  pred_norm_v4[is.na(pred_norm_v4)] <- 0
  pred_norm_v5[is.na(pred_norm_v5)] <- 0
  for (i in 1:nrow(lm1_mcmc_ref)){
    for (t in 2:n_2050){
      if (pred_norm_v1[i,t-1] >- 0){
        pred_norm_v1[i,t-1]<-pred_norm_v1[i,t-1]}
      else {pred_norm_v1[i,t-1]<-10^(-50)}
      if (pred_norm_v2[i,t-1] >- 0){
        pred_norm_v2[i,t-1]<-pred_norm_v2[i,t-1]}
      else {pred_norm_v2[i,t-1]<-10^(-50)}
      if (pred_norm_v3[i,t-1] >- 0){
        pred_norm_v3[i,t-1]<-pred_norm_v3[i,t-1]}
      else {pred_norm_v3[i,t-1]<-10^(-50)}
      if (pred_norm_v4[i,t-1] >- 0){
        pred_norm_v4[i,t-1]<-pred_norm_v4[i,t-1]}
      else {pred_norm_v4[i,t-1]<-10^(-50)}
      if (pred_norm_v5[i,t-1] >- 0){
        pred_norm_v5[i,t-1]<-pred_norm_v5[i,t-1]}
      else {pred_norm_v5[i,t-1]<-10^(-50)}
    }
  }  
  
  pred_yr1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050)
  pred_y1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_yr2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050)
  pred_y2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_yr3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050)
  pred_y3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 1:n_2050){
      pred_yr1[i,t] <- lm1_mcmc_ref[i,69+1]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,69+4]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,69+7])
      pred_yr2[i,t] <- lm1_mcmc_ref[i,69+2]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,69+5]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,69+8])
      pred_yr3[i,t] <- lm1_mcmc_ref[i,69+3]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,69+6]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,69+9])
    }
    for (t in 2:29){
      pred_y1[i,t-1] <- yo[t,j]
      pred_y2[i,t-1] <- yo[t+29*1,j]
      pred_y3[i,t-1] <- yo[t+29*2,j]
    }
    for (t in 30:n_2050){
      pred_y1[i,t-1] <-  pred_yr1[i,t]+(pred_y1[i,t-2]-pred_yr1[i,t-1])*lm1_mcmc_y[i,j]+rnorm(1,mean=0,sd=lm1_mcmc_y[i,j+120])
      pred_y2[i,t-1] <-  pred_yr2[i,t]+(pred_y2[i,t-2]-pred_yr2[i,t-1])*lm1_mcmc_y[i,j+40*1]+rnorm(1,mean=0,sd=lm1_mcmc_y[i,j+120+40*1])
      pred_y3[i,t-1] <-  pred_yr3[i,t]+(pred_y3[i,t-2]-pred_yr3[i,t-1])*lm1_mcmc_y[i,j+40*2]+rnorm(1,mean=0,sd=lm1_mcmc_y[i,j+120+40*2])
    }
  }
  pred_norm_y1<-pred_y1/(pred_y1+pred_y2+pred_y3)
  pred_norm_y2<-pred_y2/(pred_y1+pred_y2+pred_y3)
  pred_norm_y3<-pred_y3/(pred_y1+pred_y2+pred_y3)
  pred_norm_y1[is.na(pred_norm_y1)] <- 0
  pred_norm_y2[is.na(pred_norm_y2)] <- 0
  pred_norm_y3[is.na(pred_norm_y3)] <- 0
  for (i in 1:nrow(lm1_mcmc_ref)){
    for (t in 2:n_2050){
      if (pred_norm_y1[i,t-1] >- 0){
        pred_norm_y1[i,t-1]<-pred_norm_y1[i,t-1]}
      else {pred_norm_y1[i,t-1]<-10^(-50)}
      if (pred_norm_y2[i,t-1] >- 0){
        pred_norm_y2[i,t-1]<-pred_norm_y2[i,t-1]}
      else {pred_norm_y2[i,t-1]<-10^(-50)}
      if (pred_norm_y3[i,t-1] >- 0){
        pred_norm_y3[i,t-1]<-pred_norm_y3[i,t-1]}
      else {pred_norm_y3[i,t-1]<-10^(-50)}
    }
  }   
  
  pred_e_s1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s4 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s5 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s6 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v4 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v5 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  
  
  for (i in 1:nrow(lm1_mcmc_fh)){
    for (t in 2:n_2050){
      pred_e_s1[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s1[i,t-1]*pred_eih1[i,t-1]/1000
      pred_e_s2[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s2[i,t-1]*pred_eih2[i,t-1]/1000
      pred_e_s3[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s3[i,t-1]*pred_eih3[i,t-1]/1000
      pred_e_s4[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s4[i,t-1]*pred_eih4[i,t-1]/1000
      pred_e_s5[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s5[i,t-1]*pred_eih5[i,t-1]/1000
      pred_e_s6[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s6[i,t-1]*pred_eih6[i,t-1]/1000
      pred_e_v1[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v1[i,t-1]*pred_eiv1[i,t-1]/1000
      pred_e_v2[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v2[i,t-1]*pred_eiv2[i,t-1]/1000
      pred_e_v3[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v3[i,t-1]*pred_eiv3[i,t-1]/1000
      pred_e_v4[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v4[i,t-1]*pred_eiv4[i,t-1]/1000
      pred_e_v5[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v5[i,t-1]*pred_eiv5[i,t-1]/1000
      pred_e_g1[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y3[i,t-1]*pred_eig1[i,t-1]/1000
    }
  }
  pred_e1<- pred_e_s1
  pred_e2<- pred_e_s2
  pred_e3<- pred_e_s3+pred_e_v1
  pred_e4<- pred_e_s4+pred_e_v2
  pred_e5<- pred_e_s5+pred_e_v3
  pred_e6<- pred_e_s6+pred_e_v4+pred_e_g1
  pred_e7<- pred_e_v5
  
  pred_e_s<- pred_e_s1+pred_e_s2+pred_e_s3+pred_e_s4+pred_e_s5+pred_e_s6
  pred_e_v<- pred_e_v1+pred_e_v2+pred_e_v3+pred_e_v4+pred_e_v5
  pred_e_g<- pred_e_g1
  pred_e<-pred_e1+pred_e2+pred_e3+pred_e4+pred_e5+pred_e6+pred_e7
  
  for (i in 1:nrow(lm1_mcmc_fh)){
    for (t in 1:n_2050-1){
      pred_e1_euro[i,t]<-pred_e1_euro[i,t]+pred_e1[i,t]
      pred_e2_euro[i,t]<-pred_e2_euro[i,t]+pred_e2[i,t]
      pred_e3_euro[i,t]<-pred_e3_euro[i,t]+pred_e3[i,t]
      pred_e4_euro[i,t]<-pred_e4_euro[i,t]+pred_e4[i,t]
      pred_e5_euro[i,t]<-pred_e5_euro[i,t]+pred_e5[i,t]
      pred_e6_euro[i,t]<-pred_e6_euro[i,t]+pred_e6[i,t]
      pred_e7_euro[i,t]<-pred_e7_euro[i,t]+pred_e7[i,t]
      pred_e_s_euro[i,t]<-pred_e_s_euro[i,t]+pred_e_s[i,t]
      pred_e_v_euro[i,t]<-pred_e_v_euro[i,t]+pred_e_v[i,t]
      pred_e_g_euro[i,t]<-pred_e_g_euro[i,t]+pred_e_g[i,t]
      pred_e_euro[i,t]<-pred_e_euro[i,t]+pred_e[i,t]
    }
  }
  
  fac_gdp<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_p<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_r<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_y<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_s<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_f<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_f1<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_f2<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_f3<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_f4<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_f5<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_f6<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_f7<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_s1<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_s2<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_s3<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_s4<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_s5<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_s6<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  fac_s7<- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=1)
  L <- matrix(NA, nrow = nrow(lm1_mcmc_y),ncol=12)
  
  a1=55
  for (i in 1:nrow(lm1_mcmc_y)){
    L[i,1]<-(pred_e_s1[i,a1]-his_e[1,j])/(log(pred_e_s1[i,a1]/his_e[1,j]))
    L[i,2]<-(pred_e_s2[i,a1]-his_e[2,j])/(log(pred_e_s2[i,a1]/his_e[2,j]))
    L[i,3]<-(pred_e_s3[i,a1]-his_e[3,j])/(log(pred_e_s3[i,a1]/his_e[3,j]))
    L[i,4]<-(pred_e_s4[i,a1]-his_e[4,j])/(log(pred_e_s4[i,a1]/his_e[4,j]))
    L[i,5]<-(pred_e_s5[i,a1]-his_e[5,j])/(log(pred_e_s5[i,a1]/his_e[5,j]))
    L[i,6]<-(pred_e_s6[i,a1]-his_e[6,j])/(log(pred_e_s6[i,a1]/his_e[6,j]))
    
    L[i,7]<-(pred_e_v1[i,a1]-his_e[7,j])/(log(pred_e_v1[i,a1]/his_e[7,j]))
    L[i,8]<-(pred_e_v2[i,a1]-his_e[8,j])/(log(pred_e_v2[i,a1]/his_e[8,j]))
    L[i,9]<-(pred_e_v3[i,a1]-his_e[9,j])/(log(pred_e_v3[i,a1]/his_e[9,j]))
    L[i,10]<-(pred_e_v4[i,a1]-his_e[10,j])/(log(pred_e_v4[i,a1]/his_e[10,j]))
    L[i,11]<-(pred_e_v5[i,a1]-his_e[11,j])/(log(pred_e_v5[i,a1]/his_e[11,j]))
    
    L[i,12]<-(pred_e_g1[i,a1]-his_e[12,j])/(log(pred_e_g1[i,a1]/his_e[12,j]))
  }
  
  for (i in 1:nrow(lm1_mcmc_y)){
    fac_f1[i,]=L[i,1]*log(pred_eih1[i,a1]/fh[21,j])
    fac_f2[i,]=L[i,2]*log(pred_eih2[i,a1]/fh[42,j])
    fac_f3[i,]=L[i,3]*log(pred_eih3[i,a1]/fh[63,j])+L[i,7]*log(pred_eiv1[i,a1]/fv[21,j])
    fac_f4[i,]=L[i,4]*log(pred_eih4[i,a1]/fh[84,j])+L[i,8]*log(pred_eiv2[i,a1]/fv[42,j])
    fac_f5[i,]=L[i,5]*log(pred_eih5[i,a1]/fh[105,j])+L[i,9]*log(pred_eiv3[i,a1]/fv[63,j])
    fac_f6[i,]=L[i,6]*log(pred_eih6[i,a1]/fh[126,j])+L[i,10]*log(pred_eiv4[i,a1]/fv[84,j])+L[i,12]*log(pred_eig1[i,a1]/fg[21,j])
    fac_f7[i,]=L[i,11]*log(pred_eiv5[i,a1]/fv[105,j])
    
    fac_s1[i,]=L[i,1]*(log(pred_norm_s1[i,a1]/so[21,j]))
    fac_s2[i,]=L[i,2]*(log(pred_norm_s2[i,a1]/so[42,j]))
    fac_s3[i,]=L[i,3]*(log(pred_norm_s3[i,a1]/so[63,j]))+L[i,7]*(log(pred_norm_v1[i,a1]/vo[21,j]))
    fac_s4[i,]=L[i,4]*(log(pred_norm_s4[i,a1]/so[84,j]))+L[i,8]*(log(pred_norm_v2[i,a1]/vo[42,j]))
    fac_s5[i,]=L[i,5]*(log(pred_norm_s5[i,a1]/so[105,j]))+L[i,9]*(log(pred_norm_v3[i,a1]/vo[63,j]))
    fac_s6[i,]=L[i,6]*(log(pred_norm_s6[i,a1]/so[126,j]))+L[i,10]*(log(pred_norm_v4[i,a1]/vo[84,j]))
    fac_s7[i,]=L[i,11]*(log(pred_norm_v5[i,a1]/vo[105,j]))
    
    fac_f[i,]=fac_f1[i,]+fac_f2[i,]+fac_f3[i,]+fac_f4[i,]+fac_f5[i,]+fac_f6[i,]+fac_f7[i,]
    
    fac_s[i,]=fac_s1[i,]+fac_s2[i,]+fac_s3[i,]+fac_s4[i,]+fac_s5[i,]+fac_s6[i,]+fac_s7[i,]
    
    fac_y[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6])*(log(pred_norm_y1[i,a1]/yo[21,j]))+(L[i,7]++L[i,8]+L[i,9]+L[i,10]+L[i,11])*(log(pred_norm_y2[i,a1]/yo[50,j]))+(L[i,12])*(log(pred_norm_y3[i,a1]/yo[79,j]))
    
    fac_r[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6]+L[i,7]+L[i,8]+L[i,9]+L[i,10]+L[i,11]+L[i,12])*(log(pred_r[i,a1+1]/ro[21,j]))
    
    fac_gdp[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6]+L[i,7]+L[i,8]+L[i,9]+L[i,10]+L[i,11]+L[i,12])*(log(pred_gdp[i,a1+1]/gdp[21,j]))
    
    fac_p[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6]+L[i,7]+L[i,8]+L[i,9]+L[i,10]+L[i,11]+L[i,12])*(log(pred_p[i,a1+1]/p[21,j]))
  }
  
  for (i in 1: nrow(fac_f)){
    fac_f1_reg[i,1]<-fac_f1_reg[i,1]+fac_f1[i,]
    fac_f2_reg[i,1]<-fac_f2_reg[i,1]+fac_f2[i,]
    fac_f3_reg[i,1]<-fac_f3_reg[i,1]+fac_f3[i,]
    fac_f4_reg[i,1]<-fac_f4_reg[i,1]+fac_f4[i,]
    fac_f5_reg[i,1]<-fac_f5_reg[i,1]+fac_f5[i,]
    fac_f6_reg[i,1]<-fac_f6_reg[i,1]+fac_f6[i,]
    fac_f7_reg[i,1]<-fac_f7_reg[i,1]+fac_f7[i,]
    fac_s1_reg[i,1]<-fac_s1_reg[i,1]+fac_s1[i,]
    fac_s2_reg[i,1]<-fac_s2_reg[i,1]+fac_s2[i,]
    fac_s3_reg[i,1]<-fac_s3_reg[i,1]+fac_s3[i,]
    fac_s4_reg[i,1]<-fac_s4_reg[i,1]+fac_s4[i,]
    fac_s5_reg[i,1]<-fac_s5_reg[i,1]+fac_s5[i,]
    fac_s6_reg[i,1]<-fac_s6_reg[i,1]+fac_s6[i,]
    fac_s7_reg[i,1]<-fac_s7_reg[i,1]+fac_s7[i,]
    fac_f_reg[i,1]<-fac_f_reg[i,1]+fac_f[i,]
    fac_s_reg[i,1]<-fac_s_reg[i,1]+fac_s[i,]
    fac_y_reg[i,1]<-fac_y_reg[i,1]+fac_y[i,]
    fac_gdp_reg[i,1]<-fac_gdp_reg[i,1]+fac_gdp[i,]
    fac_p_reg[i,1]<-fac_p_reg[i,1]+fac_p[i,]
    fac_r_reg[i,1]<-fac_r_reg[i,1]+fac_r[i,]
  }
}
#changes of emission intensity
for (i in 1:nrow(lm1_mcmc_y)){
  for (t in 1:n_2050-1){
    pred_ei1_euro[i,t]<-pred_e1_euro[i,t]/pred_o1_euro[i,t]*1000
    pred_ei2_euro[i,t]<-pred_e2_euro[i,t]/pred_o2_euro[i,t]*1000
    pred_ei3_euro[i,t]<-pred_e3_euro[i,t]/pred_o3_euro[i,t]*1000
    pred_ei4_euro[i,t]<-pred_e4_euro[i,t]/pred_o4_euro[i,t]*1000
    pred_ei5_euro[i,t]<-pred_e5_euro[i,t]/pred_o5_euro[i,t]*1000
    pred_ei6_euro[i,t]<-pred_e6_euro[i,t]/pred_o6_euro[i,t]*1000
    pred_ei7_euro[i,t]<-pred_e7_euro[i,t]/pred_o7_euro[i,t]*1000
    pred_ei_euro[i,t]<- pred_e7_euro[i,t]/(pred_o1_euro[i,t]+pred_o2_euro[i,t]+pred_o3_euro[i,t]+pred_o4_euro[i,t]+pred_o5_euro[i,t]+pred_o6_euro[i,t]+pred_o7_euro[i,t])*1000
     }
}
cha_ei1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = 18)
cha_ei2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = 18)
cha_ei3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = 18)
cha_ei4 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = 18)
cha_ei5 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = 18)
cha_ei6 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = 18)
cha_ei7 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = 18)
for (i in 1:nrow(lm1_mcmc_y)){
  cha_ei1[i,1]=(pred_ei1_euro[i,55]-his_ei[1,1])/his_ei[1,1]*100
  cha_ei2[i,1]=(pred_ei2_euro[i,55]-his_ei[2,1])/his_ei[2,1]*100
  cha_ei3[i,1]=(pred_ei3_euro[i,55]-his_ei[3,1])/his_ei[3,1]*100
  cha_ei4[i,1]=(pred_ei4_euro[i,55]-his_ei[4,1])/his_ei[4,1]*100
  cha_ei5[i,1]=(pred_ei5_euro[i,55]-his_ei[5,1])/his_ei[5,1]*100
  cha_ei6[i,1]=(pred_ei6_euro[i,55]-his_ei[6,1])/his_ei[6,1]*100
  cha_ei7[i,1]=(pred_ei7_euro[i,55]-his_ei[7,1])/his_ei[7,1]*100
}
#predict intervals for EU emissions
e1_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e1_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e1_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e2_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e2_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e2_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e3_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e3_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e3_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e4_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e4_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e4_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e5_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e5_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e5_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e6_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e6_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e6_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e7_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e7_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e7_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_s_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_s_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_s_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_v_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_v_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_v_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_g_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_g_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_g_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
e_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))

ei1_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei1_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei1_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei2_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei2_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei2_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei3_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei3_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei3_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei4_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei4_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei4_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei5_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei5_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei5_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei6_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei6_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei6_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei7_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei7_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei7_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei_upper= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei_mean= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))
ei_lower= as.data.frame(matrix(NA, nrow =n_2050-1,ncol=17))

uncertain_mean_reg<- matrix(NA, nrow =n_2050-1,ncol=1)
for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e_euro[,t])}  
uncertain_lower_reg <- apply(pred_e_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_e_euro, MARGIN = 2, quantile, prob = 0.975)
e_upper[,1]<- uncertain_upper_reg
e_lower[,1]<- uncertain_lower_reg
e_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e1_euro[,t])}  
uncertain_lower_reg <- apply(pred_e1_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_e1_euro, MARGIN = 2, quantile, prob = 0.975)
e1_upper[,1]<- uncertain_upper_reg
e1_lower[,1]<- uncertain_lower_reg
e1_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e2_euro[,t])}  
uncertain_lower_reg <- apply(pred_e2_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_e2_euro, MARGIN = 2, quantile, prob = 0.975)
e2_upper[,1]<- uncertain_upper_reg
e2_lower[,1]<- uncertain_lower_reg
e2_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e3_euro[,t])}  
uncertain_lower_reg <- apply(pred_e3_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_e3_euro, MARGIN = 2, quantile, prob = 0.975)
e3_upper[,1]<- uncertain_upper_reg
e3_lower[,1]<- uncertain_lower_reg
e3_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e4_euro[,t])}  
uncertain_lower_reg <- apply(pred_e4_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_e4_euro, MARGIN = 2, quantile, prob = 0.975)
e4_upper[,1]<- uncertain_upper_reg
e4_lower[,1]<- uncertain_lower_reg
e4_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e5_euro[,t])}  
uncertain_lower_reg <- apply(pred_e5_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_e5_euro, MARGIN = 2, quantile, prob = 0.975)
e5_upper[,1]<- uncertain_upper_reg
e5_lower[,1]<- uncertain_lower_reg
e5_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e6_euro[,t])}  
uncertain_lower_reg <- apply(pred_e6_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_e6_euro, MARGIN = 2, quantile, prob = 0.975)
e6_upper[,1]<- uncertain_upper_reg
e6_lower[,1]<- uncertain_lower_reg
e6_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e7_euro[,t])}  
uncertain_lower_reg <- apply(pred_e7_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_e7_euro, MARGIN = 2, quantile, prob = 0.975)
e7_upper[,1]<- uncertain_upper_reg
e7_lower[,1]<- uncertain_lower_reg
e7_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e_s_euro[,t])}  
uncertain_lower_reg <- apply(pred_e_s_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_e_s_euro, MARGIN = 2, quantile, prob = 0.975)
e_s_upper[,1]<- uncertain_upper_reg
e_s_lower[,1]<- uncertain_lower_reg
e_s_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e_v_euro[,t])}  
uncertain_lower_reg <- apply(pred_e_v_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_e_v_euro, MARGIN = 2, quantile, prob = 0.975)
e_v_upper[,1]<- uncertain_upper_reg
e_v_lower[,1]<- uncertain_lower_reg
e_v_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e_g_euro[,t])}  
uncertain_lower_reg <- apply(pred_e_g_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_e_g_euro, MARGIN = 2, quantile, prob = 0.975)
e_g_upper[,1]<- uncertain_upper_reg
e_g_lower[,1]<- uncertain_lower_reg
e_g_mean[,1]<- uncertain_mean_reg


for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei1_euro[,t])}  
uncertain_lower_reg <- apply(pred_ei1_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei1_euro, MARGIN = 2, quantile, prob = 0.975)
ei1_upper[,1]<- uncertain_upper_reg
ei1_lower[,1]<- uncertain_lower_reg
ei1_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei2_euro[,t])}  
uncertain_lower_reg <- apply(pred_ei2_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei2_euro, MARGIN = 2, quantile, prob = 0.975)
ei2_upper[,1]<- uncertain_upper_reg
ei2_lower[,1]<- uncertain_lower_reg
ei2_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei3_euro[,t])}  
uncertain_lower_reg <- apply(pred_ei3_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei3_euro, MARGIN = 2, quantile, prob = 0.975)
ei3_upper[,1]<- uncertain_upper_reg
ei3_lower[,1]<- uncertain_lower_reg
ei3_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei4_euro[,t])}  
uncertain_lower_reg <- apply(pred_ei4_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei4_euro, MARGIN = 2, quantile, prob = 0.975)
ei4_upper[,1]<- uncertain_upper_reg
ei4_lower[,1]<- uncertain_lower_reg
ei4_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei5_euro[,t])}  
uncertain_lower_reg <- apply(pred_ei5_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei5_euro, MARGIN = 2, quantile, prob = 0.975)
ei5_upper[,1]<- uncertain_upper_reg
ei5_lower[,1]<- uncertain_lower_reg
ei5_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei6_euro[,t])}  
uncertain_lower_reg <- apply(pred_ei6_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei6_euro, MARGIN = 2, quantile, prob = 0.975)
ei6_upper[,1]<- uncertain_upper_reg
ei6_lower[,1]<- uncertain_lower_reg
ei6_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei7_euro[,t])}  
uncertain_lower_reg <- apply(pred_ei7_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei7_euro, MARGIN = 2, quantile, prob = 0.975)
ei7_upper[,1]<- uncertain_upper_reg
ei7_lower[,1]<- uncertain_lower_reg
ei7_mean[,1]<- uncertain_mean_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei_euro[,t])}  
uncertain_lower_reg <- apply(pred_ei_euro, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei_euro, MARGIN = 2, quantile, prob = 0.975)
ei_upper[,1]<- uncertain_upper_reg
ei_lower[,1]<- uncertain_lower_reg
ei_mean[,1]<- uncertain_mean_reg

#projection of emissions for other 37 countries
pred_e_global<- pred_e_euro
pred_e1_global<- pred_e1_euro
pred_e2_global<- pred_e2_euro
pred_e3_global<- pred_e3_euro
pred_e4_global<- pred_e4_euro
pred_e5_global<- pred_e5_euro
pred_e6_global<- pred_e6_euro
pred_e7_global<- pred_e7_euro
pred_e_s_global<- pred_e_s_euro
pred_e_v_global<- pred_e_v_euro
pred_e_g_global<- pred_e_g_euro

for (j in 25:37){
  pred_p <- matrix(NA, nrow = nrow(lm1_mcmc_p),ncol=n_2050)
  for (i in 1:nrow(lm1_mcmc_p)){
    for (t in 1:26){pred_p[i,t] <- exp(x1[t,j])}
    for (t in 27:56){
      if (pred_p[i,t-1] > 0){pred_p[i,t]<-exp(rnorm(1,mean=0,sd=sd(lm1_mcmc_p[,j]))+log(pred_p[i,t-1])+d1[t-1,j])+rnorm(1,mean=0,sd=lm1_mcmc_p[i,j+40])}
      else {pred_p[i,t]<-0}}
  }
  for (i in 1:nrow(lm1_mcmc_gdp_pre)){
    for (t in 1:n_2050){
      if (pred_p[i,t] > 0){pred_p[i,t] <- pred_p[i,t]}
      else {pred_p[i,t]<-10^(-50)}
    }
  }
  
  
  pred_gdp <- matrix(NA, nrow = nrow(lm1_mcmc_gdp_pre),ncol=n_2050)
  for (i in 1:nrow(lm1_mcmc_gdp_pre)){
    for (t in 1:24){pred_gdp[i,t] <- exp(x[t,j])}
    for (t in 25:25){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*1]}
      else {pred_gdp[i,t]<-0}}
    for (t in 26:26){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*2]}
      else {pred_gdp[i,t]<-0}}
    for (t in 27:27){
      if (pred_gdp[i,t-1] >- 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*3]}
      else {pred_gdp[i,t]<-0}}
    for (t in 28:28){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*4]}
      else {pred_gdp[i,t]<-0}}
    for (t in 29:29){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*5]}
      else {pred_gdp[i,t]<-0}}
    for (t in 30:30){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*6]}
      else {pred_gdp[i,t]<-0}}
    for (t in 31:31){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*7]}
      else {pred_gdp[i,t]<-0}}
    for (t in 32:32){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*8]}
      else {pred_gdp[i,t]<-0}}
    for (t in 33:33){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*9]}
      else {pred_gdp[i,t]<-0}}
    for (t in 34:34){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*10]}
      else {pred_gdp[i,t]<-0}}
    for (t in 35:35){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*11]}
      else {pred_gdp[i,t]<-0}}
    for (t in 36:36){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*12]}
      else {pred_gdp[i,t]<-0}}
    for (t in 37:37){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*13]}
      else {pred_gdp[i,t]<-0}}
    for (t in 38:38){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*14]}
      else {pred_gdp[i,t]<-0}}
    for (t in 39:39){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*15]}
      else {pred_gdp[i,t]<-0}}
    for (t in 40:40){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*16]}
      else {pred_gdp[i,t]<-0}}
    for (t in 41:41){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*17]}
      else {pred_gdp[i,t]<-0}}
    for (t in 42:42){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*18]}
      else {pred_gdp[i,t]<-0}}
    for (t in 43:43){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*19]}
      else {pred_gdp[i,t]<-0}}
    for (t in 44:44){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*20]}
      else {pred_gdp[i,t]<-0}}
    for (t in 45:45){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*21]}
      else {pred_gdp[i,t]<-0}}
    for (t in 46:46){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*22]}
      else {pred_gdp[i,t]<-0}}
    for (t in 47:47){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*23]}
      else {pred_gdp[i,t]<-0}}
    for (t in 48:48){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*24]}
      else {pred_gdp[i,t]<-0}}
    for (t in 49:49){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*25]}
      else {pred_gdp[i,t]<-0}}
    for (t in 50:50){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*26]}
      else {pred_gdp[i,t]<-0}}
    for (t in 51:51){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*27]}
      else {pred_gdp[i,t]<-0}}
    for (t in 52:52){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*28]}
      else {pred_gdp[i,t]<-0}}
    for (t in 53:53){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*29]}
      else {pred_gdp[i,t]<-0}}
    for (t in 54:54){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*30]}
      else {pred_gdp[i,t]<-0}}
    for (t in 55:55){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*31]}
      else {pred_gdp[i,t]<-0}}
    for (t in 56:56){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*32]}
      else {pred_gdp[i,t]<-0}}
  }
  for (i in 1:nrow(lm1_mcmc_gdp_pre)){
    for (t in 1:n_2050){
      if (pred_gdp[i,t] > 0){pred_gdp[i,t] <- pred_gdp[i,t]}
      else {pred_gdp[i,t]<-10^(-50)}
    }
  }
  
  
  pred_r <- matrix(NA, nrow = nrow(lm1_mcmc_p),ncol=n_2050)
  for (i in 1:nrow(lm1_mcmc_r)){
    for (t in 1:n_2015){pred_r[i,t] <- x2[t,j]}
    for (t in 22:56){
      if (pred_r[i,t-1] > 0){pred_r[i,t]<-rnorm(1,mean=0,sd=lm1_mcmc_r[i,j+80])+pred_r[i,t-1]+d2[t-1,j]+rnorm(1,mean=0,sd=lm1_mcmc_r[i,j+40])}
      else {pred_r[i,t]<-0}}
  }
  for (i in 1:nrow(lm1_mcmc_r)){
    for (t in 1:n_2050){
      if (pred_r[i,t] > 0){pred_r[i,t] <- pred_r[i,t]}
      else {pred_r[i,t]<-10^(-50)}
    }
  }
  
  pred_eih1 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih2 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih3 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih4 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih5 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih6 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv1 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv2 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv3 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv4 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv5 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eig1 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 2:21){
      pred_eih1[i,t-1] <- eih[t,j]
      pred_eih2[i,t-1] <- eih[t+21,j]
      pred_eih3[i,t-1] <- eih[t+21*2,j] 
      pred_eih4[i,t-1] <- eih[t+21*3,j]       
      pred_eih5[i,t-1] <- eih[t+21*4,j] 
      pred_eih6[i,t-1] <- eih[t+21*5,j] 
      
      pred_eiv1[i,t-1] <- eiv[t,j]
      pred_eiv2[i,t-1] <- eiv[t+21,j]
      pred_eiv3[i,t-1] <- eiv[t+21*2,j] 
      pred_eiv4[i,t-1] <- eiv[t+21*3,j]       
      pred_eiv5[i,t-1] <- eiv[t+21*4,j] 
      pred_eig1[i,t-1] <- eig[t,j]
    }
    for (t in 22:n_2050){
      pred_eih1[i,t-1] <-(pred_eihr1[i,t]+(pred_eih1[i,t-2]-pred_eihr1[i,t-1])*lm1_mcmc_fh[i,j+40*0])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240])
      pred_eih2[i,t-1] <-(pred_eihr2[i,t]+(pred_eih2[i,t-2]-pred_eihr2[i,t-1])*lm1_mcmc_fh[i,j+40*1])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*1])
      pred_eih3[i,t-1] <-(pred_eihr3[i,t]+(pred_eih3[i,t-2]-pred_eihr3[i,t-1])*lm1_mcmc_fh[i,j+40*2])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*2])
      pred_eih4[i,t-1] <-(pred_eihr4[i,t]+(pred_eih4[i,t-2]-pred_eihr4[i,t-1])*lm1_mcmc_fh[i,j+40*3])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*3])
      pred_eih5[i,t-1] <-(pred_eihr5[i,t]+(pred_eih5[i,t-2]-pred_eihr5[i,t-1])*lm1_mcmc_fh[i,j+40*4])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*4])
      pred_eih6[i,t-1] <-(pred_eihr6[i,t]+(pred_eih6[i,t-2]-pred_eihr6[i,t-1])*lm1_mcmc_fh[i,j+40*5])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*5])
      pred_eiv1[i,t-1] <-(pred_eivr1[i,t]+(pred_eiv1[i,t-2]-pred_eivr1[i,t-1])*lm1_mcmc_fv[i,j+40*0])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200])
      pred_eiv2[i,t-1] <-(pred_eivr2[i,t]+(pred_eiv2[i,t-2]-pred_eivr2[i,t-1])*lm1_mcmc_fv[i,j+40*1])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*1])
      pred_eiv3[i,t-1] <-(pred_eivr3[i,t]+(pred_eiv3[i,t-2]-pred_eivr3[i,t-1])*lm1_mcmc_fv[i,j+40*2])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*2])
      pred_eiv4[i,t-1] <-(pred_eivr4[i,t]+(pred_eiv4[i,t-2]-pred_eivr4[i,t-1])*lm1_mcmc_fv[i,j+40*3])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*3])
      pred_eiv5[i,t-1] <-(pred_eivr5[i,t]+(pred_eiv5[i,t-2]-pred_eivr5[i,t-1])*lm1_mcmc_fv[i,j+40*4])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*4])
      pred_eig1[i,t-1] <-(pred_eigr1[i,t]+(pred_eig1[i,t-2]-pred_eigr1[i,t-1])*lm1_mcmc_fg[i,j+40*0])+rnorm(1,mean=0,sd=lm1_mcmc_fg[i,j+40])
    }
    for (t in 2:n_2050){
      pred_eih1[i,t-1]=exp(pred_eih1[i,t-1])
      pred_eih2[i,t-1]=exp(pred_eih2[i,t-1])
      pred_eih3[i,t-1]=exp(pred_eih3[i,t-1])
      pred_eih4[i,t-1]=exp(pred_eih4[i,t-1])
      pred_eih5[i,t-1]=exp(pred_eih5[i,t-1])
      pred_eih6[i,t-1]=exp(pred_eih6[i,t-1])
      pred_eiv1[i,t-1]=exp(pred_eiv1[i,t-1])
      pred_eiv2[i,t-1]=exp(pred_eiv2[i,t-1])
      pred_eiv3[i,t-1]=exp(pred_eiv3[i,t-1])
      pred_eiv4[i,t-1]=exp(pred_eiv4[i,t-1])
      pred_eiv5[i,t-1]=exp(pred_eiv5[i,t-1])
      pred_eig1[i,t-1]=exp(pred_eig1[i,t-1])
    }
    for (t in 2:n_2050){
      if (pred_eih1[i,t-1] >- 0){
        pred_eih1[i,t-1]<-pred_eih1[i,t-1]}
      else {pred_eih1[i,t-1]<-10^(-50)}
      if (pred_eih2[i,t-1] >- 0){
        pred_eih2[i,t-1]<-pred_eih2[i,t-1]}
      else {pred_eih2[i,t-1]<-10^(-50)}
      if (pred_eih3[i,t-1] >- 0){
        pred_eih3[i,t-1]<-pred_eih3[i,t-1]}
      else {pred_eih3[i,t-1]<-10^(-50)}
      if (pred_eih4[i,t-1] >- 0){
        pred_eih4[i,t-1]<-pred_eih4[i,t-1]}
      else {pred_eih4[i,t-1]<-10^(-50)}
      if (pred_eih5[i,t-1] >- 0){
        pred_eih5[i,t-1]=pred_eih5[i,t-1]}
      else {pred_eih5[i,t-1]<-10^(-50)}
      if (pred_eih6[i,t-1] >- 0){
        pred_eih6[i,t-1]<-pred_eih6[i,t-1]}
      else {pred_eih6[i,t-1]<-10^(-50)}
      if (pred_eiv1[i,t-1] >- 0){
        pred_eiv1[i,t-1]<-pred_eiv1[i,t-1]}
      else {pred_eiv1[i,t-1]<-10^(-50)}
      if (pred_eiv2[i,t-1] >- 0){
        pred_eiv2[i,t-1]<-pred_eiv2[i,t-1]}
      else {pred_eiv2[i,t-1]<-10^(-50)}
      if (pred_eiv3[i,t-1] >- 0){
        pred_eiv3[i,t-1]<-pred_eiv3[i,t-1]}
      else {pred_eiv3[i,t-1]<-10^(-50)}
      if (pred_eiv4[i,t-1] >- 0){
        pred_eiv4[i,t-1]<-pred_eiv4[i,t-1]}
      else {pred_eiv4[i,t-1]<-10^(-50)}
      if (pred_eiv5[i,t-1] >- 0){
        pred_eiv5[i,t-1]=pred_eiv5[i,t-1]}
      else {pred_eiv5[i,t-1]<-10^(-50)}
      if (pred_eig1[i,t-1] >- 0){
        pred_eig1[i,t-1]<-pred_eig1[i,t-1]}
      else {pred_eig1[i,t-1]<-10^(-50)}
    }
  }
  
  pred_sr1 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s1 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr2 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s2 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr3 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s3 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr4 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s4 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr5 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s5 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr6 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s6 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 1:n_2050){
      pred_sr1[i,t] <- lm1_mcmc_ref[i,36+1]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+7]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+13])
      pred_sr2[i,t] <- lm1_mcmc_ref[i,36+2]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+8]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+14])
      pred_sr3[i,t] <- lm1_mcmc_ref[i,36+3]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+9]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+15])
      pred_sr4[i,t] <- lm1_mcmc_ref[i,36+4]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+10]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+16])
      pred_sr5[i,t] <- lm1_mcmc_ref[i,36+5]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+11]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+17])
      pred_sr6[i,t] <- lm1_mcmc_ref[i,36+6]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+12]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+18])
    }
    for (t in 2:21){
      pred_s1[i,t-1] <- s[t,j]
      pred_s2[i,t-1] <- s[t+21*1,j]
      pred_s3[i,t-1] <- s[t+21*2,j]
      pred_s4[i,t-1] <- s[t+21*3,j]
      pred_s5[i,t-1] <- s[t+21*4,j]
      pred_s6[i,t-1] <- s[t+21*5,j]
    }
    for (t in 22:n_2050){
      pred_s1[i,t-1] <-  pred_sr1[i,t]+(pred_s1[i,t-2]-pred_sr1[i,t-1])*lm1_mcmc_s[i,j]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240])
      pred_s2[i,t-1] <-  pred_sr2[i,t]+(pred_s2[i,t-2]-pred_sr2[i,t-1])*lm1_mcmc_s[i,j+40*1]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*1])
      pred_s3[i,t-1] <-  pred_sr3[i,t]+(pred_s3[i,t-2]-pred_sr3[i,t-1])*lm1_mcmc_s[i,j+40*2]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*2])
      pred_s4[i,t-1] <-  pred_sr4[i,t]+(pred_s4[i,t-2]-pred_sr4[i,t-1])*lm1_mcmc_s[i,j+40*3]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*3])
      pred_s5[i,t-1] <-  pred_sr5[i,t]+(pred_s5[i,t-2]-pred_sr5[i,t-1])*lm1_mcmc_s[i,j+40*4]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*4])
      pred_s6[i,t-1] <-  pred_sr6[i,t]+(pred_s6[i,t-2]-pred_sr6[i,t-1])*lm1_mcmc_s[i,j+40*5]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*5])
    }
    for (t in 2:n_2050){
      pred_s1[i,t-1]=exp(pred_s1[i,t-1])
      pred_s2[i,t-1]=exp(pred_s2[i,t-1])
      pred_s3[i,t-1]=exp(pred_s3[i,t-1])
      pred_s4[i,t-1]=exp(pred_s4[i,t-1])
      pred_s5[i,t-1]=exp(pred_s5[i,t-1])
      pred_s6[i,t-1]=exp(pred_s6[i,t-1])
    }
  }
  pred_norm_s1<-pred_s1/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s2<-pred_s2/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s3<-pred_s3/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s4<-pred_s4/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s5<-pred_s5/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s6<-pred_s6/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s1[is.na(pred_norm_s1)] <- 0
  pred_norm_s2[is.na(pred_norm_s2)] <- 0
  pred_norm_s3[is.na(pred_norm_s3)] <- 0
  pred_norm_s4[is.na(pred_norm_s4)] <- 0
  pred_norm_s5[is.na(pred_norm_s5)] <- 0
  pred_norm_s6[is.na(pred_norm_s6)] <- 0
  for (i in 1:nrow(lm1_mcmc_ref)){
    for (t in 2:n_2050){
      if (pred_norm_s1[i,t-1] > 0){
        pred_norm_s1[i,t-1]<-pred_norm_s1[i,t-1]}
      else {pred_norm_s1[i,t-1]<-10^(-50)}
      if (pred_norm_s2[i,t-1] > 0){
        pred_norm_s2[i,t-1]<-pred_norm_s2[i,t-1]}
      else {pred_norm_s2[i,t-1]<-10^(-50)}
      if (pred_norm_s3[i,t-1] >- 0){
        pred_norm_s3[i,t-1]=pred_norm_s3[i,t-1]}
      else {pred_norm_s3[i,t-1]<-10^(-50)}
      if (pred_norm_s4[i,t-1] > 0){
        pred_norm_s4[i,t-1]<-pred_norm_s4[i,t-1]}
      else {pred_norm_s4[i,t-1]<-10^(-50)}
      if (pred_norm_s5[i,t-1] > 0){
        pred_norm_s5[i,t-1]<-pred_norm_s5[i,t-1]}
      else {pred_norm_s5[i,t-1]<-10^(-50)}
      if (pred_norm_s6[i,t-1] > 0){
        pred_norm_s6[i,t-1]<-pred_norm_s6[i,t-1]}
      else {pred_norm_s6[i,t-1]<-10^(-50)}
    }
  } 
  
  pred_vr1 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v1 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr2 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v2 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr3 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v3 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr4 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v4 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr5 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v5 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 1:n_2050){
      pred_vr1[i,t] <- lm1_mcmc_ref[i,54+1]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+6]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+11])
      pred_vr2[i,t] <- lm1_mcmc_ref[i,54+2]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+7]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+12])
      pred_vr3[i,t] <- lm1_mcmc_ref[i,54+3]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+8]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+13])
      pred_vr4[i,t] <- lm1_mcmc_ref[i,54+4]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+9]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+14])
      pred_vr5[i,t] <- lm1_mcmc_ref[i,54+5]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+10]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+15])
    }
    for (t in 2:21){
      pred_v1[i,t-1] <- v[t,j]
      pred_v2[i,t-1] <- v[t+21*1,j]
      pred_v3[i,t-1] <- v[t+21*2,j]
      pred_v4[i,t-1] <- v[t+21*3,j]
      pred_v5[i,t-1] <- v[t+21*4,j]
    }
    for (t in 22:n_2050){
      pred_v1[i,t-1] <-  pred_vr1[i,t]+(pred_v1[i,t-2]-pred_vr1[i,t-1])*lm1_mcmc_v[i,j]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200])
      pred_v2[i,t-1] <-  pred_vr2[i,t]+(pred_v2[i,t-2]-pred_vr2[i,t-1])*lm1_mcmc_v[i,j+40*1]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*1])
      pred_v3[i,t-1] <-  pred_vr3[i,t]+(pred_v3[i,t-2]-pred_vr3[i,t-1])*lm1_mcmc_v[i,j+40*2]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*2])
      pred_v4[i,t-1] <-  pred_vr4[i,t]+(pred_v4[i,t-2]-pred_vr4[i,t-1])*lm1_mcmc_v[i,j+40*3]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*3])
      pred_v5[i,t-1] <-  pred_vr5[i,t]+(pred_v5[i,t-2]-pred_vr5[i,t-1])*lm1_mcmc_v[i,j+40*4]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*4])
    }
    for (t in 2:n_2050){
      pred_v1[i,t-1]=exp(pred_v1[i,t-1])
      pred_v2[i,t-1]=exp(pred_v2[i,t-1])
      pred_v3[i,t-1]=exp(pred_v3[i,t-1])
      pred_v4[i,t-1]=exp(pred_v4[i,t-1])
      pred_v5[i,t-1]=exp(pred_v5[i,t-1])
    }
  }
  pred_norm_v1<-pred_v1/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v2<-pred_v2/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v3<-pred_v3/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v4<-pred_v4/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v5<-pred_v5/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v1[is.na(pred_norm_v1)] <- 0
  pred_norm_v2[is.na(pred_norm_v2)] <- 0
  pred_norm_v3[is.na(pred_norm_v3)] <- 0
  pred_norm_v4[is.na(pred_norm_v4)] <- 0
  pred_norm_v5[is.na(pred_norm_v5)] <- 0
  for (i in 1:nrow(lm1_mcmc_ref)){
    for (t in 2:n_2050){
      if (pred_norm_v1[i,t-1] >- 0){
        pred_norm_v1[i,t-1]<-pred_norm_v1[i,t-1]}
      else {pred_norm_v1[i,t-1]<-10^(-50)}
      if (pred_norm_v2[i,t-1] >- 0){
        pred_norm_v2[i,t-1]<-pred_norm_v2[i,t-1]}
      else {pred_norm_v2[i,t-1]<-10^(-50)}
      if (pred_norm_v3[i,t-1] >- 0){
        pred_norm_v3[i,t-1]<-pred_norm_v3[i,t-1]}
      else {pred_norm_v3[i,t-1]<-10^(-50)}
      if (pred_norm_v4[i,t-1] >- 0){
        pred_norm_v4[i,t-1]<-pred_norm_v4[i,t-1]}
      else {pred_norm_v4[i,t-1]<-10^(-50)}
      if (pred_norm_v5[i,t-1] >- 0){
        pred_norm_v5[i,t-1]<-pred_norm_v5[i,t-1]}
      else {pred_norm_v5[i,t-1]<-10^(-50)}
    }
  }  
  
  pred_yr1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050)
  pred_y1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_yr2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050)
  pred_y2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_yr3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050)
  pred_y3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 1:n_2050){
      pred_yr1[i,t] <- lm1_mcmc_ref[i,69+1]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,69+4]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,69+7])
      pred_yr2[i,t] <- lm1_mcmc_ref[i,69+2]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,69+5]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,69+8])
      pred_yr3[i,t] <- lm1_mcmc_ref[i,69+3]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,69+6]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,69+9])
    }
    for (t in 2:29){
      pred_y1[i,t-1] <- yo[t,j]
      pred_y2[i,t-1] <- yo[t+29*1,j]
      pred_y3[i,t-1] <- yo[t+29*2,j]
    }
    for (t in 30:n_2050){
      pred_y1[i,t-1] <-  pred_yr1[i,t]+(pred_y1[i,t-2]-pred_yr1[i,t-1])*lm1_mcmc_y[i,j]+rnorm(1,mean=0,sd=lm1_mcmc_y[i,j+120])
      pred_y2[i,t-1] <-  pred_yr2[i,t]+(pred_y2[i,t-2]-pred_yr2[i,t-1])*lm1_mcmc_y[i,j+40*1]+rnorm(1,mean=0,sd=lm1_mcmc_y[i,j+120+40*1])
      pred_y3[i,t-1] <-  pred_yr3[i,t]+(pred_y3[i,t-2]-pred_yr3[i,t-1])*lm1_mcmc_y[i,j+40*2]+rnorm(1,mean=0,sd=lm1_mcmc_y[i,j+120+40*2])
    }
  }
  pred_norm_y1<-pred_y1/(pred_y1+pred_y2+pred_y3)
  pred_norm_y2<-pred_y2/(pred_y1+pred_y2+pred_y3)
  pred_norm_y3<-pred_y3/(pred_y1+pred_y2+pred_y3)
  pred_norm_y1[is.na(pred_norm_y1)] <- 0
  pred_norm_y2[is.na(pred_norm_y2)] <- 0
  pred_norm_y3[is.na(pred_norm_y3)] <- 0
  for (i in 1:nrow(lm1_mcmc_ref)){
    for (t in 2:n_2050){
      if (pred_norm_y1[i,t-1] >- 0){
        pred_norm_y1[i,t-1]<-pred_norm_y1[i,t-1]}
      else {pred_norm_y1[i,t-1]<-10^(-50)}
      if (pred_norm_y2[i,t-1] >- 0){
        pred_norm_y2[i,t-1]<-pred_norm_y2[i,t-1]}
      else {pred_norm_y2[i,t-1]<-10^(-50)}
      if (pred_norm_y3[i,t-1] >- 0){
        pred_norm_y3[i,t-1]<-pred_norm_y3[i,t-1]}
      else {pred_norm_y3[i,t-1]<-10^(-50)}
    }
  }   
  
  pred_e_s1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s4 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s5 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s6 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s7 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v4 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v5 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v6 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v7 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g4 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g5 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g6 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g7 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  
  for (i in 1:nrow(lm1_mcmc_fh)){
    for (t in 2:n_2050){
      pred_e_s1[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s1[i,t-1]*pred_eih1[i,t-1]/1000
      pred_e_s2[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s2[i,t-1]*pred_eih2[i,t-1]/1000
      pred_e_s3[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s3[i,t-1]*pred_eih3[i,t-1]/1000
      pred_e_s4[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s4[i,t-1]*pred_eih4[i,t-1]/1000
      pred_e_s5[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s5[i,t-1]*pred_eih5[i,t-1]/1000
      pred_e_s6[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s6[i,t-1]*pred_eih6[i,t-1]/1000
      pred_e_v1[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v1[i,t-1]*pred_eiv1[i,t-1]/1000
      pred_e_v2[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v2[i,t-1]*pred_eiv2[i,t-1]/1000
      pred_e_v3[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v3[i,t-1]*pred_eiv3[i,t-1]/1000
      pred_e_v4[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v4[i,t-1]*pred_eiv4[i,t-1]/1000
      pred_e_v5[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v5[i,t-1]*pred_eiv5[i,t-1]/1000
      pred_e_g1[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y3[i,t-1]*pred_eig1[i,t-1]/1000
    }
  }
  pred_e1<- pred_e_s1
  pred_e2<- pred_e_s2
  pred_e3<- pred_e_s3+pred_e_v1
  pred_e4<- pred_e_s4+pred_e_v2
  pred_e5<- pred_e_s5+pred_e_v3
  pred_e6<- pred_e_s6+pred_e_v4+pred_e_g1
  pred_e7<- pred_e_v5
  
  pred_e_s<- pred_e_s1+pred_e_s2+pred_e_s3+pred_e_s4+pred_e_s5+pred_e_s6
  pred_e_v<- pred_e_v1+pred_e_v2+pred_e_v3+pred_e_v4+pred_e_v5
  pred_e_g<- pred_e_g1
  pred_e<-pred_e1+pred_e2+pred_e3+pred_e4+pred_e5+pred_e6+pred_e7
  
  for (i in 1:nrow(lm1_mcmc_y)){
    for (t in 1:n_2050-1){
      pred_e1_global[i,t]<-pred_e1_global[i,t]+pred_e1[i,t]
      pred_e2_global[i,t]<-pred_e2_global[i,t]+pred_e2[i,t]
      pred_e3_global[i,t]<-pred_e3_global[i,t]+pred_e3[i,t]
      pred_e4_global[i,t]<-pred_e4_global[i,t]+pred_e4[i,t]
      pred_e5_global[i,t]<-pred_e5_global[i,t]+pred_e5[i,t]
      pred_e6_global[i,t]<-pred_e6_global[i,t]+pred_e6[i,t]
      pred_e7_global[i,t]<-pred_e7_global[i,t]+pred_e7[i,t]
      pred_e_s_global[i,t]<-pred_e_s_global[i,t]+pred_e_s[i,t]
      pred_e_v_global[i,t]<-pred_e_v_global[i,t]+pred_e_v[i,t]
      pred_e_g_global[i,t]<-pred_e_g_global[i,t]+pred_e_g[i,t]
      pred_e_global[i,t]<-pred_e_global[i,t]+pred_e[i,t]
    }
  }
  
  for (i in 1:nrow(lm1_mcmc_y)){
    cha_ei1[i,j-23]=(pred_e1[i,55]/pred_o1[i,55]*1000-his_ei[1,j-23])/his_ei[1,j-23]*100
    cha_ei2[i,j-23]=(pred_e2[i,55]/pred_o2[i,55]*1000-his_ei[2,j-23])/his_ei[2,j-23]*100
    cha_ei3[i,j-23]=(pred_e3[i,55]/pred_o3[i,55]*1000-his_ei[3,j-23])/his_ei[3,j-23]*100
    cha_ei4[i,j-23]=(pred_e4[i,55]/pred_o4[i,55]*1000-his_ei[4,j-23])/his_ei[4,j-23]*100
    cha_ei5[i,j-23]=(pred_e5[i,55]/pred_o5[i,55]*1000-his_ei[5,j-23])/his_ei[5,j-23]*100
    cha_ei6[i,j-23]=(pred_e6[i,55]/pred_o6[i,55]*1000-his_ei[6,j-23])/his_ei[6,j-23]*100
    cha_ei7[i,j-23]=(pred_e7[i,55]/pred_o7[i,55]*1000-his_ei[7,j-23])/his_ei[7,j-23]*100
  }
  
  pred_ei <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  for (i in 1:nrow(lm1_mcmc_y)){
    for (t in 1:n_2050-1){
      pred_ei[i,t]<- pred_e[i,t]/(pred_o1[i,t]+pred_o2[i,t]+pred_o3[i,t]+pred_o4[i,t]+pred_o5[i,t]+pred_o6[i,t]+pred_o7[i,t])*1000
    }
  }
# predict intervals for emissions
  uncertain_mean_reg<- matrix(NA, nrow =n_2050-1,ncol=1)
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e[,t])}  
  uncertain_lower_reg <- apply(pred_e, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e, MARGIN = 2, quantile, prob = 0.975)
  e_upper[,j-23]<- uncertain_upper_reg
  e_lower[,j-23]<- uncertain_lower_reg
  e_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e1[,t])}  
  uncertain_lower_reg <- apply(pred_e1, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e1, MARGIN = 2, quantile, prob = 0.975)
  e1_upper[,j-23]<- uncertain_upper_reg
  e1_lower[,j-23]<- uncertain_lower_reg
  e1_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e2[,t])}  
  uncertain_lower_reg <- apply(pred_e2, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e2, MARGIN = 2, quantile, prob = 0.975)
  e2_upper[,j-23]<- uncertain_upper_reg
  e2_lower[,j-23]<- uncertain_lower_reg
  e2_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e3[,t])}  
  uncertain_lower_reg <- apply(pred_e3, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e3, MARGIN = 2, quantile, prob = 0.975)
  e3_upper[,j-23]<- uncertain_upper_reg
  e3_lower[,j-23]<- uncertain_lower_reg
  e3_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e4[,t])}  
  uncertain_lower_reg <- apply(pred_e4, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e4, MARGIN = 2, quantile, prob = 0.975)
  e4_upper[,j-23]<- uncertain_upper_reg
  e4_lower[,j-23]<- uncertain_lower_reg
  e4_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e5[,t])}  
  uncertain_lower_reg <- apply(pred_e5, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e5, MARGIN = 2, quantile, prob = 0.975)
  e5_upper[,j-23]<- uncertain_upper_reg
  e5_lower[,j-23]<- uncertain_lower_reg
  e5_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e6[,t])}  
  uncertain_lower_reg <- apply(pred_e6, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e6, MARGIN = 2, quantile, prob = 0.975)
  e6_upper[,j-23]<- uncertain_upper_reg
  e6_lower[,j-23]<- uncertain_lower_reg
  e6_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e7[,t])}  
  uncertain_lower_reg <- apply(pred_e7, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e7, MARGIN = 2, quantile, prob = 0.975)
  e7_upper[,j-23]<- uncertain_upper_reg
  e7_lower[,j-23]<- uncertain_lower_reg
  e7_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e_s[,t])}  
  uncertain_lower_reg <- apply(pred_e_s, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e_s, MARGIN = 2, quantile, prob = 0.975)
  e_s_upper[,j-23]<- uncertain_upper_reg
  e_s_lower[,j-23]<- uncertain_lower_reg
  e_s_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e_v[,t])}  
  uncertain_lower_reg <- apply(pred_e_v, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e_v, MARGIN = 2, quantile, prob = 0.975)
  e_v_upper[,j-23]<- uncertain_upper_reg
  e_v_lower[,j-23]<- uncertain_lower_reg
  e_v_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e_g[,t])}  
  uncertain_lower_reg <- apply(pred_e_g, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e_g, MARGIN = 2, quantile, prob = 0.975)
  e_g_upper[,j-23]<- uncertain_upper_reg
  e_g_lower[,j-23]<- uncertain_lower_reg
  e_g_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei[,t])}  
  uncertain_lower_reg <- apply(pred_ei, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_ei, MARGIN = 2, quantile, prob = 0.975)
  ei_upper[,j-23]<- uncertain_upper_reg
  ei_lower[,j-23]<- uncertain_lower_reg
  ei_mean[,j-23]<- uncertain_mean_reg
# LMDI decompostion 
  a1=55
  for (i in 1:nrow(lm1_mcmc_y)){
    L[i,1]<-(pred_e_s1[i,a1]-his_e[1,j])/(log(pred_e_s1[i,a1]/his_e[1,j]))
    L[i,2]<-(pred_e_s2[i,a1]-his_e[2,j])/(log(pred_e_s2[i,a1]/his_e[2,j]))
    L[i,3]<-(pred_e_s3[i,a1]-his_e[3,j])/(log(pred_e_s3[i,a1]/his_e[3,j]))
    L[i,4]<-(pred_e_s4[i,a1]-his_e[4,j])/(log(pred_e_s4[i,a1]/his_e[4,j]))
    L[i,5]<-(pred_e_s5[i,a1]-his_e[5,j])/(log(pred_e_s5[i,a1]/his_e[5,j]))
    L[i,6]<-(pred_e_s6[i,a1]-his_e[6,j])/(log(pred_e_s6[i,a1]/his_e[6,j]))
    
    L[i,7]<-(pred_e_v1[i,a1]-his_e[7,j])/(log(pred_e_v1[i,a1]/his_e[7,j]))
    L[i,8]<-(pred_e_v2[i,a1]-his_e[8,j])/(log(pred_e_v2[i,a1]/his_e[8,j]))
    L[i,9]<-(pred_e_v3[i,a1]-his_e[9,j])/(log(pred_e_v3[i,a1]/his_e[9,j]))
    L[i,10]<-(pred_e_v4[i,a1]-his_e[10,j])/(log(pred_e_v4[i,a1]/his_e[10,j]))
    L[i,11]<-(pred_e_v5[i,a1]-his_e[11,j])/(log(pred_e_v5[i,a1]/his_e[11,j]))
    
    L[i,12]<-(pred_e_g1[i,a1]-his_e[12,j])/(log(pred_e_g1[i,a1]/his_e[12,j]))
  }
  
  for (i in 1:nrow(lm1_mcmc_y)){
    fac_f1[i,]=L[i,1]*log(pred_eih1[i,a1]/fh[21,j])
    fac_f2[i,]=L[i,2]*log(pred_eih2[i,a1]/fh[42,j])
    fac_f3[i,]=L[i,3]*log(pred_eih3[i,a1]/fh[63,j])+L[i,7]*log(pred_eiv1[i,a1]/fv[21,j])
    fac_f4[i,]=L[i,4]*log(pred_eih4[i,a1]/fh[84,j])+L[i,8]*log(pred_eiv2[i,a1]/fv[42,j])
    fac_f5[i,]=L[i,5]*log(pred_eih5[i,a1]/fh[105,j])+L[i,9]*log(pred_eiv3[i,a1]/fv[63,j])
    fac_f6[i,]=L[i,6]*log(pred_eih6[i,a1]/fh[126,j])+L[i,10]*log(pred_eiv4[i,a1]/fv[84,j])+L[i,12]*log(pred_eig1[i,a1]/fg[21,j])
    fac_f7[i,]=L[i,11]*log(pred_eiv5[i,a1]/fv[105,j])
    
    fac_s1[i,]=L[i,1]*(log(pred_norm_s1[i,a1]/so[21,j]))
    fac_s2[i,]=L[i,2]*(log(pred_norm_s2[i,a1]/so[42,j]))
    fac_s3[i,]=L[i,3]*(log(pred_norm_s3[i,a1]/so[63,j]))+L[i,7]*(log(pred_norm_v1[i,a1]/vo[21,j]))
    fac_s4[i,]=L[i,4]*(log(pred_norm_s4[i,a1]/so[84,j]))+L[i,8]*(log(pred_norm_v2[i,a1]/vo[42,j]))
    fac_s5[i,]=L[i,5]*(log(pred_norm_s5[i,a1]/so[105,j]))+L[i,9]*(log(pred_norm_v3[i,a1]/vo[63,j]))
    fac_s6[i,]=L[i,6]*(log(pred_norm_s6[i,a1]/so[126,j]))+L[i,10]*(log(pred_norm_v4[i,a1]/vo[84,j]))
    fac_s7[i,]=L[i,11]*(log(pred_norm_v5[i,a1]/vo[105,j]))
    
    fac_f[i,]=fac_f1[i,]+fac_f2[i,]+fac_f3[i,]+fac_f4[i,]+fac_f5[i,]+fac_f6[i,]+fac_f7[i,]
    
    fac_s[i,]=fac_s1[i,]+fac_s2[i,]+fac_s3[i,]+fac_s4[i,]+fac_s5[i,]+fac_s6[i,]+fac_s7[i,]
    
    fac_y[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6])*(log(pred_norm_y1[i,a1]/yo[21,j]))+(L[i,7]++L[i,8]+L[i,9]+L[i,10]+L[i,11])*(log(pred_norm_y2[i,a1]/yo[50,j]))+(L[i,12])*(log(pred_norm_y3[i,a1]/yo[79,j]))
    
    fac_r[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6]+L[i,7]+L[i,8]+L[i,9]+L[i,10]+L[i,11]+L[i,12])*(log(pred_r[i,a1+1]/ro[21,j]))
    
    fac_gdp[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6]+L[i,7]+L[i,8]+L[i,9]+L[i,10]+L[i,11]+L[i,12])*(log(pred_gdp[i,a1+1]/gdp[21,j]))
    
    fac_p[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6]+L[i,7]+L[i,8]+L[i,9]+L[i,10]+L[i,11]+L[i,12])*(log(pred_p[i,a1+1]/p[21,j]))
  }
  
  for (i in 1: nrow(fac_f)){
    fac_f1_reg[i,j-23]<-fac_f1[i,]
    fac_f2_reg[i,j-23]<-fac_f2[i,]
    fac_f3_reg[i,j-23]<-fac_f3[i,]
    fac_f4_reg[i,j-23]<-fac_f4[i,]
    fac_f5_reg[i,j-23]<-fac_f5[i,]
    fac_f6_reg[i,j-23]<-fac_f6[i,]
    fac_f7_reg[i,j-23]<-fac_f7[i,]
    fac_s1_reg[i,j-23]<-fac_s1[i,]
    fac_s2_reg[i,j-23]<-fac_s2[i,]
    fac_s3_reg[i,j-23]<-fac_s3[i,]
    fac_s4_reg[i,j-23]<-fac_s4[i,]
    fac_s5_reg[i,j-23]<-fac_s5[i,]
    fac_s6_reg[i,j-23]<-fac_s6[i,]
    fac_s7_reg[i,j-23]<-fac_s7[i,]
    fac_f_reg[i,j-23]<-fac_f[i,]
    fac_s_reg[i,j-23]<-fac_s[i,]
    fac_y_reg[i,j-23]<-fac_y[i,]
    fac_gdp_reg[i,j-23]<-fac_gdp[i,]
    fac_p_reg[i,j-23]<-fac_p[i,]
    fac_r_reg[i,j-23]<-fac_r[i,]
  }
}
# projection for India, Indonesia and Rest of World 
for (j in 38:40){
  pred_p <- matrix(NA, nrow = nrow(lm1_mcmc_p),ncol=n_2050)
  for (i in 1:nrow(lm1_mcmc_p)){
    for (t in 1:26){pred_p[i,t] <- exp(x1[t,j])}
    for (t in 27:56){
      if (pred_p[i,t-1] > 0){pred_p[i,t]<-exp(rnorm(1,mean=0,sd=sd(lm1_mcmc_p[,j]))+log(pred_p[i,t-1])+d1[t-1,j])+rnorm(1,mean=0,sd=lm1_mcmc_p[i,j+40])}
      else {pred_p[i,t]<-0}}
  }
  for (i in 1:nrow(lm1_mcmc_gdp_pre)){
    for (t in 1:n_2050){
      if (pred_p[i,t] > 0){pred_p[i,t] <- pred_p[i,t]}
      else {pred_p[i,t]<-10^(-50)}
    }
  }
  
  
  
  pred_gdp <- matrix(NA, nrow = nrow(lm1_mcmc_gdp_pre),ncol=n_2050)
  for (i in 1:nrow(lm1_mcmc_gdp_pre)){
    for (t in 1:24){pred_gdp[i,t] <- exp(x[t,j])}
    for (t in 25:25){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*1]}
      else {pred_gdp[i,t]<-0}}
    for (t in 26:26){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*2]}
      else {pred_gdp[i,t]<-0}}
    for (t in 27:27){
      if (pred_gdp[i,t-1] >- 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*3]}
      else {pred_gdp[i,t]<-0}}
    for (t in 28:28){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*4]}
      else {pred_gdp[i,t]<-0}}
    for (t in 29:29){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*5]}
      else {pred_gdp[i,t]<-0}}
    for (t in 30:30){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*6]}
      else {pred_gdp[i,t]<-0}}
    for (t in 31:31){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*7]}
      else {pred_gdp[i,t]<-0}}
    for (t in 32:32){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*8]}
      else {pred_gdp[i,t]<-0}}
    for (t in 33:33){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*9]}
      else {pred_gdp[i,t]<-0}}
    for (t in 34:34){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*10]}
      else {pred_gdp[i,t]<-0}}
    for (t in 35:35){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*11]}
      else {pred_gdp[i,t]<-0}}
    for (t in 36:36){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*12]}
      else {pred_gdp[i,t]<-0}}
    for (t in 37:37){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*13]}
      else {pred_gdp[i,t]<-0}}
    for (t in 38:38){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*14]}
      else {pred_gdp[i,t]<-0}}
    for (t in 39:39){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*15]}
      else {pred_gdp[i,t]<-0}}
    for (t in 40:40){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*16]}
      else {pred_gdp[i,t]<-0}}
    for (t in 41:41){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*17]}
      else {pred_gdp[i,t]<-0}}
    for (t in 42:42){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*18]}
      else {pred_gdp[i,t]<-0}}
    for (t in 43:43){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*19]}
      else {pred_gdp[i,t]<-0}}
    for (t in 44:44){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*20]}
      else {pred_gdp[i,t]<-0}}
    for (t in 45:45){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*21]}
      else {pred_gdp[i,t]<-0}}
    for (t in 46:46){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*22]}
      else {pred_gdp[i,t]<-0}}
    for (t in 47:47){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*23]}
      else {pred_gdp[i,t]<-0}}
    for (t in 48:48){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*24]}
      else {pred_gdp[i,t]<-0}}
    for (t in 49:49){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*25]}
      else {pred_gdp[i,t]<-0}}
    for (t in 50:50){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*26]}
      else {pred_gdp[i,t]<-0}}
    for (t in 51:51){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*27]}
      else {pred_gdp[i,t]<-0}}
    for (t in 52:52){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*28]}
      else {pred_gdp[i,t]<-0}}
    for (t in 53:53){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*29]}
      else {pred_gdp[i,t]<-0}}
    for (t in 54:54){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*30]}
      else {pred_gdp[i,t]<-0}}
    for (t in 55:55){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*31]}
      else {pred_gdp[i,t]<-0}}
    for (t in 56:56){
      if (pred_gdp[i,t-1] > 0){pred_gdp[i,t]<-exp(lm1_mcmc_gdp_pre[i,j]+log(pred_gdp[i,t-1])+d[t-1,j])+lm1_mcmc_gdp_pre[i,j+40*32]}
      else {pred_gdp[i,t]<-0}}
  }
  for (i in 1:nrow(lm1_mcmc_gdp_pre)){
    for (t in 1:n_2050){
      if (pred_gdp[i,t] > 0){pred_gdp[i,t] <- pred_gdp[i,t]}
      else {pred_gdp[i,t]<-10^(-50)}
    }
  }
  
  
  pred_r <- matrix(NA, nrow = nrow(lm1_mcmc_p),ncol=n_2050)
  for (i in 1:nrow(lm1_mcmc_r)){
    for (t in 1:n_2015){pred_r[i,t] <- x2[t,j]}
    for (t in 22:56){
      if (pred_r[i,t-1] > 0){pred_r[i,t]<-rnorm(1,mean=0,sd=lm1_mcmc_r[i,j+80])+pred_r[i,t-1]+d2[t-1,j]+rnorm(1,mean=0,sd=lm1_mcmc_r[i,j+40])}
      else {pred_r[i,t]<-0}}
  }
  for (i in 1:nrow(lm1_mcmc_r)){
    for (t in 1:n_2050){
      if (pred_r[i,t] > 0){pred_r[i,t] <- pred_r[i,t]}
      else {pred_r[i,t]<-10^(-50)}
    }
  }
  
  pred_eih1 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih2 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih3 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih4 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih5 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eih6 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv1 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv2 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv3 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv4 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eiv5 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  pred_eig1 <- matrix(NA, nrow = nrow(lm1_mcmc_fh), ncol=n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 2:21){
      pred_eih1[i,t-1] <- eih[t,j]
      pred_eih2[i,t-1] <- eih[t+21,j]
      pred_eih3[i,t-1] <- eih[t+21*2,j] 
      pred_eih4[i,t-1] <- eih[t+21*3,j]       
      pred_eih5[i,t-1] <- eih[t+21*4,j] 
      pred_eih6[i,t-1] <- eih[t+21*5,j] 
      
      pred_eiv1[i,t-1] <- eiv[t,j]
      pred_eiv2[i,t-1] <- eiv[t+21,j]
      pred_eiv3[i,t-1] <- eiv[t+21*2,j] 
      pred_eiv4[i,t-1] <- eiv[t+21*3,j]       
      pred_eiv5[i,t-1] <- eiv[t+21*4,j] 
      pred_eig1[i,t-1] <- eig[t,j]
    }
    for (t in 22:n_2050){
      pred_eih1[i,t-1] <-(pred_eihr1[i,t]+(pred_eih1[i,t-2]-pred_eihr1[i,t-1])*lm1_mcmc_fh[i,j+40*0])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240])
      pred_eih2[i,t-1] <-(pred_eihr2[i,t]+(pred_eih2[i,t-2]-pred_eihr2[i,t-1])*lm1_mcmc_fh[i,j+40*1])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*1])
      pred_eih3[i,t-1] <-(pred_eihr3[i,t]+(pred_eih3[i,t-2]-pred_eihr3[i,t-1])*lm1_mcmc_fh[i,j+40*2])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*2])
      pred_eih4[i,t-1] <-(pred_eihr4[i,t]+(pred_eih4[i,t-2]-pred_eihr4[i,t-1])*lm1_mcmc_fh[i,j+40*3])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*3])
      pred_eih5[i,t-1] <-(pred_eihr5[i,t]+(pred_eih5[i,t-2]-pred_eihr5[i,t-1])*lm1_mcmc_fh[i,j+40*4])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*4])
      pred_eih6[i,t-1] <-(pred_eihr6[i,t]+(pred_eih6[i,t-2]-pred_eihr6[i,t-1])*lm1_mcmc_fh[i,j+40*5])+rnorm(1,mean=0,sd=lm1_mcmc_fh[i,j+240+40*5])
      pred_eiv1[i,t-1] <-(pred_eivr1[i,t]+(pred_eiv1[i,t-2]-pred_eivr1[i,t-1])*lm1_mcmc_fv[i,j+40*0])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200])
      pred_eiv2[i,t-1] <-(pred_eivr2[i,t]+(pred_eiv2[i,t-2]-pred_eivr2[i,t-1])*lm1_mcmc_fv[i,j+40*1])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*1])
      pred_eiv3[i,t-1] <-(pred_eivr3[i,t]+(pred_eiv3[i,t-2]-pred_eivr3[i,t-1])*lm1_mcmc_fv[i,j+40*2])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*2])
      pred_eiv4[i,t-1] <-(pred_eivr4[i,t]+(pred_eiv4[i,t-2]-pred_eivr4[i,t-1])*lm1_mcmc_fv[i,j+40*3])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*3])
      pred_eiv5[i,t-1] <-(pred_eivr5[i,t]+(pred_eiv5[i,t-2]-pred_eivr5[i,t-1])*lm1_mcmc_fv[i,j+40*4])+rnorm(1,mean=0,sd=lm1_mcmc_fv[i,j+200+40*4])
      pred_eig1[i,t-1] <-(pred_eigr1[i,t]+(pred_eig1[i,t-2]-pred_eigr1[i,t-1])*lm1_mcmc_fg[i,j+40*0])+rnorm(1,mean=0,sd=lm1_mcmc_fg[i,j+40])
    }
    for (t in 2:n_2050){
      pred_eih1[i,t-1]=exp(pred_eih1[i,t-1])
      pred_eih2[i,t-1]=exp(pred_eih2[i,t-1])
      pred_eih3[i,t-1]=exp(pred_eih3[i,t-1])
      pred_eih4[i,t-1]=exp(pred_eih4[i,t-1])
      pred_eih5[i,t-1]=exp(pred_eih5[i,t-1])
      pred_eih6[i,t-1]=exp(pred_eih6[i,t-1])
      pred_eiv1[i,t-1]=exp(pred_eiv1[i,t-1])
      pred_eiv2[i,t-1]=exp(pred_eiv2[i,t-1])
      pred_eiv3[i,t-1]=exp(pred_eiv3[i,t-1])
      pred_eiv4[i,t-1]=exp(pred_eiv4[i,t-1])
      pred_eiv5[i,t-1]=exp(pred_eiv5[i,t-1])
      pred_eig1[i,t-1]=exp(pred_eig1[i,t-1])
    }
    for (t in 2:n_2050){
      if (pred_eih1[i,t-1] >- 0){
        pred_eih1[i,t-1]<-pred_eih1[i,t-1]}
      else {pred_eih1[i,t-1]<-10^(-50)}
      if (pred_eih2[i,t-1] >- 0){
        pred_eih2[i,t-1]<-pred_eih2[i,t-1]}
      else {pred_eih2[i,t-1]<-10^(-50)}
      if (pred_eih3[i,t-1] >- 0){
        pred_eih3[i,t-1]<-pred_eih3[i,t-1]}
      else {pred_eih3[i,t-1]<-10^(-50)}
      if (pred_eih4[i,t-1] >- 0){
        pred_eih4[i,t-1]<-pred_eih4[i,t-1]}
      else {pred_eih4[i,t-1]<-10^(-50)}
      if (pred_eih5[i,t-1] >- 0){
        pred_eih5[i,t-1]=pred_eih5[i,t-1]}
      else {pred_eih5[i,t-1]<-10^(-50)}
      if (pred_eih6[i,t-1] >- 0){
        pred_eih6[i,t-1]<-pred_eih6[i,t-1]}
      else {pred_eih6[i,t-1]<-10^(-50)}
      if (pred_eiv1[i,t-1] >- 0){
        pred_eiv1[i,t-1]<-pred_eiv1[i,t-1]}
      else {pred_eiv1[i,t-1]<-10^(-50)}
      if (pred_eiv2[i,t-1] >- 0){
        pred_eiv2[i,t-1]<-pred_eiv2[i,t-1]}
      else {pred_eiv2[i,t-1]<-10^(-50)}
      if (pred_eiv3[i,t-1] >- 0){
        pred_eiv3[i,t-1]<-pred_eiv3[i,t-1]}
      else {pred_eiv3[i,t-1]<-10^(-50)}
      if (pred_eiv4[i,t-1] >- 0){
        pred_eiv4[i,t-1]<-pred_eiv4[i,t-1]}
      else {pred_eiv4[i,t-1]<-10^(-50)}
      if (pred_eiv5[i,t-1] >- 0){
        pred_eiv5[i,t-1]=pred_eiv5[i,t-1]}
      else {pred_eiv5[i,t-1]<-10^(-50)}
      if (pred_eig1[i,t-1] >- 0){
        pred_eig1[i,t-1]<-pred_eig1[i,t-1]}
      else {pred_eig1[i,t-1]<-10^(-50)}
    }
  }
  
  pred_sr1 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s1 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr2 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s2 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr3 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s3 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr4 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s4 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr5 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s5 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  pred_sr6 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050)
  pred_s6 <- matrix(NA, nrow = nrow(lm1_mcmc_s), ncol = n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 1:n_2050){
      pred_sr1[i,t] <- lm1_mcmc_ref[i,36+1]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+7]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+13])
      pred_sr2[i,t] <- lm1_mcmc_ref[i,36+2]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+8]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+14])
      pred_sr3[i,t] <- lm1_mcmc_ref[i,36+3]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+9]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+15])
      pred_sr4[i,t] <- lm1_mcmc_ref[i,36+4]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+10]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+16])
      pred_sr5[i,t] <- lm1_mcmc_ref[i,36+5]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+11]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+17])
      pred_sr6[i,t] <- lm1_mcmc_ref[i,36+6]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,36+12]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,36+18])
    }
    for (t in 2:21){
      pred_s1[i,t-1] <- s[t,j]
      pred_s2[i,t-1] <- s[t+21*1,j]
      pred_s3[i,t-1] <- s[t+21*2,j]
      pred_s4[i,t-1] <- s[t+21*3,j]
      pred_s5[i,t-1] <- s[t+21*4,j]
      pred_s6[i,t-1] <- s[t+21*5,j]
    }
    for (t in 22:n_2050){
      pred_s1[i,t-1] <-  pred_sr1[i,t]+(pred_s1[i,t-2]-pred_sr1[i,t-1])*lm1_mcmc_s[i,j]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240])
      pred_s2[i,t-1] <-  pred_sr2[i,t]+(pred_s2[i,t-2]-pred_sr2[i,t-1])*lm1_mcmc_s[i,j+40*1]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*1])
      pred_s3[i,t-1] <-  pred_sr3[i,t]+(pred_s3[i,t-2]-pred_sr3[i,t-1])*lm1_mcmc_s[i,j+40*2]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*2])
      pred_s4[i,t-1] <-  pred_sr4[i,t]+(pred_s4[i,t-2]-pred_sr4[i,t-1])*lm1_mcmc_s[i,j+40*3]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*3])
      pred_s5[i,t-1] <-  pred_sr5[i,t]+(pred_s5[i,t-2]-pred_sr5[i,t-1])*lm1_mcmc_s[i,j+40*4]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*4])
      pred_s6[i,t-1] <-  pred_sr6[i,t]+(pred_s6[i,t-2]-pred_sr6[i,t-1])*lm1_mcmc_s[i,j+40*5]+rnorm(1,mean=0,sd=lm1_mcmc_s[i,j+240+40*5])
    }
    for (t in 2:n_2050){
      pred_s1[i,t-1]=exp(pred_s1[i,t-1])
      pred_s2[i,t-1]=exp(pred_s2[i,t-1])
      pred_s3[i,t-1]=exp(pred_s3[i,t-1])
      pred_s4[i,t-1]=exp(pred_s4[i,t-1])
      pred_s5[i,t-1]=exp(pred_s5[i,t-1])
      pred_s6[i,t-1]=exp(pred_s6[i,t-1])
    }
  }
  pred_norm_s1<-pred_s1/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s2<-pred_s2/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s3<-pred_s3/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s4<-pred_s4/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s5<-pred_s5/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s6<-pred_s6/(pred_s1+pred_s2+pred_s3+pred_s4+pred_s5+pred_s6)
  pred_norm_s1[is.na(pred_norm_s1)] <- 0
  pred_norm_s2[is.na(pred_norm_s2)] <- 0
  pred_norm_s3[is.na(pred_norm_s3)] <- 0
  pred_norm_s4[is.na(pred_norm_s4)] <- 0
  pred_norm_s5[is.na(pred_norm_s5)] <- 0
  pred_norm_s6[is.na(pred_norm_s6)] <- 0
  for (i in 1:nrow(lm1_mcmc_ref)){
    for (t in 2:n_2050){
      if (pred_norm_s1[i,t-1] > 0){
        pred_norm_s1[i,t-1]<-pred_norm_s1[i,t-1]}
      else {pred_norm_s1[i,t-1]<-10^(-50)}
      if (pred_norm_s2[i,t-1] > 0){
        pred_norm_s2[i,t-1]<-pred_norm_s2[i,t-1]}
      else {pred_norm_s2[i,t-1]<-10^(-50)}
      if (pred_norm_s3[i,t-1] >- 0){
        pred_norm_s3[i,t-1]=pred_norm_s3[i,t-1]}
      else {pred_norm_s3[i,t-1]<-10^(-50)}
      if (pred_norm_s4[i,t-1] > 0){
        pred_norm_s4[i,t-1]<-pred_norm_s4[i,t-1]}
      else {pred_norm_s4[i,t-1]<-10^(-50)}
      if (pred_norm_s5[i,t-1] > 0){
        pred_norm_s5[i,t-1]<-pred_norm_s5[i,t-1]}
      else {pred_norm_s5[i,t-1]<-10^(-50)}
      if (pred_norm_s6[i,t-1] > 0){
        pred_norm_s6[i,t-1]<-pred_norm_s6[i,t-1]}
      else {pred_norm_s6[i,t-1]<-10^(-50)}
    }
  } 
  
  pred_vr1 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v1 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr2 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v2 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr3 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v3 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr4 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v4 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  pred_vr5 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050)
  pred_v5 <- matrix(NA, nrow = nrow(lm1_mcmc_v), ncol = n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 1:n_2050){
      pred_vr1[i,t] <- lm1_mcmc_ref[i,54+1]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+6]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+11])
      pred_vr2[i,t] <- lm1_mcmc_ref[i,54+2]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+7]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+12])
      pred_vr3[i,t] <- lm1_mcmc_ref[i,54+3]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+8]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+13])
      pred_vr4[i,t] <- lm1_mcmc_ref[i,54+4]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+9]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+14])
      pred_vr5[i,t] <- lm1_mcmc_ref[i,54+5]+log(pred_gdp[i,t])*lm1_mcmc_ref[i,54+10]+rnorm(1,mean=0,sd=lm1_mcmc_ref[i,54+15])
    }
    for (t in 2:21){
      pred_v1[i,t-1] <- v[t,j]
      pred_v2[i,t-1] <- v[t+21*1,j]
      pred_v3[i,t-1] <- v[t+21*2,j]
      pred_v4[i,t-1] <- v[t+21*3,j]
      pred_v5[i,t-1] <- v[t+21*4,j]
    }
    for (t in 22:n_2050){
      pred_v1[i,t-1] <-  pred_vr1[i,t]+(pred_v1[i,t-2]-pred_vr1[i,t-1])*lm1_mcmc_v[i,j]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200])
      pred_v2[i,t-1] <-  pred_vr2[i,t]+(pred_v2[i,t-2]-pred_vr2[i,t-1])*lm1_mcmc_v[i,j+40*1]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*1])
      pred_v3[i,t-1] <-  pred_vr3[i,t]+(pred_v3[i,t-2]-pred_vr3[i,t-1])*lm1_mcmc_v[i,j+40*2]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*2])
      pred_v4[i,t-1] <-  pred_vr4[i,t]+(pred_v4[i,t-2]-pred_vr4[i,t-1])*lm1_mcmc_v[i,j+40*3]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*3])
      pred_v5[i,t-1] <-  pred_vr5[i,t]+(pred_v5[i,t-2]-pred_vr5[i,t-1])*lm1_mcmc_v[i,j+40*4]+rnorm(1,mean=0,sd=lm1_mcmc_v[i,j+200+40*4])
    }
    for (t in 2:n_2050){
      pred_v1[i,t-1]=exp(pred_v1[i,t-1])
      pred_v2[i,t-1]=exp(pred_v2[i,t-1])
      pred_v3[i,t-1]=exp(pred_v3[i,t-1])
      pred_v4[i,t-1]=exp(pred_v4[i,t-1])
      pred_v5[i,t-1]=exp(pred_v5[i,t-1])
    }
  }
  pred_norm_v1<-pred_v1/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v2<-pred_v2/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v3<-pred_v3/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v4<-pred_v4/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v5<-pred_v5/(pred_v1+pred_v2+pred_v3+pred_v4+pred_v5)
  pred_norm_v1[is.na(pred_norm_v1)] <- 0
  pred_norm_v2[is.na(pred_norm_v2)] <- 0
  pred_norm_v3[is.na(pred_norm_v3)] <- 0
  pred_norm_v4[is.na(pred_norm_v4)] <- 0
  pred_norm_v5[is.na(pred_norm_v5)] <- 0
  for (i in 1:nrow(lm1_mcmc_ref)){
    for (t in 2:n_2050){
      if (pred_norm_v1[i,t-1] >- 0){
        pred_norm_v1[i,t-1]<-pred_norm_v1[i,t-1]}
      else {pred_norm_v1[i,t-1]<-10^(-50)}
      if (pred_norm_v2[i,t-1] >- 0){
        pred_norm_v2[i,t-1]<-pred_norm_v2[i,t-1]}
      else {pred_norm_v2[i,t-1]<-10^(-50)}
      if (pred_norm_v3[i,t-1] >- 0){
        pred_norm_v3[i,t-1]<-pred_norm_v3[i,t-1]}
      else {pred_norm_v3[i,t-1]<-10^(-50)}
      if (pred_norm_v4[i,t-1] >- 0){
        pred_norm_v4[i,t-1]<-pred_norm_v4[i,t-1]}
      else {pred_norm_v4[i,t-1]<-10^(-50)}
      if (pred_norm_v5[i,t-1] >- 0){
        pred_norm_v5[i,t-1]<-pred_norm_v5[i,t-1]}
      else {pred_norm_v5[i,t-1]<-10^(-50)}
    }
  }  
  
  pred_yr1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050)
  pred_y1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_yr2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050)
  pred_y2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_yr3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050)
  pred_y3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  for (i in 1:nrow(lm1_mcmc_s)){
    for (t in 1:n_2050){
      pred_yr1[i,t] <- lm1_mcmc_ref1[i,1]+log(pred_gdp[i,t])*lm1_mcmc_ref1[i,4]+log(pred_gdp[i,t])*log(pred_gdp[i,t])*lm1_mcmc_ref1[i,7]+rnorm(1,mean=0,sd=lm1_mcmc_ref1[i,10])
      pred_yr2[i,t] <- lm1_mcmc_ref1[i,2]+log(pred_gdp[i,t])*lm1_mcmc_ref1[i,5]+log(pred_gdp[i,t])*log(pred_gdp[i,t])*lm1_mcmc_ref1[i,8]+rnorm(1,mean=0,sd=lm1_mcmc_ref1[i,11])
      pred_yr3[i,t] <- lm1_mcmc_ref1[i,3]+log(pred_gdp[i,t])*lm1_mcmc_ref1[i,6]+log(pred_gdp[i,t])*log(pred_gdp[i,t])*lm1_mcmc_ref1[i,9]+rnorm(1,mean=0,sd=lm1_mcmc_ref1[i,12])
    }
    for (t in 2:29){
      pred_y1[i,t-1] <- yo[t,j]
      pred_y2[i,t-1] <- yo[t+29*1,j]
      pred_y3[i,t-1] <- yo[t+29*2,j]
    }
    for (t in 30:n_2050){
      pred_y1[i,t-1] <-  pred_yr1[i,t-1]+(pred_y1[i,t-2]-pred_yr1[i,t-2])*lm1_mcmc_y[i,j]+rnorm(1,mean=0,sd=lm1_mcmc_y[i,j+120])
      pred_y2[i,t-1] <-  pred_yr2[i,t-1]+(pred_y2[i,t-2]-pred_yr2[i,t-2])*lm1_mcmc_y[i,j+40*1]+rnorm(1,mean=0,sd=lm1_mcmc_y[i,j+120+40*1])
      pred_y3[i,t-1] <-  pred_yr3[i,t-1]+(pred_y3[i,t-2]-pred_yr3[i,t-2])*lm1_mcmc_y[i,j+40*2]+rnorm(1,mean=0,sd=lm1_mcmc_y[i,j+120+40*2])
    }
  }
  pred_norm_y1<-pred_y1/(pred_y1+pred_y2+pred_y3)
  pred_norm_y2<-pred_y2/(pred_y1+pred_y2+pred_y3)
  pred_norm_y3<-pred_y3/(pred_y1+pred_y2+pred_y3)
  pred_norm_y1[is.na(pred_norm_y1)] <- 0
  pred_norm_y2[is.na(pred_norm_y2)] <- 0
  pred_norm_y3[is.na(pred_norm_y3)] <- 0
  for (i in 1:nrow(lm1_mcmc_ref)){
    for (t in 2:n_2050){
      if (pred_norm_y1[i,t-1] >- 0){
        pred_norm_y1[i,t-1]<-pred_norm_y1[i,t-1]}
      else {pred_norm_y1[i,t-1]<-10^(-50)}
      if (pred_norm_y2[i,t-1] >- 0){
        pred_norm_y2[i,t-1]<-pred_norm_y2[i,t-1]}
      else {pred_norm_y2[i,t-1]<-10^(-50)}
      if (pred_norm_y3[i,t-1] >- 0){
        pred_norm_y3[i,t-1]<-pred_norm_y3[i,t-1]}
      else {pred_norm_y3[i,t-1]<-10^(-50)}
    }
  }  
  
  pred_e_s1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s4 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s5 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s6 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_s7 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v4 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v5 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v6 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_v7 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g1 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g2 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g3 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g4 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g5 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g6 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  pred_e_g7 <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  
  for (i in 1:nrow(lm1_mcmc_fh)){
    for (t in 2:n_2050){
      pred_e_s1[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s1[i,t-1]*pred_eih1[i,t-1]/1000
      pred_e_s2[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s2[i,t-1]*pred_eih2[i,t-1]/1000
      pred_e_s3[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s3[i,t-1]*pred_eih3[i,t-1]/1000
      pred_e_s4[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s4[i,t-1]*pred_eih4[i,t-1]/1000
      pred_e_s5[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s5[i,t-1]*pred_eih5[i,t-1]/1000
      pred_e_s6[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y1[i,t-1]*pred_norm_s6[i,t-1]*pred_eih6[i,t-1]/1000
      pred_e_v1[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v1[i,t-1]*pred_eiv1[i,t-1]/1000
      pred_e_v2[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v2[i,t-1]*pred_eiv2[i,t-1]/1000
      pred_e_v3[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v3[i,t-1]*pred_eiv3[i,t-1]/1000
      pred_e_v4[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v4[i,t-1]*pred_eiv4[i,t-1]/1000
      pred_e_v5[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y2[i,t-1]*pred_norm_v5[i,t-1]*pred_eiv5[i,t-1]/1000
      pred_e_g1[i,t-1]<-pred_p[i,t]*pred_gdp[i,t]*pred_r[i,t]*pred_norm_y3[i,t-1]*pred_eig1[i,t-1]/1000
    }
  }
  pred_e1<- pred_e_s1
  pred_e2<- pred_e_s2
  pred_e3<- pred_e_s3+pred_e_v1
  pred_e4<- pred_e_s4+pred_e_v2
  pred_e5<- pred_e_s5+pred_e_v3
  pred_e6<- pred_e_s6+pred_e_v4+pred_e_g1
  pred_e7<- pred_e_v5
  
  pred_e_s<- pred_e_s1+pred_e_s2+pred_e_s3+pred_e_s4+pred_e_s5+pred_e_s6
  pred_e_v<- pred_e_v1+pred_e_v2+pred_e_v3+pred_e_v4+pred_e_v5
  pred_e_g<- pred_e_g1
  pred_e<-pred_e1+pred_e2+pred_e3+pred_e4+pred_e5+pred_e6+pred_e7
  
  for (i in 1:nrow(lm1_mcmc_y)){
    for (t in 1:n_2050-1){
      pred_e1_global[i,t]<-pred_e1_global[i,t]+pred_e1[i,t]
      pred_e2_global[i,t]<-pred_e2_global[i,t]+pred_e2[i,t]
      pred_e3_global[i,t]<-pred_e3_global[i,t]+pred_e3[i,t]
      pred_e4_global[i,t]<-pred_e4_global[i,t]+pred_e4[i,t]
      pred_e5_global[i,t]<-pred_e5_global[i,t]+pred_e5[i,t]
      pred_e6_global[i,t]<-pred_e6_global[i,t]+pred_e6[i,t]
      pred_e7_global[i,t]<-pred_e7_global[i,t]+pred_e7[i,t]
      pred_e_s_global[i,t]<-pred_e_s_global[i,t]+pred_e_s[i,t]
      pred_e_v_global[i,t]<-pred_e_v_global[i,t]+pred_e_v[i,t]
      pred_e_g_global[i,t]<-pred_e_g_global[i,t]+pred_e_g[i,t]
      pred_e_global[i,t]<-pred_e_global[i,t]+pred_e[i,t]
    }
  }
  
  for (i in 1:nrow(lm1_mcmc_y)){
    cha_ei1[i,j-23]=(pred_e1[i,55]/pred_o1[i,55]*1000-his_ei[1,j-23])/his_ei[1,j-23]*100
    cha_ei2[i,j-23]=(pred_e2[i,55]/pred_o2[i,55]*1000-his_ei[2,j-23])/his_ei[2,j-23]*100
    cha_ei3[i,j-23]=(pred_e3[i,55]/pred_o3[i,55]*1000-his_ei[3,j-23])/his_ei[3,j-23]*100
    cha_ei4[i,j-23]=(pred_e4[i,55]/pred_o4[i,55]*1000-his_ei[4,j-23])/his_ei[4,j-23]*100
    cha_ei5[i,j-23]=(pred_e5[i,55]/pred_o5[i,55]*1000-his_ei[5,j-23])/his_ei[5,j-23]*100
    cha_ei6[i,j-23]=(pred_e6[i,55]/pred_o6[i,55]*1000-his_ei[6,j-23])/his_ei[6,j-23]*100
    cha_ei7[i,j-23]=(pred_e7[i,55]/pred_o7[i,55]*1000-his_ei[7,j-23])/his_ei[7,j-23]*100
  }
  
  pred_ei <- matrix(NA, nrow = nrow(lm1_mcmc_y), ncol = n_2050-1)
  for (i in 1:nrow(lm1_mcmc_y)){
    for (t in 1:n_2050-1){
      pred_ei[i,t]<- pred_e[i,t]/(pred_o1[i,t]+pred_o2[i,t]+pred_o3[i,t]+pred_o4[i,t]+pred_o5[i,t]+pred_o6[i,t]+pred_o7[i,t])*1000
    }
  }
  
  uncertain_mean_reg<- matrix(NA, nrow =n_2050-1,ncol=1)
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e[,t])}  
  uncertain_lower_reg <- apply(pred_e, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e, MARGIN = 2, quantile, prob = 0.975)
  e_upper[,j-23]<- uncertain_upper_reg
  e_lower[,j-23]<- uncertain_lower_reg
  e_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e1[,t])}  
  uncertain_lower_reg <- apply(pred_e1, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e1, MARGIN = 2, quantile, prob = 0.975)
  e1_upper[,j-23]<- uncertain_upper_reg
  e1_lower[,j-23]<- uncertain_lower_reg
  e1_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e2[,t])}  
  uncertain_lower_reg <- apply(pred_e2, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e2, MARGIN = 2, quantile, prob = 0.975)
  e2_upper[,j-23]<- uncertain_upper_reg
  e2_lower[,j-23]<- uncertain_lower_reg
  e2_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e3[,t])}  
  uncertain_lower_reg <- apply(pred_e3, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e3, MARGIN = 2, quantile, prob = 0.975)
  e3_upper[,j-23]<- uncertain_upper_reg
  e3_lower[,j-23]<- uncertain_lower_reg
  e3_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e4[,t])}  
  uncertain_lower_reg <- apply(pred_e4, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e4, MARGIN = 2, quantile, prob = 0.975)
  e4_upper[,j-23]<- uncertain_upper_reg
  e4_lower[,j-23]<- uncertain_lower_reg
  e4_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e5[,t])}  
  uncertain_lower_reg <- apply(pred_e5, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e5, MARGIN = 2, quantile, prob = 0.975)
  e5_upper[,j-23]<- uncertain_upper_reg
  e5_lower[,j-23]<- uncertain_lower_reg
  e5_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e6[,t])}  
  uncertain_lower_reg <- apply(pred_e6, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e6, MARGIN = 2, quantile, prob = 0.975)
  e6_upper[,j-23]<- uncertain_upper_reg
  e6_lower[,j-23]<- uncertain_lower_reg
  e6_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e7[,t])}  
  uncertain_lower_reg <- apply(pred_e7, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e7, MARGIN = 2, quantile, prob = 0.975)
  e7_upper[,j-23]<- uncertain_upper_reg
  e7_lower[,j-23]<- uncertain_lower_reg
  e7_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e_s[,t])}  
  uncertain_lower_reg <- apply(pred_e_s, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e_s, MARGIN = 2, quantile, prob = 0.975)
  e_s_upper[,j-23]<- uncertain_upper_reg
  e_s_lower[,j-23]<- uncertain_lower_reg
  e_s_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e_v[,t])}  
  uncertain_lower_reg <- apply(pred_e_v, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e_v, MARGIN = 2, quantile, prob = 0.975)
  e_v_upper[,j-23]<- uncertain_upper_reg
  e_v_lower[,j-23]<- uncertain_lower_reg
  e_v_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_e_g[,t])}  
  uncertain_lower_reg <- apply(pred_e_g, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_e_g, MARGIN = 2, quantile, prob = 0.975)
  e_g_upper[,j-23]<- uncertain_upper_reg
  e_g_lower[,j-23]<- uncertain_lower_reg
  e_g_mean[,j-23]<- uncertain_mean_reg
  
  for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei[,t])}  
  uncertain_lower_reg <- apply(pred_ei, MARGIN = 2, quantile, prob = 0.025)
  uncertain_upper_reg <- apply(pred_ei, MARGIN = 2, quantile, prob = 0.975)
  ei_upper[,j-23]<- uncertain_upper_reg
  ei_lower[,j-23]<- uncertain_lower_reg
  ei_mean[,j-23]<- uncertain_mean_reg

  a1=55
  for (i in 1:nrow(lm1_mcmc_y)){
    L[i,1]<-(pred_e_s1[i,a1]-his_e[1,j])/(log(pred_e_s1[i,a1]/his_e[1,j]))
    L[i,2]<-(pred_e_s2[i,a1]-his_e[2,j])/(log(pred_e_s2[i,a1]/his_e[2,j]))
    L[i,3]<-(pred_e_s3[i,a1]-his_e[3,j])/(log(pred_e_s3[i,a1]/his_e[3,j]))
    L[i,4]<-(pred_e_s4[i,a1]-his_e[4,j])/(log(pred_e_s4[i,a1]/his_e[4,j]))
    L[i,5]<-(pred_e_s5[i,a1]-his_e[5,j])/(log(pred_e_s5[i,a1]/his_e[5,j]))
    L[i,6]<-(pred_e_s6[i,a1]-his_e[6,j])/(log(pred_e_s6[i,a1]/his_e[6,j]))
    
    L[i,7]<-(pred_e_v1[i,a1]-his_e[7,j])/(log(pred_e_v1[i,a1]/his_e[7,j]))
    L[i,8]<-(pred_e_v2[i,a1]-his_e[8,j])/(log(pred_e_v2[i,a1]/his_e[8,j]))
    L[i,9]<-(pred_e_v3[i,a1]-his_e[9,j])/(log(pred_e_v3[i,a1]/his_e[9,j]))
    L[i,10]<-(pred_e_v4[i,a1]-his_e[10,j])/(log(pred_e_v4[i,a1]/his_e[10,j]))
    L[i,11]<-(pred_e_v5[i,a1]-his_e[11,j])/(log(pred_e_v5[i,a1]/his_e[11,j]))
    
    L[i,12]<-(pred_e_g1[i,a1]-his_e[12,j])/(log(pred_e_g1[i,a1]/his_e[12,j]))
  }
  
  for (i in 1:nrow(lm1_mcmc_y)){
    fac_f1[i,]=L[i,1]*log(pred_eih1[i,a1]/fh[21,j])
    fac_f2[i,]=L[i,2]*log(pred_eih2[i,a1]/fh[42,j])
    fac_f3[i,]=L[i,3]*log(pred_eih3[i,a1]/fh[63,j])+L[i,7]*log(pred_eiv1[i,a1]/fv[21,j])
    fac_f4[i,]=L[i,4]*log(pred_eih4[i,a1]/fh[84,j])+L[i,8]*log(pred_eiv2[i,a1]/fv[42,j])
    fac_f5[i,]=L[i,5]*log(pred_eih5[i,a1]/fh[105,j])+L[i,9]*log(pred_eiv3[i,a1]/fv[63,j])
    fac_f6[i,]=L[i,6]*log(pred_eih6[i,a1]/fh[126,j])+L[i,10]*log(pred_eiv4[i,a1]/fv[84,j])+L[i,12]*log(pred_eig1[i,a1]/fg[21,j])
    fac_f7[i,]=L[i,11]*log(pred_eiv5[i,a1]/fv[105,j])
    
    fac_s1[i,]=L[i,1]*(log(pred_norm_s1[i,a1]/so[21,j]))
    fac_s2[i,]=L[i,2]*(log(pred_norm_s2[i,a1]/so[42,j]))
    fac_s3[i,]=L[i,3]*(log(pred_norm_s3[i,a1]/so[63,j]))+L[i,7]*(log(pred_norm_v1[i,a1]/vo[21,j]))
    fac_s4[i,]=L[i,4]*(log(pred_norm_s4[i,a1]/so[84,j]))+L[i,8]*(log(pred_norm_v2[i,a1]/vo[42,j]))
    fac_s5[i,]=L[i,5]*(log(pred_norm_s5[i,a1]/so[105,j]))+L[i,9]*(log(pred_norm_v3[i,a1]/vo[63,j]))
    fac_s6[i,]=L[i,6]*(log(pred_norm_s6[i,a1]/so[126,j]))+L[i,10]*(log(pred_norm_v4[i,a1]/vo[84,j]))
    fac_s7[i,]=L[i,11]*(log(pred_norm_v5[i,a1]/vo[105,j]))
    
    fac_f[i,]=fac_f1[i,]+fac_f2[i,]+fac_f3[i,]+fac_f4[i,]+fac_f5[i,]+fac_f6[i,]+fac_f7[i,]
    
    fac_s[i,]=fac_s1[i,]+fac_s2[i,]+fac_s3[i,]+fac_s4[i,]+fac_s5[i,]+fac_s6[i,]+fac_s7[i,]
    
    fac_y[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6])*(log(pred_norm_y1[i,a1]/yo[21,j]))+(L[i,7]++L[i,8]+L[i,9]+L[i,10]+L[i,11])*(log(pred_norm_y2[i,a1]/yo[50,j]))+(L[i,12])*(log(pred_norm_y3[i,a1]/yo[79,j]))
    
    fac_r[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6]+L[i,7]+L[i,8]+L[i,9]+L[i,10]+L[i,11]+L[i,12])*(log(pred_r[i,a1+1]/ro[21,j]))
    
    fac_gdp[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6]+L[i,7]+L[i,8]+L[i,9]+L[i,10]+L[i,11]+L[i,12])*(log(pred_gdp[i,a1+1]/gdp[21,j]))
    
    fac_p[i,]=(L[i,1]+L[i,2]+L[i,3]+L[i,4]+L[i,5]+L[i,6]+L[i,7]+L[i,8]+L[i,9]+L[i,10]+L[i,11]+L[i,12])*(log(pred_p[i,a1+1]/p[21,j]))
  }
  
  for (i in 1: nrow(fac_f)){
    fac_f1_reg[i,j-23]<-fac_f1[i,]
    fac_f2_reg[i,j-23]<-fac_f2[i,]
    fac_f3_reg[i,j-23]<-fac_f3[i,]
    fac_f4_reg[i,j-23]<-fac_f4[i,]
    fac_f5_reg[i,j-23]<-fac_f5[i,]
    fac_f6_reg[i,j-23]<-fac_f6[i,]
    fac_f7_reg[i,j-23]<-fac_f7[i,]
    fac_s1_reg[i,j-23]<-fac_s1[i,]
    fac_s2_reg[i,j-23]<-fac_s2[i,]
    fac_s3_reg[i,j-23]<-fac_s3[i,]
    fac_s4_reg[i,j-23]<-fac_s4[i,]
    fac_s5_reg[i,j-23]<-fac_s5[i,]
    fac_s6_reg[i,j-23]<-fac_s6[i,]
    fac_s7_reg[i,j-23]<-fac_s7[i,]
    fac_f_reg[i,j-23]<-fac_f[i,]
    fac_s_reg[i,j-23]<-fac_s[i,]
    fac_y_reg[i,j-23]<-fac_y[i,]
    fac_gdp_reg[i,j-23]<-fac_gdp[i,]
    fac_p_reg[i,j-23]<-fac_p[i,]
    fac_r_reg[i,j-23]<-fac_r[i,]
  }
}
#changes of ei for other countries and global
pred_ei1_global<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei2_global<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei3_global<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei4_global<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei5_global<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei6_global<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei7_global<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
pred_ei_global<- matrix(0, nrow = nrow(lm1_mcmc_p),ncol=n_2050-1)
for (i in 1:nrow(lm1_mcmc_y)){
  for (t in 1:n_2050-1){
    pred_ei1_global[i,t]<-pred_e1_global[i,t]/pred_o1_global[i,t]*1000
    pred_ei2_global[i,t]<-pred_e2_global[i,t]/pred_o2_global[i,t]*1000
    pred_ei3_global[i,t]<-pred_e3_global[i,t]/pred_o3_global[i,t]*1000
    pred_ei4_global[i,t]<-pred_e4_global[i,t]/pred_o4_global[i,t]*1000
    pred_ei5_global[i,t]<-pred_e5_global[i,t]/pred_o5_global[i,t]*1000
    pred_ei6_global[i,t]<-pred_e6_global[i,t]/pred_o6_global[i,t]*1000
    pred_ei7_global[i,t]<-pred_e7_global[i,t]/pred_o7_global[i,t]*1000
    pred_ei_global[i,t]<-pred_e_global[i,t]/(pred_o1_global[i,t]+pred_o2_global[i,t]+pred_o3_global[i,t]+pred_o4_global[i,t]+pred_o5_global[i,t]+pred_o6_global[i,t]+pred_o7_global[i,t])*1000
  }
}

for (i in 1:nrow(lm1_mcmc_y)){
  cha_ei1[i,18]=(pred_ei1_global[i,55]-his_ei[1,18])/his_ei[1,18]*100
  cha_ei2[i,18]=(pred_ei2_global[i,55]-his_ei[2,18])/his_ei[2,18]*100
  cha_ei3[i,18]=(pred_ei3_global[i,55]-his_ei[3,18])/his_ei[3,18]*100
  cha_ei4[i,18]=(pred_ei4_global[i,55]-his_ei[4,18])/his_ei[4,18]*100
  cha_ei5[i,18]=(pred_ei5_global[i,55]-his_ei[5,18])/his_ei[5,18]*100
  cha_ei6[i,18]=(pred_ei6_global[i,55]-his_ei[6,18])/his_ei[6,18]*100
  cha_ei7[i,18]=(pred_ei7_global[i,55]-his_ei[7,18])/his_ei[7,18]*100
}

ei_global<- matrix(NA, nrow = 55, ncol =24)

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei1_global[,t])}  
uncertain_lower_reg <- apply(pred_ei1_global, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei1_global, MARGIN = 2, quantile, prob = 0.975)
ei_global[,1]<- uncertain_lower_reg
ei_global[,8]<- uncertain_mean_reg
ei_global[,15]<- uncertain_upper_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei2_global[,t])}  
uncertain_lower_reg <- apply(pred_ei2_global, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei2_global, MARGIN = 2, quantile, prob = 0.975)
ei_global[,2]<- uncertain_lower_reg
ei_global[,9]<- uncertain_mean_reg
ei_global[,16]<- uncertain_upper_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei3_global[,t])}  
uncertain_lower_reg <- apply(pred_ei3_global, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei3_global, MARGIN = 2, quantile, prob = 0.975)
ei_global[,3]<- uncertain_lower_reg
ei_global[,10]<- uncertain_mean_reg
ei_global[,17]<- uncertain_upper_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei4_global[,t])}  
uncertain_lower_reg <- apply(pred_ei4_global, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei4_global, MARGIN = 2, quantile, prob = 0.975)
ei_global[,4]<- uncertain_lower_reg
ei_global[,11]<- uncertain_mean_reg
ei_global[,18]<- uncertain_upper_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei5_global[,t])}  
uncertain_lower_reg <- apply(pred_ei5_global, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei5_global, MARGIN = 2, quantile, prob = 0.975)
ei_global[,5]<- uncertain_lower_reg
ei_global[,12]<- uncertain_mean_reg
ei_global[,19]<- uncertain_upper_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei6_global[,t])}  
uncertain_lower_reg <- apply(pred_ei6_global, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei6_global, MARGIN = 2, quantile, prob = 0.975)
ei_global[,6]<- uncertain_lower_reg
ei_global[,13]<- uncertain_mean_reg
ei_global[,20]<- uncertain_upper_reg


for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei7_global[,t])}  
uncertain_lower_reg <- apply(pred_ei7_global, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei7_global, MARGIN = 2, quantile, prob = 0.975)
ei_global[,7]<- uncertain_lower_reg
ei_global[,14]<- uncertain_mean_reg
ei_global[,21]<- uncertain_upper_reg

for (t in 1:n_2050-1){uncertain_mean_reg[t,]<-median(pred_ei_global[,t])}  
uncertain_lower_reg <- apply(pred_ei_global, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(pred_ei_global, MARGIN = 2, quantile, prob = 0.975)
ei_global[,22]<- uncertain_lower_reg
ei_global[,23]<- uncertain_mean_reg
ei_global[,24]<- uncertain_upper_reg

#LMDI for world
for (i in 1: nrow(fac_f)){
  fac_f1_reg[i,18]<-sum(fac_f1_reg[i,1:17])
  fac_f2_reg[i,18]<-sum(fac_f2_reg[i,1:17])
  fac_f3_reg[i,18]<-sum(fac_f3_reg[i,1:17])
  fac_f4_reg[i,18]<-sum(fac_f4_reg[i,1:17])
  fac_f5_reg[i,18]<-sum(fac_f5_reg[i,1:17])
  fac_f6_reg[i,18]<-sum(fac_f6_reg[i,1:17])
  fac_f7_reg[i,18]<-sum(fac_f7_reg[i,1:17])
  
  fac_s1_reg[i,18]<-sum(fac_s1_reg[i,1:17])
  fac_s2_reg[i,18]<-sum(fac_s2_reg[i,1:17])
  fac_s3_reg[i,18]<-sum(fac_s3_reg[i,1:17])
  fac_s4_reg[i,18]<-sum(fac_s4_reg[i,1:17])
  fac_s5_reg[i,18]<-sum(fac_s5_reg[i,1:17])
  fac_s6_reg[i,18]<-sum(fac_s6_reg[i,1:17])
  fac_s7_reg[i,18]<-sum(fac_s7_reg[i,1:17])
  
  fac_f_reg[i,18]<-sum(fac_f_reg[i,1:17])
  fac_s_reg[i,18]<-sum(fac_s_reg[i,1:17])
  fac_y_reg[i,18]<-sum(fac_y_reg[i,1:17])
  fac_gdp_reg[i,18]<-sum(fac_gdp_reg[i,1:17])
  fac_p_reg[i,18]<-sum(fac_p_reg[i,1:17])
  fac_r_reg[i,18]<-sum(fac_r_reg[i,1:17])
}


e_global<- matrix(NA, nrow = n_2050-1,ncol=33)
e_global[,1] <- apply(pred_e_global, MARGIN = 2, quantile, prob = 0.025)
for (t in 1:n_2050-1){e_global[t,2] <- median(pred_e_global[,t])}
e_global[,3] <- apply(pred_e_global, MARGIN = 2, quantile, prob = 0.975)

e_global[,4] <- apply(pred_e_s_global, MARGIN = 2, quantile, prob = 0.025)
e_global[,5] <- apply(pred_e_v_global, MARGIN = 2, quantile, prob = 0.025)
e_global[,6] <- apply(pred_e_g_global, MARGIN = 2, quantile, prob = 0.025)

for (t in 1:n_2050-1){e_global[t,7] <- median(pred_e_s_global[,t])}
for (t in 1:n_2050-1){e_global[t,8] <- median(pred_e_v_global[,t])}
for (t in 1:n_2050-1){e_global[t,9] <- median(pred_e_g_global[,t])}

e_global[,10] <- apply(pred_e_s_global, MARGIN = 2, quantile, prob = 0.975)
e_global[,11] <- apply(pred_e_v_global, MARGIN = 2, quantile, prob = 0.975)
e_global[,12] <- apply(pred_e_g_global, MARGIN = 2, quantile, prob = 0.975)

e_global[,13] <- apply(pred_e1_global, MARGIN = 2, quantile, prob = 0.025)
e_global[,14] <- apply(pred_e2_global, MARGIN = 2, quantile, prob = 0.025)
e_global[,15] <- apply(pred_e3_global, MARGIN = 2, quantile, prob = 0.025)
e_global[,16] <- apply(pred_e4_global, MARGIN = 2, quantile, prob = 0.025)
e_global[,17] <- apply(pred_e5_global, MARGIN = 2, quantile, prob = 0.025)
e_global[,18] <- apply(pred_e6_global, MARGIN = 2, quantile, prob = 0.025)
e_global[,19] <- apply(pred_e7_global, MARGIN = 2, quantile, prob = 0.025)
for (t in 1:n_2050-1){e_global[t,20] <- median(pred_e1_global[,t])}
for (t in 1:n_2050-1){e_global[t,21] <- median(pred_e2_global[,t])}
for (t in 1:n_2050-1){e_global[t,22] <- median(pred_e3_global[,t])}
for (t in 1:n_2050-1){e_global[t,23] <- median(pred_e4_global[,t])}
for (t in 1:n_2050-1){e_global[t,24] <- median(pred_e5_global[,t])}
for (t in 1:n_2050-1){e_global[t,25] <- median(pred_e6_global[,t])}
for (t in 1:n_2050-1){e_global[t,26] <- median(pred_e7_global[,t])}
e_global[,27] <- apply(pred_e1_global, MARGIN = 2, quantile, prob = 0.975)
e_global[,28] <- apply(pred_e2_global, MARGIN = 2, quantile, prob = 0.975)
e_global[,29] <- apply(pred_e3_global, MARGIN = 2, quantile, prob = 0.975)
e_global[,30] <- apply(pred_e4_global, MARGIN = 2, quantile, prob = 0.975)
e_global[,31] <- apply(pred_e5_global, MARGIN = 2, quantile, prob = 0.975)
e_global[,32] <- apply(pred_e6_global, MARGIN = 2, quantile, prob = 0.975)
e_global[,33] <- apply(pred_e7_global, MARGIN = 2, quantile, prob = 0.975)

fac_reg<- matrix(NA, nrow = 60,ncol=18)
fac_reg[1,] <- apply(fac_f1_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[2,] <- apply(fac_f2_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[3,] <- apply(fac_f3_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[4,] <- apply(fac_f4_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[5,] <- apply(fac_f5_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[6,] <- apply(fac_f6_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[7,] <- apply(fac_f7_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[8,] <- apply(fac_s1_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[9,] <- apply(fac_s2_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[10,] <- apply(fac_s3_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[11,] <- apply(fac_s4_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[12,] <- apply(fac_s5_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[13,] <- apply(fac_s6_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[14,] <- apply(fac_s7_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[15,] <- apply(fac_f_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[16,] <- apply(fac_s_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[17,] <- apply(fac_y_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[18,] <- apply(fac_r_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[19,] <- apply(fac_gdp_reg, MARGIN = 2, quantile, prob = 0.025)
fac_reg[20,] <- apply(fac_p_reg, MARGIN = 2, quantile, prob = 0.025)

for (i in 1:18){fac_reg[21,i] <- median(fac_f1_reg[,i])}
for (i in 1:18){fac_reg[22,i] <- median(fac_f2_reg[,i])}
for (i in 1:18){fac_reg[23,i] <- median(fac_f3_reg[,i])}
for (i in 1:18){fac_reg[24,i] <- median(fac_f4_reg[,i])}
for (i in 1:18){fac_reg[25,i] <- median(fac_f5_reg[,i])}
for (i in 1:18){fac_reg[26,i] <- median(fac_f6_reg[,i])}
for (i in 1:18){fac_reg[27,i] <- median(fac_f7_reg[,i])}
for (i in 1:18){fac_reg[28,i] <- median(fac_s1_reg[,i])}
for (i in 1:18){fac_reg[29,i] <- median(fac_s2_reg[,i])}
for (i in 1:18){fac_reg[30,i] <- median(fac_s3_reg[,i])}
for (i in 1:18){fac_reg[31,i] <- median(fac_s4_reg[,i])}
for (i in 1:18){fac_reg[32,i] <- median(fac_s5_reg[,i])}
for (i in 1:18){fac_reg[33,i] <- median(fac_s6_reg[,i])}
for (i in 1:18){fac_reg[34,i] <- median(fac_s7_reg[,i])}
for (i in 1:18){fac_reg[35,i] <- median(fac_f_reg[,i])}
for (i in 1:18){fac_reg[36,i] <- median(fac_s_reg[,i])}
for (i in 1:18){fac_reg[37,i] <- median(fac_y_reg[,i])}
for (i in 1:18){fac_reg[38,i] <- median(fac_r_reg[,i])}
for (i in 1:18){fac_reg[39,i] <- median(fac_gdp_reg[,i])}
for (i in 1:18){fac_reg[40,i] <- median(fac_p_reg[,i])}

fac_reg[41,] <- apply(fac_f1_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[42,] <- apply(fac_f2_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[43,] <- apply(fac_f3_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[44,] <- apply(fac_f4_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[45,] <- apply(fac_f5_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[46,] <- apply(fac_f6_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[47,] <- apply(fac_f7_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[48,] <- apply(fac_s1_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[49,] <- apply(fac_s2_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[50,] <- apply(fac_s3_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[51,] <- apply(fac_s4_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[52,] <- apply(fac_s5_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[53,] <- apply(fac_s6_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[54,] <- apply(fac_s7_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[55,] <- apply(fac_f_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[56,] <- apply(fac_s_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[57,] <- apply(fac_y_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[58,] <- apply(fac_r_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[59,] <- apply(fac_gdp_reg, MARGIN = 2, quantile, prob = 0.975)
fac_reg[60,] <- apply(fac_p_reg, MARGIN = 2, quantile, prob = 0.975)

#prediction interval for change of emission intensity 
uncertain_mean_reg<- matrix(NA, nrow =18,ncol=1)
uncertain_lower_reg<- matrix(NA, nrow =18,ncol=1)
uncertain_upper_reg<- matrix(NA, nrow =18,ncol=1)
cha_ei<- matrix(NA, nrow =18,ncol=21)
for (t in 1:18){uncertain_mean_reg[t,]<-median(cha_ei1[,t])}  
uncertain_lower_reg <- apply(cha_ei1, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(cha_ei1, MARGIN = 2, quantile, prob = 0.975)
cha_ei[,1]<- uncertain_lower_reg
cha_ei[,8]<- uncertain_mean_reg
cha_ei[,15]<- uncertain_upper_reg

for (t in 1:18){uncertain_mean_reg[t,]<-median(cha_ei2[,t])}  
uncertain_lower_reg <- apply(cha_ei2, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(cha_ei2, MARGIN = 2, quantile, prob = 0.975)
cha_ei[,2]<- uncertain_lower_reg
cha_ei[,9]<- uncertain_mean_reg
cha_ei[,16]<- uncertain_upper_reg

for (t in 1:18){uncertain_mean_reg[t,]<-median(cha_ei3[,t])}  
uncertain_lower_reg <- apply(cha_ei3, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(cha_ei3, MARGIN = 2, quantile, prob = 0.975)
cha_ei[,3]<- uncertain_lower_reg
cha_ei[,10]<- uncertain_mean_reg
cha_ei[,17]<- uncertain_upper_reg

for (t in 1:18){uncertain_mean_reg[t,]<-median(cha_ei4[,t])}  
uncertain_lower_reg <- apply(cha_ei4, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(cha_ei4, MARGIN = 2, quantile, prob = 0.975)
cha_ei[,4]<- uncertain_lower_reg
cha_ei[,11]<- uncertain_mean_reg
cha_ei[,18]<- uncertain_upper_reg

for (t in 1:18){uncertain_mean_reg[t,]<-median(cha_ei5[,t])}  
uncertain_lower_reg <- apply(cha_ei5, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(cha_ei5, MARGIN = 2, quantile, prob = 0.975)
cha_ei[,5]<- uncertain_lower_reg
cha_ei[,12]<- uncertain_mean_reg
cha_ei[,19]<- uncertain_upper_reg

for (t in 1:18){uncertain_mean_reg[t,]<-median(cha_ei6[,t])}  
uncertain_lower_reg <- apply(cha_ei6, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(cha_ei6, MARGIN = 2, quantile, prob = 0.975)
cha_ei[,6]<- uncertain_lower_reg
cha_ei[,13]<- uncertain_mean_reg
cha_ei[,20]<- uncertain_upper_reg

for (t in 1:18){uncertain_mean_reg[t,]<-median(cha_ei7[,t])}  
uncertain_lower_reg <- apply(cha_ei7, MARGIN = 2, quantile, prob = 0.025)
uncertain_upper_reg <- apply(cha_ei7, MARGIN = 2, quantile, prob = 0.975)
cha_ei[,7]<- uncertain_lower_reg
cha_ei[,14]<- uncertain_mean_reg
cha_ei[,21]<- uncertain_upper_reg