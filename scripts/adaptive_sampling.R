library(here)
library(readxl)
library(moments)
# nc: size of cluster (mSy)
hh_acs<- function(y, nc){
  if(length(y) != length(nc)) stop('y != cluster size')
  return(mean(y/nc))
}  

ht_acs<- function(N, n, y, nc){
  if(length(y) != length(nc)) stop('y != cluster size')
  k<- length(y)
  ak<- 1 - (choose(N-nc,n)/choose(N,n))
  return((1/N) * sum(y/ak))
}

rfp1<- function(x, y, nc, a, b, eta, theta){ #w1= alpha_s w2= beta_s
  
  o1<- 1 + cvy^2 + cvx^2 * lambda1 * (lambda1-2*deltaS) 
  o2<- mu_y^(-2) * (1 + lambda1^2 * ((2*gam-1) / 2) * cvx^2)
  o3<- mu_y^(-1) * (1 + lambda1 * (gam-1) * ((gam/2)+deltaS-1) * cvx^2)
  o4<- 1 + (lambda1/8) * (3*lambda1 - 4*deltaS) * cvx^2
  o5<- mu_y^(-1) * (1 + lambda1^2 * ((4*gam^2 - 8*gam + 3)/8) * cvx^2)
  
  w1<<- (o2*o4-o3*o5)/(o1*o2-o3^2)
  w2<<- (o1*o5-o3*o4)/(o1*o2-o3^2)
  
  mrfp1<- (w1 * hh_acs(y,nc) + w2*((a*hh_acs(x,nc) + b)/(a*mu_x + b))^eta) * exp((a*(mu_x-hh_acs(x,nc)))/(a*(mu_x+hh_acs(x,nc))+2*b))
  return(mrfp1)
}

rfp2<- function(N, n, x, y, nc, a, b, eta, theta){ #w1= alpha_s w2= beta_s
  
  q1<- 1 + cvy^2 - 4*eta*lambda2*rS*cvx*cvy + eta*(2*eta+1)*lambda2^2*cvx^2
  q2<- mu_y^(-2) * (1+lambda2^2*gam*(2*gam-1)*cvx^2)
  q3<- mu_y^(-1) * (1+lambda2^2 * (((eta+gam)^2)/2 - (1/8))*cvx^2 - lambda2*(eta+gam-(1/2))*rS*cvx*cvy)
  q4<- 1 + eta*lambda2 * (((eta+1)/4) * lambda2*cvx^2 - rS*cvx*cvy)
  q5<- mu_y^(-1) * (1+lambda2^2*((4*gam^2 - 1)/8)*cvx^2)
  
  w1<<- (q2*q4-q3*q5)/(q1*q2-q3^2)
  w2<<- (q1*q5-q3*q4)/(q1*q2-q3^2)
  
  mrfp2<- (w1 * ht_acs(N,n,y,nc) + w2*exp(theta*((a*(mu_x-ht_acs(N,n,x,nc)))/(a*(mu_x+ht_acs(N,n,x,nc))+2*b)))) * ((a*mu_x+b)/(a*ht_acs(N,n,x,nc)+b))^eta
  return(mrfp2)
}

y_srs<- function(N, n, y, x, t){
  mu<- mean(y)
  pu<- mean(x)
  r<- cor(y,x)
  ind<- 1:length(y)
  my<- mx<- beta<- rep(NA, t)
  
  for (i in 1:t) {
    smp<- sample(ind, n)
    smp_y<- y[smp]
    my[i]<- mean(smp_y)
    smp_x<- x[smp]
    mx[i]<- mean(smp_x)
    beta[i]<- lm(smp_y ~ smp_x)$coeff[[2]]
  }
  # estimate
  y_r<- (my/mx)*pu
  y_R<- mean(na.omit(y_r))
  y_p<- (my*mx)/pu
  y_P<- mean(na.omit(y_p))
  y_lr<- my + beta *(pu- mx)
  y_LR<- mean(na.omit(y_lr))
  y_re<- my * exp((pu-mx)/(pu+mx))
  y_RE<- mean(na.omit(y_re))
  y_s1e<- (my + beta * (pu-mx)) * exp((pu-mx)/(pu+mx))
  y_S1e<- mean(na.omit(y_s1e))
  # Bias
  bias_y<- mean(na.omit(my-mu))
  bias_R<- mean(na.omit(y_r-mu))
  bias_P<- mean(na.omit(y_p-mu))
  bias_LR<- mean(na.omit(y_lr-mu))
  bias_RE<- mean(na.omit(y_re-mu))
  bias_S1e<- mean(na.omit(y_s1e-mu))
  # MSE
  mse_y<- mean(na.omit((my-mu)^2))
  mse_R<- mean(na.omit((y_r-mu)^2))
  mse_P<- mean(na.omit((y_p-mu)^2))
  mse_LR<- mean(na.omit((y_lr-mu)^2))
  mse_RE<- mean(na.omit((y_re-mu)^2))
  mse_S1e<- mean(na.omit((y_s1e-mu)^2))
  
  return(list(Population = c(mu=mu, P=pu, r=r),
              Estimate = c(y_bar=mean(my), y_R=y_R, y_P=y_P, y_LR=y_LR, y_RE=y_RE, y_S1e=y_S1e),
              Bias = c(bias_y=bias_y, bias_R=bias_R, bias_P=bias_P, bias_lr=bias_LR, bias_RE=bias_RE, bias_S1e=bias_S1e),
              MSE = c(MSE_y=mse_y, MSE_R=mse_R, MSE_P=mse_P, MSE_lr=mse_LR, MSE_RE=mse_RE, MSE_S1e=mse_S1e)
  ))
}


## read data
data<- read_xlsx(here('data/data_covid.xlsx'))

## define population parameters
N<- nrow(data)
y<- data$y
x<- data$x
mu_y<- mean(data$y)
mu_x<- mean(data$x)
cvy<- sd(data$ySy)/mean(data$ySy)
cvx<- sd(data$xSy)/mean(data$xSy)
rS<- cor(data$ySy, data$xSy)
beta1<- lm(data$ySy ~ data$xSy)$coeff[[2]] 
h<- skewness(data$xSy)
k<- kurtosis(data$xSy)
deltaS<- rS*cvy/cvx
a1<- cvx; b1<- N*mu_x 
a2<- cvx; b2<- N*mu_x
lambda1<- (a1*mu_x)/(a1*mu_x+b1)
lambda2<- (a2*mu_x)/(a2*mu_x+b2)
## sample size
n<- c(2, 5, 8, 12, 28, 42, 88)

ind<- 1:nrow(data)
# variables: estimate
en<- rep(NA, length(n))
hh_y<- ht_y<- rep(NA, length(n))
dchh<- dcht<- rep(NA, length(n))
chlr_hh<- chd_hh<- rep(NA, length(n))
yer3<- rep(NA, length(n))
y_bar<- rep(NA, length(n))
mrfp1<- mrfp2<- rep(NA, length(n))
# variables: Bias
bias_hh<- bias_ht<- rep(NA, length(n))
bias_srs<- rep(NA, length(n))
bias_dchh<- bias_dcht<- rep(NA, length(n))
bias_chlr_hh<- bias_chd_hh<- rep(NA, length(n))
bias_yer3<- rep(NA, length(n))
bias_mrfp1<- bias_mrfp2<- rep(NA, length(n))
# variables: MSE
mse_hh<- mse_ht<- rep(NA, length(n))
mse_srs<- rep(NA, length(n))
mse_dchh<- mse_dcht<- rep(NA, length(n))
mse_chlr_hh<- mse_chd_hh<- rep(NA, length(n))
mse_yer3<- rep(NA, length(n))
mse_mrfp1<- mse_mrfp2<- rep(NA, length(n))
# variables PRE
pre_srs<- rep(NA, length(n))
pre_hh<- pre_ht<-  rep(NA, length(n))
pre_dchh<- pre_dcht<- rep(NA, length(n))
pre_chlr_hh<- pre_chd_hh<- rep(NA, length(n))
pre_yer3<- rep(NA, length(n))
pre_mrfp1<- pre_mrfp2<-  rep(NA, length(n))

# define loop
for(j in 1:length(n)){
  hhy<- hty<-c()
  hhx<- htx<- c()
  mdc_hh<- mdc_ht<- c()
  mch_hh<- mch_d<- c()
  yer3_hh<- c()
  mrfp1_hh<- mrfp2_ht<- c()
  ns<- c()
  y_wor<- c()
  gam<- (1/n[j] - 1/N) 
  # y
  for(i in 1:10000){
    smp<- sample(ind, n[j])
    smp_y<- data$ySy[smp]
    smp_x<- data$xSy[smp]
    mc_y<- data$mSy[smp]
    nc_y<- data$nSy[smp]
    cy<- data$c_y[smp]
    ys<- smp_y[!duplicated(cy)]
    xs<- smp_x[!duplicated(cy)]
    ms<- mc_y[!duplicated(cy)]
    fns<- nc_y[!duplicated(cy)]
    hhy[i]<- hh_acs(smp_y, mc_y)
    hty[i]<- ht_acs(N, n[j], ys, ms)
    hhx[i]<- hh_acs(smp_x, mc_y)
    htx[i]<- ht_acs(N, n[j], xs, ms)
    mdc_hh[i]<- (hhy[i]/hhx[i])*mu_x
    mdc_ht[i]<- (hty[i]/htx[i])*mu_x
    mch_hh[i]<- hhy[i]+beta1*(mu_x-hhx[i])
    mch_d[i]<- hhy[i]+(mu_x-hhx[i])
    yer3_hh[i]<- hhy[i]*((mu_x*k+h)/(hhx[i]*k+h))
    mrfp1_hh[i]<- rfp1(y=smp_y, x=smp_x, nc=mc_y, a=a1, b=b1, eta=-1, theta=0)
    mrfp2_ht[i]<- rfp2(N=N, n=n[j], x=xs, y=ys, nc=ms, a=a2, b=b2, eta=1, theta=-1)
    ns[i]<- sum(fns)
    y_wor[i]<- mean(sample(data$y, ns[i]))
  }
  # Estimate
  en[j]<- mean(ns)
  hh_y[j]<- mean(hhy)
  ht_y[j]<- mean(hty)
  y_bar[j]<- mean(y_wor)
  dchh[j]<- mean(na.omit(mdc_hh))
  dcht[j]<- mean(na.omit(mdc_ht))
  chlr_hh[j]<- mean(na.omit(mch_hh))
  chd_hh[j]<- mean(na.omit(mch_d))
  yer3[j]<- mean(yer3_hh)
  mrfp1[j]<- mean(mrfp1_hh)
  mrfp2[j]<- mean(mrfp2_ht)
  # Bias
  bias_hh[j]<- mean(hhy - mu_y)
  bias_ht[j]<- mean(hty - mu_y)
  bias_srs[j]<- mean(y_wor - mu_y)
  bias_dchh[j]<- mean(na.omit(mdc_hh - mu_y))
  bias_dcht[j]<- mean(na.omit(mdc_ht - mu_y))
  bias_chlr_hh[j]<- mean(na.omit(mch_hh - mu_y))
  bias_chd_hh[j]<- mean(na.omit(mch_d - mu_y))
  bias_yer3[j]<- mean(yer3_hh-mu_y)
  bias_mrfp1[j]<- mean(mrfp1_hh-mu_y)
  bias_mrfp2[j]<- mean(mrfp2_ht-mu_y)
  # MSE
  mse_hh[j]<- mean((hhy - mu_y)^2)
  mse_ht[j]<- mean((hty - mu_y)^2)
  mse_srs[j]<- mean((y_wor - mu_y)^2)
  mse_dchh[j]<- mean(na.omit((mdc_hh - mu_y)^2))
  mse_dcht[j]<- mean(na.omit((mdc_ht - mu_y)^2))
  mse_chlr_hh[j]<- mean(na.omit((mch_hh - mu_y)^2))
  mse_chd_hh[j]<- mean(na.omit((mch_d - mu_y)^2))
  mse_yer3[j]<- mean((yer3_hh-mu_y)^2)
  mse_mrfp1[j]<- mean((mrfp1_hh-mu_y)^2)
  mse_mrfp2[j]<- mean((mrfp2_ht-mu_y)^2)
}

# y_srs
n1<- ceiling(en)
for(i in 1:length(n1)){
  assign(paste0('srs', i), y_srs(N=N, n=n1[i], y=data$y, x=data$x, t=10000))
}

srs_est<- rbind(srs1$Estimate, srs2$Estimate, srs3$Estimate, srs4$Estimate, srs5$Estimate, srs6$Estimate, srs7$Estimate)
srs_bias<- rbind(srs1$Bias, srs2$Bias, srs3$Bias, srs4$Bias, srs5$Bias, srs6$Bias, srs7$Bias)
srs_mse<- rbind(srs1$MSE, srs2$MSE, srs3$MSE, srs4$MSE, srs5$MSE, srs6$MSE, srs7$MSE)
srs_pre<- round((srs_mse[,'MSE_y']/srs_mse)*100, 1)
colnames(srs_pre)<- c('PRE_y', 'PRE_R', 'PRE_P', 'PRE_lr', 'PRE_RE', 'PRE_S1e')

msey<- srs_mse[,'MSE_y']
# PRE
for(j in 1:length(n)) {
  pre_srs[j]<- round(msey[j]/msey[j]*100, 1)
  pre_hh[j]<- round((msey[j]/mse_hh[j])*100, 1)
  pre_ht[j]<- round((msey[j]/mse_ht[j])*100, 1)
  pre_dchh[j]<- round((msey[j]/mse_dchh[j])*100, 1)
  pre_dcht[j]<- round((msey[j]/mse_dcht[j])*100, 1)
  pre_chlr_hh[j]<- round((msey[j]/mse_chlr_hh[j])*100, 1)
  pre_chd_hh[j]<- round((msey[j]/mse_chd_hh[j])*100, 1)
  pre_yer3[j]<- round((msey[j]/mse_yer3[j])*100, 1)
  pre_mrfp1[j]<- round((msey[j]/mse_mrfp1[j])*100, 1)
  pre_mrfp2[j]<- round((msey[j]/mse_mrfp2[j])*100, 1)
}


# result
result<- list(Estimate=(cbind(n0=n, 'E(n)'=n1, Ybar=srs_est[,'y_bar'], HHy=hh_y , HTy=ht_y, 
                              YR=srs_est[,'y_R'], YP=srs_est[,'y_P'], 
                              YLR=srs_est[,'y_LR'], YRE=srs_est[,'y_RE'], YS1e=srs_est[,'y_S1e'],
                              DC_HH=dchh, DC_HT=dcht, CHlr_HH=chlr_hh, CHd_HH=chd_hh, 
                              YER3=yer3, RFp1=mrfp1, RFp2=mrfp2)),
              Bias=cbind(n0=n, 'E(n)'=n1, SRS=srs_bias[,'bias_y'], HH=bias_hh, HT=bias_ht, 
                         R=srs_bias[,'bias_R'], P=srs_bias[,'bias_P'], LR=srs_bias[,'bias_lr'], RE=srs_bias[,'bias_RE'], S1e=srs_bias[,'bias_S1e'],
                         DC_HH=bias_dchh, DC_HT=bias_dcht, CHlr_HH=bias_chlr_hh, CHd_HH=bias_chd_hh, 
                         YER3=bias_yer3, RFp1=bias_mrfp1, RFp2=bias_mrfp2),
              MSE=cbind(n0=n, 'E(n)'=n1, SRS=srs_mse[,'MSE_y'], HH=mse_hh, HT=mse_ht, 
                        R=srs_mse[,'MSE_R'], P=srs_mse[,'MSE_P'], LR=srs_mse[,'MSE_lr'], RE=srs_mse[,'MSE_RE'], S1e=srs_mse[,'MSE_S1e'],
                        DC_HH=mse_dchh, DC_HT=mse_dcht, CHlr_HH=mse_chlr_hh, CHd_HH=mse_chd_hh, 
                        YER3=mse_yer3, RFp1=mse_mrfp1, RFp2=mse_mrfp2),
              PRE=cbind(n0=n, 'E(n)'=n1, SRS=pre_srs, HH=pre_hh, HT=pre_ht, 
                        R=srs_pre[,'PRE_R'], P=srs_pre[,'PRE_P'], LR=srs_pre[,'PRE_lr'], RE=srs_pre[,'PRE_RE'], S1e=srs_pre[,'PRE_S1e'],
                        DC_HH=pre_dchh, DC_HT=pre_dcht, CHlr_HH=pre_chlr_hh, CHd_HH=pre_chd_hh, 
                        YER3=pre_yer3, RFp1=pre_mrfp1, RFp2=pre_mrfp2)
)
result
estimate<- data.frame(result$Estimate)
bias<- data.frame(result$Bias)
mse<- data.frame(result$MSE)
pre<- data.frame(result$PRE)
result1<- list(estimate, bias, mse, pre)
