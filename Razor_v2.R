require(Kendall)
require(hydroTSM)
require(lubridate)
source('juli_ET.R')
source('stage1.R')
source('euler.R')


########################## Long term rainfall
input_file="Tirumungulam_rainfall_modified_v2.csv"
input=read.csv(input_file,header=TRUE)
in_date=as.Date(input[,1],"%Y-%m-%d")
x=zoo(input[,2],in_date)
t=window(x,start=as.Date(input[1,1]))

########################## 2013 field collected rainfall
input_file="Hourly_rain.csv"
input2=read.csv(input_file,header=TRUE)
date_seq=seq(as.POSIXct("2013-09-26 00:00"),by="hour",length.out = 5324)

t2=zoo(input2[7],date_seq)

t_daily=subdaily2daily(t2,FUN=sum)

#############################
##Box 1: Fast flow box
#############################
maxS1=0.062
ini_S1=0.0*maxS1
tw=0.159

#############################
##Box 2: Slow flow box
#############################
maxS2=14.184
ini_S2=0.2*maxS2
tc=1.538
tu=187.987

Se=0.9*maxS2

#############################
##Variables
#############################
Q1f=vector()
ET1=vector()
ET2=vector()
Qw=vector()
Q2u=vector()
S1=vector()
S2=vector()

#############################
##FE parameters
#############################
a=0 #days
b=2000#length(x) #days
h=0.01
time=(b-a)/h

#############################
##Stretch Rainfall
#############################
tp=coredata(x[1:b])#*1000
test=rep(tp,each=1/h)/(1/h)

PET=rep(4.1,1/h)/(1/h)
PET2=rep(PET,b)

i=1

for(i in 1:(time)) {
  # cur_date=index(x[i])
  # cur_month=month(cur_date)
  # cur_jET=juli_ET(cur_month)*h
  
  cur_PET=PET2[i]
  
  cur_P=test[i]
  
  if (i>1){
    cur_S1=S1[i-1]
  } else {cur_S1=ini_S1}
  
  Q1f[i]=(cur_S1-maxS1)/h
  if (Q1f[i]<0){Q1f[i]=0}
  Qw[i]=cur_S1/tw
  ET1[i]=cur_PET*(cur_S1/maxS1)
  
  S1[i]=cur_S1+h*(cur_P-Q1f[i]-Qw[i]-ET1[i])
  
  if (i>1){
    cur_S2=S2[i-1]
  } else {cur_S2=ini_S2}
  
  Q2u[i]=cur_S2/tu
  
  ET2[i]=cur_PET*(cur_S2/Se)
  
  S2[i]=cur_S2+h*(Qw[i]-Q2u[i]-ET2[i])
  
}

tp1=matrix(Q1f,nrow=1/h,ncol = b,byrow = FALSE)
daily_Q1f=apply(tp1, 2, FUN=sum)
plot(daily_Q1f,type = 'l')

tp2=matrix(S1,nrow=1/h,ncol = b,byrow = FALSE)
daily_S1=apply(tp2, 2, FUN=sum)
#plot(daily_S1,type = 'l')

tp3=matrix(S2,nrow=1/h,ncol = b,byrow = FALSE)
daily_S2=apply(tp3, 2, FUN=mean)
plot(daily_S2,type = 'l')

plot(tp,type='l')
