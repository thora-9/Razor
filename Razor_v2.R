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

#############################
##Box 1: Fast flow box
#############################
maxS1=0.069
ini_S1=0.0*maxS1
tw=0.189

#############################
##Box 2: Fast flow box
#############################
maxS2=326.358
ini_S2=0.3*maxS2
tc=1.538

maxSe=159.756

#############################
##Variables
#############################
Q1f=vector()
AET_S1=vector()
Qw=vector()
S1=vector()
S2=vector()

#############################
##FE parameters
#############################
a=0 #days
b=1000 #days
h=0.01
time=(b-a)/h

#############################
##Stretch Rainfall
#############################
tp=coredata(x[1:b])
test=rep(tp,each=1/h)/(1/h)

PET=rep(4.1,1/h)/(1/h)
PET2=rep(PET,b)

i=1

for(i in 1:(time-1)) {
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
  AET_S1[i]=cur_PET*(cur_S1/maxS1)
  
  S1[i]=cur_S1+h*(cur_P-Q1f[i]-Qw[i]-AET_S1[i])
  i=i+1
}
  