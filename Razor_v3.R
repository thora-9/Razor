require(Kendall)
require(hydroTSM)
require(lubridate)
source('juli_ET.R')
source('juli_ET_v2.R')
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

t_daily=subdaily2daily(t2,FUN=sum)*1000
#t_daily[1:length(t_daily)]=1

#############################
##Box 1: Fast flow box
#############################
maxS1=3
ini_S1=0.0*maxS1
tw=0.08

#############################
##Box 2: Slow flow box
#############################
maxS2=300
ini_S2=0.0*maxS2
tc=5
tu=100

Se=maxS2

ini_Sc1=0

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
Qf=vector()
Qu=vector()
Sc1=vector()
Sc2=vector()

#Tank variables
inflow_f1=vector()
inflow_s1=vector()
t1_area=vector()
t1_area0=0
t1_vol=vector()
t1_stage=vector()
t1_area=vector()
t1_spill=vector()
t1_sluice=vector()
t1_GW=vector()
t1_ET=vector()
t1_const=cbind(5e6,3.595,30,276405)
colnames(t1_const)=c("max_catch","weir_height","spill_len","max_volume")#Units c(m2,meter,meter,m3)

#############################
##FE parameters
#############################
a=0 #days
b=length(t_daily)#2000#length(x) #days
h=0.01
time=(b-a)/h

#############################
##Stretch Rainfall
#############################
tp=coredata(t_daily[1:b])#*1000
#tp=coredata(x[1:b])#*1000
test=rep(tp,each=1/h)/(1/h)
dates=rep(index(t_daily),each=1/h)
months=month(dates)


PET=rep(9.1,1/h)/(1/h)
PET2=rep(PET,b)

i=100

for(i in 1:(time)) {
  #print(i)
  #cur_date=dates[i]
  cur_month=months[i]
  # cur_date=index(x[i])
  # cur_month=month(cur_date)
  # cur_jET=juli_ET(cur_month)*h
  cur_PET=juli_ET(cur_month)*h
  
  cur_P=test[i]
#Bucket 1  
  if (i>1){
    cur_S1=S1[i-1]
  } else {cur_S1=ini_S1}
  
  Q1f[i]=(cur_S1-maxS1)
  if (Q1f[i]<0){Q1f[i]=0}
  Qw[i]=(cur_S1/tw)*h #units of tw are days, so reduce that to timestep of FE
  ET1[i]=cur_PET*(cur_S1/maxS1)*0
  
  S1[i]=cur_S1+(cur_P-Q1f[i]-Qw[i]-ET1[i])
#Bucker 2  
  if (i>1){
    cur_S2=S2[i-1]
  } else {cur_S2=ini_S2}
  
  Q2u[i]=(cur_S2/tu)*h #units of tu are days, so reduce that to timestep of FE
  
  ET2[i]=cur_PET*(cur_S2/Se)
  
  S2[i]=cur_S2+(Qw[i]-Q2u[i]-ET2[i])
# #Bucket 3  
#   if (i>1){
#     cur_Sc1=Sc1[i-1]
#   } else {cur_Sc1=ini_Sc1}
#   
#   Qf[i]=(cur_Sc1/tc)*h
#   
#   Sc1[i]=cur_Sc1+h*(Q1f[i]-Qf[i])
#   
#i=i+1  
  #print(i)

#Start the tank water balance every after every (1/h) steps  
  if (i%%(1/h)==0){
    j=i*h
    inflow_f1[j]=sum(Q1f[(i-(1/h)):i])
    inflow_s1[j]=sum(Q2u[(i-(1/h)):i])
    inflow1[j]=inflow_f1[j]+inflow_s1[j]
    
    if(j==1){
      t1_area_cur=t1_area0
    } else {t1_area_cur=t1_area[j-1]}
    #Rainfall directly on the tank in that timestep
    rain_t1=t1_area_cur*t_daily[j]*(1/1000) #in meters
    
    t1_vol[j]=rain_t1+((inflow1[j])*(t1_const[1]-t1_area_cur)*(1/1000))
    
    #Stage-Volume relationship from Mike data
    t1_stage[j]=(t1_vol[j]/22914)^(1/1.9461)
    #Stage Area 
    t1_area[j]=42942*(t1_stage[j])^1.0993
    
    # ####Spillway taken from Sakti paper modified to take absolute value of difference
    # if(tank1_const[2]-stage1<0){
    #   form_spill1=1.7*tank1_const[3]*(abs(tank1_const[2]-stage1))^(3/2)*(24*60*60)
    # } else{form_spill1=0}
    
    ####ET from tanks
    t1_ET[j]=t1_area[j]*juli_ET(cur_month)*(1/1000) #m3/day
    
    ####Compare the tank capacity and current volume of water in tank.
    vol_diff1=t1_vol[j]-tank1_const[4]
    if (vol_diff1>=0){
      t1_spill[j]=vol_diff1
    } else{t1_spill[j]=0}
    
    
    ####Sluice Outflow
    Qo1a = ((t1_stage[j]-.785)*5.1903)*86.4 #86.4 Converst L/s to m3/d
    Qo1b = (((t1_stage[j]-1.185)*9.6768)+((t1_stage[j]-1.185)*4.9196))*86.4 #86.4 Converst L/s to m3/d
    
    if (Qo1a<0){Qo1a = 0}
    if (Qo1b<0) {Qo1b = 0}
    
    if (t1_stage[j]<0.785){
     t1_sluice[j] = 0} else if(0.785<t1_stage[j] & t1_stage[j]<1.185) {
     t1_sluice[j] = Qo1a} else if(t1_stage[j]>1.185){
     t1_sluice[j] = (Qo1a + Qo1b)}
    
    ###Spillage from upstream tank
    spill_add1=0
    
    ####GW exchange- Mike paper
    GW_loss1=8.6*t1_stage[j]-6.5 #mm/day
    t1_GW[j]=(GW_loss1/1000)*t1_area[j]#m3/day
    
    #Total Storage change
    t1_vol[j]=t1_vol[j]-(t1_ET[j]+t1_sluice[j]+t1_spill[j]+t1_GW[j])
    
    #Stage-Volume relationship from Mike data
    t1_stage[j]=(t1_vol[j]/22914)^(1/1.9461)
    #Stage Area 
    t1_area[j]=42942*(t1_stage[j])^1.0993
    
  }

}


i=1
tp1=matrix(Q1f,nrow=1/h,ncol = b,byrow = FALSE)
daily_Q1f=apply(tp1, 2, FUN=sum)


tp2=matrix(ET2,nrow=1/h,ncol = b,byrow = FALSE)
daily_ET2=apply(tp2, 2, FUN=sum)
#plot(daily_S1,type = 'l')

tp3=matrix(S2,nrow=1/h,ncol = b,byrow = FALSE)
daily_S2=apply(tp3, 2, FUN=mean)

tp4=matrix(Qw,nrow=1/h,ncol = b,byrow = FALSE)
daily_Qw=apply(tp4, 2, FUN=mean)

par(mfrow=c(2,2))

plot(Qr,type='l')
plot(daily_S2,type = 'l')
#plot(daily_S1,type = 'l')

plot(daily_Q1f,type = 'l')
plot(t_daily,type='l')


#sum(Q1f)
sum(daily_Q1f)
sum(Qr,na.rm = TRUE)

