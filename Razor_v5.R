require(Kendall)
require(hydroTSM)
require(lubridate)
require(tibble)
source('juli_ET.R')
source('juli_ET_v2.R')
source('stage1.R')
source('euler.R')
source('SCS_curve.R')
source('LU_calendar.R')


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
#t_daily[1:length(t_daily)]=10

#############################
##Box 1: Runoff generation
#############################
LU_details=tribble(
              ~LU,~Irrigation,~Area,~plant_month,~LI,~LD,~LM,~LL,~KCL,~KCLP,~KCM,~KCL,~RD,~CN,
              'Rice',NA,50,10,20,30,30,25,1,1.15,1.20,0.9,300,NA,
              'Kharif_Millet',NA,20,10,20,35,40,30,0.3,0.3,1.2,0.5,300,70,
              'Kharif_Fallow',NA,0,0,0,0,0,0,0,0,0,0,0,58,
              'Juliflora',NA,0,0,0,0,0,0,0,0,0,0,0,66,
              'Cotton',NA,0,0,0,0,0,0,0,0,0,0,0,66,
              'Rabi_Fallow',NA,40,0,0,0,0,0,0,0,0,0,0,66,
              'Summer_Fallow',NA,100,0,0,0,0,0,0,0,0,0,0,66)


LU_calen=LU_calendar(LU_details)

HRU_catch1=tribble(
            ~Percent,~Calendar,
            50,c(1,0,0,0,0,0,0,0,0,1,1,1),
            50,c(0,0,0,0,0,0,0,0,0,0,0,0))

CN_List=vector(mode = 'list',length = 3)
names(CN_List) <- c('0','2','3')
CN_List[1:3]=c(70,80,90)
#LU codes: Fallow=0; Rice=1; Juliflora=2

AMC=0

#Using curve number to get the maximum soil water retention
HRU1_CN=70
HRU2_CN=70 #
HRU3_CN=70 #
#This basically converts the orignal curve number formula to mm
S=(25400/CN)-254

#############################
##Box 2: Soil Moisture box
#############################
RD=300
theta_wp=RD*0.1725
theta_fc=RD*0.2825
theta_sat=RD*0.415
ini_S2=theta_wp

#For now, anything greater than theta_fc is routed straight to s3

#############################
##Box 3: Groundwater box
#############################
maxS3=600
ini_S3=500

#############################
##Variables
#############################
Q1f=vector()
ET1=vector()
ET2=vector()
Qw=vector()
Qr=vector()
Qx=vector()
S1=vector()
S2=vector()
S3=vector()
Qf=vector()
Qu=vector()
Sc1=vector()
Sc2=vector()
runoff=matrix(,nrow = nrow(HRU_catch1),ncol=length(t_daily))

#Tank variables
inflow_f1=vector()
inflow_s1=vector()
t1_area=vector()
t1_area0=0
t1_vol=vector()
t1_vol0=0
t1_stage=vector()
t1_area=vector()
t1_spill=vector()
t1_sluice=vector()
t1_GW=vector()
t1_ET=vector()
t1_all=data.frame(matrix(ncol = 8, nrow = 0))
t1_const=as.data.frame(cbind(5e6,3.595,30,276405))
colnames(t1_const)=c("max_catch","weir_height","spill_len","max_volume")#Units c(m2,meter,meter,m3)


#############################
#Initialize
#############################
samay=length(t_daily)
i=1

for(i in 1:(time)) {
  
  cur_date=t_daily[i]
  cur_month=months[i]
  cur_P=test[i]
  
  j=1
  #Runoff Generation
  #Estimate the antecendant conditions of the catchment
  if(i>5){
    rain_5=sum(t_daily[(i-5):(i-1)])*1000
  } else {rain_5=30}
  
  for (j in 1:nrow(HRU_catch1)){
    cur_LU=HRU_catch1[[j,2]][cur_month]
    if(cur_LU!=1){
      runoff[j,i]=SCS_curve(cur_LU,cur_P,rain_5,CN_List)
    } else {runoff[j,i]=0}
  }
  
  #Estimating the Soil Moisture Balance for each HRU
  
  
  
  if (cur_P==0){}
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
  
  ET2[i]=cur_PET*(cur_S2/Se)
  S2[i]=cur_S2+(Qw[i])#-ET2[i])
  if (S2[i]>theta_fc){
    Qr[i]=S2[i]-theta_fc
  } else {Qr[i]=0}
  
  #Bucket 3
  if (i>1){
    cur_S3=S3[i-1]
  } else {cur_S3=ini_S3}
  
  Qx[i]=ET2[i]
  
  S3[i]=cur_S3+(Qr[i]-Qx[i])
  
  #i=i+1
  #print(i)
  
  #Start the tank water balance every after every (1/h) steps  
  if (i%%(1/h)==0){
    j=i*h
    inflow_f1[j]=sum(Q1f[(i-(1/h)):i])
    inflow_s1[j]=0#sum(Q2u[(i-(1/h)):i])*0
    inflow1[j]=inflow_f1[j]+inflow_s1[j]
    
    if(j==1){
      t1_area_cur=t1_area0
    } else {t1_area_cur=t1_area[j-1]}
    
    if(j==1){
      t1_vol_cur=t1_vol0
    } else {t1_vol_cur=t1_vol[j-1]}
    #Rainfall directly on the tank in that timestep
    rain_t1=t1_area_cur*t_daily[j]*(1/1000) #in meters
    
    inflow_vol=rain_t1+((inflow1[j])*(t1_const[1]-t1_area_cur)*(1/1000))
    t1_vol[j]=t1_vol_cur+rain_t1+inflow_vol#((inflow1[j])*(t1_const[1]-t1_area_cur)*(1/1000))
    
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
    
    cur_all=cbind(t1_stage[j],t1_area[j],t1_vol[j],inflow1[j],coredata(inflow_vol),t1_ET[j],t1_sluice[j],t1_spill[j],t1_GW[j])
    t1_all=rbind(t1_all,cur_all)
  }
  
}




