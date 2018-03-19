require(Kendall)
require(hydroTSM)
require(lubridate)
source('juli_ET.R')


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
S1=vector()
S1[1]=0.30*maxS1
tw=0.189

#############################
##Box 2: Fast flow box
#############################
maxS2=326.358
S2=vector()
S2[1]=0.3*maxS2
tc=1.538

maxSe=159.756

#############################
##Variables
#############################
runoff_S1=vector()
AET_S1=vector()

time_period=365*3
i=2
for(i in 1:time_period){
  cur_date=index(x[i])
  cur_month=month(cur_date)
  cur_jET=juli_ET(cur_month)
  
  cur_P=20#coredata(x[i])
  
  if (S1[i]-maxS1>0){
    runoff_S1[i]=S1[i]-maxS1
  } else if (S1[i]-maxS1<0){
    runoff_S1[i]=0
  }
  
  temp_S1=S1[i]-runoff_S1[i]
  
  infil_S1=S1[i]/tw
  
  ET_S1=cur_jET*(S1[i]/maxS1)
  
  AET_S1[i]=ET_S1
  
  change_S1=cur_P-runoff_S1[i]-infil_S1-AET_S1
  
  if (S1[i]+change_S1>0){
    S1[i+1]=S1[i]-change_S1
  } else if (S1[i]+change_S1<0){
    S1[i+1]=0
  }

}

