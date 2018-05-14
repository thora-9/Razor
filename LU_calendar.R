LU_calendar <- function (LU_details){
  cal_out=vector()
  a=1
  for (a in 1:nrow(LU_details)){
    cur_LU=LU_details[[a,'LU']]
    plant_date=LU_details[[a,'plant_month']]
    if(plant_date>0){
    grow_start=as.Date(paste(plant_date,"-15",sep = ""),'%m-%d')
    doy=yday(grow_start)
    grow_end=doy+LU_details$LI[a]+LU_details$LD[a]+LU_details$LM[a]+LU_details$LL[a]
    #Yearly calendar
    cal1=rep(0,doy)
    cal2=rep(1,abs(doy-grow_end))
    if((365-grow_end)>0){
      cal3=rep(0,365-grow_end)
    } else {
      tmp1=grow_end-365
      cal1[1:tmp1]=1
      cal2=cal2[1:(length(cal2)-tmp1)]
      cal3=NULL}
    cal_all=c(cal1,cal2,cal3)
    } else {cal_all=rep(0,365)}
    
    #Transpose individual calendars
    cal_out=cbind(cal_out,cal_all)
    colnames(cal_out)[a]=cur_LU
  } 
  return(cal_out)
  
}
