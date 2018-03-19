stage1<-function(maxS1,S1,tw,cur_jET,delta_t,P){
  Q1f=(S1-maxS1)/delta_t
  Qw=S1/tw
  ET1=cur_jET*(S1/maxS1)
  
  del_S1=P-Q1f-Qw-ET1
  
  return(del_S1)
}

