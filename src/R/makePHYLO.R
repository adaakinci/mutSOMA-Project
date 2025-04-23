Sample_info_csv <- "/Users/adaakinci/Desktop/sample_info.csv"
sample_info <- read.csv(Sample_info_csv, header=T, sep=",")

makePHYLO<-function(tall, pedigree, sample.info = sample_info){

  t13<-330
  t14<-330
  s13<-sample_info[which(sample_info[,"Stem"] == 13),]
  s13[,"Branchpoint_date"]<-t13- s13[,"Branchpoint_date"]
  s14<-sample_info[which(sample_info[,"Stem"] == 14),]
  s14[,"Branchpoint_date"]<-t14 - s14[,"Branchpoint_date"]
  sample_info<-rbind(s13, s14)
  
  time1<-NULL
  time2<-NULL
  time0<-NULL
  delta.t<-NULL
  
  for (a in 1:nrow(pedigree))
  {
    
    temp0<-pedigree[a,]
    pair1<-sample_info[which(as.character(sample_info[,1]) == as.character(temp0[,1])),]
    pair2<-sample_info[which(as.character(sample_info[,1]) == as.character(temp0[,2])),]
    
    if (pair1[,"Stem"] == pair2[,"Stem"])
    {
      t1<-pair1[,2] + pair1[,3]
      t2<-pair2[,2] + pair2[,3]
      t0<-min(pair1[,3], pair2[,3])
      time1[a]<-t1
      time2[a]<-t2
      time0[a]<-t0
      
    }
    
    if (pair1[,"Stem"] != pair2[,"Stem"])
    {
      t1<-pair1[,2]+pair1[,3]
      t2<-pair2[,2]+pair2[,3]
      t0<-min(pair1[,3], pair2[,3])
      time1[a]<-t1
      time2[a]<-t2
      time0[a]<-0
    }
    
    
  }
  
  delta.t<-time1 + time2 - 2*time0
  pdata<-cbind(pedigree, time1, time2, time0, delta.t) 
  
  pslim<-cbind(pdata[,"time0"],  pdata[,"time1"], pdata[,"time2"], pdata[,"D.value"])
  colnames(pslim)<-c("time0", "time1", "time2", "D.value")
  
  out<-list(pslim, pdata)
  print(out)
  
} # End of function        





