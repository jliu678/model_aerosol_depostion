getdP<- function(relP,P11,nNiches) {
  
  Prow1<- c(P11,nNiches-P11)
  crP<-f1(Prow1,log(c(1,relP)/c(c(1,relP)%*%Prow1)),CFU,LsA,LsD,Ln1,Sumini1full,logCircle2)
  
  Prow1d<- c(P11+100L,nNiches-P11-100L)
  dP11<-f1(Prow1d,log(c(1,relP)/c(c(1,relP)%*%Prow1d)),CFU,LsA,LsD,Ln1,Sumini1full,logCircle2)/crP-1
  
  relPd<-relP+0.1
  drelP<-f1(Prow1,log(c(1,relPd)/c(c(1,relPd)%*%Prow1)),CFU,LsA,LsD,Ln1,Sumini1full,logCircle2)/crP-1
  
  list(P11=dP11,relP=drelP,crP=crP)
}
GD<-function(relP=10,P11=2e5L,nNiches=1e6L,nItr=1e3,mrelP=0.5,mP11=0.5,
             k=round(c(1:2)*c(83,17)/c(c(1:2)%*%c(83,17)),3),
             CFU=3e3L){
  #1>mP11>0 the percentage of last dP11 contribute to new dP11
  #1>mrelP>0 the percentage of last drelP contribute to new drelP  
  
  precomputed<-pre_computed(k,CFU)
  precomputed_names<-c('Ln1','LsA','LsD','Sumini1full','logCircle2')
  for (i in (1:length(precomputed))){
    assign(precomputed_names[i],precomputed[i])
  }
  
  m<-matrix(ncol = 5, nrow= nItr)
  colnames(m)<-c("relP","P11","P", "drelP","dP11")
  dP<-getdP(relP,P11,nNiches)
  if (is.na(dP$P11)) {stop("current P is 0, try set a better initializaiton")}
  if (abs(P11*dP$P11)<=1000 && abs(P11*dP$P11)>=100){
    dP11<-P11*dP$P11}
  else if (abs(P11*dP$P11)>1000) {dP11<-sign(dP$P11)*1000L}
  else {dP11<-sign(dP$P11)*100L}
  
  if (abs(relP*dP$relP)<=1 && abs(relP*dP$relP)>=0.1){
    drelP<-relP*dP$relP}
  else if (abs(relP*dP$relP)>1) {drelP<-  sign(dP$relP)}
  else {drelP<-  sign(dP$relP)*0.1}
  
  m[1,]<-c(relP,P11,dP$crP,drelP,dP11)
  relP<-relP+drelP
  P11<-round(P11+dP11)
  
  for (i in 2:nItr) {
    tm<-proc.time()
    dP<-getdP(relP,P11,nNiches)
    if (abs(P11*dP$P11)<=1000 && abs(P11*dP$P11)>=100){
      dP11<-mP11*dP11+(1-mP11)*P11*dP$P11}
    else if (abs(P11*dP$P11)>1000) {dP11<-mP11*dP11+(1-mP11)*sign(dP$P11)*1000L}
    else {dP11<-mP11*dP11+(1-mP11)*sign(dP$P11)*100L}
    
    if (abs(relP*dP$relP)<=1 && abs(relP*dP$relP)>=0.1){
      drelP<-mrelP*drelP+(1-mrelP)*relP*dP$relP}
    else if (abs(relP*dP$relP)>1) {drelP<-mrelP*drelP+(1-mrelP)*sign(dP$relP)}
    else {drelP<- mrelP*drelP+(1-mrelP)*sign(dP$relP)*0.1}
    
    m[i,]<-c(relP,P11,dP$crP,drelP,dP11)
    
    relP<-relP+drelP
    P11<-round(P11+dP11)
    cat("iteration",i,"done",(proc.time()-tm)[3],"s\n")
  }
  return(m)
}

#1. set up CFU in "whatGivesMaxP.R" and source it====
#2. initial P setup=====
nNiches<-2e5L

relP<-500
P11<-9.5e4L
# P<- rbind(number=Prow1,prob=c(1,relP)/c(c(1,relP)%*%Prow1))

#3. GD,observe,optimize and store results===== 
m<-GD(relP=relP,P11=P11,nNiches=nNiches,nItr=500)

m<-GD(relP=m[300,1],P11=m[50,2],nNiches=3e5L,nItr=50)
# 500.3432 294832.5 3.055627e-19  0.7285322   599.7825

m<-GD(relP=m[50,1],P11=m[50,2],nNiches=4e5L,nItr=50)
# 658.6710 395298.9 2.056486e-19 0.3704280  -600.0000

m<-GD(relP=m[50,1],P11=m[50,2],nNiches=5e5L,nItr=50)
# 761.2309 494782.4 1.117586e-19 0.9174192  333.3332

m<-GD(relP=m[50,1],P11=m[50,2],nNiches=6e5L,nItr=50)
# 899.5130 594915.7 1.224963e-19  0.891580435  600.0000

m<-GD(relP=m[50,1],P11=m[50,2],nNiches=7e5L,nItr=1000)
# relP           P11             P         drelP          dP11 
# 1.458354e+03  6.949697e+05  8.742549e-19  3.074121e-01 -6.000000e+02 

m<-GD(relP=m[1000,1],P11=m[1000,2],nNiches=8e5L,nItr=500)
# 1701.588 794538.2 9.787099e-19 0.5187465  333.3333

m<-GD(relP=m[500,1],P11=m[500,2],nNiches=9e5L,nItr=500)
# 1936.127 894647.1 1.181415e-18 0.4178746  333.3333

for (i in seq.int(2e5L,24e5L,by=1e5))
{m<-GD(relP=m[500,1],P11=i-5e3L,nNiches=i,nItr=500)
  saveRDS(m,file = paste0("nNiches",i,"Itr500.rds"))
}

m<-GD(relP=277,P11=1e5L,nNiches=2e5L,nItr=1000)
saveRDS(m,file = paste0("nNiches",nNiches,"Itr",nItr,".rds"))

#observations====
#over iteration will keep increasing relP but not P
#comparing "nNiches2e5Itr1e3.rds" and "nNiches2e5,Succeeding24e5.rds"-->
#   600.9509 194322.3 4.007808e-18 0.31014295  600.0000   
#   720.6564 194544.8 5.607676e-18 0.10197788 -600
