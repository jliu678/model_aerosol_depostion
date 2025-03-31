#====
#main function

# N is a matrix whose element Nij is the number of niches with Pj (the
# probability that a pathogen will be deposit into the niche) and occupied
# by i (an integer number, any of 1:1:3) pathogens
# N=[N11 N12
#   N21 N22
#   N31 N32]
# CFU should be able to be calculated from N, but for speed we directly input CFU

f1<- function(Prow1,logProw2,CFU,LsA,LsD,Ln1,Sumini1full,logCircle2){
  # f1 gets the probability given a fixed P;
  
  # 1. generate all possible N compatible with P and nAlveoli/pathogen distribution
  # in the paper
  l<- (LsA<=Prow1[1]) & (LsD<=Prow1[2])
  Ln<-Ln1[,l]
  Sumini1<-Sumini1full[l]
  # Sumini2<-CFU-Sumini1
  # generally Ln is vector, but if ncol(Ln) can equal to ncol(Ln1), it would be matrix
  
  # 2. get logCircle3 of each N
  logCircle3<-logProw2[1]*Sumini1+logProw2[2]*(CFU-Sumini1)+logCmbn(sum(Prow1),Prow1[1])
  
  
  # 3. logCircle1 of each N
  logCircle1<-getlogCircle1(Ln,Prow1)
  
  # tm<-logCircle2+logCircle3+logCircle1
  # summary(tm)
  # summary(exp(tm))
  sum(exp(logCircle2+logCircle3+logCircle1))
  # why sum here?
}
#====
getlogCircle1<-function(Ln,Prow1){
  lLn<-length(Ln)
  Ln2sbstr<-Ln
  Ln2sbstr[3:4,]<-0L
  Ln2sbstr<-c(0L,0L,Ln2sbstr[-((lLn-1):lLn)])
  #all components of logCircle1
  alllogCircle1components<-logCmbn(Prow1-Ln2sbstr,Ln)
  #get logCircle1 of each N
  #speed on mac: aggregation<data.table<built matrix+colsums
  # system.time({for (i in 1:1e2){
  colSums(matrix(alllogCircle1components,nrow=4),na.rm = T)
  # }})
}
#====
#necessary functions
ramanujan <- function(n){
  n*log(n) - n + log(n*(1 + 4*n*(1+2*n)))/6 + log(pi)/2
}

logCmbn<- function(n,k){
  ramanujan(n) - ramanujan(k) - ramanujan(n-k)
}

logPermute<- function(n,k){
  ramanujan(n) - ramanujan(n-k)
}


