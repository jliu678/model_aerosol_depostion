#given CFU,seek a matrix P to maximize f1(P)======
# P =[nP1,np2
#     P1 ,p2]
# rowsum(P)=[nAlveoli;]
# P[1,] %*% P[2,]=1

# N is a matrix whose element Nij is the number of niches with Pj (the
# probability that a pathogen will be deposit into the niche) and occupied
# by i (an integer number, any of 1:1:3) pathogens
# N=[N11 N12
#   N21 N22
#   N31 N32]
# for N used in downstream, dim(N) is c(2,2) 
# N=[N11 N12
#   N21 N22]
# CFU should be able to be calculated from N, but for speed we directly input CFU

#** k=the portion of mtb that goes to a niche of one, two or three-mtb occupied niche**
#** as reported in Fig1D **
# k<-round(c(1:3)*c(81.25,15.625,3.125)/c(c(1:3)%*%c(81.25,15.625,3.125)),3)
# k<-round(c(1:2)*c(83,17)/c(c(1:2)%*%c(83,17)),3)

source('Funs.R')

pre_computed<-function(k=round(c(1:2)*c(83,17)/c(c(1:2)%*%c(83,17)),3),
               CFU=3e3L){
  #** k*CFU is the number of mtb that  goes to a niche of one, two or three-mtb occupied niche**
  #** SUMni is the number of niches that are occupied by one, two or three-mtb **
  SUMni<-as.integer(round(1/(1:length(k*CFU))*k*CFU))
  
  #** below pre-compute matrices or vectors to select N matrices that obeys P matrix **
  #** Ln1 is argument of f1(). **
  #** Ln1 is a 4 X I matrix, I is the number of all possible N constrained by CFU and paper's Fig1D **
  #**     Rows 1 to 4 of Ln1 is possible values of {N11}, {N12}, {N21}, {N22} respectively **
  #** 
  
  Ln0<-reshape2::melt(matrix(nrow = SUMni[2]+1, ncol =  SUMni[1]+1))[,c("Var2","Var1")]-1L
  Ln1<-cbind(Ln0[,1],SUMni[1]-Ln0[,1],Ln0[,2],SUMni[2]-Ln0[,2])
  Ln1<-t(Ln1)#0.13s/run
  
  #** LsA is argument of f1(). **
  #** LsA= Ln1[1,] + Ln1[3,], i.e. N11+N21 the total number of niches that have p1 chance to be occupied**
  LsA<-Ln1[1,] + Ln1[3,]
  
  #** LsD is argument of f1(). **
  #** LsD= Ln1[2,] + Ln1[4,], i.e. N12+N22 the total number of niches that have p2 chance to be occupied**
  LsD<- sum(SUMni)-LsA
  # identical(Ln1[2,] + Ln1[4,],LsD)
  
  #logCirCle3 component,Sumini1, and Sumini2=CFU-Sumini1 
  Sumini1full<-LsA+(0L:SUMni[2]) 
  
  #logCircle2 is independent of N P, and is the same for different N P as long as CFU remains the same
  logCircle2<- logCmbn(CFU,SUMni[1])+
    sum(log(seq.int(CFU-SUMni[1]-1L,1L,by = -2L)))+logPermute(SUMni[2],SUMni[2]-1)
  
  list(Ln1,LsA,LsD,Sumini1full,logCircle2)
}

