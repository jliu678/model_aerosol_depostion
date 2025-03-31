#faster Monte_Carlo

#uniform
Monte_Carlo_1<- function (CFU, nAlveoli, nSim=1e6,extraLabel=""){
  t1<-integer(3);
  names(t1)<-c("nNicheOccupiedBy1Bug","nNicheOccupiedBy2Bug",
               "nNicheOccupiedByMore")
  for (i in 1:nSim) {
    t<-table(sample(1L: nAlveoli, CFU, replace=T))
    t1[1]<-t1[1]+as.integer(sum(t==1L))
    t1[2]<-t1[2]+as.integer(sum(t==2L))
    t1[3]<-t1[3]+as.integer(sum(t>2L))
  }
  saveRDS(t,paste0("CFU", CFU, ",nAlveoli", nAlveoli,".rds"))
}
for (CFU in seq.int(2e3,4e3,1e2)) {
  for (nAlveoli in seq.int(1e5,2.4e6,1e5)){
    tm<-proc.time()
    Monte_Carlo_1(CFU,nAlveoli,extraLabel="_1",)
    proc.time()-tm 
  }
}


#using the above optimized relP and P11
Monte_Carlo_4<- function (CFU, nAlveoli,P11,relP,MaxnBugCounted=9L,folder,extraLabel="",nSim=1e6){
  t<-integer(MaxnBugCounted*2L+2L);
  c1<-paste0("nNicheOccupiedBy",c(1L:MaxnBugCounted,"more"),"Bug")
  names(t)<-paste0(rep(c("P11_","P12_"),each=(MaxnBugCounted+1L)),c1)
  
  p<-c(rep(x=1,P11),rep(x= as.integer(relP),nAlveoli-P11))
  itr<-(0L:10L)*nSim/10L
  for (j in 1L:10L) {
    for (i in (itr[j]+1L):itr[j+1L]) {
      # t2<-sample.int(nAlveoli,size = CFU, replace = T)
      t1<-table(sample.int(nAlveoli,prob =p ,size = CFU, replace = T))
      l<-names(t1)<=P11
      
      l1<-lapply(1L:MaxnBugCounted,function(m) t1==m)
      l1[[MaxnBugCounted+1L]]<-t1>MaxnBugCounted
      
      for (k in 1L:(MaxnBugCounted+1L)){
        t[k]<-t[k]+as.integer(sum(l & l1[[k]]))
        t[k+MaxnBugCounted+1L]<-t[k+MaxnBugCounted+1L]+as.integer(sum((!l) & l1[[k]]))
      }
    }
    saveRDS(t/(nSim/10)/j,paste0(folder,"CFU", CFU, ",nAlveoli", nAlveoli,extraLabel,j,".rds")) 
  }
}
for (i in seq.int(from = 0L,to = 10L,by = 1L)){
  bootstraping4(CFU=2e3L*2L^(i),nAlveoli=12e5L,P11=1196347L,
                relP=2.640163e3,folder="BootstrapMaxp/",extraLabel=paste0("optP11relP_nNiche2bugLagerThan15_"),nSim = 1e3L)
  
}