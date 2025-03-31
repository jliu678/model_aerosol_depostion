#batch read rds====
BatchReadRDS<-function(folder="./GD,cfu4e3/",fileofInt="",
                       strBefore="_relP",strAfter="\\.rds",str2rm=""){
  lf<-list.files(path = folder, pattern = NULL, all.files = FALSE,
                 full.names = F, recursive = FALSE,
                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  lf<-lf[!dir.exists(list.files(path = folder,full.names = T))]
  
  lf<-lf[grep(pattern =fileofInt,x=lf)]
  
  mCFU4e3<-lapply(lf,function(i) readRDS(paste0(folder,i)))
  
  # library(qdap)

  names(mCFU4e3)<-sapply(lf,function(i){
    # unname(genXtract(i, strBefore, strAfter))
    # substr(i,regexpr(pattern =strBefore,i)+nchar(strBefore),regexpr(pattern =strAfter,i)-1)
    gsub(pattern = str2rm,replacement="",x = i)
    })
  
 mCFU4e3 <- mCFU4e3[order(names(mCFU4e3))]
}
