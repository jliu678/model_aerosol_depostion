#succeeding from 'bootstrap,OptimizeRelpP11.R'

#data prep====

#1.uniform====
source("BatchReadRDS.r")
t<-BatchReadRDS(folder="./BootstrapRandomDeposit/",str2rm = "uniform|\\.rds",
                fileofInt = "_uniform[0-9]")
t1<-t[[1]]
for (i in 2:length(t)){
  t1<-rbind(t1,t[[i]])
}
rownames(t1)<-names(t)
rm(t,i)
t1<-data.frame(t1)
t1$CFU<- as.numeric(gsub(pattern ="CFU|,nAlveoli[0-9].*" ,x = rownames(t1),replacement = ""))
t1$nAlveoli<- as.numeric(gsub(pattern =".*,nAlveoli|_[0-9].*" ,x = rownames(t1),replacement = ""))
t1$Possible1Bug<-t1$CFU-3*t1$nNicheOccupiedByMore-2*t1$nNicheOccupiedBy2Bug
t1$nNicheOccupiedBy1Bug[is.na(t1$nNicheOccupiedBy1Bug)]<-t1$Possible1Bug[is.na(t1$nNicheOccupiedBy1Bug)]
t1$totalnNiche<-rowSums(t1[,1:3])
summary(t1$totalnNiche)
t1<-t1[,-ncol(t1)]
saveRDS(t1,"uniform,bugDistribution.rds")



#2.model====
source("BatchReadRDS.r")
t<-BatchReadRDS(folder="./BootstrapMaxp/",str2rm = ",nAlveoli1200000optP11relP|\\.rds",
                fileofInt = ",nAlveoli1200000optP11relP_[0-9]")
# t<-BatchReadRDS(folder="./BootstrapMaxp/",str2rm = ",nAlveoli1200000optP11relP_nNiche2bugLagerThan15|\\.rds",
#                 fileofInt = ",nAlveoli1200000optP11relP_nNiche2bugLagerThan15_")
t1<-t[[1]]
for (i in 2:length(t)){
  t1<-rbind(t1,t[[i]])
}
rownames(t1)<-names(t)
rm(t,i)
t1<-data.frame(t1)
t1$totalnNiche<-rowSums(t1)
summary(t1$totalnNiche)
t1$CFU<- gsub(pattern ="CFU|_[0-9]{1,2}" ,x = rownames(t1),replacement = "")
saveRDS(t1,"CFUvary,P11relPopt,bugDistribution.rds")
# saveRDS(t1,"CFUvary,P11relPopt_nNiche2bugLagerThan15,bugDistribution.rds")



#Demonstration====

#1. initial deposition with uniformed distribution (varied CFU and nNiches)
#   and with the best-modeled distribution to show successful modeling

#2. initial deposition with varying CFU to show max niches before neutrophil recruitment
#   by demonstrating bug distribution with the best model and the one would give most niches 
#   among all good models


#1.plot "uniform"====
t2<-readRDS("uniform,bugDistribution.rds")
colnames(t2)<-gsub(pattern = "nNicheOccupiedBy|Bug",replacement = "",x=colnames(t2))
colnames(t2)<-gsub(pattern = "More",replacement = "3+",x=colnames(t2))
t2[,7:9]<-t2[,1:3]/t2$totalnNiche*100
t2$total<- rep(100,nrow(t2))
library(tidyverse)
t2a<- t2 %>% group_by(CFU,nAlveoli) %>% summarise_if(is.numeric, list(m=mean,std=function(i) {sqrt(var(i))}))
t2a<-reshape2::melt(data = t2a, id.vars = c("CFU","nAlveoli"), variable.name = "nBugOrTotalNiche",value.name = "PercentageOrNumber")
t2a<-cbind(t2a[1:(nrow(t2a)/2),],t2a[(nrow(t2a)/2+1):nrow(t2a),])
colnames(t2a)[ncol(t2a)]<-"stdofPercentageOrNumber"
t2a<-t2a[,!duplicated(colnames(t2a))]
t2a$nBugOrTotalNiche<-factor(gsub(x = t2a$nBugOrTotalNiche,pattern = "_m|\\.1",replacement = ""))
t2t<-t2a[is.element(t2a$nBugOrTotalNiche,c("totalnNiche","total")),]
t2t<-t2t[1:(nrow(t2t)/2),]
t2a<-t2a[!is.element(t2a$nBugOrTotalNiche,c("totalnNiche","total")),]
t2a<-t2a[(nrow(t2a)/2+1):nrow(t2a),]
colnames(t2a)[4:5]<-c("percentage","stdofPercentage")
t2a<-plyr::rbind.fill(t2a,data.frame(stdofPercentage=rep("paper",3),
                                     CFU=rep(0,3),nAlveoli=rep(0,3),
                                     nBugOrTotalNiche=unique(t2a$nBugOrTotalNiche),
                                     percentage=100*c(0.8125,0.159375,1-(0.8125+0.159375))))
t2a$id<-factor(paste(t2a$CFU,t2a$nAlveoli,sep = "_"))


t2a<-t2a %>% mutate(id=factor(id,levels=(levels(id)[c(1,3,2,4,6,5,7,9,8,10)])))%>% 
  mutate(nBugOrTotalNiche=factor(nBugOrTotalNiche,levels=(levels(nBugOrTotalNiche)[1:3])))

library(ggplot2)

cl1<-c("#edf2ee","green","#f76d6d")
cl2<-c("#000000", "#ffd6d6","#e06e6e","#cc0202",
      "#c9d1ff","#667aed","#001bb5","#c4ffc5","#67eb69","#00c903")
# lb1<-c(expression(AM~containing~one~italic(Mtb)),
#        expression(AM~containing~two~italic(Mtb)),
#        expression(AM~containing~>two~italic(Mtb)))
lb1<-c("AM containing one Mtb    ",
        "AM containing two Mtb    ",
        "AM containing >two Mtb    ")

lb2<- paste0(rep(paste0(c("2K ","3K ","4K "),"Mtb,"),each=3),
            rep(paste0(c("0.1M ","1.2M ","2.4M "),"Niches"),3))



p<-ggplot(data = t2a,aes(x=id,y=percentage,group=nBugOrTotalNiche,fill=nBugOrTotalNiche))+
  scale_fill_manual(values=cl1,labels=lb1)+
  # scale_color_manual(values=cl)+
  geom_col()+
  scale_y_continuous(expand= c(0,0),limits = c(0,100.3), breaks=seq(0, 100, 20),name = "% of all infected AM")+
  scale_x_discrete(labels=c("paper reported",lb2))+
  coord_flip(clip = "off", expand = FALSE) +
  guides(fill=guide_legend(nrow=1))+
  geom_hline(yintercept=seq.int(20,80,20),linetype =3,size = 0.2,color="grey")+
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(face="bold",size=14,color = 'black'),
        axis.text.y = element_text(face="bold",size=12,color = c('black',rep("blue",9))),
        axis.ticks.y.left = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face="bold",size=14),
        legend.position = c(0.51,1.07),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0.7,"cm"), 
        legend.title=element_blank(),
        legend.text = element_text(face="bold",size=13),
        panel.background = element_blank(),
        # panel.border = element_rect(linetype = "dashed", fill = NA),
        plot.title = element_text(size = 18, hjust = 0, face = "bold", colour = "black", vjust = 1),
        plot.subtitle=element_text(size = 18, hjust = 0, face = "bold", colour = "black"),
        plot.caption = element_text(size = 16, hjust = 0.5, face = "italic", colour = "black"),
        plot.background = element_blank(),
        plot.margin = margin(1.4, 0.4, 0, 0, "cm"))
p

p+guides(fill=guide_legend(nrow=3,ncol=1))+theme(plot.margin = margin(1.8, 0.4, 0, 0, "cm"),
                                                 legend.position = c(0.51,1.1),
                                                 legend.key = element_rect(fill = NA, color = NA))
#using export to check legend location etc.
ggsave("1.png",dpi = 800)

saveRDS(p,"uniform,boostrapVisulization.rds")



#2.plot "model"====
t2<-readRDS("CFUvary,P11relPopt,bugDistribution.rds")
 # t2<-readRDS("CFUvary,P11relPopt_nNiche2bugLagerThan15,bugDistribution.rds")

t2[,1:10]<-t2[,1:10]+t2[,11:20]
t2<-t2[,-c(11:20)]
colnames(t2)<-gsub(pattern = "P11_nNicheOccupiedBy|Bug",replacement = "",x=colnames(t2))
colnames(t2)<-gsub(pattern = "more",replacement = "10+",x=colnames(t2))
t2$CFU<-as.integer(t2$CFU)
t2<-t2[order(t2$CFU),]
t2[,13:22]<-t2[,1:10]/t2$totalnNiche
t2$total<- rep(1,nrow(t2))



library(tidyverse)
t2a<- t2 %>% group_by(CFU) %>% summarise_if(is.numeric, mean)

library(reshape2)
t2a<-melt(data = t2a, id.vars = "CFU", variable.name = "nBugOrTotalNiche",value.name = "PercentageOrNumber")
t2a<-cbind(t2a[1:(nrow(t2a)/2),],t2a[(nrow(t2a)/2+1):nrow(t2a),])

if (identical(t2a[,1],t2a[,4])) {t2a<- t2a[,-c(4:5)]}
colnames(t2a)[3:4]<-c("number","percentage")
t2a$nBugOrTotalNiche[t2a$nBugOrTotalNiche=="totalnNiche"]<-"total"
t2a$percentage<-t2a$percentage*100
t2a$number<-round(t2a$number)
t2b<-t2a[(nrow(t2a)-length(unique(t2a$CFU))+1):nrow(t2a),] 
rownames(t2b)<-t2b[,1]


#2a.plot "modeled" fig1D day1====
#and compare with paper
t2c<-t2a[t2a$CFU==2000 & t2a$nBugOrTotalNiche!='total',]
t2c<-rbind(t2c,matrix(nrow=nrow(t2c),ncol=ncol(t2c),dimnames = list(paste0(rownames(t2c),"paper"), colnames(t2c))))
t2c$nBugOrTotalNiche[(nrow(t2c)/2+1):nrow(t2c)]<-t2c$nBugOrTotalNiche[1:(nrow(t2c)/2)]
t2c$percentage[(nrow(t2c)/2+1):nrow(t2c)]<-100*c(0.8125,0.159375,0.01875,0,0.007375,rep(0,5))
t2c$number<-rep(c('modeled','paper reported'),each=nrow(t2c)/2)
t2c<- t2c %>% arrange(desc(number)) %>% # First sort the dataframe, but NOT the factor levels
  mutate(number=factor(number, levels=unique(number))) # Then update the factor levels


library(ggplot2)
ggplot(data = t2c,aes(x=nBugOrTotalNiche,y=percentage,group= number,fill=number))+
  scale_fill_manual(values=c("#000000", "#808080"))+
  geom_col(position = position_dodge2(width = 0.7))+
  scale_y_continuous(expand= c(0,0),limits = c(0,100.3), breaks=seq(0, 100, 20))+
  # guides(fill=guide_legend(label.theme = element_text(face="bold",size=14)))+
  guides(fill = guide_legend(nrow = 2,ncol = 1))+
  labs(x ="number of Mtb per AM", y = "% of all infected AMs")+
  theme(axis.line = element_line(linetype = 1,size = 1),
        axis.text = element_text(face="bold",size=14,colour = "#000000"),
        axis.ticks = element_line(linetype = 1,size = 0.7,colour = "#000000"),
        axis.title = element_text(face="bold",size=14),
        legend.position = c(0.84,0.94),
        legend.title=element_blank(),
        legend.text = element_text(face="bold",size=14),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(0.2, 0, 0, 0, "cm"))
#to be uniformed,specify same size and dpi
ggsave("1.png",dpi = 800,width = 8.22125, height = 4.59625,units = "in")



#2b.animated visulizaiton=====
# initial deposition with varying CFU to show max niches before neutrophil recruitment

#attach funs and package
Rcpp::cppFunction('NumericVector rcpp_clip( NumericVector x, double a, double b){
    return clamp( a, x, b ) ;
}')
generateC<-function(nCFUwithoutNeutr=integer(1),nCFUgiveNeutr=integer(1)){
  color1<-viridis(begin = 0.34,n = 100,option = "D")[as.integer(seq(from=1,by=99/nCFUwithoutNeutr,length.out = nCFUwithoutNeutr))]
  
  color2<-viridis(begin = 0.18,end = 0.58,n = 100,option = "C")
  
  color2a<-switch(rcpp_clip(nCFUgiveNeutr,1,3),
         '#c70000',
         c('#ff0000','#c70000'),
         c(color2[as.integer(seq(from=100,by=99/(2-nCFUgiveNeutr),length.out = nCFUgiveNeutr-2))],'#ff0000','#c70000'))
  
  c(color1,color2a)
}
library(gganimate)
library(gifski)


levels(factor(t2a$CFU))
# nCFUwithoutNeutr<-which(levels(factor(t2a$CFU))==27000)#27000 is the largest CFU that dosen't give neutrophil recruitment
nCFUwithoutNeutr<-which(levels(factor(t2a$CFU))==23000)

nCFUgiveNeutr<-nlevels(factor(t2a$CFU))-nCFUwithoutNeutr

cl<-generateC(nCFUwithoutNeutr,nCFUgiveNeutr)

lb<-tibble(l=c(rep("",nCFUwithoutNeutr),rep("Neutrophils recruited!",nCFUgiveNeutr)),CFU=levels(factor(t2a$CFU)))


staticplot <- ggplot(t2a %>% filter(nBugOrTotalNiche != 'total'), 
                     aes(nBugOrTotalNiche, group = CFU, fill = as.factor(CFU), 
                         colour = as.factor(CFU)))+
  scale_fill_manual(values=cl)+scale_color_manual(values=cl)+
  geom_tile(aes(y = percentage/2, height = percentage, width = 0.7), alpha = 0.8) +
  # geom_errorbar(aes(ymin=PercentageOrNumber-sd, ymax=PercentageOrNumber+sd), width=.1)+
  # geom_text(aes(y = 0, label = paste(country, " ")), vjust = 0.2, hjust = 1, size = 8) +
  geom_text(aes(y = percentage, label = round(number), hjust = 0.5,angle=0,vjust = -0.14), size = 5.7) +
  guides(colour = FALSE, fill = FALSE) +
  scale_y_continuous(expand= c(0,0),limits = c(0,104))+
  labs(x ="number of Mtb per AM",y = "% of all infected AMs")+
  theme(axis.line = element_line(linetype = 1,size = 1),
        axis.text = element_text(face="bold",size=16,colour = "#000000"),
        axis.ticks = element_line(linetype = 1,size = 0.7,colour = "#000000"),
        axis.title = element_text(face="bold",size=16),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 24, hjust = 0.5, face = "bold", colour = "red",vjust = -4),
        plot.subtitle=element_text(size = 18, hjust = 0, face = "bold", colour = "black"),
        plot.caption = element_text(size = 16, hjust = 0.5, face = "italic", colour = "black"),
        plot.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"))
staticplot


anim = staticplot + 
  transition_states(CFU, transition_length = 1, state_length = 2)+
  labs(subtitle = "Infection dose: {closest_state} CFU\n{t2b$number[t2b$CFU==closest_state]} AMs initially infected",
       title = "{lb$l[lb$CFU==closest_state]}",
       caption = "")

animate(anim, nframes = 100, fps = 7,  width = 1100, height = 600, 
        renderer = gifski_renderer("DynamicBar.gif"))

animate(anim, nframes = 100, fps = 7,  width = 1100, height = 600, 
        renderer = av_renderer("DynamicBar.mp4"))
