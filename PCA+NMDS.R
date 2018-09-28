library("vegan")
library("ggbiplot")

####### 输入数据 #######
otu = read.table("otu_table.xls",sep="\t",head=T,row.names=1)
data = t(otu)
group = read.table("group4.txt",sep="\t",row.names= 1)
sample = rownames(group)
class = as.factor(group$V2)

######### PCA #########
data.pca <- prcomp(data, scale. = TRUE)
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.title.x =element_text(size=14), 
          axis.title.y=element_text(size=14),
          axis.title = element_text(color='black',vjust=0.1),
          legend.title=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'))}

ggbiplot(data.pca, obs.scale = 1, var.scale = 1,groups = class, 
         ellipse = T, circle = F,var.axes = F,ellipse.prob = 0.68) +
  scale_color_discrete(name = '') +
  theme_zg() 


########### 基于bray-curtis的差异分析 ############
distance.bray<-vegdist(data,method = 'bray')
anosim.result<-anosim(distance.bray,class,permutations = 999)
adonis.result = adonis(data~class,data = group,permutations = 999,method="bray")

################ NMDS ########################
sol <- metaMDS(data)
NMDS = data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2])

ggplot(data = NMDS, aes(MDS1,MDS2))+theme_bw()+
  theme(panel.grid=element_blank())+
  geom_point(aes(shape= group$V2, color = group$V2))+
  geom_text(label=sample,size=3,hjust=0.5, vjust=0)+
  theme(legend.key=element_rect(colour="white",fill="white",size=3,linetype="dashed"))+
  theme(legend.title=element_blank())+
  geom_polygon(aes(fill=group$V2,alpha=0.8),show.legend=FALSE)+
  annotate("text",x=min(NMDS$MDS1),y=min(NMDS$MDS2),hjust=0,vjust=0,label=paste("Stress:",sol$stress))+
  theme_zg()

#################### 输出 ###################
NMDS_out = data.frame(sample,NMDS$MDS1,NMDS$MDS2)
write.table (NMDS_out, file ="NMDS.txt", sep ="\t", row.names =F, col.names =TRUE, quote =FALSE)
sink("anosim.txt")
anosim.result
sink("adonis.txt")
adonis.result
sink()
