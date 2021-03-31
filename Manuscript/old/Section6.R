require(ggplot2)
require(foreach)
source('client-side/code/Manuscript/ggplot.style.R')
load('client-side/output/TME.between.primary.and.metastatic.R.output/TME.between.primary.and.metastatic.RData')

####### Fig 6a, show how we compare immune cell abundance between MET and PRI, taking purity into consideration
draw.df <- M.vs.P.Basal.rs$data$`CD8+ T cells`
ggplot(draw.df,aes(x=purity,y=immune.score ))          + geom_point(size=8.5,color='black') + ggplot.style + geom_smooth( method='loess',se=FALSE,color='red',lwd=5.5) 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section6/basal.CD8.immune.score.vs.purity.pdf',width = 20,height=20)
cor(draw.df$purity,draw.df$immune.score,method='spearman')

ggplot(draw.df,aes(x=purity,y=adjusted.immune.score )) + geom_point(size=8.5,color='black') + ggplot.style + geom_smooth( method='loess',se=FALSE,color='red',lwd=5.5) 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section6/basal.CD8.adjusted.immune.score.vs.purity.pdf',width = 20,height=20)
cor(draw.df$purity,draw.df$adjusted.immune.score,method='spearman')




##### Fig 6b, Basal subtype #################
data.list <- M.vs.P.Basal.rs$data
pooled.df <- foreach(e = names(data.list),.combine='rbind') %do% {
    df <- data.list[[e]]  
    df$cell.type  <-  e
    df
  
}

ggplot(pooled.df,aes(x=cell.type,y=adjusted.immune.score,fill= class )) +
                               geom_boxplot(lwd=1.5,outlier.shape = NA) + 
                   geom_point(position=position_jitterdodge(),size=3.0) + 
            scale_fill_manual(values = c('MET500'='red','TCGA'='grey')) + 
    theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('') 

ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/basal.TME.pdf',width = 20,height=20)



##### Fig 6c, LuminalB subtype #################
data.list <- M.vs.P.LumB.rs$data
pooled.df <- foreach(e = names(data.list),.combine='rbind') %do% {
  df <- data.list[[e]]  
  df$cell.type  <-  e
  df
  
}

ggplot(pooled.df,aes(x=cell.type,y=adjusted.immune.score,fill= class )) +
  geom_boxplot(lwd=1.5,outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(),size=3.0) + 
  scale_fill_manual(values = c('MET500'='red','TCGA'='grey')) + 
  theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('') 

ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/lumb.TME.pdf',width = 20,height=20)


##### Fig 6d, Her2 subtype #################
data.list <- M.vs.P.Her2.rs$data
pooled.df <- foreach(e = names(data.list),.combine='rbind') %do% {
  df <- data.list[[e]]  
  df$cell.type  <-  e
  df
  
}

ggplot(pooled.df,aes(x=cell.type,y=adjusted.immune.score,fill= class )) +
  geom_boxplot(lwd=1.5,outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(),size=3.0) + 
  scale_fill_manual(values = c('MET500'='red','TCGA'='grey')) + 
  theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('')  + ylim(-3.5,3.5)

ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/her2.TME.pdf',width = 20,height=20)


####################################### Trash code ####################################### 
# #+ geom_point(size=3.5,color='black') + ggplot.style + geom_smooth( method='loess') + geom_point(data=draw.df[draw.df$class == 'MET500',],aes(x=purity,y=adjusted.immune.score),color='red',size=3.5)
# 
# 
# rs <- M.vs.P.Basal.rs$rs
# rs <- rs[order(rs$p.value),]
# rs[rs$fdr < 0.05,]
# 
# draw.df <- M.vs.P.Basal.rs$data$`B cells`
# draw.df$class <- factor(draw.df$class,levels = c('TCGA','MET500'))
# ggplot(draw.df) + geom_boxplot(aes(x=class,y=adjusted.immune.score),lwd=1.0,outlier.shape = NA)+ ggplot.style + geom_jitter(aes(x=class,y=adjusted.immune.score),size=3.5,width=0.05) 
# 
# 
# draw.df <- M.vs.P.Basal.rs$data$`CD8+ T cells`
# draw.df$class <- factor(draw.df$class,levels = c('TCGA','MET500'))
# ggplot(draw.df) + geom_boxplot(aes(x=class,y=adjusted.immune.score),lwd=1.5,outlier.shape = NA)+ ggplot.style + geom_jitter(aes(x=class,y=adjusted.immune.score),size=3.5,width=0.05) 


