a=1:10
#维度
dim(a)=c(2,5)

a=data.frame(a)
#列名 行名
names(a)
rownames(a)

rownames(a)=a[,2]
pheatmap::pheatmap(a)

a[2,2]='4'

class(a)
View(a)
str(a)
pheatmap::pheatmap(a)

#清空变量
rm(list = ls())
rm(a)
ls()
b=as.data.frame(a)
b[2,2]='4'
save(b,file = 'binput.Rdata')

rm(list=ls())

load(file = 'binput.Rdata')
pheatmap::pheatmap(b)
str(b)
b[,2]=as.integer(b[,2])
str(b)
pheatmap::pheatmap(b)

d=options()
str(d[[2]])
str(d$add.smooth)

lapply(d, length)
e = d$connectionObserver
f = d[[2]]
s = unlist(lapply(d,length))
View(s)
class(s)
ss=as.numeric(s)
View(ss)
unlist(x = ss,use.names = F)

b=c(1,100,4)

grep('4',b[,1])

gse428=read.csv('GSE42872_eprsData.csv')
grep('4$',gse428[,1])



name <- c("张三","李四","王五","小明","张华","李然","马涛","魏然")
chinese <- c(88,55,56,89,58,65,75,56)
english <- c(89,48,57,78,29,68,89,64)
cj_data <- data.frame(name,chinese,english)
name_data <- data.frame(name)

name_data=name_data[c(-2,-4),]
name_data=as.data.frame(name_data)
ss = merge(name_data,cj_data,x = name_data$name_data,y=cj_data$name)
??merge

write.table(GSE5949_eprsData,'gse5949.csv')
s1=read.table('gse5949.csv')

rowMeans(GSE5949_eprsData)
head(rowMeans(GSE5949_eprsData))

row1=GSE5949_eprsData[1,]
row1s=as.numeric(row1)
class(row1s)
View(row1s)
s2=c(1,2,3)
View(s2)

apply(GSE5949_eprsData,1,function(x){
  mean(x)
})

for (i in 1:nrow(GSE5949_eprsData)) {
  print(mean(as.numeric(GSE5949_eprsData[i,])))
}


GSE42872_eprsData <- read.csv("D:/rwork/GEOR01/GSE42872_eprsData.csv")
GSE42872_eprsData=read.csv(file = 'GSE42872_eprsData.csv',row.names = 1)
pheatmap::pheatmap(GSE42872_eprsData[1:50,])

sd=apply(GSE42872_eprsData,1,sd)
sdTop50=GSE42872_eprsData[names(sort(sd,decreasing = T)[1:50]),]
pheatmap::pheatmap(sdTop50)


downGSE('GSE151180')
GSEGSE151180_eprsData=read.csv(file = 'GSE151180_eprsData.csv',row.names = 1)
pheatmap::pheatmap(GSEGSE151180_eprsData[1:50,])
GSEGSE151180_eprsData_noNa=na.omit(GSEGSE151180_eprsData)
pheatmap::pheatmap(GSEGSE151180_eprsData_noNa[1:50,])
names(sort(apply(GSEGSE151180_eprsData_noNa,1,sd),decreasing = T)[1:50])
pheatmap::pheatmap(GSEGSE151180_eprsData_noNa[names(sort(apply(GSEGSE151180_eprsData_noNa,1,sd),decreasing = T)[1:50]),])


gene151180=read.table('MergeExpro_contrib1-GPL21575.txt',row.names = 1,header = T,sep = '\t')
pheatmap::pheatmap(gene151180[1:30,])
pheatmap::pheatmap(gene151180[names(sort(apply(gene151180,1,sd),decreasing = T)[1:30]),])

a1=rnorm(100)
dim(a1)=c(5,20)
library(pheatmap)
pheatmap(a1)
a2=rnorm(100)+2
dim(a2)=c(5,20)
pheatmap(a2)

pheatmap(cbind(a1,a2))
pheatmap(cbind(a1,a2),cluster_cols = F,cluster_rows = F)

b=as.data.frame(cbind(a1,a2))
names(b)=c(paste('a1',1:ncol(a1),sep = '_'),paste('a2',1:ncol(a2),sep = '_'))
pheatmap(b,cluster_cols = F)

pheatmap(gene151180[1:50,])
pheatmap(gene151180[1:50,], scale = "row", clustering_distance_rows = "correlation",cluster_cols = F)

gene151180_name=names(sort(apply(gene151180,1,sd),decreasing = T)[1:50])
pheatmap(gene151180[gene151180_name,])
pheatmap(gene151180[gene151180_name,],scale = 'row',clustering_distance_rows = 'correlation',cluster_cols = F)

pheatmap(gene151180)
pheatmap(gene151180, scale = "row", clustering_distance_rows = "correlation",cluster_cols = F)
pheatmap(gene151180, scale = "row",cluster_cols = F)

ttest151180=read.table('ttest151180.txt',header = T,sep = '\t',row.names = 1)
ttest1511801=subset(ttest151180,p.value<0.05,)

ttest1511803=subset(ttest151180,logFC<1)

ttest1511802=subset(ttest151180,logFC<1&logFC>-1)

ttest1511804=subset(ttest151180,logFC<1&logFC>-1&ttest151180$p.value<0.01)


nrow(gene151180[rownames(ttest1511801),])
pheatmap(gene151180[rownames(ttest1511804),][1:100,], scale = "row",cluster_cols = F)
pheatmap(gene151180[rownames(ttest1511804),][1:100,], scale = "row")
pheatmap(gene151180[rownames(ttest1511804),][1:100,])
pheatmap(gene151180[rownames(ttest1511804),][1:100,],cluster_cols = F)

pheatmap(gene151180[rownames(ttest1511804),],scale = 'row',cluster_cols = F)

b=as.data.frame(cbind(a1,a2))
names(b)=c(paste('a1',1:ncol(a1),sep = '_'),paste('a2',1:ncol(a2),sep = '_'))
pheatmap(b)
temp=data.frame(group=c(rep('a1',20),rep('a2',20)))
rownames(temp)=colnames(b)
pheatmap(b,annotation_col = temp)

gene151180Log=t(scale(t(log2(gene151180+1))))
pheatmap(gene151180Log[1:100,])











gene151180Log=log2(gene151180)
gene151180s=gene151180
gene151180s[1:4,1:4]
dim(gene151180s)

cg=names(tail(sort(apply(gene151180Log,1,sd)),1000))
pheatmap(gene151180Log[cg,],show_rownames = F,show_colnames = F)
n=t(scale(t(gene151180Log[cg,])))
n[n>2]=2
n[n< -2]=-2
n[1:4,1:4]
pheatmap(n,show_rownames = F,show_colnames = F)

#ID转换
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#加载失败解决方式
options(connectionObserver = NULL) #加载org.Hs.eg.db失败时的解决方法
options(stringsAsFactors = F)
suppressMessages(library(org.Hs.eg.db))

g2s=toTable(org.Hs.egSYMBOL)
g2e=toTable(org.Hs.egENSEMBL)
ensembl=read.table('ensembl.txt')
library(stringr)
ensembl_id=str_split(ensembl$V1,'[.]',simplify = T)[,1]
ensembl$ensembl_id=ensembl_id
b=merge(ensembl,g2e,by='ensembl_id',all.x=T)
d=merge(b,g2s,by='gene_id',all.x=T)

table(d$ensembl_id)
table(table(d$ensembl_id)>1)
d[table(d$ensembl_id)>1,]
table(d$ensembl_id)[table(d$ensembl_id)>1]
#按某列进行排序
d1=d[order(d$V1),]
d2=d[!duplicated(d$V1),]
#两个表格按某列进行匹配排序
d3=d2[match(ensembl$V1,d2$V1),]

#################生存分析
rm(list = ls())

a=read.table('LGG_93663_50_50.csv',header = T,sep=',',fill = T)
colnames(a)
dat=a
#
options(stringsAsFactors = F)
BiocManager::install('ggstatsplot')
library(ggstatsplot)
ggbetweenstats(data = dat,x = Group,y=Expression)
library(ggplot2)
#BiocManager::install('survminer')
install.packages("survminer")
install.packages("survival")
install.packages("vctrs")
library(ggpubr)
library(vctrs)
library(survminer)
library(survival)
table(dat$Status)
dat$Status=ifelse(dat$Status=='Dead',1,0)
sfit=survfit(Surv(Days,Status)~Group,data = dat)
sfit
summary(sfit)
ggsurvplot(sfit,conf.int = F,pval = T)
ggsurvplot(sfit,palette = c('#E7B800','#2E9FDF'),risk.table = T,pval = T,conf.int = T,xlab='Time in months',ggtheme = theme_light())




