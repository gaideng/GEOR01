#GEO  
library(githubinstall)
# 22
BiocManager::install('ArrayExpress')
library(GEOquery)

downGSE = function(studyID,destDir='.'){
  library(GEOquery)
  eSet = getGEO(studyID,destdir = destDir,getGPL = F)
  epresData = exprs(eSet[[1]])
  pData = pData(eSet[[1]])
  write.csv(epresData,paste0(studyID,'_eprsData.csv'))
  write.csv(pData,paste0(studyID,'_metadata.csv'))
  return(eSet)
}
eSet = downGSE('GSE42872')
eRaad = read.csv('GSE42872_eprsData.csv')
rownames(eRaad) = eRaad[,1]
eRaad = eRaad[,-1]
eRaad2 = exprs(eSet[[1]])

#GEO ID 转换  GSE42872
BiocManager::install('hugene10sttranscriptcluster.db')
### 加载db失败
options(connectionObserver = NULL)
library(hugene10sttranscriptcluster.db)
?hugene10sttranscriptcluster.db
ls("package:hugene10sttranscriptcluster.db")
ids = toTable(hugene10sttranscriptclusterSYMBOL)
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol))))

exprSet42872 = exprs(eSet[[1]])
table(rownames(exprSet42872) %in% ids$probe_id)
dim(exprSet42872)
#过滤出有基因对应的探针行
exprSet42872s = exprSet42872[rownames(exprSet42872) %in% ids$probe_id,]
dim(exprSet42872s)
#顺序
ids2 = ids[match(rownames(exprSet42872s),ids$probe_id),]
head(ids2)
tmp = by(exprSet42872s,
         ids$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])


eSet5949 = downGSE('GSE5949')
exprsData1 = exprs(eSet5949[[1]])
BiocManager::install('hgu95av2.db')
library(hgu95av2.db)















