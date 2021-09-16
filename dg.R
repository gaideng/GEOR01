##R语言去除重复数据,过滤空行
###data 原始数据
###isNA 是否过滤空值 error：row.name是否有遗漏值 isNa==0代表不过滤;isNa==1代表过滤
###isNA missing values in 'row.names' are not allowed 
###移除重复数据 error:data are assentially constant
removeRepData = function(data,isNa=1){
  ##读取原始数据,读取之前先过滤是否有空值,有空值会导致无法读取
  if(isNa==1){
    befor = read.table(data,header = T,sep = '\t')
    befor=na.omit(befor)
    row.names(befor) = befor[,1]
    befor=befor[,-1]
  }else{
    befor = read.table(data,header = T,sep = '\t',row.names = 1)
  }
  #定义重复数据向量容器
  repRow=vector(mode="numeric",length=0)
  #循环所有原始数据获取所有列数据相同的行
  coln = ncol(befor)
  for (i in 1:nrow(befor)) {
    a = 0
    for (j in 1:coln) {
      if(j == 1){
        a=befor[i,j]
      }
      if(a != befor[i,j]){
        break
      }
      if(j == coln){
        repRow=append(repRow,-i)
      }
    }
  }
  
  #移除重复行 先判断是否有重复行,没有则直接赋值
  if(length(repRow)==0){
    after=befor
  }else{
    after=befor[repRow,]
  }
  #写出文件
  write.table(after,paste("New",data,sep = "-"),row.names = T,sep = '\t',quote = F,col.names = T)
}

removeRepData("MergeExpro_contrib1-GPL341.txt")