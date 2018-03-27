library(org.Hs.eg.db)
library(annotate)

TopNDisease = function(N,data){
  result = {}
  for(i in 1:ncol(data)){
    temp = data[,i]
    names(temp)=rownames(data)
    tmp = rev(sort(temp))
    symbol = lookUp(colnames(data)[i], 'org.Hs.eg', 'SYMBOL')
    for(j in 1:N){
      result = rbind(result,c(symbol,names(tmp)[j]))
    }
  }
  return(result)
}

setwd("/path/to/Prediction")

read_data = read.csv("eqHydraPredictionResults.txt",header=F)
disease = read.csv("index_disease.txt",header=F)
gene = read.csv("index_gene.txt",header=F)
rownames(read_data)=disease[,2]
colnames(read_data)=gene[,2]
number_of_disease_related_gene = 4838
top1 = TopNDisease(1,read_data[,-c(1:number_of_disease_related_gene)])
top5 = TopNDisease(5,read_data[,-c(1:number_of_disease_related_gene)])

write.table(top1,"eqhydraTop1PredictionPathwayIntegrated.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(top5,"eqhydraTop5PredictionPathwayIntegrated.txt",row.names=F,col.names=F,sep="\t",quote=F)
