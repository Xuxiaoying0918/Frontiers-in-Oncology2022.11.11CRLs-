######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")


#引用包
library(limma)
library(ggplot2)
library(ggpubr)

lncName="AP000695.2"     #目标lncRNA的名称
expFile="lncRNA.txt"     #表达数据文件
setwd("C:\\biowolf\\lncRNA\\14.diff")      #设置工作目录

#读取表达数据文件,并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=t(data[lncName,,drop=F])

#获取正常样品和肿瘤样品的数目
group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       #正常组样品数目
treatNum=length(group[group==0])     #肿瘤组样品数目
Type=c(rep(1,conNum), rep(2,treatNum))

#样品分组
exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("lncRNA", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")

#输出目标lncRNA的表达量
outTab=exp
colnames(outTab)=c(lncName, "Type")
outTab=cbind(ID=row.names(outTab), outTab)
write.table(outTab, file="singleLncExp.txt", sep="\t", quote=F, row.names=F)

#设置比较组
group=levels(factor(exp$Type))
exp$Type=factor(exp$Type, levels=group)
comp=combn(group,2)
exp$lncRNA=log2(exp$lncRNA+1)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制箱线图
boxplot=ggboxplot(exp, x="Type", y="lncRNA", color="Type",
		          xlab="",
		          ylab=paste0(lncName, " expression"),
		          legend.title="Type",
		          palette = c("blue","red"),
		          add = "jitter")+ 
	stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

#输出图形
pdf(file=paste0(lncName,".diff.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

