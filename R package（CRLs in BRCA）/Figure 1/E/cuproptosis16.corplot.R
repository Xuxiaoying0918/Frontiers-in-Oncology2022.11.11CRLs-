######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("ggExtra")


#引用包
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)

riskFile="risk.all.txt"                #风险文件
cuprExpFile="cuproptosisExp.txt"       #基因表达文件
lncExpFile="cuproptosisLncExp.txt"     #lncRNA表达文件
setwd("C:\\biowolf\\cuproptosis\\16.corplot")     #设置工作目录

#读取铜死亡基因的表达文件,并对数据进行处理
rt1=read.table(cuprExpFile, header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
cuproptosis=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
cuproptosis=avereps(cuproptosis)
cuproptosis=cuproptosis[rowMeans(cuproptosis)>0.1,]

#删掉正常样品
group=sapply(strsplit(colnames(cuproptosis),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
cuproptosis=t(cuproptosis[,group==0])

#读取lncRNA的表达文件,并对数据进行处理
rt=read.table(lncExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
lncRNA=data[,group==0]

#读取风险文件,提取模型lncRNA的表达量
riskRT=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
lncRNA=t(lncRNA[colnames(riskRT)[3:(ncol(riskRT)-2)],])

#相关性分析
outTab=data.frame()
for(lncrna in colnames(lncRNA)){
	for(gene in colnames(cuproptosis)){
		x=as.numeric(lncRNA[,lncrna])
		y=as.numeric(cuproptosis[,gene])
		corT=cor.test(x, y)
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(cuproptosis=gene, lncrna=lncrna, cor, text, pvalue))
	}
}

#绘制相关性热图
outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=7, height=5.6)
ggplot(outTab, aes(lncrna, cuproptosis)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() +    #去掉背景
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),   #x轴字体
	      axis.text.y = element_text(size = 12, face = "bold")) +       #y轴字体
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #设置图例
	scale_x_discrete(position = "bottom")      #定义X轴名称显示的位置
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

