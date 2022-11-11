######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("ggpubr")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


#引用包
library(limma)
library(ggpubr)
gene="CD274"          #基因的标准名字
showName="PD-L1"      #图形里面显示的基因名称
setwd("D:\\biowolf\\irgICI\\29.riskGene")      #设置工作目录

#定义绘制图形的函数
riskGene=function(riskFile=null, expFile=null, boxFile=null){
	#读取表达数据文件
	rt=read.table(expFile, header=T, sep="\t", check.names=F)
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	
	#删掉正常样品
	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	group=gsub("2", "1", group)
	data=data[,group==0]
	
	#提取目标基因表达量
	data=rbind(data, gene=data[gene,])
	exp=t(data[c("gene",gene),])
	row.names(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\",  row.names(exp))
	exp=avereps(exp)
	
	#读取风险数据文件
	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	
	#合并数据
	sameSample=intersect(row.names(exp), row.names(risk))
	exp=exp[sameSample,]
	exp[exp>quantile(exp,0.975)]=quantile(exp,0.975)
	risk=risk[sameSample,]
	data=cbind(as.data.frame(exp), as.data.frame(risk))
	
	#设置比较组
	data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
	group=levels(factor(data$risk))
	data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
	comp=combn(group,2)
	my_comparisons=list()
	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
	#绘制boxplot
	boxplot=ggboxplot(data, x="risk", y="gene", color="risk",
			          xlab="",
			          ylab=paste(showName, "expression"),
			          legend.title="",
			          palette = c("blue", "red"),
			          add = "jitter")+ 
		    stat_compare_means(comparisons = my_comparisons)
	
	#输出图片
	pdf(file=boxFile, width=5, height=4.5)
	print(boxplot)
	dev.off()
}

riskGene(riskFile="risk.all.txt", expFile="normalize.txt", boxFile=paste0(gene,".pdf"))


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
