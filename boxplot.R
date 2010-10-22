

#SMM  GD vs MC boxplot 

tsmm = read.table("tab_smm005.txt",header=F,sep="\t",as.is=T)
names (tsmm) = c("l","Size","PCC")
tsmm_mc = read.table("tab_smm_mc_005.txt",header=F,sep="\t",as.is=T)
names (tsmm_mc) = c("l","Size","PCC")

tsmm_sub = tsmm[tsmm$l==0.01,]
tsmm_mc_sub = tsmm_mc[tsmm_mc$l==0.005,]

pdf("boxplot.pdf",width=10,height=10)
boxplot(tsmm_mc$PCC,tsmm$PCC)
dev.off()
