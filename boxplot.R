

#SMM  GD vs MC boxplot 

tsmm = read.table("tab_smm005.txt",header=F,sep="\t",as.is=T)
names (tsmm) = c("l","Size","PCC")
tsmm_mc = read.table("tab_smm_mc_005.txt",header=F,sep="\t",as.is=T)
names (tsmm_mc) = c("l","Size","PCC")

tsmm_sub = tsmm[tsmm$l==0.02,]
tsmm_mc_sub = tsmm_mc[tsmm_mc$l==0.005,]

pdf("boxplot.pdf",width=10,height=10)
boxplot(tsmm_mc_sub$PCC,tsmm_sub$PCC,widh=0.1,boxwex=0.5,col=c("white","grey"),outline=F,names=c("Monte Carlo","Gradient decent"),ylab="PCC")
dev.off()

tsmm_1 = tsmm[tsmm$l==0,]
tsmm_2 = tsmm[tsmm$l==0.005,]
tsmm_3 = tsmm[tsmm$l==0.01,]
tsmm_4 = tsmm[tsmm$l==0.02,]
tsmm_5 = tsmm[tsmm$l==0.03,]
tsmm_6 = tsmm[tsmm$l==0.04,]
tsmm_7 = tsmm[tsmm$l==0.05,]
tsmm_8 = tsmm[tsmm$l==0.1,]

tab_gd = cbind(tsmm_1$PCC,tsmm_2$PCC,tsmm_3$PCC,tsmm_4$PCC,tsmm_5$PCC,tsmm_6$PCC,tsmm_7$PCC,tsmm_8$PCC)

pdf("choice_lambda.pdf",width=10,height=10)
boxplot(tab_gd,widh=0.1,boxwex=0.5,col="grey",outline=F,names=c("0.0","0.005","0.01","0.02","0.03","0.04","0.05","0.1"),ylab="PCC",xlab="lambda",main="SMM gradient decent")
dev.off()

tsmm_mc_1 = tsmm_mc[tsmm_mc$l==0,]
tsmm_mc_2 = tsmm_mc[tsmm_mc$l==0.005,]
tsmm_mc_3 = tsmm_mc[tsmm_mc$l==0.01,]
tsmm_mc_4 = tsmm_mc[tsmm_mc$l==0.02,]
tsmm_mc_5 = tsmm_mc[tsmm_mc$l==0.03,]
tsmm_mc_6 = tsmm_mc[tsmm_mc$l==0.04,]
tsmm_mc_7 = tsmm_mc[tsmm_mc$l==0.05,]
tsmm_mc_8 = tsmm_mc[tsmm_mc$l==0.1,]

tab_mc = cbind(tsmm_mc_1$PCC,tsmm_mc_2$PCC,tsmm_mc_3$PCC,tsmm_mc_4$PCC,tsmm_mc_5$PCC,tsmm_mc_6$PCC,tsmm_mc_7$PCC,tsmm_mc_8$PCC)

pdf("choice_lambda_mc.pdf",width=10,height=10)
boxplot(tab_mc,widh=0.1,boxwex=0.5,col="grey",outline=F,names=c("0.0","0.005","0.01","0.02","0.03","0.04","0.05","0.1"),ylab="PCC",xlab="lambda",main="SMM Monte Carlo")
dev.off()
