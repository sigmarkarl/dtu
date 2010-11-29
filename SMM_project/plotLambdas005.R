

#SMM

tsmm = read.table("tab_smm005.txt",header=F,sep="\t",as.is=T)
names (tsmm) = c("l","Size","PCC")

max_smm=max(tsmm[,"PCC"])
min_smm=min(tsmm[,"PCC"])

#max_smm=1
#min_smm=0

tsmm = tsmm[order(tsmm[,"Size"]),]

pdf("smm_l005_ppc_size.pdf",width=10,height=10)

plot (x=c(0,0.1),y=c(min_smm,max_smm),xlab="l",ylab="PCC",type="n",main="Pearson correlation in SMM Gradient Decent")
#,main="Pearson correlation trend in SMM on lambda value changing")
nlambda = 8
col=0
smm_small=tsmm[tsmm$Size<100,]
smm_normal=tsmm[tsmm$Size>200 & tsmm$Size<500,]
smm_big=tsmm[tsmm$Size>2000,]

c = 1
while (c < (nrow(smm_small))) {
	c1 = c + (nlambda -1 )
#	lines(spline(smm_small[c:c1,c(1,3)],n=100),type="l",col="red")
	lines(smm_small[c:c1,c(1,3)],type="l",col="red")
	lines(smm_small[c:c1,c(1,3)],type="p",pch=2,col="red")
	c   = c+ nlambda
#	col = col+1
} 
c = 1
while (c < (nrow(smm_normal))) {
	c1 = c + (nlambda -1 )
#	lines(spline(smm_normal[c:c1,c(1,3)],n=100),type="l",col="green")
	lines(smm_normal[c:c1,c(1,3)],type="l",col="green")
	lines(smm_normal[c:c1,c(1,3)],type="p",pch=1,col="green")
	c   = c+ nlambda
#	col = col+1
} 
c = 1
while (c < (nrow(smm_big))) {
	c1 = c + (nlambda -1 )
#	lines(spline(smm_big[c:c1,c(1,3)],n=100),type="l",col="black")
	lines(smm_big[c:c1,c(1,3)],type="l",col="black")
	lines(smm_big[c:c1,c(1,3)],type="p",pch=5,col="black")
	c   = c+ nlambda
#	col = col+1
} 

legend(x=0.06,y=0.45,c("N > 2000","500 > N > 200 "," N < 150"),lty=c(1,1,1),col=c("black","green","red"),pch=c(5,1,2),title="Size of samples",box.lty=0)
dev.off()

#SMM Monte Carlo

tsmm = read.table("tab_smm_mc_005.txt",header=F,sep="\t",as.is=T)
names (tsmm) = c("l","Size","PCC")

max_smm=max(tsmm[,"PCC"])
min_smm=min(tsmm[,"PCC"])

#max_smm=1
#min_smm=0

tsmm = tsmm[order(tsmm[,"Size"]),]

pdf("smm_mc_l005_ppc_size.pdf",width=10,height=10)

plot (x=c(0,0.1),y=c(min_smm,max_smm),xlab="l",ylab="PCC",type="n",main="Pearson correlation in SMM Monte Carlo")
#main="Pearson correlation trend in SMM on lambda value changing")
nlambda = 8
col=0
smm_small=tsmm[tsmm$Size<100,]
smm_normal=tsmm[tsmm$Size>200 & tsmm$Size<500,]
smm_big=tsmm[tsmm$Size>2000,]

c = 1
while (c < (nrow(smm_small))) {
	c1 = c + (nlambda -1 )
#	lines(spline(smm_small[c:c1,c(1,3)],n=100),type="l",col="red")
	lines(smm_small[c:c1,c(1,3)],type="l",col="red")
	lines(smm_small[c:c1,c(1,3)],type="p",pch=2,col="red")
	c   = c+ nlambda
#	col = col+1
} 
c = 1
while (c < (nrow(smm_normal))) {
	c1 = c + (nlambda -1 )
#	lines(spline(smm_normal[c:c1,c(1,3)],n=100),type="l",col="green")
	lines(smm_normal[c:c1,c(1,3)],type="l",col="green")
	lines(smm_normal[c:c1,c(1,3)],type="p",pch=1,col="green")
	c   = c+ nlambda
#	col = col+1
} 
c = 1
while (c < (nrow(smm_big))) {
	c1 = c + (nlambda -1 )
#	lines(spline(smm_big[c:c1,c(1,3)],n=100),type="l",col="black")
	lines(smm_big[c:c1,c(1,3)],type="l",col="black")
	lines(smm_big[c:c1,c(1,3)],type="p",pch=5,col="black")
	c   = c+ nlambda
#	col = col+1
} 

legend(x=0.06,y=0.45,c("N > 2000","500 > N > 200 "," N < 150"),lty=c(1,1,1),col=c("black","green","red"),pch=c(5,1,2),title="Size of samples",box.lty=0)
dev.off()
