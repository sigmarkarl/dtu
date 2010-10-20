tsmm = read.table("tab_smm.txt",header=F,sep="\t",as.is=T)
tsmm_mc = read.table("tab_smm_mc.txt",header=F,sep="\t",as.is=T)
names (tsmm) = c("l","Size","PCC")
names (tsmm_mc) = c("l","Size","PCC")

#SMM
max_smm=max(tsmm[,"PCC"])
min_smm=min(tsmm[,"PCC"])

#max_smm=1
#min_smm=0

tsmm = tsmm[order(tsmm[,"Size"]),]

plot (x=c(0,0.1),y=c(min_smm,max_smm),xlab="l",ylab="PCC",type="n")
nlambda = 6
col=0
smm_small=tsmm[tsmm$Size<100,]
smm_normal=tsmm[tsmm$Size>200 & tsmm$Size<500,]
smm_big=tsmm[tsmm$Size>2000,]

c = 1
while (c < (nrow(smm_small))) {
	c1 = c + (nlambda -1 )
	lines(spline(smm_small[c:c1,c(1,3)],n=100),type="l",col="red")
	c   = c+6
#	col = col+1
} 
c = 1
while (c < (nrow(smm_normal))) {
	c1 = c + (nlambda -1 )
	lines(spline(smm_normal[c:c1,c(1,3)],n=100),type="l",col="green")
	c   = c+6
#	col = col+1
} 
c = 1
while (c < (nrow(smm_big))) {
	c1 = c + (nlambda -1 )
	lines(spline(smm_big[c:c1,c(1,3)],n=100),type="l",col="black")
	c   = c+6
#	col = col+1
} 
