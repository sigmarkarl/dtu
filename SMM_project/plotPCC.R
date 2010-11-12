pccC= read.table("pccC",header=F,sep=" ",as.is=T)
pccF= read.table("pccf",header=F,sep=" ",as.is=T)
errC= read.table("errC",header=F,sep=" ",as.is=T)
errF= read.table("errF",header=F,sep=" ",as.is=T)


brp = read.table("barplot.txt",header=T,sep="\t")
vetbp= as.numeric(as.matrix(brp[1,]))

pdf("barplot.pdf",height=10,width=10)
barplot(vetbp,col=c("blue","red","blue","red","blue","red","blue","red","blue","red","cyan"),ylim=c(0,1),xlim=c(0,17),space=c(0,0,1,0,1,0,1,0,1,0,1),names=c("Test 1","","Test 2","","Test 3","","Test 4","","Test 5","","Evaluation"),)
legend(x=8,y=1,lty=c(1,1,1),lwd=c(5,5,5),col=c("blue","red","cyan"),legend=c("Training","Evaluation of Trainig","Concatenated Evaluation"))
dev.off()

