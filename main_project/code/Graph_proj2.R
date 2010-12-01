PSSM <- read.table("tab_pssm",header=F,sep="\t",as.is=T)
names(PSSM) <- c("Param","Allele","Sample","Size","PCC")

NN <- read.table("tab_nn",header=F,sep="\t",as.is=T)
names(NN) <- c("Param","Allele","Sample","Size","PCC")




plot_algo <- function (tab=NULL) 
{
	library(akima) ## Used to interpolate the points...(can be better)
	min_pcc <- min(tab$PCC)
	max_pcc <- max(tab$PCC)
	min_size <- min(tab$Size)
	max_size <- max(tab$Size)

	tab_p <- unique(tab$Param)
	tab_lp <- length(tab_p)

	plot(x=c(min_size,max_size),y=c(min_pcc,max_pcc),type="n",xlab="Size",ylab="PPC")

	legend(x=500,y=0,legend=tab_p,lty=c(1:tab_lp),col=c(1:tab_lp))

	for (i in 1:tab_lp) {

		tab_part <- tab[tab$"Param"== tab_p[i],]
		tab_part <- tab_part[order(tab_part$Size),]
		lines(aspline(x=tab_part$Size,y=tab_part$PCC),col=i,lty=i)
	}
}
