PSSM <- read.table("tab_pssm",header=F,sep="\t",as.is=T)
names(PSSM) <- c("Param","Allele","Sample","Size","PCC")

NN <- read.table("tab_nn_new",header=F,sep="\t",as.is=T)
NN = (na.exclude(NN))
names(NN) <- c("Param","Allele","Sample","Size","PCC")

plot_algo <- function (tab=NULL,exp=FALSE,approx=FALSE,xbyname=FALSE) 
{

	if (exp == TRUE) {
		tab$PCC = exp(tab$PCC)
	} 

	min_pcc <- min(tab$PCC)
	max_pcc <- max(tab$PCC)
	min_size <- min(tab$Size)
	max_size <- max(tab$Size)
	tab_s <- unique(tab$Sample)
	tab_p <- unique(tab$Param)
	tab_lp <- length(tab_p)
	if (xbyname==TRUE) {
		min_size = 0
		max_size = length(unique(tab$Allele))
	}	
	plot(x=c(min_size,max_size),y=c(0.4,max_pcc),type="n",xlab="Size",ylab="PPC")
	legend(x="bottomright",legend=tab_p,lty=c(1:tab_lp),col=c(1:tab_lp),bty="n")
	
	for (i in 1:tab_lp) {

		tab_part <- tab[tab$"Param"== tab_p[i],]
		tab_part <- tab_part[order(tab_part$Size),]
		xvals = tab_part$Size
		tab_lines = matrix(nrow=35,ncol=2,0)
		if (xbyname==TRUE) {
			xvals = c(1:max_size)
		}
#		if (approx==TRUE) {
#			library(akima) ## Used to interpolate the points...(can be better)
#			tab_lines = aspline(x=xvals,y=tab_part$PCC)
#		} else {
		for (s in 1:length(tab_s)) {
			tab_lines = tab_lines + cbind(x=xvals,y=tab_part$PCC[tab_part$Sample==tab_s[s]])
#		}
		}
		lines(tab_lines/length(tab_s),col=i,lty=i)

	} 
}

boxplot_algo <- function (tab=NULL) 
{

	min_pcc <- min(tab$PCC)
	max_pcc <- max(tab$PCC)
	alleles <- unique(tab$Allele)
	alleles <- cbind(Allele = alleles, nulls="") 
	samp <- unique(tab$Sample)
	tab_p <- unique(tab$Param)
	tab_lp <- length(tab_p)
	
	for (i in 1:tab_lp) {

		tab_part <- tab[tab$"Param"== tab_p[i],]
		tab_part <- tab_part[order(tab_part$Size),]
		bxpobj <- as.numeric(tab_part$PCC[tab_part$Sample==samp[1]])
		if (length(samp) > 1) {
			barpos <- c(1:tab_lp)
			barpos[1]=0
			for (s in 2:tab_lp) {
			barpos[s] <- barpos[s-1]+1/tab_lp
			}
			barpos - median(barpos)
			bwex = 1/tab_lp*3/4
			for (l in 2:length(samp)) {
				bxpobj <- cbind(bxpobj,as.numeric(tab_part$PCC[tab_part$Sample==samp[l]]))
			}
			bxcolor=i
			if (i == 1) {
				names_x = unique(tab_part$Allele)
				bxcolor = 0 	
				boxplot(t(bxpobj),boxwex = bwex,at = 1:nrow(bxpobj) + barpos[i],col=bxcolor, names = names_x,outline=F,ylab="PCC")
			

## to reduce the names size on x is possible to add the parameter to boxplot like : cex.axis=0.2
			} else {
				boxplot(t(bxpobj),boxwex = bwex,at = 1:nrow(bxpobj) + barpos[i],col=bxcolor,add=T, names = alleles[,"nulls"],outline=F)
			}
		}
	}
	#return(bxpobj)
	legend(x="bottomright",unique(tab$Param),fill=c(0,2:tab_lp),bty="n")
}

