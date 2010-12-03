PSSM <- read.table("tab_pssm",header=F,sep="\t",as.is=T)
names(PSSM) <- c("Param","Allele","Sample","Size","PCC")

NN <- read.table("~/dtu/main_project/tabs/tab_nn",header=F,sep="\t",as.is=T)
names(NN) <- c("Param","Allele","Sample","Size","PCC")

plot_algo <- function (tab=NULL,exp=FALSE,approx=FALSE) 
{

	if (exp == TRUE) {
		tab$PCC = exp(tab$PCC)
	} 

	min_pcc <- min(tab$PCC)
	max_pcc <- max(tab$PCC)
	min_size <- min(tab$Size)
	max_size <- max(tab$Size)

	tab_p <- unique(tab$Param)
	tab_lp <- length(tab_p)
	
	plot(x=c(min_size,max_size),y=c(min_pcc,max_pcc),type="n",xlab="Size",ylab="PPC")
	legend(x="bottomright",legend=tab_p,lty=c(1:tab_lp),col=c(1:tab_lp),bty="n")
	
	for (i in 1:tab_lp) {

		tab_part <- tab[tab$"Param"== tab_p[i],]
		tab_part <- tab_part[order(tab_part$Size),]
		
		if (approx==TRUE) {
			library(akima) ## Used to interpolate the points...(can be better)
			tab_lines = aspline(x=tab_part$Size,y=tab_part$PCC)
		} else {
			tab_lines = cbind(x=tab_part$Size,y=tab_part$PCC)
		}
		lines(tab_lines,col=i,lty=i)
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
			if (i == 1) {
				boxplot(t(bxpobj),boxwex = bwex,at = 1:nrow(bxpobj) + barpos[i],col=i, names = alleles[,"Allele"],outline=F,ylab="PCC")
			

## to reduce the names size on x is possible to add the parameter to boxplot like : cex.axis=0.2
			} else {
				boxplot(t(bxpobj),boxwex = bwex,at = 1:nrow(bxpobj) + barpos[i],col=i,add=T, names = alleles[,"nulls"],outline=F)
			}
		}
	}
	#return(bxpobj)
	legend(x="bottomright",unique(tab$Param),fill=c(1:tab_lp),bty="n")
}

