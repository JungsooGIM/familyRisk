#	Title: R pakcage for evaluating disease risk (conditional mean of being affected, or CM) using family history
#	Developed by: Jungsoo Gim (iedenkim@gmail.com)

#' A function evaluting a conditional mean (CM) of being affected of the disease using family history
#' @param  fam	Input file in FAM format of PLINK. Target individuals and their relatives should be all included in the FAM file.
#' @param  h	Heritability of the disease
#' @param  prev	Prevalence of the disease (should be specified according to the population)
#' @param  target	Target individuals whose conditional mean being evaluated. This should be a index vector in the samples in the FAM file
#' @details Please see the reference paper.
#' @return  \item{res_rp_I}{a matrix including IID and CM} 
#' @export
#' @author Jungsoo Gim
#' @references Jungsoo Gim, 2017, Genetics, "Incorporating family disease history in risk prediction models with large-scale genetic data substantially reduces unexplained variation."
#' @examples
#' (you example code, each line start with #')
cal_rp <- function(fam,h,prev,target){

	library(tmvtnorm)
	library(kinship2)

	target_fam <- fam[target,]
	famid <- unique(as.character(target_fam[,1]))
	n <- length(famid)

	fam.info <- function(i){
		a <- which(as.character(fam[,1])==famid[i])
	}

	resI <- lapply(1:n,fam.info)

	res_rp <- lapply(1:n, cal_rp_fam, fam=fam, resI=resI, target=target, h=h, prev=prev)
	res_rp_I <- do.call(rbind,res_rp)
	res_rp_return <- data.frame(ID = res_rp_I[,1], ConditionalMean = as.numeric(res_rp_I[,2]))
	return(res_rp_return)
}


pack_res <- function(ind,L,h,prev,V,fam,resI,i){

	famV <- h*V+(1-h)*diag(length(L))

	k <- which(L==2)
	kk <- which(L==1)
	t <- qnorm(prev,0,1,lower.tail=F)

	aa <- rep(-Inf,length(L)); aa[k] <- t
	bb <- rep(Inf,length(L)); bb[kk] <- t
	aa[ind] <- -Inf; bb[ind] <- Inf

	para <- mtmvnorm(mean=rep(0,length(L)),sigma=famV,lower=aa,upper=bb,doComputeVariance=FALSE)
	mu <- para$tmean[ind]
	var <- para$tvar[ind,ind]
	prob <- pnorm((t-mu)/sqrt(var),lower.tail=F)

	return(matrix(c(fam[resI[[i]][ind],2],mu),length(ind),2))
}

cal_rp_fam <- function(i,fam,resI,target,h,prev){

	print(paste("Evaluating conditional mean for ",i,"th family ... ",sep=""))

	new_fam <- fam[resI[[i]],]
	if(nrow(new_fam)==0){
		stop("Check the number of IDs")
	}
	else if(nrow(new_fam)==1){
		RP <- c(new_fam[,2], 0)
	}
	else{
		ped <- with(new_fam,pedigree(id=V2,dadid=V3,momid=V4,sex=V5,missid='0'))
		V <- 2*as.matrix(kinship(ped))[!new_fam[,6]==-9,!new_fam[,6]==-9]
		# V <- 2*as.matrix(kinship(ped))
		ind <- which(resI[[i]] %in% target)
		L <- new_fam[!new_fam[,6]==-9,6]
		# L <- new_fam[,6]
		RP <- pack_res(ind,L,h,prev,V,fam,resI,i)		
	}
	return(RP)
}
