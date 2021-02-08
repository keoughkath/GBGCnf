

classify.har <- function(HAR,SIM){

	likes.from.fit <- function(fit){
	#-------------------------------
		#- fix fitting issues; needs to go in SGE code...

		tmp = fit[-1,1] - fit[1,1]
    tmp = round(tmp,2)
		if(any(is.na(tmp))) {return(rep(NA,5))}
		if(tmp[1] < min(tmp[2],tmp[3])) { tmp[1] =  min(tmp[2],tmp[3]) }
		#- throw out messed up fits
		if(tmp[1] < 0) { return(rep(NA,5))}
		if(tmp[4] < 0) { return(rep(NA,5))}
		if(tmp[5] < 0) { return(rep(NA,5))}
		tmp
	}

	#- need to get classification boundaries
	#========================================

	#- extract likelihood scores for simulations under the NULL
	#----------------------------------------------------------
	l.null.s = matrix(unlist(lapply(SIM$S[[1]]$fit, likes.from.fit)),ncol=5,byrow=TRUE)
	colnames(l.null.s) = c("s","s.relax","s.pos","bgc","s.bgc")
	l.null.b = matrix(unlist(lapply(SIM$B[[1]]$fit, likes.from.fit)),ncol=5,byrow=TRUE)
	colnames(l.null.b) = c("s","s.relax","s.pos","bgc","s.bgc")
	l.null   = rbind(l.null.s, l.null.b)
	na.null  = sum(is.na(rowSums(l.null)))

	alpha = 0.05
	da.0  = quantile(l.null[,1],1-alpha, na.rm=TRUE)
	db.0  = quantile(l.null[,4],1-alpha, na.rm=TRUE)
	dab.0 = quantile(l.null[,5],1-alpha, na.rm=TRUE)

	#- extract likelihood scores for simulations under gBGC and selection
	#--------------------------------------------------------------------
	apfu <- function(x) {
	#--------------------
		likes = matrix(unlist(lapply(x$fit, likes.from.fit)),ncol=5,byrow=TRUE)
		colnames(likes) = c("s","s.relax","s.pos","bgc","s.bgc")
		return(list(likes=likes,val=x$val))
	}
	l.gbgc = lapply(SIM$B[-1],apfu )
	l.sel  = lapply(SIM$S[-1],apfu )
	na.s = sum(unlist(lapply(l.sel, function(x) sum(is.na(rowSums(x$likes))))))
	na.b = sum(unlist(lapply(l.gbgc,function(x) sum(is.na(rowSums(x$likes))))))

	da.bs  = unlist(lapply(l.gbgc, function(x) quantile(x$likes[,1]-x$likes[,4],1-alpha, na.rm=TRUE)))
	dab.bs = unlist(lapply(l.gbgc, function(x) quantile(x$likes[,5]-x$likes[,4],1-alpha, na.rm=TRUE)))

	da.b  = max(c(da.bs,0))
	dab.b = max(c(dab.bs,0))

	db.as  = unlist(lapply(l.sel, function(x) quantile(x$likes[,4]-x$likes[,1],1-alpha, na.rm=TRUE)))
        dab.as = unlist(lapply(l.sel, function(x) quantile(x$likes[,5]-x$likes[,1],1-alpha, na.rm=TRUE)))

	db.a = max(c(db.as,0))
	dab.a = max(c(dab.as,0))

	#- Assign the class of the HAR
	#=============================

	assign.class <- function(la,lb,lab,db.0, da.0, dab.0, da.b, db.a, dab.b, dab.a){
	#-------------------------------------------------------------------------------

		cuts  = c(db.0, da.0, dab.0,  da.b, db.a , dab.b , dab.a )
		vals  = c(lb  ,   la,   lab, la-lb, lb-la, lab-lb, lab-la)
		comps = vals > cuts
		names(comps) = c("b>0","a>0","ab>0","a>b","b>a","ab>b","ab>a")

		#- the easy cases, only one model rejects null
		#---------------------------------------------
		if(all(comps[1:3] == c(FALSE,FALSE,FALSE))) return(c(0,comps))   #- null
		if(all(comps[1:3] == c(TRUE ,FALSE,FALSE))) {
			if(comps[5] == TRUE  ) 		    return(c(22,comps))  #- bgc + high-conf
			if(comps[5] == FALSE )		    return(c(20,comps))  #- bgc
		}
		if(all(comps[1:3] == c(FALSE,TRUE ,FALSE))) return(c(1,comps))   #- sel
		if(all(comps[1:3] == c(FALSE,FALSE,TRUE ))) return(c(3,comps))   #- bgc & sel

		#- two models reject null
		#------------------------

		#- bgc and sel reject null
		if(all(comps[1:3] == c(TRUE ,TRUE ,FALSE))){
			if(all(comps[4:5] == c(TRUE ,FALSE))) return(c(1,comps))   #- sel
			if(all(comps[4:5] == c(FALSE,TRUE ))) return(c(22,comps))  #- bgc + high-conf
			if(all(comps[4:5] == c(FALSE,FALSE))) return(c(20,comps))  #- bgc
			#if(comps[4:5] == c(TRUE ,FALSE)) cannot be
		}

		#- sel and sel&bgc reject null
		if(all(comps[1:3] == c(FALSE ,TRUE ,TRUE ))){
			if(comps[7] == TRUE ) return(c(3,comps)) #- bgc & sel
			if(comps[7] == FALSE) return(c(1,comps)) #- sel
                }

		#- bgc and sel&bgc reject null
		if(all(comps[1:3] == c(TRUE ,FALSE,TRUE ))){
                        if(comps[6] == TRUE ) return(c(3,comps)) #- bgc & sel
                        if(comps[6] == FALSE) {
				if(comps[5] == TRUE  ) return(c(22,comps))
				if(comps[5] == FALSE ) return(c(20,comps))
			}
                }

		#- all three models reject null
		#------------------------------
		if(all(comps[1:3] == c(TRUE ,TRUE ,TRUE ))){
			#- sel rejects bgc
			if(comps[4] == TRUE){
				if(comps[7] == TRUE ) return(c(3,comps)) #- bgc & sel
				if(comps[7] == FALSE) return(c(1,comps)) #- sel
			}
			#- bgc rejects sel
			if(comps[5] == TRUE){
				if(comps[6] == TRUE ) return(c(3,comps))  #- bgc & sel
				if(comps[6] == FALSE) return(c(22,comps)) #- bgc + high-conf
			}
			#- no bgc vs sel decision
			if(all(comps[4:5] == c(FALSE,FALSE))){
				if(all(comps[6:7] == c(TRUE ,TRUE ))) return(c(3,comps))
				if(all(comps[6:7] == c(TRUE ,FALSE))) return(c(3,comps))  #- cons
				if(all(comps[6:7] == c(FALSE,TRUE ))) return(c(20,comps)) #- cons
				if(all(comps[6:7] == c(FALSE,FALSE))) return(c(20,comps)) #- cons
			}
		}

		return(NA)
	}


	l0    = HAR$fit$summary[1,1]
	la    = HAR$fit$summary[2,1]- l0
	la.r  = HAR$fit$summary[3,1]- l0
	la.p  = HAR$fit$summary[4,1]- l0
	lb    = HAR$fit$summary[5,1]- l0
	lab   = HAR$fit$summary[6,1]- l0

	class = assign.class(la,lb,lab,db.0,da.0,dab.0,da.b,db.a,dab.b,dab.a)

	#- refine the selection class
	#============================

	SG = abs(HAR$fit$summary[1,3])
	sel.neu.inds = which(unlist(lapply(SIM$S , function(x) x$val <= SG)))

	apfu <- function(x){
	#-------------------
		#- slightly different, need index here
		quantile( l.sel[[x]]$likes[,3] -  l.sel[[x]]$likes[,2], .95,na.rm=TRUE)
	}
	dsps = sapply(sel.neu.inds,apfu)
	dsp  = max(dsps)
	if(class[1] ==1){
		if(la.p - la.r > dsp){
			class[1] = 11
		}
		else{
			class[1] = 10
		}
	}

	#- parameter estimates
	#=====================
	if(class[1]==0){
		est=HAR$fit$summary[1,]
	}
	if(class[1]==22 | class[1] == 20){
		est=HAR$fit$summary[5,]
	}
	if(class[1]==3){
		est=HAR$fit$summary[6,]
	}
	if(class[1]==10){
		est=HAR$fit$summary[3,]
	}
	if(class[1]==11){
		est=HAR$fit$summary[4,]
	}

  return(list(class=class,est=est))

}
