
fit.har	= function(	mods				,
			ali				,
			n.bound  = "sel[-10,100]"	,
			a.bound  = "sel[-10,100]"	,
			b.bound  = "sel[-10,100]"	,
			ab.bound = "sel[-10,100]"	,
			fitBG    = FALSE		){
#=========================================================

	if(fitBG){
                no.opt.str = c("ratematrix","branches")
        } else{
                no.opt.str = c("ratematrix","branches","backgd")
        }


	#- fit the null model
	#====================
	res.null        = phyloFit(     ali,
                                        init.mod=mods$null,
                                        no.opt=no.opt.str,
                                        quiet=TRUE,
                                        bound=n.bound)
        bl.null         = summary.tree(res.null$tree)[ which(summary.tree(res.null$tree)[,"name"]=="hg38")  ,"tparent"]
        bl.null         = bl.null*(diag(res.null$rate.matrix) %*% res.null$backgd)*(-1)

	#- fit the selection model
	#=========================
	if(mods$sel$selection==0 & mods$sel$ls.model$selection == 0){
		mods$sel$selection = res.null$selection * .95
		mods$sel$rate.matrix = apply.bgc.sel(mods$sel$rate.matrix,s=mods$sel$selection)
		mods$sel$ls.model$rate.matrix = mods$sel$rate.matrix
	}
	res.s           = phyloFit(     ali,
                                        init.mod=mods$sel,
                                        no.opt=no.opt.str,
                                        quiet=TRUE,
                                        bound=a.bound)
        bl.s            = summary.tree(res.s$tree)[ which(summary.tree(res.s$tree)[,"name"]=="hg38")  ,"tparent"]
        bl.s            = bl.s*(diag(res.s$ls.model$rate.matrix) %*% res.s$backgd)*(-1)

	#- relaxation of constraint
	#--------------------------
	if(mods$sel.r$ls.model[[1]]$selection== 0){
                mods$sel.r$ls.model[[1]]$selection = res.null$selection
                mods$sel.r$ls.model[[1]]$rate.matrix = apply.bgc.sel(mods$sel.r$ls.model[[1]]$rate.matrix,s=res.null$selection)
		mods$sel.r$ls.model[[2]]$selection = res.null$selection
                mods$sel.r$ls.model[[2]]$rate.matrix = apply.bgc.sel(mods$sel.r$ls.model[[2]]$rate.matrix,s=res.null$selection)
        }
	res.s.r         = phyloFit(     ali,
                                        init.mod=mods$sel.r,
                                        no.opt=no.opt.str,
                                        quiet=TRUE)
        #- need to make sure that we have no smaller-than-other human brach
        #- and fudge it if we do
        if(res.s.r$ls.model[[2]]$selection < res.s.r$ls.model[[1]]$selection){
                res.s.r$ls.model[[2]]$selection   = res.s.r$ls.model[[1]]$selection
                res.s.r$ls.model[[2]]$rate.matrix = apply.bgc.sel(     res.s.r$rate.matrix,
                                                                        s=res.s.r$ls.model[[2]]$selection)
                res.s.r$likelihood = likelihood.msa(ali,res.s.r)
        }
	bl.s.r          = summary.tree(res.s.r$tree)[ which(summary.tree(res.s.r$tree)[,"name"]=="hg38")  ,"tparent"]
        bl.s.r          = bl.s.r*(diag(res.s.r$ls.model[[2]]$rate.matrix) %*% res.s.r$backgd)*(-1)

	#- positive selection
	#--------------------
	if(mods$sel.p$ls.model[[1]]$selection== 0){
                mods$sel.p$ls.model[[1]]$selection = res.null$selection
                mods$sel.p$ls.model[[1]]$rate.matrix = apply.bgc.sel(mods$sel.p$ls.model[[1]]$rate.matrix,s=res.null$selection)
	}
	res.s.p         = phyloFit(     ali,
                                        init.mod=mods$sel.p,
                                        no.opt=no.opt.str,
                                        quiet=TRUE)
        bl.s.p          = summary.tree(res.s.p$tree)[ which(summary.tree(res.s.p$tree)[,"name"]=="hg38")  ,"tparent"]
        bl.s.p          = bl.s.p*(diag(res.s.p$ls.model[[2]]$rate.matrix) %*% res.s.p$backgd)*(-1)

	#- fit the gBGC only model
	#=========================

	 res.b           = phyloFit(    ali,
                                        init.mod=mods$bgc,
                                        no.opt=no.opt.str,
                                        quiet=TRUE,
                                        bound=b.bound)
        bl.b            = summary.tree(res.b$tree)[ which(summary.tree(res.b$tree)[,"name"]=="hg38")  ,"tparent"]
        bl.b            = bl.b*(diag(res.b$ls.model$rate.matrix) %*% res.b$backgd)*(-1)



        res.sb          = phyloFit(     ali,
                                        init.mod=mods$both,
                                        no.opt=no.opt.str,
                                        quiet=TRUE,
                                        bound=ab.bound)
        bl.sb           = summary.tree(res.sb$tree)[ which(summary.tree(res.sb$tree)[,"name"]=="hg38")  ,"tparent"]
        bl.sb           = bl.sb*(diag(res.sb$ls.model$rate.matrix) %*% res.sb$backgd)*(-1)

	#- make the summary
	#==================

	res = matrix(0,nrow=6,ncol=5)

        colnames(res) = c(      "likelihood"    ,
                                "hg-branchlen"  ,
                                "shat2"         ,
                                "shat1"         ,
                                "bhat"          )

        rownames(res) = c(      "null"          ,
                                "sel"           ,
                                "sel.r"         ,
                                "sel.p"         ,
                                "bgc"           ,
                                "sel.bgc"       )

        res[1,1] = res.null$likelihood
        res[1,2] = bl.null
        res[1,3] = res.null$selection
        res[1,4] = res.null$selection

        res[2,1] = res.s$likelihood
        res[2,2] = bl.s
        res[2,3] = res.s$selection
        res[2,4] = res.s$selection + res.s$ls.model$selection

        res[3,1] = res.s.r$likelihood
        res[3,2] = bl.s.r
        res[3,3] = res.s.r$ls.model[[1]]$selection
        res[3,4] = res.s.r$ls.model[[2]]$selection

        res[4,1] = res.s.p$likelihood
        res[4,2] = bl.s.p
        res[4,3] = res.s.p$ls.model[[1]]$selection
        res[4,4] = res.s.p$ls.model[[2]]$selection

        res[5,1] = res.b$likelihood
        res[5,2] = bl.b
        res[5,3] = res.b$selection
        res[5,4] = res.b$selection #+ res.b$ls.model$selection
        res[5,5] = res.b$ls.model$bgc

        res[6,1] = res.sb$likelihood
        res[6,2] = bl.sb
        res[6,3] = res.sb$selection
        res[6,4] = res.sb$selection + res.sb$ls.model$selection
        res[6,5] = res.sb$ls.model$bgc

  	MODS = list(    null    = res.null,
                        sel     = res.s,
                        sel.r   = res.s.r,
                        sel.p   = res.s.p,
                        bgc     = res.b,
                        both    = res.sb        )
        return(list(summary=res, mods=MODS))


}



makeRanges <- function(summary){

        makeDel <- function(x){
                del              = abs(x) * 0.5
                if(del < .1) del = 0.1
                del
        }

	#- bounds for global selection (for fitting)
        tmp         = round(summary[1,3],1)
        del         = makeDel(tmp)
        n.bound     = paste("sel[",tmp-del ,",",tmp+del,"]",sep="")
        tmp         = round(summary[2,3],1)
        del         = makeDel(tmp)
        a.bound     = paste("sel[",tmp-del ,",",tmp+del,"]",sep="")
        tmp         = round(summary[5,3],1)
        del         = makeDel(tmp)
        b.bound     = paste("sel[",tmp-del ,",",tmp+del,"]",sep="")
        tmp         = round(summary[6,3],1)
        del         = makeDel(tmp)
        ab.bound    = paste("sel[",tmp-del ,",",tmp+del,"]",sep="")


	#- selection, bgc and both model ranges (for model creation)
        tmp         	= round(summary[2,4] - summary[2,3],1)
        del         	= makeDel(tmp)
        a.shg.range 	= paste("sel[",max(0,tmp-del),",",tmp+del,"]",sep="")
        tmp         	= round(summary[5,5],1)
        del         	= makeDel(tmp)
        b.bhg.range     = paste("bgc[",max(0,tmp-del) ,",",tmp+del,"]",sep="")
        tmp1        	= round(summary[6,4] - summary[6,3],1)
        del1        	= makeDel(tmp1)
        tmp2        	= round(summary[6,5],1)
        del2        	= makeDel(tmp2)
        ab.sbhg.range   = paste("sel[",max(0,tmp1-del1),",",tmp1+del1,"]",sep="")
        ab.sbhg.range   = paste(ab.sbhg.range, paste("bgc[",max(0,tmp2-del2),",",tmp2+del2,"]",sep=""), sep=",")

	#- the relaxation and pos selection ranges
        tmp         	= round(summary[3,3],1)
        del         	= makeDel(tmp)
        ar.sg.range   	= paste("sel[",max(-10,tmp-del),",",min(0,tmp+del),"]",sep="")
        tmp         	= round(summary[3,4],1)
        del         	= makeDel(tmp)
        ar.shg.range	= paste("sel[",max(-10,tmp-del),",",min(0,tmp+del),"]",sep="")

        tmp         	= round(summary[4,3],1)
        del         	= makeDel(tmp)
        ap.sg.range   	= paste("sel[",max(-10,tmp-del),",",min(0,tmp+del),"]",sep="")
        tmp         	= round(summary[4,4],1)
        del         	= makeDel(tmp)
        ap.shg.range   	= paste("sel[",max(0,tmp-del),",",min(100,tmp+del),"]",sep="")

	 return(rbind(	n.bound,
                        a.bound,
                        b.bound,
                        ab.bound,
                        a.shg.range,
                        b.bhg.range,
                        ab.sbhg.range,
                        ar.sg.range,
                        ar.shg.range,
                        ap.sg.range,
                        ap.shg.range))
}






fm2 <- function(        mod,
                        msa,
                        fitBG=FALSE){
#======================================

	tmp.mod = models.har(mod)
	tmp = fit.har(tmp.mod,msa,fitBG=fitBG)

        initvals = list(        n.sg 	 = tmp$summary[1,3] * (1+ 0.05*(runif(1)-.5)),
                                a.sg   	 = tmp$summary[2,3] * (1+ 0.05*(runif(1)-.5))   ,
                                a.shg    = (tmp$summary[2,4]- tmp$summary[2,3]) * (1+ 0.05*(runif(1)-.5)) ,
                                b.sg     = tmp$summary[5,3] * (1+ 0.05*(runif(1)-.5)),
                                b.bhg    = tmp$summary[5,5] * (1+ 0.05*(runif(1)-.5)),
                                ar.sg    = tmp$summary[3,3] * (1+ 0.05*(runif(1)-.5)),
                                ar.shg   = tmp$summary[3,4] * (1+ 0.05*(runif(1)-.5)),
                                ap.sg    = tmp$summary[4,3] * (1+ 0.05*(runif(1)-.5)),
                                ap.shg   = tmp$summary[4,4] * (1+ 0.05*(runif(1)-.5)),
                                ab.sg  	 = tmp$summary[6,3] * (1+ 0.05*(runif(1)-.5)),
                                ab.shg   = (tmp$summary[6,4]- tmp$summary[6,3]) * (1+ 0.05*(runif(1)-.5)),
                                ab.bgc   = tmp$summary[6,5]  * (1+ 0.05*(runif(1)-.5)))

        rr      = makeRanges(tmp$summary)

	tmp.mod = models.har(	mod,
				a.shg.range     = rr["a.shg.range",]               ,
                        	b.bhg.range     = rr["b.bhg.range",]               ,
                        	ab.sbhg.range   = rr["ab.sbhg.range",]    ,
                        	ar.sg.range     = rr["ar.sg.range",]               ,
                        	ar.shg.range    = rr["ar.shg.range",]               ,
                        	ap.sg.range     = rr["ap.sg.range",]               ,
                        	ap.shg.range    = rr["ap.shg.range",]     	,
                        	initvals        = initvals			)

	res = fit.har(	tmp.mod,
			msa,
			n.bound = rr["n.bound",],
			a.bound = rr["a.bound",],
			ab.bound = rr["ab.bound",],
			fitBG=fitBG)

        return(tmp)

}


fit.neutral = function( mod,
                        ali){
#============================

        #- scale-only
        res.neutral.s          = phyloFit(      ali                             ,
                                                init.mod=mod                    ,
                                                scale.only=TRUE                 ,
                                                no.opt=c("ratematrix")          ,
                                                quiet=TRUE                      )
        #- scale subtree
        res.neutral.t           = phyloFit(     ali                             ,
                                                init.mod=mod                    ,
                                                scale.only=TRUE                 ,
                                                scale.subtree="hg38:loss"       ,
                                                no.opt="ratematrix"             ,
                                                quiet=TRUE                      )

        #- scale tree and re-estimate Q
        res.neutral.q           = phyloFit(     ali                             ,
                                                init.mod=mod                    ,
                                                scale.only=TRUE                 ,
                                                no.opt=NULL                     ,
                                                quiet=TRUE)

        return(list(    neutral.s = res.neutral.s,
                        neutral.t = res.neutral.t,
                        neutral.q = res.neutral.q))

}


fitMixBgrid <- function(        modB.mix,
                                ali){

        #pts = seq(0,1,len=11)
        #grid = expand.grid(pts,pts)



        afu = function(x,modB.mix,ali){
                mod = modB.mix
                mod$rate.weights = c(x[1],1-x[1])
                mod$ls.model   = NULL
                mod     = add.ls.mod(  mod,
                                                        branch="hg38",
                                                        category=2,
                                                        separate.params="bgc[0,200]", bgc=x[2])
                #- wait for melissa to clamp rate.weights
                #bgc                  = phyloFit(ali, init.mod=modB.mix,
                #                               no.opt=c( "ratematrix","branches"))
                return(likelihood.msa(ali,mod))
        }

        res = apply(grid, 1, afu)
        contour(matrix(res,ncol=11,byrow=T))

}


fitMixGBGC2 <- function(        modB,
                                modAB,
                                ali){

        minB  = list(likelihood=-Inf)
        minAB = list(likelihood=-Inf)

        bmin = modB$ls.model$bgc
        bmax = max(20,bmin*2)
        Bs   = seq(bmin,bmax,len=11)

        #- fit the models with fixed B
        bfu = function(x){
                mod = modB
                mod$nratecats           = 2
                mod$rate.consts         = c(1,1)
                mod$rate.weights        = c(.5,.5)
                mod$ls.model      = NULL
                # these rate.consts should always be 1, this means that the
                mod                = add.ls.mod(  mod,
                                                   branch="hg38",
                                                   category=2,
                                                   separate.params="bgc", bgc=x,const.params="bgc")
                bgc                  = phyloFit(ali, init.mod=mod,
                                                no.opt=c( "ratematrix","branches","backgd"))
                return(bgc)
        }

        tmp   = lapply(Bs,bfu)
        ind   = which.max(unlist(lapply(tmp,function(x) x$likelihood)))
        modB.mix = tmp[[ind]]
        modB.mix$ls.model$defn = "hg38#1:bgc[0,1000]"
        bgc                  = phyloFit(ali, init.mod=modB.mix,
                                                no.opt=c( "ratematrix","branches","backgd"))

        #- backgd breaks the thing -> waiting for melissa
        modAB.mix = bgc
        modAB.mix$ls.model = NULL
        modAB.mix   = add.ls.mod(  modAB.mix, branch="hg38", category=1,
                                 separate.params="sel[0,1000]", sel=0)
        modAB.mix  =  add.ls.mod( modAB.mix, branch="hg38", category=2,
                                 separate.params=c("bgc[0,1000]", "sel#1"), bgc=modB.mix$ls.model$bgc,sel=0)

        both                  = phyloFit(ali, init.mod=modAB.mix,
                                        no.opt=c( "ratematrix","branches","backgd"),quiet=TRUE)

        minB    = bgc
        minAB   = both

        facs = 1+(runif(10)-1/2)*.4
        for(fac in facs){
                mod.b                   = modB
                mod.b$ls.model                 = NULL
                mod.b$nratecats                 = 2
                # these rate.consts should always be 1, this means that the
                #rate matrix is scaled the same in both categories
                mod.b$rate.consts       = c(1,1)
                # rate.weights gives initial frequencies for the two categories
                mod.b$rate.weights      = c(1-fac/max(facs),fac/max(facs)) #- all weight on the gBGC model
                # add gBGC to the human branch in the 2nd category only.  It seems to be necessary
                # to set the gBGC parameter to something > 0 or else it often gets stuck at 0
                mod.b                   = add.ls.mod(mod.b, branch="hg38", category=2, separate.params="bgc[1,100]", bgc=fac*modB$ls.model$bgc)
                bgc                     = phyloFit(ali, init.mod=mod.b, no.opt=c( "ratematrix","branches","backgd"),quiet=TRUE)

                if(bgc$likelihood > minB$likelihood) {
                        minB = bgc
                }

                mod.ab                  = modAB
                mod.ab$nratecats        = 2
                mod.ab$rate.consts       = c(1,1)
                mod.ab$rate.weights      = c(1-fac/max(facs),fac/max(facs))
                mod.ab$ls.model        = NULL
                mod.ab                  = add.ls.mod(  mod.ab, branch="hg38", category=1,
                                                separate.params="sel[0,1000]", sel=fac*modAB$ls.model$selection)
                # The sel#1 notation in the following command is currently undocumented, but implies
                # that the value is shared with selection from category 1
                mod.ab                  =  add.ls.mod( mod.ab, branch="hg38", category=2,
                                                separate.params=c("bgc[1,1000]", "sel#1"), bgc=fac*modAB$ls.model$bgc,sel=fac*modAB$ls.model$selection)
                both                     = phyloFit(ali, init.mod=mod.ab, no.opt=c( "ratematrix","branches","backgd"),quiet=TRUE)

                if(both$likelihood > minAB$likelihood) {
                        minAB = both
                }
        }

        return(list(bgc=minB, both=minAB))
}










fitMixGBGC2 <- function(        modB,
                                modAB,
                                ali){

        minB  = list(likelihood=-Inf)
        minAB = list(likelihood=-Inf)

        bmin = modB$ls.model$bgc
        bmax = max(20,bmin*2)
        Bs   = seq(bmin,bmax,len=11)

        #- fit the models with fixed B
        bfu = function(x){
                mod = modB
                mod$nratecats           = 2
                mod$rate.consts         = c(1,1)
                mod$rate.weights        = c(.5,.5)
                mod$ls.model      = NULL
                # these rate.consts should always be 1, this means that the
                mod                = add.ls.mod(  mod,
                                                   branch="hg38",
                                                   category=2,
                                                   separate.params="bgc", bgc=x,const.params="bgc")
                bgc                  = phyloFit(ali, init.mod=mod,
                                                no.opt=c( "ratematrix","branches","backgd"))
                return(bgc)
        }

        tmp   = lapply(Bs,bfu)
        ind   = which.max(unlist(lapply(tmp,function(x) x$likelihood)))
        modB.mix = tmp[[ind]]
        modB.mix$ls.model$defn = "hg38#1:bgc[0,1000]"
        bgc                  = phyloFit(ali, init.mod=modB.mix,
                                                no.opt=c( "ratematrix","branches","backgd"))

        #- backgd breaks the thing -> waiting for melissa
        modAB.mix = bgc
        modAB.mix$ls.model = NULL
        modAB.mix   = add.ls.mod(  modAB.mix, branch="hg38", category=1,
                                 separate.params="sel[0,1000]", sel=0)
        modAB.mix  =  add.ls.mod( modAB.mix, branch="hg38", category=2,
                                 separate.params=c("bgc[0,1000]", "sel#1"), bgc=modB.mix$ls.model$bgc,sel=0)

        both                  = phyloFit(ali, init.mod=modAB.mix,
                                        no.opt=c( "ratematrix","branches","backgd"),quiet=TRUE)

        minB    = bgc
        minAB   = both

        facs = 1+(runif(10)-1/2)*.4
        for(fac in facs){
                mod.b                   = modB
                mod.b$ls.model                 = NULL
                mod.b$nratecats                 = 2
                # these rate.consts should always be 1, this means that the
                #rate matrix is scaled the same in both categories
                mod.b$rate.consts       = c(1,1)
                # rate.weights gives initial frequencies for the two categories
                mod.b$rate.weights      = c(1-fac/max(facs),fac/max(facs)) #- all weight on the gBGC model
                # add gBGC to the human branch in the 2nd category only.  It seems to be necessary
                # to set the gBGC parameter to something > 0 or else it often gets stuck at 0
                mod.b                   = add.ls.mod(mod.b, branch="hg38", category=2, separate.params="bgc[1,100]", bgc=fac*modB$ls.model$bgc)
                bgc                     = phyloFit(ali, init.mod=mod.b, no.opt=c( "ratematrix","branches","backgd"),quiet=TRUE)

                if(bgc$likelihood > minB$likelihood) {
                        minB = bgc
                }

                mod.ab                  = modAB
                mod.ab$nratecats        = 2
                mod.ab$rate.consts       = c(1,1)
                mod.ab$rate.weights      = c(1-fac/max(facs),fac/max(facs))
                mod.ab$ls.model        = NULL
                mod.ab                  = add.ls.mod(  mod.ab, branch="hg38", category=1,
                                                separate.params="sel[0,1000]", sel=fac*modAB$ls.model$selection)
                # The sel#1 notation in the following command is currently undocumented, but implies
                # that the value is shared with selection from category 1
                mod.ab                  =  add.ls.mod( mod.ab, branch="hg38", category=2,
                                                separate.params=c("bgc[1,1000]", "sel#1"), bgc=fac*modAB$ls.model$bgc,sel=fac*modAB$ls.model$selection)
                both                     = phyloFit(ali, init.mod=mod.ab, no.opt=c( "ratematrix","branches","backgd"),quiet=TRUE)

                if(both$likelihood > minAB$likelihood) {
                        minAB = both
                }
        }

        return(list(bgc=minB, both=minAB))
}






generate.null.model = function(init.mod){
#========================================

        mod             = init.mod
        mod$selection   = 0
return(mod)
}

generate.a.model = function(    init.mod,
                                null.fit,
                                hg.sel = 0){
#=========================================

        mod                     = init.mod
        mod$selection           = null.fit$selection
	mod$rate.matrix         = apply.bgc.sel(mod$rate.matrix,s=mod$selection)
        mod                     = add.ls.mod(mod, branch="hg38", separate.params="sel[,1000]",sel=hg.sel)
        return(mod)
}


generate.b.model = function(    init.mod,                                null.fit,
                                IN.rate.consts=c(1,1),
                                IN.rate.weights=c(.5,.5),
                                IN.bgc=0){
#=========================================================

        # fit a model with bgc on the human branch in some fraction of sites        # This has 3 free parameters: scale, gBGC, and the fraction
        mod                     = init.mod
        mod$selection           = null.fit$selection
        mod$rate.matrix         = apply.bgc.sel(mod$rate.matrix,s=mod$selection)
        mod$ls.model           = NULL
        mod$nratecats           = 2
        # these rate.consts should always be 1, this means that the
        #rate matrix is scaled the same in both categories
        mod$rate.consts         = IN.rate.consts
        # rate.weights gives initial frequencies for the two categories
        mod$rate.weights        = IN.rate.weights
        # add gBGC to the human branch in the 2nd category only.  It seems to be necessary
        # to set the gBGC parameter to something > 0 or else it often gets stuck at 0
        mod                     = add.ls.mod(mod, branch="hg38", category=2, separate.params="bgc[0,1000]", bgc=IN.bgc)
        return(mod)
}


generate.ab.model = function(   init.mod ,
                                a.fit ,
                                b.fit ,
                                IN.rate.consts=c(1,1)){
#======================================================

        if(a.fit$likelihood > b.fit$likelihood){
                ref.model       = a.fit
                IN.rate.weights = c(1,0)
                IN.sel          = a.fit$ls.model$selection
        }
        else{
                ref.model       = b.fit
                IN.rate.weights = b.fit$rate.weights
                IN.sel          = 0
        }

        mod                     = init.mod
        mod$ls.model           = NULL
        mod$nratecats           = 2
        mod$selection           = ref.model$selection
        mod$rate.consts         = IN.rate.consts
        mod$rate.weights        = IN.rate.weights
        mod$backgd              = ref.model$backgd
        mod$rate.matrix         = apply.bgc.sel(mod$rate.matrix,s=mod$selection)

        mod                     = add.ls.mod(  mod, branch="hg38", category=1,
                                                separate.params="sel[0,1000]", sel=IN.sel)
        # The sel#1 notation in the following command is currently undocumented, but implies
        # that the value is shared with selection from category 1
        mod                     = add.ls.mod(  mod, branch="hg38", category=2,
                                                separate.params=c("bgc[0,1000]", "sel#1"), bgc=b.fit$ls.model$bgc,sel=IN.sel)


}
