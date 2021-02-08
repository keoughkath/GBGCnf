
models.har = function(	init.mod,
                        a.shg.range    	= "sel[0,100]"               ,
                        b.bhg.range    	= "bgc[0,100]"               ,
                        ab.sbhg.range	= "bgc[0,100],sel[0,100]"    ,
                        ar.sg.range  	= "sel[-10,0]"               ,
                        ar.shg.range  	= "sel[-10,0]"               ,
                        ap.sg.range	= "sel[-10,0]"               ,
                        ap.shg.range  	= "sel[0,100]"               ,
                        initvals   	= NULL			     ){
#======================================================================

        #- null model
        #------------

        mod.null                = init.mod
        mod.null$selection      = 0
        if(! is.null(initvals$n.sg)){
                mod.null$selection = initvals$n.sg
                mod.null$rate.matrix = apply.bgc.sel(mod.null$rate.matrix,s=mod.null$selection)
        }

        #- gBGC model
        #------------

        mod.b            = init.mod
        mod.b$selection  = 0
        #mod.b            = add.ls.mod(mod.b,"Homo_sapiens",sep=b.range,selection=0,const.params="sel")
        mod.b            = add.ls.mod(mod.b,"Homo_sapiens",sep=b.bhg.range, bgc=0)
        if(! is.null(initvals$b.sg)){
                mod.b$selection   = initvals$b.sg
                mod.b$rate.matrix = apply.bgc.sel(mod.b$rate.matrix,s=mod.b$selection)
        }
        if(! is.null(initvals$b.bhg)){
                mod.b$alt.model$bgc          = initvals$b.bhg
                mod.b$alt.model$rate.matrix  = apply.bgc.sel(init.mod$rate.matrix,b=mod.b$alt.model$bgc)
        }

	#- selection model
        #-----------------

        mod.s            = init.mod
        mod.s$selection  = 0
        mod.s            = add.ls.mod(mod.s,"Homo_sapiens",sep=a.shg.range)
        if(! is.null(initvals$a.sg)){
                mod.s$selection   = initvals$a.sg
                mod.s$rate.matrix = apply.bgc.sel(mod.s$rate.matrix,s=mod.s$selection)
        }
        if(! is.null(initvals$a.shg)){
                mod.s$alt.model$selection    = initvals$a.shg
                mod.s$alt.model$rate.matrix  = apply.bgc.sel(init.mod$rate.matrix,s=mod.s$selection + mod.s$alt.model$selection)
        }

	#- relaxation and positive selection
        #-----------------------------------

        # apply a label to all the branches besides Homo_sapiens
        mod   = init.mod
        mod$tree        <- name.ancestors(mod$tree)
        backgd.branches <- summary.tree(mod$tree)$name
        backgd.branches <- backgd.branches[backgd.branches != "Homo_sapiens" & !is.na(backgd.branches)]
        backgd.branches <- as.character(backgd.branches)
        mod$tree        <- label.branches(mod$tree, backgd.branches, "backgd")
        mod$tree        <- label.branches(mod$tree, "Homo_sapiens", "Homo_sapiens")

        # relaxation of constraint model
        mod.s.r <- mod
        mod.s.r <- add.ls.mod(mod.s.r, label="backgd", separate.params=ar.sg.range)
        mod.s.r <- add.ls.mod(mod.s.r, label="Homo_sapiens"  , separate.params=ar.shg.range)
        # positive selection model
        mod.s.p <- mod
        mod.s.p <- add.ls.mod(mod.s.p, label="backgd", separate.params=ap.sg.range)
        mod.s.p <- add.ls.mod(mod.s.p, label="Homo_sapiens"  , separate.params=ap.shg.range)
        rm(mod)

        if(!is.null(initvals$ar.sg)){
                mod.s.r$alt.model[[1]]$selection = initvals$ar.sg
                mod.s.r$alt.model[[1]]$rate.matrix = apply.bgc.sel(     mod.s.r$alt.model[[1]]$rate.matrix,
                                                                        s=mod.s.r$alt.model[[1]]$selection)
        }
        if(!is.null(initvals$ar.shg)){
                mod.s.r$alt.model[[2]]$selection = initvals$ar.shg
                mod.s.r$alt.model[[2]]$rate.matrix = apply.bgc.sel(     mod.s.r$alt.model[[2]]$rate.matrix,
                                                                        s=mod.s.r$alt.model[[2]]$selection)
        }
        if(!is.null(initvals$ap.sg)){
                mod.s.p$alt.model[[1]]$selection = initvals$ap.sg
                mod.s.p$alt.model[[1]]$rate.matrix = apply.bgc.sel(     mod.s.p$alt.model[[1]]$rate.matrix,
                                                                        s=mod.s.p$alt.model[[1]]$selection)
        }
        if(!is.null(initvals$ap.shg)){
                mod.s.p$alt.model[[2]]$selection = initvals$ap.shg
                mod.s.p$alt.model[[2]]$rate.matrix = apply.bgc.sel(     mod.s.p$alt.model[[2]]$rate.matrix,
                                                                        s=mod.s.p$alt.model[[2]]$selection)
        }

	#- selection and gBGC
        #--------------------

        mod.sb            = init.mod
        mod.sb$selection  = 0
        mod.sb            = add.ls.mod(mod.sb,"Homo_sapiens",sep=ab.sbhg.range)

        if(!is.null(initvals$ab.sg)){
                mod.sb$selection   = initvals$ab.sg
                mod.sb$rate.matrix = apply.bgc.sel( mod.sb$rate.matrix, s = mod.sb$selection)

                if((!is.null(initvals$ab.shg)) | (!is.null(initvals$ab.bgc))){
                        c.sel = initvals$ab.shg
                        c.bgc = initvals$ab.bgc
                        if(is.null(c.bgc)) c.bgc = 0
                        if(is.null(c.sel)) c.sel = 0
                        mod.sb$alt.model$selection   = c.sel
                        mod.sb$alt.model$bgc         = c.bgc
                        mod.sb$alt.model$rate.matrix =  apply.bgc.sel(  mod.sb$alt.model$rate.matrix,
                                                                s= c.sel + mod.sb$selection, b= c.bgc)
                        rm(c.bgc)
                        rm(c.sel)
                } else {

                        mod.sb$alt.model$selection   = 0
                        mod.sb$alt.model$bgc         = 0
                        mod.sb$alt.model$rate.matrix = mod.sb$rate.matrix
                }
        }

	return(list(	null  = mod.null,
			sel   = mod.s,
			sel.r = mod.s.r,
			sel.p = mod.s.p,
			bgc   = mod.b,
			both  = mod.sb ))

}
