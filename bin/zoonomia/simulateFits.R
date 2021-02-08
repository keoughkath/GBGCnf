



#- SIMULATIONS for LRTs
#======================

nsim = 5

simulate_fits <- function(init.mod,har_gc,anc_gc,har.msa,har_fit_summary,mod.neu){
#=================================================================================

    #- the gap pattern for the alignment
    #------------------------------------
    har.map 	= make.map(har.msa)


    #- get a tree model for simulation; alignment without human
    #----------------------------------------------------------
    tmp.mod  	            = init.mod
    tmp.mod  	            = mod.backgd.tm(tmp.mod, gc=har_gc)
    har.msa.nohg 	      	= har.msa
    har.msa.nohg[[1]][1] 	= gsub("[ACGT]","N",har.msa.nohg["Homo_sapiens"]$seq)
    res.scale 	        	= phyloFit(	har.msa.nohg,
    				                        	init.mod=tmp.mod,
    				                        	scale.only=TRUE,
    			                        		no.opt=c("backgd", "ratematrix"))

    mod.sim               = res.scale
    mod.sim               = mod.backgd.tm(mod.sim, gc=anc_gc)


    #- make a grid for S and B
    #-------------------------
    smax = ceiling(max(c(har_fit_summary[,3],har_fit_summary[,4])))
    bmax = ceiling(max(har_fit_summary[,5]))

    #- estimate the boost to neutral rate (approximately)
    s.neu = (-1)*har_fit_summary["sel.r","shat2"]

    #-  minimum for s.max and b.max
    if(smax < ceiling(s.neu)*5){
    	smax = ceiling(s.neu)*5
    }
    if(bmax < round(ceiling(s.neu)*5*1.5)){
    	bmax = round(ceiling(s.neu)*5*1.5)
    }

    ses = sort( c(seq(-s.neu,smax,len=51),0));
    ses = ses - ses[1] #- since the model we're simulating from has no sel parameter
    bes = seq(0,bmax,len=51)


    #- SIMULATE
    #----------

    SIM       = list(S = vector(len=52,mode="list"), B = vector(len=51,mode="list"))

    for(i in 1:52){

      print(i)

    	#- transform generating model
    	#============================
    	cur.s = ses[i]
    	tmp.mod = add.ls.mod(mod.sim,"Homo_sapiens",selection=cur.s)
    	SIM$S[[i]]$mod = tmp.mod
    	SIM$S[[i]]$val = cur.s
    	msas = vector(len=NSIM,mode="list")
    	rslt = vector(len=NSIM,mode="list")
    	for(j in 1:NSIM){
    		msas[[j]] <- simulate.msa(tmp.mod, nsim=dim(har.msa)[2])
    		msas[[j]] <- gappify(msas[[j]],har.map)
        sim.mods  <- models.har(mod.neu)
        sim.fit   <- fit.har(sim.mods,msas[[j]],fitBG=TRUE)
    		rslt[[j]] <- sim.fit$summary
    	}
    	SIM$S[[i]]$dat = msas
    	SIM$S[[i]]$fit = rslt
    }

    for(i in 1:51){

      print(i)

      #- transform generating model
    	#============================
    	cur.b = bes[i]
    	tmp.mod = add.ls.mod(mod.sim,"Homo_sapiens",bgc=cur.b)
    	SIM$B[[i]]$mod = tmp.mod
    	SIM$B[[i]]$val = cur.b
    	msas = vector(len=NSIM,mode="list")
    	rslt = vector(len=NSIM,mode="list")
    	for(j in 1:NSIM){
    		msas[[j]] <- simulate.msa(tmp.mod, nsim=dim(har.msa)[2])
    		msas[[j]] <- gappify(msas[[j]],har.map)
        sim.mods  <- models.har(mod.neu)
        sim.fit   <- fit.har(sim.mods,msas[[j]],fitBG=TRUE)
    		rslt[[j]] <- sim.fit$summary
    	}
    	SIM$B[[i]]$dat = msas
    	SIM$B[[i]]$fit = rslt
    }

    return(SIM)
  }
