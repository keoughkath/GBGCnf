

parsimony.func = function(x){
#============================

        if(x[1] == x[2]) return(x[1])
        if(x[1] == x[3]) return(x[1])
        if(x[2] == x[3]) return(x[2])
        return("N")
}


maskCpG <- function(ali,names=c("hg38","panTro4","rheMac3")){
#============================================================

	inds = sapply(names, function(x) which(ali$names == x))
	seqs = sapply(inds,  function(x) toupper(ali[[1]][x]))	 #- keep order, all UC
	SEQS = sapply(seqs,  function(x) strsplit(x,split=""))
	SEQS = matrix(unlist(SEQS),byrow=TRUE,nrow=3)
	anc  = apply(SEQS,2,parsimony.func)	
	len  = dim(SEQS)[2]	

	#- informative sites
	inf  = SEQS[1,] != SEQS[2,]

	cls  = rep(NA, len)
	if( any(SEQS[,2]=="G")){
		 if(anc[1] == "C" & anc[2] == "G") {
			cls[1] = 2
		} else {
			cls[1] = 3
		}
	} else {
		cls[1] = 1
	}

	for(i in 2:(len-1)){
	
		#- potential CpG
		if(any(SEQS[,i-1] == "C") | any(SEQS[,i+1]=="G")){
			if(anc[i] == "C" & anc[i+1] == "G"){
				cls[i] = 2
				next
			} 
			if(anc[i] == "G" & anc[i-1] == "C"){
                                cls[i] = 2
                                next
                        }
			cls[i] = 3

		} else{
			cls[i] = 1
		}
	}

	if( any(SEQS[,len-1]=="C")){
                 if(anc[len-1] == "C" & anc[len] == "G") {
                        cls[len] = 2
                } else {
                        cls[len] = 3
                }
        } else {
                cls[len] = 1
        }
	
	msk.ind  = (cls != 1) & (inf)

	LSEQS    = sapply(ali[[1]],function(x)strsplit(x,split=""))
	LMSEQS   = lapply(LSEQS,function(x) {x[msk.ind] = "N"; x }) 
	ali[[1]] = sapply(1:length(LSEQS),function(i) paste(LMSEQS[[i]],collapse=""))
	return(ali)
}



countSub <- function(ali,names=c("hg38","panTro4","rheMac3")){
#=============================================================

        inds = sapply(names, function(x) which(ali$names == x))
        seqs = sapply(inds,  function(x) toupper(ali[[1]][x]))   #- keep order, all UC
        SEQS = sapply(seqs,  function(x) strsplit(x,split=""))
        SEQS = matrix(unlist(SEQS),byrow=TRUE,nrow=3)
        anc  = apply(SEQS,2,parsimony.func)
	table(unlist(apply(cbind(anc,SEQS[1,]),1,function(x) if(x[1] != x[2]) paste(x[1],"->",x[2],sep="") )))
}

getAnc <- function(ali,names=c("hg38","panTro4","rheMac3")){
#===========================================================

        inds = sapply(names, function(x) which(ali$names == x))
        seqs = sapply(inds,  function(x) toupper(ali[[1]][x]))   #- keep order, all UC
        SEQS = sapply(seqs,  function(x) strsplit(x,split=""))
        SEQS = matrix(unlist(SEQS),byrow=TRUE,nrow=3)
        anc  = apply(SEQS,2,parsimony.func)
	anc
}


get.subs <- function(ali){
#=========================
        anc     = getAnc(ali)
        hg      = toupper(ali[[1]][which(ali$names=="hg38")])
        hg      = strsplit(hg, split="")[[1]]
        pt      = toupper(ali[[1]][which(ali$names=="panTro4")])
        pt      = strsplit(pt, split="")[[1]]
        ind1    = anc %in% c("A","C","G","T")
        ind2    = hg %in% c("A","C","G","T")

        #anc     = anc[ind1&ind2]
        #hg      = hg[ind1&ind2]

        sbs.t     = anc!=hg
        sbs      = as.numeric(sbs.t)
        sbs[sbs==1] = "?"
       	sbs[( ind1 & (! ind2) ) | ( ind2 & (! ind1))] = "indel"
	ind      = ( (anc %in% c("A","T")) & (hg %in% c("C","G")) ) & sbs.t
        sbs[ind] = "ws"
        ind      = ( (anc %in% c("C","G")) & (hg %in% c("A","T")) ) & sbs.t
        sbs[ind] = "sw"
        ind      = ( (anc %in% c("C","G")) & (hg %in% c("C","G")) ) & sbs.t
        sbs[ind] = "ss"
        ind      = ( (anc %in% c("A","T")) & (hg %in% c("A","T")) ) & sbs.t
        sbs[ind] = "ww"

        return(sbs)
}





