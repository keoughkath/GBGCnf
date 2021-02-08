

readMap <- function(filename){
#=============================

        map        = vector(2,mode="list")
        names(map) = c("mask","mlt")

        garbl = scan(filename, what=character(1),quiet=TRUE)
        specs = garbl[3]
        specs = strsplit(specs,split=",")[[1]]
        for( j in 4:length(garbl)){
                if(j%%2 == 0) {
                        tmp = strsplit(garbl[j],split="")[[1]]
                        tmp = as.logical(as.numeric(tmp))
                        map$mask = cbind(map$mask,tmp)
                }else{
                        map$mlt = c(map$mlt,as.numeric(garbl[j]))
                } 
        }
        rownames(map$mask) = specs
        map
}

make.map <- function(ali){
#=========================	

	SEQS = sapply(ali[[1]],function(x) strsplit(x,split=""))
	SEQS = matrix(unlist(SEQS),nrow=dim(ali)[1],byrow=TRUE)
	inds = matrix(SEQS %in% c("A","a","C","c","G","g","T","t"),nrow=dim(ali)[1])
	cols = apply(inds,2,function(x) paste(as.integer(x),collapse="",sep=""))
	
	map     = list(mask=NULL,mlt=NULL)
	map$mlt = as.vector(table(cols))
	names(map$mlt)=NULL
	cols = sapply(names(table(cols)),function(x) strsplit(x,split=""))
	map$msk = matrix(unlist(lapply(cols,function(x) x != "1")),nrow=dim(ali)[1])
	rownames(map$msk) = names(ali)
	return(map)

}

gappify <- function(ali,map){
#============================

        map$msk        	= map$msk[ali$names,,drop=FALSE]
        seqs            = ali[[1]]
        tmp             = sapply(seqs,function(x) strsplit(x,split=""))
        seqs            = matrix(unlist(tmp), byrow=T, ncol=dim.msa(ali)[2],nrow=dim.msa(ali)[1])
        prob            = map$mlt/sum(map$mlt)
        tmp             = rmultinom(1,sum(map$mlt),prob)
        types           = which(tmp > 0)
        mults           = tmp[tmp>0]
        off             = 1

        for( i in 1:length(types)){
                seqs[map$msk[,types[i]],off:(off+mults[i]-1)] = "N"
                off = off+mults[i]
        }

        seqs     = apply(seqs,1,function(x) paste(x,collapse=""))
        ali[[1]] = seqs
        return(ali)
}





