
#==============================
#- get hg-pt-ancestral GC-freqs
#==============================

get_anc_gc <- function(har.msa){
#===============================

  seqs = rbind(	strsplit(har.msa["hg38"]$seq   ,split="")[[1]],
	strsplit(har.msa["panTro4"]$seq,split="")[[1]],
	strsplit(har.msa["rheMac3"]$seq,split="")[[1]])

  anc.seq      = apply(seqs,2,parsimony.func)
  inds         = anc.seq %in% c("A","C","G","T")
  anc.gcfreq   = sum((table(anc.seq[inds]) / sum(table(anc.seq[inds])))[c("C","G")])
  anc.gcfreq   = round(anc.gcfreq,3)

  seqs = rbind(   strsplit(har.msk.msa["hg38"]$seq   ,split="")[[1]],
                  strsplit(har.msk.msa["panTro4"]$seq,split="")[[1]],
                  strsplit(har.msk.msa["rheMac3"]$seq,split="")[[1]])

  anc.msk.seq  = apply(seqs,2,parsimony.func)
  inds         = anc.msk.seq %in% c("A","C","G","T")
  anc.m.gcfreq = sum((table(anc.msk.seq[inds]) / sum(table(anc.msk.seq[inds])))[c("C","G")])
  anc.m.gcfreq = round(anc.m.gcfreq,3)

  res = c(anc.gcfreq,anc.m.gcfreq)
  names(res) = c("unmasked","masked")
  return(res)

}

#====================================
#- get the flanking sequence GC-freqs
#====================================

get_flank_gc <- function(har.flank.msa){
#=======================================

  freqs 		= base.freq.msa(har.flank.msa)
  bind 		= which(as.character(freqs[,"states"]) %in% c("A","C","G","T"))
  freqs 		= freqs[bind,]
  freqs[,"freq"] 	= freqs[,"freq"] / sum(freqs[,"freq"])
  cind 		= which(as.character(freqs[,"states"]) == "C")
  gind 		= which(as.character(freqs[,"states"]) == "G")
  gcfreq.flnk	= sum(freqs[c(cind,gind),"freq"])
  gcfreq.flnk	= round(gcfreq.flnk,3)
  rm(cind,gind,freqs,bind)

  return(gcfreq.flnk)
}

#==========================
#- get the HAR-ali GC-freqs
#==========================

get_har_gc <- function(har.msa){
#===============================

  freqs 		= base.freq.msa(har.msa)
  bind 		= which(as.character(freqs[,"states"]) %in% c("A","C","G","T"))
  freqs 		= freqs[bind,]
  freqs[,"freq"] 	= freqs[,"freq"] / sum(freqs[,"freq"])
  cind 		= which(as.character(freqs[,"states"]) == "C")
  gind 		= which(as.character(freqs[,"states"]) == "G")
  gcfreq		= sum(freqs[c(cind,gind),"freq"])
  gcfreq		= round(gcfreq,3)
  rm(cind,gind,freqs,bind)

  return(gcfreq)
}
