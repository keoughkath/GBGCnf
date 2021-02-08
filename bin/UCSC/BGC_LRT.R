
library(rphast)

#- BRAND NEW CHEF
#================

#- GC-content stuff
source("./getGCfunctions.R")
source("./maskCpG.R")
#- model fitting stuff
source("./hars_models.R")
source("./hars_fit.R")
#- for gap patterns
source("./gappify_msa.R")
#- for simulating fits
source("./simulateFits.R")
#- for classifying the HAR
source("./classify_har.R")


#- INPUTS:

NSIM                = 5 # orig was 55
modfile			        = "./indat/28way_ssrev_4d.mod"
har.msa.file		    = "./indat/HAR_001_ali_umskd.ss"
har.flank.msa.file	= "./indat/har_1_flnk.ss"
#har.map.file		    = "./indat/HAR_001_ali_umskd.map"
ali.orgs.file		    = "/pollard/home/kathleen/projects/enhancer_hijacking/dat/analyses_dat/call_hars_hg38/species_names/original_HAR_species.lst"#"./indat/moi_orgs.txt"
resdir = "."


#- the species we keep
orgs	      	  = scan(ali.orgs.file,what="character")
#- the initial model
init.mod      	= read.tm(modfile)
init.mod$tree   = prune.tree(init.mod$tree,orgs,all.but=TRUE)
#- the HAR sequence
har.msa       	= read.msa(har.msa.file,seqnames=orgs)
har.msk.msa     = maskCpG(har.msa)
#- the flanking sequence
har.flank.msa 	= read.msa(har.flank.msa.file,seqnames=orgs)
#- outfile
outfile         = paste(resdir,"/",basename(har.msa.file),".rdat",sep="")


#- GC-CONTENT(S)
#===============

#- ancestral
anc_gc   = get_anc_gc(har.msa)
#- flanking
flank_gc = get_flank_gc(har.flank.msa)
#- har
har_gc   = get_har_gc(har.msa)

#- MODEL FITTING
#===============

#- calibrate initial model to local GC and neutral rate:
tmp.mod       = init.mod
tmp.mod       = mod.backgd.tm(tmp.mod, gc=flank_gc)
res.neutral		= phyloFit(	har.flank.msa,
				                  init.mod=tmp.mod,
				                  scale.only=TRUE,
				                  no.opt=c("backgd","ratematrix"))

#- generate and fit models.
mods = models.har(res.neutral)
fits = fit.har(mods,har.msa,fitBG=TRUE)
HAR  = list(mod.neu = res.neutral,mod.fit=mods,fit=fits)

#- SIMULATE FITS
#===============

SIM  = simulate_fits(init.mod,har_gc,anc_gc["unmasked"],har.msa,HAR$fit$summary,HAR$mod.neu)

#- ASSIGN CLASS
#==============

class = classify.har(HAR,SIM)
print(class)
lapply(class, function(x) write.table( data.frame(x), "test_class.csv"  , append= T, sep=',' ))
# write.table(class,file="test_class.csv", sep=",")
