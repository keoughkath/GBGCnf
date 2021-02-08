library(rphast)

#- BRAND NEW CHEF
#================


#- INPUTS:

# sample command: Rscript /pollard/home/kathleen/projects/enhancer_hijacking/src/selection_analysis/bgc_lrt/BGC_LRT_new.R /pollard/data/functional_genomics/100way_multiz/maf/hqspecies_hars/neutral_models_hqspecies/autosomes.mod /pollard/home/kathleen/projects/enhancer_hijacking/dat/analyses_dat/call_hars_hg38/selection_analysis/hqspecies/MSAs/har10898.ss /pollard/home/kathleen/projects/enhancer_hijacking/dat/analyses_dat/call_hars_hg38/selection_analysis/hqspecies/MSAs_flanking/har10898_flank.ss /pollard/home/kathleen/projects/enhancer_hijacking/dat/analyses_dat/call_hars_hg38/species_names/hq_species.lst /pollard/home/kathleen/projects/enhancer_hijacking/dat/analyses_dat/call_hars_hg38/selection_analysis/hqspecies/pos_sel_out/har10898.csv

args = commandArgs(trailingOnly=TRUE)

NSIM                = 55 #int(args[1])# orig was 55
modfile			        = args[1]#"./indat/28way_ssrev_4d.mod"
har.msa.file		    = args[2]#"./indat/HAR_001_ali_umskd.ss"
har.flank.msa.file	= args[3]#"./indat/har_1_flnk.ss"
ali.orgs.file		    = args[4]#"./indat/moi_orgs.txt"
outfile         = args[5] # must end in .csv
setwd(args[6])


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

print(args)

# print(NSIM)


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
lapply(class, function(x) write.table( data.frame(x), outfile  , append= T, sep=',' ))