VERSION = "0.0.0"

// Only using for autosomal HARs at the moment.

// parse inputs

initHarFile = file(params.har_path)
dogSyntenyFile = file(params.dog_synteny_file)
mouseSyntenyFile = file(params.mouse_synteny_file)
rhesusSyntenyFile = file(params.rhesus_synteny_file)
arFilterFile = file(params.ar_filters_path)

/*
* split HARs by chromosome
*/

process splitHars {
	tag "Splitting up HARs by region"

	// publishDir params.outdir, mode: "copy", overwrite: false

	errorStrategy 'retry'
	maxRetries 3

	input: 
	file(har_bed) from initHarFile

	output:
	file('har_beds/*.bed') into harsByChunk

	script: 
	"""
	mkdir -p har_beds
	python ${baseDir}/bin/splitHars.py ${har_bed} ${params.maf_path} har_beds
	"""
	
}

harsByChunk.into { harsByChunk1; harsByChunk2 }

/*
* extract MSAs per HAR
*/

process extractHarMsa {
	tag "Extracting MSA for HAR ${chunk}"

	// publishDir params.outdir, mode: "copy", overwrite: false

	errorStrategy 'retry'
	maxRetries 3

	input: 
	file(har_bed) from harsByChunk1.flatten()

	output:
	file('har_msas/*.maf') into harMafs

	script:
	chunk = har_bed.baseName
	"""
	mafsInRegion ${har_bed} -outDir har_msas ${params.maf_path}/${chunk}.maf.gz
	"""

}


/*
* filter HAR flanking regions
*/

process getHarFlankingRegions {
	tag "Filter HAR ${chunk} flanking regions"

	// publishDir params.outdir, mode: "copy", overwrite: false

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(har_bed) from harsByChunk2.flatten()
	file(ar_filters) from arFilterFile

	output: 
	file('har_flanking_beds/*.bed') into harFlankBeds

	script:
	chunk = har_bed.baseName
	"""
	mkdir -p har_flanking_beds
	bedtools slop -i ${har_bed} -b ${params.flank_size} -g ${params.genome_file} | bedtools intersect -a - -b ${dogSyntenyFile} | bedtools intersect -a - -b ${mouseSyntenyFile} | bedtools intersect -a - -b ${rhesusSyntenyFile} | bedtools subtract -a - -b ${ar_filters} | bedtools sort -i - > har_flanking_beds/${chunk}.bed
	"""
}

/*
* extract MSAs for regions flanking HARs
*/

process extractHarFlanks {
	tag "Extracting MSA for HAR ${chunk} flanking region"

	// publishDir params.outdir, mode: "copy", overwrite: false

	errorStrategy 'retry'
	maxRetries 3

	input: 
	file(har_flank_bed) from harFlankBeds.flatten()

	output:
	file('har_flanking_mafs/*.maf') into harFlankMafs

	script:
	chunk = har_flank_bed.baseName
	"""
	mkdir -p har_flanking_mafs_prelim
	mkdir -p har_flanking_mafs

	vals=\$(cut -f 4 ${har_flank_bed} | sort | uniq)

	for val in \${vals[@]}; do grep \${val}'\$' ${har_flank_bed} > \${val}.bed; done
	for val in \${vals[@]}; do mafsInRegion \${val}.bed har_flanking_mafs_prelim/\${val}.prelim.maf ${params.maf_path}/${chunk}.maf.gz; grep -v 'null' har_flanking_mafs_prelim/\${val}.prelim.maf > har_flanking_mafs/\${val}.maf; done
	"""

}

harMafs.into { harMafsOne; harMafsTwo }
harFlankMafs.into { harFlankMafsOne; harFlankMafsTwo }

/*
* convert HAR MAFs to SS files
*/

process convertHARmafToSS {
	tag "Convert HAR MAFs to SS files"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(har_maf) from harMafsOne.flatten()

	output: 
	file('har_ss/*.ss') into harSSfiles

	script:
	har_name = har_maf.simpleName
	"""
	mkdir -p har_ss
	${params.phast_path}/msa_view ${har_maf} --in-format MAF --out-format SS > har_ss/${har_name}.ss
	"""
}

/*
* convert HAR-flank MAFs to SS files
*/

process convertHARflankMafToSS {
	tag "Convert HAR-flank MAFs to SS files"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(har_flank_maf) from harFlankMafsOne.flatten()

	output: 
	file('har_flank_ss/*.ss') into harFlankSSfiles

	script:
	har_name = har_flank_maf.simpleName
	"""
	mkdir -p har_flank_ss
	${params.phast_path}/msa_view ${har_flank_maf} --in-format MAF --out-format SS > har_flank_ss/${har_name}.flank.ss
	"""
}

// match HARs to HAR flank regions

harWithName = harSSfiles.flatten().map { 
	file -> tuple(file.simpleName, file) 
}

harFlankWithName = harFlankSSfiles.flatten().map { 
	file -> tuple(file.simpleName, file) 
}

harWithFlank = harWithName.combine(harFlankWithName, by: 0)


/*
* match species from HAR with flank
*/

process getMatchingSpecies {
	tag "Getting intersection of species HAR vs. flank"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	set val(har), file(har_ss), file(har_flank_ss) from harWithFlank

	output:
	file('species_matches/*.txt') into speciesMatches

	script:
	"""
	mkdir -p species_matches
	python ${baseDir}/bin/match_species.py ${har_ss} ${har_flank_ss} ${har}
	"""
}

harMafWithName = harMafsTwo.flatten().map { 
	file -> tuple(file.simpleName, file) 
}

harMafFlankWithName = harFlankMafsTwo.flatten().map { 
	file -> tuple(file.simpleName, file) 
}

harSpeciesWithName = speciesMatches.flatten().map { 
	file -> tuple(file.simpleName, file) 
}

harSpeciesWithName.into { harSpeciesWithNameOne ; harSpeciesWithNameTwo }

harMafWithSpecies = harMafWithName.combine(harSpeciesWithNameOne, by: 0)

harMafFlankWithSpecies = harMafFlankWithName.combine(harSpeciesWithNameTwo, by: 0)

process filterMatchingSpeciesHars {
	tag "Filtering HAR MAFs for matched species"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	set val(har), file(har_maf), file(har_species) from harMafWithSpecies

	output:
	file('har*_species_matched.maf') into harMafSpeciesMatched

	script:
	"""
	mafSpeciesSubset ${har_maf} ${har_species} har${har}_species_matched.maf
	"""
}

process filterMatchingSpeciesFlanks {
	tag "Filtering HAR-flanking MAFs for matched species"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	set val(har), file(har_flank_maf), file(har_species) from harMafFlankWithSpecies

	output:
	file('har*_flank_species_matched.maf') into harFlankMafSpeciesMatched

	script:
	"""
	mafSpeciesSubset ${har_flank_maf} ${har_species} har${har}_flank_species_matched.maf
	"""
}

process convertHarMSmafSS {
	tag "Convert HAR species-matched MAFs to SS files"

	publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(har_maf) from harMafSpeciesMatched

	output: 
	set val(har_name), file('har_species_matched_ss/*.ss') into harSSSpeciesMatched

	script:
	har_name = har_maf.simpleName.split('_')[0].split('har')[1]
	"""
	mkdir -p har_species_matched_ss
	${params.phast_path}/msa_view  -i MAF -o SS --unordered-ss ${har_maf} > har_species_matched_ss/${har_name}.ss
	"""
}

process convertHarFlankMSmafSS {
	tag "Convert HAR flank species-matched MAFs to SS files"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(har_flank_maf) from harFlankMafSpeciesMatched

	output: 
	set val(har_name), file('har_flank_species_matched_ss/*.ss') into harFlankSSSpeciesMatched

	script:
	har_name = har_flank_maf.simpleName.split('_')[0].split('har')[1]
	"""
	mkdir -p har_flank_species_matched_ss
	${params.phast_path}/msa_view  -i MAF -o SS --unordered-ss ${har_flank_maf} > har_flank_species_matched_ss/${har_name}_flank.ss
	"""
}

// match HARs to HAR flank regions

harWithFlankSpeciesMatched = harSSSpeciesMatched.combine(harFlankSSSpeciesMatched, by: 0)

/*
* perform GBGC-oriented selection analysis
*/

process selectionAnalysisAutosomes {
	tag "Performing selection analysis"

	publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'ignore'
	maxRetries 3

	input:
	set val(har), file(har_ss), file(har_flank_ss) from harWithFlankSpeciesMatched

	output:
	file('selection_results/*.csv')

	script:
	"""
	mkdir -p selection_results 
	cp ${baseDir}/bin/${params.dataset}/* .
	touch /wynton/home/pollard/kathleen/miniconda3/envs/gbgc/lib/R/etc/ldpaths
	Rscript BGC_LRT_new.R ${params.neutral_model_autosomes} ${har_ss} ${har_flank_ss} ${params.species_list} selection_results/HAR${har}.csv
	"""
}


