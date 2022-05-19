import pandas as pd
import os

configfile: "config.yaml"

OUTPUTDIR = config['run_info']["output"]
CONDITION_IDS = list(config["data"].keys())
BLACKLIST = config['run_info']["blacklist"]
GTF = config['run_info']["gtf"]

samples = (
    pd.read_table(config["samplefile"])
    .set_index("Condition", drop=False)
    .sort_index
)

bam_dir = "results/mapped_clean/"
peak_dir = "results/peak_calling/"



def bams_condition(wildcards):
        files = samples.loc[wildcards.conditon]
        bams = files["Filebase"].to_list()]
        bams= [bam_dir + x + ".clean.bam" for x in bams]
        return bams

# here i kinda need to figuire out how to remove certain chromosomes from the list, as I do some cleaning later
rule get_fasta_chroms:
	input:
		"data/index/hg19.p13.plusMT.no_alt_analysis_set.fa"
	output:
		txt = os.path.join(OUTPUTDIR, "flatfiles", "chromsizes.txt"),
		bed = os.path.join(OUTPUTDIR, "flatfiles", "chromsizes.bed")
	shell:
		"samtools faidx {input} --fai-idx {output.txt};"
		"awk '{{ print $1\"\t\"0\"\t\"$2 }}' {output.txt} > {output.bed}"


rule conditionbam:
	input:
		bams_condition
	output:
		bam= "results/merged_bams/{condition}_merged.bam",
		bai= "results/merged_bams/{condition}_merged_sorted.bam.bai"
	threads:
		10
	message:
		"Joining individual bamfiles from condition {wildcards.condition}"
    conda:
        "envs/samtools.yaml"
	run:
		output_dir= "results/merged_bams/"

		merge_bams= []
		for bam in input:
			prefix= os.path.basename(bam)
			temp_prefix= output_dir + prefix
			merge_bams.append(os.path.join(output_dir, prefix + ".sorted"))
			shell("samtools sort {0} -o {1} -@ {2} -T {3}".format(bam, merge_bams[-1], threads, temp_prefix))

		if len(merge_bams) > 1:
			shell("samtools merge -@ {threads} {output.bam} " + " ".join(merge_bams))
			shell("samtools index {output.bam}")
		else:
			shell("samtools sort -o {output.bam} {input}")
			shell("samtools index {output.bam}")


rule coverage_bigwig:
	input:
		bam= rules.conditionbam.output.bam,
		chromsizes= rules.get_fasta_chroms.output.txt
	output:
		bedgraph= temp(os.path.join(OUTPUTDIR, "coverage", "{condition}_coverage.bg")),
		bigwig= os.path.join(OUTPUTDIR, "coverage", "{condition}_coverage.bw")
	params:
		"-T " + os.path.join(OUTPUTDIR, "coverage")  # temporary dir for sorting
	conda:
		os.path.join(environments_dir, "coverage.yaml")
	shell:
		"bedtools genomecov -ibam {input.bam} -bg | sort -k1,1 -k2,2n {params} > {output.bedgraph};"
		"bedGraphToBigWig {output.bedgraph} {input.chromsizes} {output.bigwig}"


gsizes = {"human": "hs",
			"mouse": "mm",
			"zebrafish": 1369631918}  # https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html

def condition_for_bam(wildcards):
        file = samples.loc[{wildcards.condition]].set_index("Filebase").loc[{wildcards.sample}][filebase]
        file= "results/mapped_clean/" + file + ".bam"
        return file


rule macs:
	input:
        bam= condition_for_bam,
	output:
		macs= os.path.join(OUTPUTDIR, "peak_calling", "{condition}", "{sample_id}_peaks.broadPeak"),
		raw= os.path.join(OUTPUTDIR, "peak_calling", "{condition}", "{sample_id}_raw.bed")
	log:
		os.path.join(OUTPUTDIR, "logs", "{condition}_{sample_id}_peak_calling.log")
	message:
		"Running macs2 with .bam-file: {input}"
	conda:
		os.path.join(environments_dir, "macs.yaml")
	params:
		"--name {sample_id}",
		"--outdir " + os.path.join(OUTPUTDIR, "peak_calling", "{condition}"),
		"--gsize " + str(gsizes[config["run_info"]["organism"]]),
		config.get("macs", "--nomodel --shift -100 --extsize 200 --broad"),
	shell:
		"macs2 callpeak -t {input} {params} &> {log}; "
		"cp {output.macs} {output.raw}; "


def peaks_condition(wildcards):
        files = samples.loc[wildcards.conditon]
        peaks = files["Filebase"].to_list()]
        peaks= [peak_dir + x + ".clean.bam" for x in peaks]
        return peaks


rule process_peaks:
	input:
		peaks = peaks_condition,
		blacklist = BLACKLIST,
		whitelist = rules.get_fasta_chroms.output.bed
	output:
		peaks = os.path.join(OUTPUTDIR, "peak_calling", "{condition}_union.bed")
	message: "Processing peaks from condition {wildcards.condition}"
	conda:
		"envs/tools.yaml"
	shell:
		"cat {input.peaks} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -d 5 | "
		"bedtools subtract -a - -b {input.blacklist} -A | "
		"bedtools intersect -a - -b {input.whitelist} -wa | "
		"awk '$1 !~ /[M]/' | "  # exclude chromosome M
		# add condition name to each peak
		"awk '{{print $s\"\\t{wildcards.condition}\"}}' > {output.peaks}"

# Union peaks across all conditions
rule merge_condition_peaks:
	input:
		[os.path.join(OUTPUTDIR, "peak_calling", condition + "_union.bed") for condition in CONDITION_IDS]
	output:
		temp(os.path.join(OUTPUTDIR, "peak_calling", "all_merged.tmp"))
	message:
		"Merging peaks across conditions"
	conda:
		os.path.join(environments_dir, "tools.yaml")
	shell:
		"cat {input} | sort -k1,1 -k2,2n | bedtools merge -d 5 -c 4 -o distinct > {output}"

# Get correct sorting of peak_names
rule sort_peak_names:
	input:
		rules.merge_condition_peaks.output
	output:
		peaks = os.path.join(OUTPUTDIR, "peak_calling", "all_merged.bed")
	run:
		out = open(output[0], "w")
		with open(input[0]) as f:
			for line in f:
				columns = line.rstrip().split("\t")

				# Sort according to condition names
				peak_ids = columns[3].split(",")
				columns[3] = ",".join(sorted(peak_ids, key= lambda x: CONDITION_IDS.index(x)))

				out.write("\t".join(columns) + "\n")
		out.close()

#peaks so far are called, merged and _merged_sorted

#Config for uropa annotation
rule uropa_config:
	input:
		bed = rules.sort_peak_names.output.peaks, 	#os.path.join(OUTPUTDIR, "peak_calling", "all_merged.bed"),
		gtf = GTF
	output:
		config = os.path.join(OUTPUTDIR, "peak_annotation", "all_merged_annotated.config")
	run:

		import json
		config = {"queries":[
					{"feature":"gene", "feature.anchor":"start", "distance":[10000,1000], "filter_attribute":"gene_biotype", "attribute_values":"protein_coding", "name":"protein_coding_promoter"},
					{"feature":"gene", "distance":1, "filter_attribute":"gene_biotype", "attribute_values":"protein_coding", "internals":0.1, "name":"protein_coding_internal"},
					{"feature":"gene", "feature.anchor":"start", "distance":[10000,1000], "name":"any_promoter"},
					{"feature":"gene", "distance":1, "internals":0.1, "name":"any_internal"},
					{"feature":"gene", "distance":[50000, 50000], "name":"distal_enhancer"},
					],
				"show_attributes":["gene_biotype", "gene_id", "gene_name"],
				"priority":"True"
				}

		config["gtf"] = input.gtf
		config["bed"] = input.bed

		string_config = json.dumps(config, indent=4)

		config_file = open(output[0], "w")
		config_file.write(string_config)
		config_file.close()


# Peak annotation
# peaks per condition or across conditions, dependent on run info output
rule uropa:
	input:
		config = rules.uropa_config.output.config 	#os.path.join(OUTPUTDIR, "peak_annotation", "all_merged_annotated.config")
	output:
		finalhits = os.path.join(OUTPUTDIR, "peak_annotation", "all_merged_annotated_finalhits.txt"),
		finalhits_sub = os.path.join(OUTPUTDIR, "peak_annotation", "all_merged_annotated_finalhits_sub.txt"),
		peaks = os.path.join(OUTPUTDIR, "peak_annotation", "all_merged_annotated.bed"),
		header = os.path.join(OUTPUTDIR, "peak_annotation", "all_merged_annotated_header.txt"),
	threads:
		99
	log:
		os.path.join(OUTPUTDIR, "logs", "uropa.log")
	params:
		prefix = os.path.join(OUTPUTDIR, "peak_annotation", "all_merged_annotated")
	conda:
		os.path.join(environments_dir, "uropa.yaml")
	shell:
		"uropa --input {input.config} --prefix {params.prefix} --threads {threads} &> {log}; "
		"cut -f 1-4,7-13,16-19 {output.finalhits} > {output.finalhits_sub}; " #Get a subset of columns
		"head -n 1 {output.finalhits_sub} > {output.header};"  #header
		"tail -n +2 {output.finalhits_sub} > {output.peaks}"   #bedlines

#create header for "peaks" given via run_infow
rule create_peaks_header:
	input:
		config["run_info"]["peaks"]
	output:
		os.path.join(OUTPUTDIR, "peak_annotation", "peaks_header.txt")
	run:

		#Assumes bed6-format + additional columns
		columns = open(input[0]).readline().split("\t")
		n = len(columns)

		header = ["chr", "start", "stop", "name", "score", "strand"] + ["column_{0}".format(i) for i in range(6+1,n+1)]

		#Write header to file
		s = "\t".join(header)
		f = open(output[0], "w")
		f.write(s)
		f.close()
