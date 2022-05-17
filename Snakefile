import panda as pd

samples = (
    pd.read_table(config["samplefile"])
    .set_index("Condition", drop=False)
    .sort_index
)


def bams_condition(wildcards):
        files = samples.loc[wildcards.conditon]
        bams = files["Filebase"].to_list()]
        bams = [bam_dir + x + ".clean.bam" for x in bams]
        return bams

rule merge_bams:
	input:
		bams_condition
	output:
		temp("results/merged_bams/{condition}_merged.bam")
	threads:
		8
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools merge {output} {input} - -threads 7

rule conditionbam:
	input:
		bams_condition
	output:
		bam = "results/merged_bams/{condition}_merged.bam",
		bai = "results/merged_bams/{condition}_merged_sorted.bam.bai"
	threads:
		10
	message:
		"Joining individual bamfiles from condition {wildcards.condition}"
	run:
		output_dir = "results/merged_bams/"

		merge_bams = []
		for bam in input:
			prefix = os.path.basename(bam)
			temp_prefix = output_dir + prefix
			merge_bams.append(os.path.join(output_dir, prefix + ".sorted"))
			shell("samtools sort {0} -o {1} -@ {2} -T {3}".format(bam, merge_bams[-1], threads, temp_prefix))

		if len(merge_bams) > 1:
			shell("samtools merge -@ {threads} {output.bam} " + " ".join(merge_bams))
			shell("samtools index {output.bam}")
		else:
			shell("samtools sort -o {output.bam} {input}")
			shell("samtools index {output.bam}")



rule sort_merged:
	input:
		"results/merged_bams/{condition}_merged.bam"
	output:
		"results/merged_bams/{condition}_merged_sorted.bam"
	threads:
		4
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools sort {input} -@ {threads} -o {output}"


rule index_merged:
    input:
		"results/merged_bams/{condition}_merged_sorted.bam"
    output:
		"results/merged_bams/{condition}_merged_sorted.bam.bai"
    conda:
		"envs/samtools.yaml"
    threads:
    		4
    shell:
    		"samtools index {input} -@ {threads}"


rule get_bigwigs:
	input:
		"results_merged_bams/{condition}_merged_sorted.bam"
	output:
		"results/bigwigs/{condition}.bw
	conda:
		"envs/deeptools.yaml"
	shell:
		"bamCoverage -b {input} -o {output} -of bigwig --effectiveGenomeSize 2864785220 -e –ignoreForNormalization chrX chrM"


rule call_peaks_merged:
	input:
		"results/called_peaks_merged/{condition}.narrowPeak
	output:


# rule matrix_compute:
#	input:
#		bw="results/bigwigs/{condition}.bw"
#		regions=""
