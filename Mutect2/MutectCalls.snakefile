import os

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input:
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.tbi",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}_f1r2.tar.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.stats",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        #expand("results/{base_file_name}/gathered_unfiltered.vcf.gz",base_file_name=config["base_file_name"]),
        #expand("results/{base_file_name}/mutect_merged.stats", base_file_name = config["base_file_name"]),
        #expand("results/{base_file_name}/read_orientation_model.tar.gz", base_file_name = config["base_file_name"]),
        #expand("results/{base_file_name}/tumor_{base_file_name}_pileup.tab",base_file_name=config["base_file_name"]),
        #expand("results/{base_file_name}/normal_{base_file_name}_pileup.tab", base_file_name = config["base_file_name"]),
        #expand("results/{base_file_name}/{base_file_name}_tum_segments.tab", base_file_name = config["base_file_name"]),
        #expand("results/{base_file_name}/{base_file_name}_matched_contamination.tab", base_file_name = config["base_file_name"]),
        #expand("results/{base_file_name}/{base_file_name}_f1r2_filtered_somatic_vcf.gz", base_file_name = config["base_file_name"])

rule Mutect2:
	input:
		tumor_file = config["samples"],
		normal_file = config["normals"]
	output:
		vcf = expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
		tbi = expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.tbi",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
		tar = expand("results/{base_file_name}/unfiltered_{chromosomes}_f1r2.tar.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
		stats = expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.stats",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
	params:
		reference_genome = config["reference_genome"],
		mutect2_germline_resource = config["germline_resource"],
		gatk = config["gatk_path"],
		panel_of_normals = config["panel_of_normals"],
		chromosomes = config["chromosomes"]
		#normals = lambda wildcards: config["samples"][wildcards.tumors][1]
	log:
		expand("logs/mutect2/{base_file_name}_{chromosomes}_mutect2.txt", base_file_name = config["base_file_name"],chromosomes = config["chromosomes"])
	shell:
		"""
		all_tumor_inputs=`cat {input.tumor_file}`
		all_normal_inputs="-I `cat {input.normal_file}`"
		
		({params.gatk} Mutect2 \
		-reference {params.reference_genome} \
		-I bam.list \
		-normal B_TRCC_18_Normal \
		-intervals {params.chromosomes} \
		--germline-resource {params.mutect2_germline_resource} \
		--f1r2-tar-gz {output.tar} \
		--panel-of-normals {params.panel_of_normals} \
		-output {output.vcf}) 2> {log}"""
     
