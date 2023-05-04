import os

configfile: "config/samples.yaml"
configfile: "config/config.yaml" 

rule all:
    input:
        #expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        #expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.tbi",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        #expand("results/{base_file_name}/unfiltered_{chromosomes}_f1r2.tar.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        #expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.stats",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        #expand("results/{base_file_name}/gathered_unfiltered.vcf.gz",base_file_name=config["base_file_name"]),
        expand("results/{base_file_name}/mutect_merged.stats", base_file_name = config["base_file_name"]),
        #expand("results/{base_file_name}/read_orientation_model.tar.gz", base_file_name = config["base_file_name"])

        
rule MergeMutectStats:
     output:
        protected("results/{tumors}/mutect_merged.stats")
     params:
        gatk = config["gatk_path"]
     log:
        "logs/MergeMutectStats/{tumors}_merge_mutect_stats.txt"
     shell:
        "
	all_stat_inputs=`for chromosome in {chromosomes}; do
        printf -- "-stats results/{tumors}/unfiltered_${chromosome}.vcf.gz.stats "; done`

	({params.gatk} MergeMutectStats \
        $all_stat_inputs \
        -O {output}) 2> {log}"
