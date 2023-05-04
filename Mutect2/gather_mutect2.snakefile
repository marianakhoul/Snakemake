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


rule merge_mutect_stats:
    input:
        chr1_stats = "results/{tumors}/unfiltered_chr1.vcf.gz.stats",
        chr2_stats = "results/{tumors}/unfiltered_chr2.vcf.gz.stats",
        chr3_stats = "results/{tumors}/unfiltered_chr3.vcf.gz.stats",
        chr4_stats = "results/{tumors}/unfiltered_chr4.vcf.gz.stats",
        chr5_stats = "results/{tumors}/unfiltered_chr5.vcf.gz.stats",
        chr6_stats = "results/{tumors}/unfiltered_chr6.vcf.gz.stats",
        chr7_stats = "results/{tumors}/unfiltered_chr7.vcf.gz.stats",
        chr8_stats = "results/{tumors}/unfiltered_chr8.vcf.gz.stats",
        chr9_stats = "results/{tumors}/unfiltered_chr9.vcf.gz.stats",
        chr10_stats = "results/{tumors}/unfiltered_chr10.vcf.gz.stats",
        chr11_stats = "results/{tumors}/unfiltered_chr11.vcf.gz.stats",
        chr12_stats = "results/{tumors}/unfiltered_chr12.vcf.gz.stats",
        chr13_stats = "results/{tumors}/unfiltered_chr13.vcf.gz.stats",
        chr14_stats = "results/{tumors}/unfiltered_chr14.vcf.gz.stats",
        chr15_stats = "results/{tumors}/unfiltered_chr15.vcf.gz.stats",
        chr16_stats = "results/{tumors}/unfiltered_chr16.vcf.gz.stats",
        chr17_stats = "results/{tumors}/unfiltered_chr17.vcf.gz.stats",
        chr18_stats = "results/{tumors}/unfiltered_chr18.vcf.gz.stats",
        chr19_stats = "results/{tumors}/unfiltered_chr19.vcf.gz.stats",
        chr20_stats = "results/{tumors}/unfiltered_chr20.vcf.gz.stats",
        chr21_stats = "results/{tumors}/unfiltered_chr21.vcf.gz.stats",
        chr22_stats = "results/{tumors}/unfiltered_chr22.vcf.gz.stats",
        chrX_stats = "results/{tumors}/unfiltered_chrX.vcf.gz.stats",
        chrY_stats = "results/{tumors}/unfiltered_chrY.vcf.gz.stats"
    output:
        protected("results/{tumors}/mutect_merged.stats")
    params:
        gatk = config["gatk_path"]
    log:
        "logs/merge_mutect_stats/{tumors}_merge_mutect_stats.txt"
    shell:
        "({params.gatk} MergeMutectStats \
        -stats {input.chr1_stats} \
        -stats {input.chr2_stats} \
        -stats {input.chr3_stats} \
        -stats {input.chr4_stats} \
        -stats {input.chr5_stats} \
        -stats {input.chr6_stats} \
        -stats {input.chr7_stats} \
        -stats {input.chr8_stats} \
        -stats {input.chr9_stats} \
        -stats {input.chr10_stats} \
        -stats {input.chr11_stats} \
        -stats {input.chr12_stats} \
        -stats {input.chr13_stats} \
        -stats {input.chr14_stats} \
        -stats {input.chr15_stats} \
        -stats {input.chr16_stats} \
        -stats {input.chr17_stats} \
        -stats {input.chr18_stats} \
        -stats {input.chr19_stats} \
        -stats {input.chr20_stats} \
        -stats {input.chr21_stats} \
        -stats {input.chr22_stats} \
        -stats {input.chrX_stats} \
        -stats {input.chrY_stats} \
        -O {output}) 2> {log}"
