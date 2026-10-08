import pandas as pd
from pathlib import Path
import numpy as np

include: "rules/cw_clonal.smk"
include: "rules/cw_mixture.smk"
include: "rules/st1_samples.smk"

NAND_REF = "/storage/home/hcoda1/1/jshih41/r-alind6-0/data/all_ncbi_blast_genomes/Nand_genome.fna"
GENOME_DIR = "/storage/home/hcoda1/1/jshih41/r-alind6-0/data/all_ncbi_blast_genomes"
GENOMES = ["GCA_054643315.1_ASM5464331v1_genomic" ,"GCA_027406485.1_ASM2740648v1_genomic" ,"GCA_027406365.1_ASM2740636v1_genomic" ,"GCA_935753195.1_ERR7745403_bin.73_genomic" ,"GCA_027406405.1_ASM2740640v1_genomic" ,"GCA_935752615.1_ERR7738269_bin.84_genomic" ,"GCA_935749445.1_ERR7738585_bin.252_genomic","GCA_935749065.1_ERR7738185_bin.87_genomic" ,"GCA_935751615.1_ERR7745412_bin.2_genomic", "CW_nextpolished_rm_masked_nobact_hasgenes", "Nand_genome"]

ST1_SAMPLES = ["SRR36118145", "SRR20078973", "SRR20078977", "ERR7745403", "SRR20078971", "ERR7738269", "ERR7738585", "ERR7738185", "ERR7745412", "ERR1137517", "ERR14044664", "ERR7745291", "ERR14040741", "ERR7738582", "ERR1305896", "ERR14025145", "ERR7745755", "ERR13966020", "ERR4335103", "ERR14026386", "ERR14020803", "ERR14035067", "ERR7738560", "ERR14039449", "ERR4341940", "ERR14022091"]

rule target:
    input:
        "final_results/cw_clonal/filt_CW.bam",
        "final_results/cw_mixture/filt_cw_mixture.bam",
        expand("final_results/st1_samples/filt_{st1_sample}.bam", st1_sample = ST1_SAMPLES),
        "final_results/cw_clonal/cw_clonal_output.vcf.gz",
        "final_results/cw_mixture/cw_mixture_output.vcf.gz",
        "final_results/st1_samples/updated_filt_merged_vcf.gz",
        "final_results/fastani_all.tsv"


rule index_ref_cw:
    input:
        ref = NAND_REF
    output:
        idx_ref = NAND_REF + ".mmi"
    resources:
        cpus_per_task = 4,
        mem_mb = 8000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        "minimap2 -d {output.idx_ref} {input.ref}"


rule qc:
    input:
        filt_aln_cw = "final_results/cw_clonal/filt_CW.bam",
        filt_aln_cw_mixture = "final_results/cw_mixture/filt_cw_mixture.bam"
    output:
        multiqc = "multiqc.html"
    resources:
        cpus_per_task = 24,
        mem_mb = 64000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        samtools stats {output.filt_aln_cw} > {output.aln_cw_stats}
        samtools flagstat {output.filt_aln_cw} > {output.aln_cw_flag}
        samtools coverage {output.filt_aln_cw} > {output.aln_cw_cov}

        samtools stats {output.filt_aln_cw_mixture} > {output.aln_cw_mixture_stats}
        samtools flagstat {output.filt_aln_cw_mixture} > {output.aln_cw_mixture_flag}
        samtools coverage {output.filt_aln_cw_mixture} > {output.aln_cw_mixture_cov}
        
        """


rule create_fastani_inputs:
    output:
        ql="final_results/query_list.txt",
        rl="final_results/ref_list.txt"
    resources:
        cpus_per_task=4,
        mem_mb=8000,
        runtime="1h",
        partition="short",
        slurm_account="gts-alind6"
    run:

        full_genomes = [GENOME_DIR + "/" +  item + ".fna" for item in GENOMES]
        genomes = '\n'.join(full_genomes)
        
        os.makedirs("final_results", exist_ok=True)
        with open(output.ql, 'w') as f:
            f.write(genomes)
        with open(output.rl, 'w') as f:
            f.write(genomes)

rule fastani_many_to_many:
    input:
        ql="final_results/query_list.txt",
        rl="final_results/ref_list.txt"
    output:
        "final_results/fastani_all.tsv"
    resources:
        cpus_per_task=24,
        mem_mb=64000,
        runtime="4h",
        partition="short",
        slurm_account="gts-alind6"
    shell:
        """
        fastANI --ql {input.ql} --rl {input.rl} \
                -o {output} -t 24
        """
