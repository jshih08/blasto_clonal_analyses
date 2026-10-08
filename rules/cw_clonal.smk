NAND_REF = "/storage/home/hcoda1/1/jshih41/r-alind6-0/data/all_ncbi_blast_genomes/Nand_genome.fna"

rule aln_cw_clonal:
    input:
        ref=NAND_REF,
        ref_idx = NAND_REF + ".mmi",
        fwd="/storage/home/hcoda1/1/jshih41/r-alind6-0/data/cw/cw_illumina/CW_S10_R1_001.fastq.gz",
        rev="/storage/home/hcoda1/1/jshih41/r-alind6-0/data/cw/cw_illumina/CW_S10_R2_001.fastq.gz",
    output:
        aln_cw = "final_results/cw_clonal/CW.bam",
    resources:
        cpus_per_task = 16,
        mem_mb = 32000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        minimap2 -ax sr -t {resources.cpus_per_task} {input.ref_idx} {input.fwd} {input.rev} | \
        samtools view -b -q 30 -F 0x4 -f 0x2 -@ {resources.cpus_per_task} | \
        samtools sort -@ {resources.cpus_per_task} -o {output.aln_cw}
        samtools index {output.aln_cw}
        """


rule filt_aln_cw_clonal:
    input:
        aln_cw = "final_results/cw_clonal/CW.bam",
    output:
        trunc_bed = "final_results/cw_clonal/trunc_CW_coverage.bed",
        filt_aln_cw = "final_results/cw_clonal/filt_CW.bam",
    resources:
        cpus_per_task = 16,
        mem_mb = 32000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        samtools view -H {input.aln_cw} | grep '^@SQ' | sed 's/.*SN://;s/\tLN:/\t/' |
        awk 'BEGIN{{OFS="\t"}} {{if ($2 > 4000) print $1, 2000, $2-2000}}' > {output.trunc_bed}

        bedtools intersect -a {input.aln_cw} -b {output.trunc_bed} |
        samtools sort -@ {resources.cpus_per_task} -o {output.filt_aln_cw}
        samtools index {output.filt_aln_cw}
        """

rule call_cw_clonal_variants:
    input:
        ref=NAND_REF,
        ref_idx = NAND_REF + ".mmi",
        filt_aln_cw = "final_results/cw_clonal/filt_CW.bam"
    output:
        cw_clonal_dv_vcf = "final_results/cw_clonal/cw_clonal_output.vcf.gz"
    resources:
        cpus_per_task = 20,
        mem_mb = 96000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        apptainer run \
          --bind .:/dv_run \
          ~/shared/containers/deepvariant_1.10.0.sif \
          /opt/deepvariant/bin/run_deepvariant \
          --model_type=WGS \
          --ref={input.ref} \
          --reads={input.filt_aln_cw} \
          --output_vcf={output.cw_clonal_dv_vcf} \
          --num_shards={resources.cpus_per_task}
        """