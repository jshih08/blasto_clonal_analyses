NAND_REF = "/storage/home/hcoda1/1/jshih41/r-alind6-0/data/all_ncbi_blast_genomes/Nand_genome.fna"

rule aln_cw_mixture:
    input:
        cw_mixture_raw = "/storage/home/hcoda1/1/jshih41/scratch/cw_isolate_long_reads/XLIND_20250314_R84050_PL20051-001_1-1-B01_bc2022-bc2022.hifi_reads.fastq.gz",
        ref=NAND_REF,
        ref_idx = NAND_REF + ".mmi"
    output:
        aln_cw_mixture = "final_results/cw_mixture/CW_mixture.bam",
    resources:
        cpus_per_task = 16,
        mem_mb = 32000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        minimap2 -ax map-pb -t {resources.cpus_per_task} {input.ref_idx} {input.cw_mixture_raw} |
        samtools view -q 30 -b -F 0x4 -@ {resources.cpus_per_task} | \
        samtools sort -@ {resources.cpus_per_task} -o {output.aln_cw_mixture}
        samtools index {output.aln_cw_mixture}
        """

rule filt_aln_cw_mixture:
    input:
        aln_cw_mixture = "final_results/cw_mixture/CW_mixture.bam"
    output:
        trunc_bed = "final_results/cw_clonal/trunc_CW_mixture_coverage.bed",
        filt_aln_cw_mixture = "final_results/cw_mixture/filt_cw_mixture.bam",
    resources:
        cpus_per_task = 16,
        mem_mb = 32000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        samtools view -H {input.aln_cw_mixture} | grep '^@SQ' | sed 's/.*SN://;s/\tLN:/\t/' |
        awk 'BEGIN{{OFS="\t"}} {{if ($2 > 4000) print $1, 2000, $2-2000}}' > {output.trunc_bed}

        bedtools intersect -a {input.aln_cw_mixture} -b {output.trunc_bed} |
        samtools sort -@ {resources.cpus_per_task} -o {output.filt_aln_cw_mixture}
        samtools index {output.filt_aln_cw_mixture}

        """

rule call_cw_mixture_variants:
    input:
        ref=NAND_REF,
        ref_idx = NAND_REF + ".mmi",
        filt_aln_cw_mixture = "final_results/cw_mixture/filt_cw_mixture.bam"
    output:
        cw_mixture_dv_vcf = "final_results/cw_mixture/cw_mixture_output.vcf.gz"
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
          --model_type=PACBIO \
          --ref={input.ref} \
          --reads={input.filt_aln_cw_mixture} \
          --output_vcf={output.cw_mixture_dv_vcf} \
          --num_shards={resources.cpus_per_task}
        """