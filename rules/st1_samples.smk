NAND_REF = "/storage/home/hcoda1/1/jshih41/r-alind6-0/data/all_ncbi_blast_genomes/Nand_genome.fna"

ST1_SAMPLES = ["SRR36118145", "SRR20078973", "SRR20078977", "ERR7745403", "SRR20078971", "ERR7738269", "ERR7738585", "ERR7738185", "ERR7745412", "ERR1137517", "ERR14044664", "ERR7745291", "ERR14040741", "ERR7738582", "ERR1305896", "ERR14025145", "ERR7745755", "ERR13966020", "ERR4335103", "ERR14026386", "ERR14020803", "ERR14035067", "ERR7738560", "ERR14039449", "ERR4341940", "ERR14022091"]


rule fastp_st1_samples:
    input:
        fwd="/storage/home/hcoda1/1/jshih41/scratch/blast_sra/{st1_sample}_1.fastq.gz",
        rev="/storage/home/hcoda1/1/jshih41/scratch/blast_sra/{st1_sample}_2.fastq.gz",
    output:
        fwd_trimmed="final_results/trimmed/{st1_sample}_1.fastq.gz",
        rev_trimmed="final_results/trimmed/{st1_sample}_2.fastq.gz",
        html="final_results/qc/fastp/{st1_sample}.html",
        json="final_results/qc/fastp/{st1_sample}.json"
    params:
        extra="--detect_adapter_for_pe"
    resources:
        cpus_per_task = 8,
        mem_mb = 16000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        fastp \
            -i {input.fwd} -I {input.rev} \
            -o {output.fwd_trimmed} -O {output.rev_trimmed} \
            --html {output.html} \
            --json {output.json} \
            --thread {resources.cpus_per_task} \
            {params.extra}
        """

rule aln_st1_samples:
    input:
        ref=NAND_REF,
        ref_idx = NAND_REF + ".mmi",
        fwd_trimmed="final_results/trimmed/{st1_sample}_1.fastq.gz",
        rev_trimmed="final_results/trimmed/{st1_sample}_2.fastq.gz"
    output:
        aln_st1_sample = "final_results/st1_samples/{st1_sample}.bam",
    resources:
        cpus_per_task = 8,
        mem_mb = 24000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        minimap2 -ax sr -t {resources.cpus_per_task} {input.ref_idx} {input.fwd_trimmed} {input.rev_trimmed} | \
        samtools view -b -q 30 -F 0x4 -f 0x2 -@ {resources.cpus_per_task} | \
        samtools sort -@ {resources.cpus_per_task} -o {output.aln_st1_sample}
        samtools index {output.aln_st1_sample}
        """

rule filt_aln_st1_samples:
    input:
        aln_st1_sample = "final_results/st1_samples/{st1_sample}.bam"
    output:
        trunc_bed = "final_results/st1_samples/{st1_sample}_trunc_regions.bed",
        filt_aln_st1_sample = "final_results/st1_samples/filt_{st1_sample}.bam"
    resources:
        cpus_per_task = 8,
        mem_mb = 16000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        samtools view -H {input.aln_st1_sample} | grep '^@SQ' | sed 's/.*SN://;s/\tLN:/\t/' |
        awk 'BEGIN{{OFS="\t"}} {{if ($2 > 4000) print $1, 2000, $2-2000}}' > {output.trunc_bed}

        bedtools intersect -a {input.aln_st1_sample} -b {output.trunc_bed} |
        samtools sort -@ {resources.cpus_per_task} -o {output.filt_aln_st1_sample}
        samtools index {output.filt_aln_st1_sample}
        """


rule call_st1_sample_variants:
    input:
        ref=NAND_REF,
        ref_idx = NAND_REF + ".mmi",
        filt_aln_st1_sample = "final_results/st1_samples/filt_{st1_sample}.bam"
    output:
        st1_sample_dv_vcf = "final_results/st1_samples/{st1_sample}_output.vcf.gz"
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
          --reads={input.filt_aln_st1_sample} \
          --output_vcf={output.st1_sample_dv_vcf} \
          --sample_name={wildcards.st1_sample} \
          --num_shards={resources.cpus_per_task}
        """

rule merge_st1_vcf:
    input:
        all_st1_sample_dv_vcf = expand("final_results/st1_samples/{st1_sample}_output.vcf.gz", st1_sample = ST1_SAMPLES)
    output:
        merged_vcf = "final_results/st1_samples/merged_vcf.gz",
        filt_merged_vcf = "final_results/st1_samples/updated_filt_merged_vcf.gz"
    resources:
        cpus_per_task = 8,
        mem_mb = 32000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        bcftools merge -Oz -o {output.merged_vcf} {input.all_st1_sample_dv_vcf} --force-samples

        bcftools filter \
        -i 'QUAL>3 && DP>3' \
        {output.merged_vcf} \
        -Oz -o {output.filt_merged_vcf}
        """