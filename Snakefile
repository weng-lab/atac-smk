import os
import glob
import sys

configfile: "config/config.yml"
workdir: config["workdir"]

os.makedirs(config.get("picard_tmpdir", ".picard"), exist_ok=True)

include: "workflow/rules/common.smk"
include: "workflow/rules/fastp.smk"
include: "workflow/rules/bowtie2_build.smk"
include: "workflow/rules/bowtie2_align.smk"
include: "workflow/rules/filter.smk"
include: "workflow/rules/picard.smk"
include: "workflow/rules/bam_to_tagalign.smk"
include: "workflow/rules/macs3_signal.smk"
include: "workflow/rules/bedgraphtobigwig.smk"
include: "workflow/rules/macs3_callpeak.smk"
include: "workflow/rules/frag_len.smk"
include: "workflow/rules/tss_enrichment.smk"
include: "workflow/rules/frip_all.smk"
include: "workflow/rules/pseudoreplicated_peaks.smk"
include: "workflow/rules/multiqc.smk"
include: "workflow/rules/chrombpnet.smk"

# wildcard_constraints:
#     q="^(?:0*(?:\.\d+)?|1(\.0*)?)$"

onerror:
    shell("rm -r results/chrombpnet")

rule all:
    input:
        # fastp (replaces fastqc-cutadapt-fastqc sequence)
        expand(
          [
            f"{RESULTS_DIR}/fastp/{{sample}}_R1.fastq.gz",
            f"{RESULTS_DIR}/fastp/{{sample}}_R2.fastq.gz",
          ],
          sample=SAMPLES
        ),
       
        # # Alignment
        # bowtie2_build
        expand(
          f"{RESULTS_DIR}/bowtie2_build/{{genome}}/.continue",
          genome=GENOMES
        ),
        # bowtie2_align
        expand(
            f"{RESULTS_DIR}/bowtie2_align/{{sample}}-{{genome}}.bam",
            sample=SAMPLES,
            genome=GENOMES 
        ),

        # Filtering
        expand(
            f"{RESULTS_DIR}/filter/{{sample}}-{{genome}}.bam",
            sample=SAMPLES,
            genome=GENOMES,
        ),

        # MarkDuplicates
        expand(
            f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}.bam",
            sample=SAMPLES,
            genome=GENOMES,
        ),
        
        # BEDPE & tagAlign
        expand(
            f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.bedpe",
            sample=SAMPLES,
            genome=GENOMES,
        ),

        # MACS3 signal bedGraph
        expand(
            f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}_treat_pileup.bdg",
            sample=SAMPLES,
            genome=GENOMES,
        ),

        # bigtools bedgraphtobigwig
        expand(
            f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}.bigWig",
            sample=SAMPLES,
            genome=GENOMES,
        ),

        # fold change signal
        expand(
            f"{RESULTS_DIR}/macs3_fold_change_signal/{{sample}}-{{genome}}_FE.bw",
            sample=SAMPLES,
            genome=GENOMES,
        ),

        # pvalue signal
        expand(
            f"{RESULTS_DIR}/macs3_pvalue_signal/{{sample}}-{{genome}}_ppois.bw",
            sample=SAMPLES,
            genome=GENOMES,
        ),

        # MACS3 callpeak
        expand(
            f"{RESULTS_DIR}/macs3_callpeak/{{q}}/{{sample}}.{{genome}}_peaks.narrowPeak",
            sample=SAMPLES,
            genome=GENOMES,
            q=QVALS
        ),

        # QC: fragment length
        expand(
            [
                f"{RESULTS_DIR}/frag_len/{{sample}}-{{genome}}.png",
                f"{RESULTS_DIR}/frag_len/{{sample}}-{{genome}}.txt"
            ],
            sample=SAMPLES,
            genome=GENOMES
        ),

        # QC: TSS enrichment
        expand(
            [
                f"{RESULTS_DIR}/tss_enrichment/{{sample}}-{{genome}}.png",
                f"{RESULTS_DIR}/tss_enrichment/{{sample}}-{{genome}}.txt"
            ],
            sample=SAMPLES,
            genome=GENOMES
        ),

        # QC: FRiP
        expand(
            f"{RESULTS_DIR}/frip_all/{{q}}/{{sample}}-{{genome}}.txt",
            sample=SAMPLES,
            genome=GENOMES,
            q=QVALS
        ),

        # Pseudoreplicates & IDR
        expand(
            [
                f"{RESULTS_DIR}/pseudoreps/{{sample}}-{{genome}}-pseudorep1.tagalign",
                f"{RESULTS_DIR}/pseudoreps/{{sample}}-{{genome}}-pseudorep2.tagalign",
            ],
            sample=SAMPLES,
            genome=GENOMES
        ),
        expand(
            f"{RESULTS_DIR}/macs3_pseudoreps/{{sample}}-{{genome}}-{{pseudorep}}_peaks.narrowPeak",
            sample=SAMPLES,
            genome=GENOMES,
            pseudorep=["pseudorep1","pseudorep2"]
        ),
        expand(
            f"{RESULTS_DIR}/idr/{{sample}}-{{genome}}.IDR.bed",
            genome=GENOMES,
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS_DIR}/overlap_peaks/{{sample}}-{{genome}}.overlap_peaks.bed",
            genome=GENOMES,
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS_DIR}/overlap_peaks_jill/{{q}}/{{sample}}-{{genome}}.overlap_peaks.bed",
            q=QVALS,
            sample=SAMPLES,
            genome=GENOMES
        ),
        expand(
            [
                f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_bias.continue",
                f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_prep_nonpeaks.continue"
            ],
            sample=SAMPLES,
            genome=GENOMES
        ),

        # MultQC
        "results/multiqc/multiqc_report.html"

    








        
