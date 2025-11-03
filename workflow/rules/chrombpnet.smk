rule chrombpnet_prep_nonpeaks:
    input: 
        peaks = f"{RESULTS_DIR}/macs3_callpeak/{config['qvals'][0]}/{{sample}}.{{genome}}_peaks.narrowPeak",
        chromsizes = rules.bowtie2_build.output.chromsizes,
    output: f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_prep_nonpeaks.continue"
    log: 
        f"{config['logdir']}/chrombpnet/prep_nonpeaks/{{sample}}.{{genome}}.log"
    params:
        fasta = lambda wildcards: config['genomes'].get(wildcards.genome),
        config_json = config["input_files"]['chrombpnet_config_json'],
        blacklist = config["input_files"]['blacklist']
    singularity: 'docker://kundajelab/chrombpnet:latest'
    shell:
        '''
        exec >> {log} 2>&1

        chrombpnet prep nonpeaks \
            -g {params.fasta} \
            -p {input.peaks} \
            -c {input.chromsizes} \
            -fl {params.config_json} \
            -br {params.blacklist} \
            -o results/chrombpnet/{wildcards.sample}.{wildcards.genome}
        touch {output}
        '''

rule chrombpnet_bias:
    input: 
        prep_nonpeaks_singal = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_prep_nonpeaks.continue",
        chromsizes = rules.bowtie2_build.output.chromsizes,
        peaks = f"{RESULTS_DIR}/macs3_callpeak/{config['qvals'][0]}/{{sample}}.{{genome}}_peaks.narrowPeak",
        tagalign = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.tagalign.gz"
    output: f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_bias.continue"
    log: 
        f"{config['logdir']}/chrombpnet/bias/{{sample}}.{{genome}}.log"
    params:
        fasta = lambda wildcards: config['genomes'].get(wildcards.genome),
        config_json = config["input_files"]['chrombpnet_config_json'],
        bias_threshold_factor = config['bias_threshold_factor']
    singularity: 'docker://kundajelab/chrombpnet:latest'
    shell:
        '''
        exec >> {log} 2>&1

        OUTPUT_DIR=results/chrombpnet/{wildcards.sample}.{wildcards.genome}/bias_model

        rm -r $OUTPUT_DIR
        chrombpnet bias pipeline \
            -d ATAC -itag {input.tagalign} {params.fasta} \
            -c {input.chromsizes} \
            -p {input.peaks} \
            -n results/chrombpnet/{wildcards.sample}.{wildcards.genome}.negatives.bed \
            -fl {params.config_json} \
            -b {params.bias_threshold_factor} \
            -o $OUTPUT_DIR
        touch {output}
        '''
 
