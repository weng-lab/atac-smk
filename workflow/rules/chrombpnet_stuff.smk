rule chrombpnet_prep_nonpeaks:
    input: "results/macs3_callpeak/{sample}-{genome}-.05/{sample}-{genome}-.05_peaks.narrowPeak"
    output: "results/chrombpnet/{sample}-{genome}/chrombpnet_prep_nonpeaks.ok"
    threads: 10
    resources: 
        slurm_partition='gpu', runtime=1440, mem_mb='96G', cpus_per_task=10, slurm_extra="--gres=gpu:1"
    log: "logs/chrombpnet/prep_nonpeaks/{sample}-{genome}.log"
    singularity: 'docker://kundajelab/chrombpnet:latest'
    shell:
        '''
        (
        chrombpnet prep nonpeaks -g resources/{wildcards.genome}.fa -p {input} -c resources/{wildcards.genome}.sizes.txt -fl resources/fold_0.json -br resources/{wildcards.genome}.blacklist.bed.gz -o results/chrombpnet/{wildcards.sample}-{wildcards.genome}
        touch results/chrombpnet/{wildcards.sample}-{wildcards.genome}/chrombpnet_prep_nonpeaks.ok
        ) &> {log}
        '''


rule chrombpnet_bias:
    input: 
        prep_nonpeaks_ok = "results/chrombpnet/{sample}-{genome}/chrombpnet_prep_nonpeaks.ok",
        peaks= "results/macs3_callpeak/{sample}-{genome}-.05/{sample}-{genome}-.05_peaks.narrowPeak",
        tagalign = "results/bedpe_tagalign/{sample}-{genome}.tagalign.gz"
        
    output: "results/chrombpnet/{sample}-{genome}/chrombpnet_bias.ok"
    threads: 10
    resources: 
        slurm_partition='gpu', runtime=1440, mem_mb='96G', cpus_per_task=10, slurm_extra="--gres=gpu:1"
    log: "logs/chrombpnet/bias/{sample}-{genome}.log"
    singularity: 'docker://kundajelab/chrombpnet:latest'
    shell:
        '''
        (
        rm -rf results/chrombpnet/{wildcards.sample}-{wildcards.genome}/bias_model
        chrombpnet bias pipeline -d ATAC -itag {input.tagalign} -g resources/{wildcards.genome}.fa -c resources/{wildcards.genome}.sizes.txt -p {input.peaks} -n results/chrombpnet/{wildcards.sample}-{wildcards.genome}_negatives.bed -fl resources/fold_0.json -b .5 -o results/chrombpnet/{wildcards.sample}-{wildcards.genome}/bias_model
        touch {output}
        ) &> {log}
        '''
 
