rule macs3_signal:
    input:
        tagalign = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.tagalign.gz"
    output:
        bedgraph = f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}_treat_pileup.bdg"
    container: "docker://clarity001/atac-smk:latest"
    log: 
        f"{config['logdir']}/macs3_signal/{{sample}}-{{genome}}.log"
    params:
        prefix = lambda wildcards: "-".join(wildcards),
    shell:
        """
        exec >> {log} 2>&1
        macs3 callpeak -f BED -t {input.tagalign} -n {params.prefix} --outdir {RESULTS_DIR}/macs3_signal/ --nomodel --shift -75 --keep-dup all --extsize 150 -g hs -q .01 -B
        """

