rule tss_enrichment:
    input: 
        bam = f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}.bam"
    output:
        png = f"{RESULTS_DIR}/tss_enrichment/{{sample}}-{{genome}}.png", 
        metrics = f"{RESULTS_DIR}/tss_enrichment/{{sample}}-{{genome}}.txt"
    params:
        tss_reference_bed = config['input_files']['tss_reference_bed']
    container: "docker://clarity001/atac-smk:latest"
    log: 
        f"{config['logdir']}/tss_enrichment/{{sample}}-{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/rules/scripts/tss_enrichment.py {input.bam} {params.tss_reference_bed} {output.png} {output.metrics}
        """

