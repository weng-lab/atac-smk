rule multiqc:
    input:
        macs3_output = expand(
            f"{RESULTS_DIR}/macs3_callpeak/{{q}}/{{sample}}-{{genome}}_peaks.narrowPeak",
            sample=SAMPLES,
            genome=GENOMES,
            q=QVALS
        ),    
    output:
        datadir = directory(f"{RESULTS_DIR}/multiqc/multiqc_data"),
        html_report = f"{RESULTS_DIR}/multiqc/multiqc_report.html",
    params: 
        multiqc_conf = f"--config {config["input_files"]['multiqc_config']}" if config['input_files'].get('multiqc_config') else ""
    container: config['container']
    log: f"{config['logdir']}/mutliqc.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Running multiqc"

        multiqc {params.multiqc_conf} {RESULTS_DIR} --outdir {RESULTS_DIR}/multiqc

        echo "$(date): Finished multiqc"
        """
