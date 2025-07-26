rule frip_all:
    input: 
        peaks = f"{RESULTS_DIR}/macs3_callpeak/{{q}}/{{sample}}-{{genome}}_peaks.narrowPeak", 
        reads = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.tagalign.gz"
    output:
        frip_report = f"{RESULTS_DIR}/frip_all/{{q}}/{{sample}}-{{genome}}.txt"
    container: config['container']
    log: 
        f"{config['logdir']}/frip_all/{{sample}}-{{genome}}-{{q}}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/rules/scripts/calculate_frip.py {input.reads} {input.peaks} {output.frip_report}
        """ 
