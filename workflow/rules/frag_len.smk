rule frag_len:
    input:
        bam = f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}.bam"
    output: 
        png = f"{RESULTS_DIR}/frag_len/{{sample}}-{{genome}}.png", 
        metrics = f"{RESULTS_DIR}/frag_len/{{sample}}-{{genome}}.txt"
    container: "docker://clarity001/atac-smk:latest"
    log: 
        f"{config['logdir']}/frag_len/{{sample}}-{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/rules/scripts/plot_fragment_length_distr.py {input.bam} {resources.threads} {output.png} {output.metrics}
        """

