rule picard:
    input: 
        bam = f"{RESULTS_DIR}/filter/{{sample}}-{{genome}}.bam"
    output: 
        deduped_bam = f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}.bam", 
        picard_report = f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}.txt"
    container: "docker://clarity001/atac-smk:latest"
    log: 
        f"{config['logdir']}/picard/{{sample}}-{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        java -jar /usr/bin/picard.jar MarkDuplicates \
            -INPUT {input.bam} \
            -OUTPUT {output.deduped_bam} \
            -METRICS_FILE {output.picard_report} \
            -REMOVE_DUPLICATES TRUE \
            -ASSUME_SORT_ORDER coordinate \
            -USE_JDK_DEFLATER TRUE \
            -USE_JDK_INFLATER TRUE \
            -VALIDATION_STRINGENCY LENIENT
        samtools quickcheck -v {output.deduped_bam}
        samtools index -@ {resources.threads} {output.deduped_bam} {output.deduped_bam}.bai
        samtools flagstat -@ {resources.threads} --output-fmt "json" {output.deduped_bam} > {output.deduped_bam}.flagstat.json
        """
