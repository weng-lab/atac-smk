def get_exclude_flag(wildcards):
    if int(config['k']) == -1 or int(config['k']) > 1:
        return("1548")
    else:
        return("1804")

rule filter:
    input: 
        f"{RESULTS_DIR}/bowtie2_align/{{sample}}-{{genome}}.bam"
    output: 
        f"{RESULTS_DIR}/filter/{{sample}}-{{genome}}.bam"
    log: 
        f"{config['logdir']}/filter/{{sample}}-{{genome}}.log"
    params:
        exclude_flag = lambda wildcards: get_exclude_flag(wildcards),
        prefix = lambda wildcards: "-".join(wildcards)
    shell:
        """
        exec >> {log} 2>&1
        echo {params.prefix}
        samtools view -@ {resources.threads} -F {params.exclude_flag} -f 2 -q 30 -u {input} | \
        samtools sort -@ {resources.threads} -n /dev/stdin -o /tmp/{params.prefix}.tmp.nmsrt.bam -T {params.prefix}

        samtools view -@ {resources.threads} -h /tmp/{params.prefix}.tmp.nmsrt.bam |
        samtools fixmate -@ {resources.threads} -r /dev/stdin /tmp/{params.prefix}.tmp.fixmate.bam
        
        samtools view -@ {resources.threads} -F {params.exclude_flag} -f 2 -u /tmp/{params.prefix}.tmp.fixmate.bam | \
            samtools sort -@ {resources.threads} /dev/stdin -o /tmp/{params.prefix}.tmp.coordsort.bam 

        samtools view -@ {resources.threads} -h /tmp/{params.prefix}.tmp.coordsort.bam | grep -v "chrM" | samtools view -@ {resources.threads} -b > {output}


        samtools quickcheck -v {output}
        samtools index -@ {resources.threads} {output} {output}.bai
        samtools flagstat -@ {resources.threads} {output} > {output}.flagstat

        rm /tmp/{params.prefix}.tmp.nmsrt.bam
        rm /tmp/{params.prefix}.tmp.fixmate.bam
        rm /tmp/{params.prefix}.tmp.coordsort.bam
        """

