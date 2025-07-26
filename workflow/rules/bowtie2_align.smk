def get_k():
    if int(config['k']) == 1:
        return ""
    elif int(config['k']) == -1:
        return "--all"
    else:
        return f"-k {int(config['k'])}"

rule bowtie2_align:
    input:
        index_dir = rules.bowtie2_build.output.index_dir,
        index = multiext(
            f"{RESULTS_DIR}/bowtie2_build/{{genome}}/{{genome}}.fa",
            ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l",
            ".rev.1.bt2l", ".rev.2.bt2l"
        ),
        R1=f"{RESULTS_DIR}/processed/{{sample}}/{{sample}}_R1.fastq.gz",
        R2=f"{RESULTS_DIR}/processed/{{sample}}/{{sample}}_R2.fastq.gz"
    output:
        f"{RESULTS_DIR}/bowtie2_align/{{sample}}-{{genome}}.bam"
    params:
        k=get_k()
    container: "docker://clarity001/atac-smk:latest"
    log:
        f"{config['logdir']}/bowtie2_align/{{sample}}-{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "Running bowtie2 align"
        bowtie2 -X 2000 --mm {params.k} \
            --threads {resources.threads} \
            -x {input.index_dir}/{wildcards.genome}.fa \
            --rg-id {wildcards.sample}-{wildcards.genome} \
            --rg SM:{wildcards.sample}-{wildcards.genome} \
            -1 {input.R1} \
            -2 {input.R2} | \
            samtools view -@ {resources.threads} -1 -S /dev/stdin > /tmp/{wildcards.sample}-{wildcards.genome}.tmp.bam
        samtools sort -@ {resources.threads} -T {wildcards.sample}-{wildcards.genome} -o {output} /tmp/{wildcards.sample}-{wildcards.genome}.tmp.bam
        samtools flagstat -@ {resources.threads} --output-fmt "json" {output} > {output}.flagstat.json
        rm /tmp/{wildcards.sample}-{wildcards.genome}.tmp.bam
        """

