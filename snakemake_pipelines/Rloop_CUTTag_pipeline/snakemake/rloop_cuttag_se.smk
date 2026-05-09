rule fastp_se:
    input:
        fq1=fq1_of
    output:
        fq1=temp("clean_data/{sample}.clean.fq.gz"),
        html="qc/fastp/{sample}.html",
        json="qc/fastp/{sample}.json"
    log:
        "logs/{sample}.fastp.log"
    threads:
        config["threads"].get("fastp", 8)
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p clean_data qc/fastp logs
        fastp \
          -i {input.fq1} \
          -o {output.fq1} \
          -w {threads} \
          -h {output.html} -j {output.json} \
          &> {log}
        """

rule bowtie2_se:
    input:
        fq1="clean_data/{sample}.clean.fq.gz"
    output:
        bam=temp("aligned_data/{sample}.sorted.bam"),
        bai=temp("aligned_data/{sample}.sorted.bam.bai")
    params:
        genome=config["genome_index"],
        extra=config.get("bowtie2_extra_se", "--very-sensitive")
    log:
        "logs/{sample}.bowtie2.log"
    threads:
        config["threads"].get("bowtie2", 16)
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p aligned_data logs
        bowtie2 -p {threads} {params.extra} -x {params.genome} \
          -U {input.fq1} 2> {log} \
          | samtools sort -@ {threads} -O bam -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """
