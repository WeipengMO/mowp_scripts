'''
Author       : windz
Date         : 2022-03-30 19:18:32
LastEditTime : 2022-04-02 20:21:40
LastEditors  : windz
FilePath     : nanonome_data_parse.smk
'''

utildir = config['utildir']
codedir = config['codedir']

rule mtsv_to_mbed:
    input:
        "mcall/{sample}.{mod}.tsv"
    output:
        "mbed/{sample}.{mod}.meth.bed.gz"
    params:
        genome=config['genome'],
        call_threshold=lambda wildcards: 1.5 if wildcards.mod == 'cpg' else 1.0
    shell:
        "{utildir}/mtsv2bedGraph.py -i {input} -q {wildcards.mod} -c {params.call_threshold} --nome -g {params.genome} | "
        "sort -k1,1 -k2,2n | bgzip "
        "> {output} && "
        "tabix -p bed {output}"


rule mbed_to_mfreq:
    input:
        "mbed/{sample}.{mod}.meth.bed.gz"
    output:
        mfreq="mfreq/{sample}.{mod}.mfreq.txt.gz",
        tabix="mfreq/{sample}.{mod}.mfreq.txt.gz.tbi",
        log="mfreq/{sample}.{mod}.mfreq.log"
    shell:
        "python -u {utildir}/parseMethylbed.py frequency -v "
        "-i {input} -m {wildcards.mod} 2> {output.log} | "
        "bgzip > {output.mfreq} && "
        "tabix -b 2 -e 2 {output.mfreq}"


#isac methylbam code
rule methylbam:
    input:
        bam = 'aligned_data/{sample}.sorted.bam',
        cpg = 'mbed/{sample}.cpg.meth.bed.gz',
        gpc = 'mbed/{sample}.gpc.meth.bed.gz'
    output:
        "methylbam/{sample}.meth.bam"
    threads: 16
    run:
        ref=config['genome']
        shell("{utildir}/convert_bam_for_methylation.py -t {threads} --windowsize 1000000 --verbose -b {input.bam}"+
              " -c {input.cpg} -g {input.gpc} -f {ref} | samtools sort -@ {threads} -o {output}")
        shell("samtools index -@ {threads} {output}")


# 这个好像会过滤掉一些?
rule mfreq_to_bigwig:
    input:
        "mfreq/{sample}.{mod}.mfreq.txt.gz"
    output:
        methwig=temp("bigwig/{sample}.{mod}.methylation.wig"),
        covwig=temp("bigwig/{sample}.{mod}.methcoverage.wig"),
        log="bigwig/{sample}.{mod}.wig.log"
    shell:
        "python {codedir}/scripts/makeWig.py -v -i {input} "
        "-o {output.methwig} -c {output.covwig} &> {output.log}"


rule wig_to_bigwig:
    input:
        "bigwig/{sample}.{mod}.{type}.wig",
        config['genomesize']
    output:
        "bigwig/{sample}.{mod}.{type}.bw"
    shell:
        "wigToBigWig {input} {output}"