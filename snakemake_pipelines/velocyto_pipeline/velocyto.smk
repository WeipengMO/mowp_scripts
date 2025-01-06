configfile: 'config.yaml'


rule all:
    input:
        expand('{sample}/{sample}/velocyto/{sample}.loom', sample=config['samples'])


rule velocyto:
    input:
        sample_name= '{sample}/{sample}',  # /path/to/cellranger-runs/sample_name
    output:
        loom='{sample}/{sample}/velocyto/{sample}.loom',
    params:
        gtf = config['gtf'],
    threads: 24
    conda:
        "scvelo"
    shell:
        """
velocyto run10x -@ {threads} {input.sample_name} {params.gtf}
        """