import os


pwd = os.getcwd()

configfile: "config.yaml"


samples = [
    d for d in os.listdir() if os.path.isdir(os.path.join(d)) 
    and not d.startswith('.') 
    and not d.endswith('logs')]

rule all:
    input:
        expand('{sample}/{sample}/outs/filtered_feature_bc_matrix.h5', sample=samples)

rule cellranger_count:
    input:
        '{sample}'
    output:
        '{sample}/{sample}/outs/filtered_feature_bc_matrix.h5'
    params:
        transcriptome=config['transcriptome'],
        fastqs=pwd + '/' + '{sample}'
    threads:
        24
    log:
        "logs/cellranger_count_{sample}.log"
    shell:
        '''
cd {wildcards.sample}

if [ -d {wildcards.sample} ]; then
  rm -rf {wildcards.sample}
fi

cellranger count \
    --id={wildcards.sample} \
    --transcriptome={params.transcriptome} \
    --fastqs={params.fastqs} \
    --sample={wildcards.sample} \
    --localcores={threads} \
    --localmem=64 &> ../{log}
        '''