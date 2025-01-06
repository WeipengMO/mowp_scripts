import os


pwd = os.getcwd()

configfile: "config.yaml"


samples = [
    d for d in os.listdir() if os.path.isdir(os.path.join(d)) 
    and not d.startswith('.') 
    and not d.endswith('logs')]

rule all:
    input:
        expand('{sample}/{sample}/outs/feature_slice.h5', sample=samples)

rule spaceranger_count:
    input:
        '{sample}'
    output:
        '{sample}/{sample}/outs/feature_slice.h5'
    params:
        transcriptome=config['transcriptome'],
        fastqs='./',
        probeset=config['probe_set'],
        image='{sample}'+'.jpg',
        cytaimage='{sample}'+'.tif',
        slide=lambda wildcards: config['sample_info'][wildcards.sample]['slide'],
        area=lambda wildcards: config['sample_info'][wildcards.sample]['area'],
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

spaceranger count \
    --id={wildcards.sample} \
    --transcriptome={params.transcriptome} \
    --probe-set={params.probeset} \
    --fastqs={params.fastqs} \
    --image={params.image} \
    --cytaimage={params.cytaimage} \
    --slide={params.slide} \
    --area={params.area} \
    --reorient-images=true \
    --localcores={threads} \
    --localmem=128 \
    --create-bam=false &> ../{log}
        '''