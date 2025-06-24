sample_table = config['sample_table']

treatment_vs_control = {}
with open(sample_table, 'r') as f:
    for line in f:
        if line[0] == '#':
            continue
        treatment = line.rstrip()

def macs2_callpeak_input(wildcard):
    return {
        'treatment': expand(
            'aligned_data/{treatment}.sorted.rmdup.bam',
            treatment=wildcard.treatment),
        # No control input since it's without control
    }

rule macs2_callpeak:
    input:
        unpack(macs2_callpeak_input)
    output:
        treat='macs2_result/{treatment}_peaks.narrowPeak'
    threads: 1
    params:
        name='{treatment}',
        out_dir='macs2_result/',
        gsize=config['gsize']
    log:
        'logs/{treatment}.log'
    conda: 'chipseq'
    shell:
        '''
macs2 callpeak -t {input.treatment} -f BAM -g {params.gsize} -n {params.name} -B --SPMR -q 0.01 --outdir {params.out_dir} &> {log}
        '''

rule bamCompare:
    input:
        unpack(macs2_callpeak_input)
    output:
        'bw_compare/{treatment}.compare.bw'
    threads: 32
    shell:
        '''
bamCompare -b1 {input.treatment} -o {output} --binSize 10 --skipNAs --centerReads --scaleFactorsMethod SES -p {threads} &> {log}
        '''

import yaml
try:
    replicate_data = config['replicate']
    with open(replicate_data, 'r') as f:
        replicate_data = yaml.load(f, Loader=yaml.Loader)
except KeyError:
    print('Please provide replicate data in the config file if you want to run the replicate_intersect rule.')
    replicate_data = None


def replicate_intersect_input(wildcard):
    wildcard = str(wildcard)
    return expand(
        'macs2_result/{replicate}_peaks.narrowPeak',
        replicate=replicate_data[wildcard])

# snakemake --use-conda -pj 120 -np replicate_intersect/{co,single}_intersect.bed
rule replicate_intersect:
    input:
        unpack(replicate_intersect_input)
    output:
        'replicate_intersect/{replicate}_intersect.bed'
    shell:
        '''
bedtools intersect -a {input[0]} -b {input[1]} -f 0.50 -r | cut -f1-4 > {output} &> {log}
        '''