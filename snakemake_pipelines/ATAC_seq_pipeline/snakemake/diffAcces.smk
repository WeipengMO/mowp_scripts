from glob import glob

# snakemake --use-conda -pj 120 --rerun-triggers mtime -s snakemake/diffAcces.smk diffAcces/featureCounts.txt

configfile:
    "config.yml"

suffix = config['suffix']
sep = config['sep']

sample = [os.path.basename(fn.split(sep)[0]) for fn in glob(f'raw_data/*{sep}{suffix}')]

rule merge_peak:
    input:
        expand("genrich_result/{sample}_peaks.narrowPeak", sample=sample)
    output:
        bed="diffAcces/peaks_merged.bed",
        saf="diffAcces/peaks_merged.saf"
    shell:
        '''
cat {input} | sortBed | mergeBed > {output.bed}
awk -v OFS="\t" '{{print "peak_" NR, $0, "."}}' {output.bed} > {output.saf}
        '''


rule run_featureCounts:
    input:
        saf="diffAcces/peaks_merged.saf",
        bam=expand("aligned_data/{sample}.sorted.rmdup.bam", sample=sample)
    output:
        counts = "diffAcces/featureCounts.txt"
    threads: 32
    conda:
        "chipseq"
    log:
        "logs/featureCounts.log"
    shell:
        '''
featureCounts -T {threads} -p -F SAF -a {input.saf} --fracOverlap 0.2 -o merged.genrich.counts {input.bam} &> {log}
sed 's/aligned_data\///g' merged.genrich.counts | awk '(NR>1)' > {output.counts}
rm merged.genrich.counts
        '''

