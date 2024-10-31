
rule convert_homer_peak:
    input:
        "diffAcces/diffAcces.{sample}.fdr_05.fc_-1.annotation.tsv"
    output:
        "homer/{sample}.peaks_homer.txt",
    shell:
        '''
awk 'NR>1 {{print $6, $1, $2, $3, "+"}}' OFS='\t' {input} > {output}
        '''