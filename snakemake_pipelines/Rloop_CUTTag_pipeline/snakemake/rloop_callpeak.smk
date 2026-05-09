rule filter_bam:
    input:
        bam="aligned_data/{sample}.sorted.bam",
        bai="aligned_data/{sample}.sorted.bam.bai"
    output:
        bam=temp("aligned_data/{sample}.qfiltered.bam"),
        bai=temp("aligned_data/{sample}.qfiltered.bam.bai")
    params:
        mapq=config.get("mapq", 30),
        blacklist=config.get("blacklist", "")
    log:
        "logs/{sample}.filter_bam.log"
    threads:
        config["threads"].get("samtools", 8)
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p aligned_data logs
        samtools view -@ {threads} -b -q {params.mapq} {input.bam} > aligned_data/{wildcards.sample}.mapq.bam
        if [ -n "{params.blacklist}" ] && [ -s "{params.blacklist}" ]; then
            bedtools intersect -v -abam aligned_data/{wildcards.sample}.mapq.bam -b {params.blacklist} > {output.bam}
        else
            mv aligned_data/{wildcards.sample}.mapq.bam {output.bam}
        fi
        samtools index -@ {threads} {output.bam}
        """

rule final_bam:
    input:
        bam="aligned_data/{sample}.qfiltered.bam"
    output:
        bam="aligned_data/{sample}.final.bam",
        bai="aligned_data/{sample}.final.bam.bai",
        metrics="qc/{sample}.markdup.metrics.txt"
    params:
        picard=config["picard_path"],
        remove_duplicates=lambda wc: str(config.get("remove_duplicates", True)).lower()
    log:
        "logs/{sample}.MarkDuplicates.log"
    threads:
        config["threads"].get("samtools", 8)
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p qc logs aligned_data
        if [ "{params.remove_duplicates}" = "true" ]; then
            java -jar {params.picard} MarkDuplicates \
              REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.01 \
              I={input.bam} O={output.bam} M={output.metrics} &> {log}
        else
            cp {input.bam} {output.bam}
            echo "remove_duplicates=false; MarkDuplicates skipped" > {output.metrics}
        fi
        samtools index -@ {threads} {output.bam}
        """

rule flagstat:
    input:
        "aligned_data/{sample}.final.bam"
    output:
        "qc/{sample}.flagstat.txt"
    threads:
        config["threads"].get("samtools", 4)
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p qc
        samtools flagstat -@ {threads} {input} > {output}
        """

rule bamCoverage_CPM:
    input:
        "aligned_data/{sample}.final.bam"
    output:
        "bw_files/{sample}.CPM.bw"
    params:
        bin_size=config.get("bin_size", 10)
    log:
        "logs/{sample}.bamCoverage.log"
    threads:
        config["threads"].get("deeptools", 12)
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p bw_files logs
        bamCoverage \
          --bam {input} \
          -o {output} \
          --binSize {params.bin_size} \
          --normalizeUsing CPM \
          --skipNonCoveredRegions \
          --numberOfProcessors {threads} \
          &> {log}
        """

rule deeptools_profile:
    input:
        bw="bw_files/{sample}.CPM.bw"
    output:
        matrix=temp("deeptools_profile/{sample}.matrix.gz"),
        profile="deeptools_profile/{sample}.scale.png",
        heatmap="deeptools_profile/{sample}.heatmap.png"
    params:
        bed=config.get("regions_bed", "")
    log:
        "logs/{sample}.computeMatrix.log"
    threads:
        config["threads"].get("deeptools", 12)
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p deeptools_profile logs
        if [ -n "{params.bed}" ] && [ -s "{params.bed}" ]; then
            computeMatrix scale-regions \
              -b 1000 -a 1000 \
              -R {params.bed} -S {input.bw} \
              --skipZeros -o {output.matrix} -p {threads} &> {log}
            plotProfile -m {output.matrix} -out {output.profile} &>> {log}
            plotHeatmap -m {output.matrix} -out {output.heatmap} &>> {log}
        else
            echo "regions_bed not provided; skip computeMatrix" > {log}
            touch {output.matrix} {output.profile} {output.heatmap}
        fi
        """

rule macs2_nocontrol:
    input:
        treatment="aligned_data/{sample}.final.bam"
    output:
        "macs2/nocontrol/{sample}/{sample}_peaks.narrowPeak"
    params:
        outdir="macs2/nocontrol/{sample}",
        name="{sample}",
        gsize=config["gsize"],
        q=config["macs2"].get("qvalue", 0.01),
        fmt=lambda wc: config["macs2"].get("format_pe", "BAMPE") if config["mode"] == "pe" else "BAM",
        extra=config["macs2"].get("extra", "--keep-dup all")
    log:
        "logs/{sample}.macs2_nocontrol.log"
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p {params.outdir} logs
        macs2 callpeak \
          -t {input.treatment} \
          -f {params.fmt} \
          -g {params.gsize} \
          -n {params.name} \
          -B --SPMR \
          -q {params.q} \
          {params.extra} \
          --outdir {params.outdir} \
          &> {log}
        """

rule macs2_vs_background:
    input:
        treatment="aligned_data/{sample}.final.bam",
        control=bg_bam_of
    output:
        "macs2/vs_background/{sample}/{sample}_peaks.narrowPeak"
    params:
        outdir="macs2/vs_background/{sample}",
        name="{sample}",
        gsize=config["gsize"],
        q=config["macs2"].get("qvalue", 0.01),
        fmt=lambda wc: config["macs2"].get("format_pe", "BAMPE") if config["mode"] == "pe" else "BAM",
        extra=config["macs2"].get("extra", "--keep-dup all")
    log:
        "logs/{sample}.macs2_vs_background.log"
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p {params.outdir} logs
        macs2 callpeak \
          -t {input.treatment} \
          -c {input.control} \
          -f {params.fmt} \
          -g {params.gsize} \
          -n {params.name} \
          -B --SPMR \
          -q {params.q} \
          {params.extra} \
          --outdir {params.outdir} \
          &> {log}
        """

rule macs2_rnaseh_sensitive:
    input:
        treatment="aligned_data/{sample}.final.bam",
        control=rnaseh_bam_of
    output:
        "macs2/rnaseh_sensitive/{sample}/{sample}_peaks.narrowPeak"
    params:
        outdir="macs2/rnaseh_sensitive/{sample}",
        name="{sample}",
        gsize=config["gsize"],
        q=config["macs2"].get("qvalue", 0.01),
        fmt=lambda wc: config["macs2"].get("format_pe", "BAMPE") if config["mode"] == "pe" else "BAM",
        extra=config["macs2"].get("extra", "--keep-dup all")
    log:
        "logs/{sample}.macs2_rnaseh_sensitive.log"
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p {params.outdir} logs
        macs2 callpeak \
          -t {input.treatment} \
          -c {input.control} \
          -f {params.fmt} \
          -g {params.gsize} \
          -n {params.name} \
          -B --SPMR \
          -q {params.q} \
          {params.extra} \
          --outdir {params.outdir} \
          &> {log}
        """

rule filter_high_confidence_peaks:
    input:
        peaks="macs2/nocontrol/{sample}/{sample}_peaks.narrowPeak",
        signal="aligned_data/{sample}.final.bam",
        background=bg_bam_of,
        rnaseh=rnaseh_bam_of
    output:
        bed="high_confidence_peaks/{sample}.highconf.bed",
        metrics="high_confidence_peaks/{sample}.highconf.metrics.tsv"
    params:
        min_log2fc_bg=config["filter"].get("min_log2fc_vs_background", 1.0),
        min_log2fc_rnaseh=config["filter"].get("min_log2fc_vs_rnaseh", 1.0),
        min_signal_cpm=config["filter"].get("min_signal_cpm", 0.0),
        pseudocount=config["filter"].get("pseudocount", 1.0)
    log:
        "logs/{sample}.filter_high_confidence.log"
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p high_confidence_peaks logs
        python scripts/filter_high_confidence_peaks.py \
          --peaks {input.peaks} \
          --signal-bam {input.signal} \
          --background-bam {input.background} \
          --rnaseh-bam {input.rnaseh} \
          --out-bed {output.bed} \
          --out-metrics {output.metrics} \
          --min-log2fc-background {params.min_log2fc_bg} \
          --min-log2fc-rnaseh {params.min_log2fc_rnaseh} \
          --min-signal-cpm {params.min_signal_cpm} \
          --pseudocount {params.pseudocount} \
          &> {log}
        """

def consensus_inputs(wildcards):
    samples = replicate_data.get(str(wildcards.group), [])
    return expand("high_confidence_peaks/{sample}.highconf.bed", sample=samples)

rule consensus_replicates:
    input:
        consensus_inputs
    output:
        "consensus_peaks/{group}.reproducible.bed"
    params:
        min_reps=config.get("min_reps_for_consensus", 2)
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p consensus_peaks
        bedtools multiinter -i {input} \
          | awk -v min_reps={params.min_reps} 'BEGIN{{OFS="\t"}} $4>=min_reps {{print $1,$2,$3}}' \
          | sort -k1,1 -k2,2n \
          | bedtools merge -i - \
          > {output}
        """

rule highconf_universe:
    input:
        expand("high_confidence_peaks/{sample}.highconf.bed", sample=HIGHCONF_SAMPLES)
    output:
        "consensus_peaks/Rloop_highconf_universe.bed"
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p consensus_peaks
        cat {input} \
          | cut -f1-3 \
          | sort -k1,1 -k2,2n \
          | bedtools merge -i - \
          > {output}
        """

rule count_universe:
    input:
        bed="consensus_peaks/Rloop_highconf_universe.bed",
        bams=expand("aligned_data/{sample}.final.bam", sample=PEAK_SAMPLES)
    output:
        "counts/Rloop_highconf_universe_counts.tsv"
    params:
        samples=",".join(PEAK_SAMPLES)
    conda:
        "envs/rloop_cuttag.yml"
    shell:
        """
        mkdir -p counts
        bedtools multicov -bams {input.bams} -bed {input.bed} \
          | python scripts/add_multicov_header.py \
                --samples {params.samples} \
                --out {output}
        """
