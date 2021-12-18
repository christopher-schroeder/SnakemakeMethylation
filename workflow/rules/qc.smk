rule bamqc:
    input:
        bam="results/mapped/{sample}.bam",
        ref="resources/genome.fasta",
    output:
        "results/qc/bamqc/{sample}/qualimapReport.html"
    params:
        outdir="results/qc/bamqc/{sample}"
    log:
        "logs/bwameth/{sample}.log"
    conda:
        "../envs/bamqc.yaml"
    threads: 4
    shell:
        "(qualimap --java-mem-size=30G bamqc -bam {input.bam} -c -nw 400 -hm 3 -ip -nt {threads} -outdir {params.outdir}) > {log} 2>&1"


rule chh_bias:
    input:
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
        ref="resources/genome.fasta",
    output:
        bed="results/qc/chh_bias/{sample}.chh_bias.bed"
    params:
        outdir="results/qc/chh_bias/{sample}"
    log:
        "logs/chh_bias/{sample}.log"
    conda:
        "../envs/methyldackel.yaml"
    threads: 4
    shell:
        "(MethylDackel mbias -@ {threads} --CHH --noCpG --noSVG {input.ref} {input.bam} {params.outdir} > {output}) 2> {log}"


rule conversion_rate:
    input:
        bed="results/qc/chh_bias/{sample}.chh_bias.bed"
    output:
        tsv="results/qc/conversion_rate/{sample}.tsv"
    log:
        "logs/conversion_rate/{sample}.log"
    shell:
        """(awk '{{if(NR>1) {{M+=$4; UM+=$5}}}}END{{printf("{wildcards.sample}\\t%f\\n", 100*(1.0-M/(M+UM)))}}' {input.bed} > {output.tsv}) 2> {log}"""


