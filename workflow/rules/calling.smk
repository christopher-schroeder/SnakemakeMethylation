rule methyldackel:
    input:
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
        ref="resources/genome.fasta",
    output:
        "results/methylation/{sample}_CpG.bedGraph",
    params:
        prefix="results/methylation/{sample}"
    log:
        "logs/methyldackel/{sample}.log"
    conda:
        "../envs/methyldackel.yaml"
    threads: 4
    shell:
        "(MethylDackel extract {input.ref} {input.bam} -@ {threads} --mergeContext -o {params.prefix}) 2> {log}"
