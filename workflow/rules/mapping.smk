rule map_reads:
    input:
        r1="results/merged/{sample}_R1.fastq.gz",
        r2="results/merged/{sample}_R2.fastq.gz",
        idx=rules.bwameth_index.output,
        ref="resources/genome.fasta",
    output:
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai"
    log:
        "logs/bwameth/{sample}.log"
    conda:
        "../envs/bwameth.yaml"
    params:
        extra=get_read_group,
    threads: 100
    shell:
        "(bwameth.py {params.extra} -t {threads} --reference {input.ref} {input.r1} {input.r2} | samblaster -M | sambamba view -S -f bam /dev/stdin | sambamba sort /dev/stdin -t {threads} -m 100G -o {output.bam}) 2> {log}"

# rule mark_duplicates:
#     input:
#         "results/mapped/{sample}.sorted.bam"
#     output:
#         bam=temp("results/dedup/{sample}.sorted.bam"),
#         metrics="results/qc/dedup/{sample}.metrics.txt"
#     log:
#         "logs/picard/dedup/{sample}.log"
#     params:
#         config["params"]["picard"]["MarkDuplicates"]
#     wrapper:
#         "0.59.2/bio/picard/markduplicates"

