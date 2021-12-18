rule get_genome:
    output:
        "resources/genome.fasta"
    log:
        "logs/get-genome.log"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    cache: True
    wrapper:
        "0.59.2/bio/reference/ensembl-sequence"


rule bwameth_index:
    input:
        ref="resources/genome.fasta"
    output:
        multiext("resources/genome.fasta", ".bwameth.c2t", ".bwameth.c2t.amb", ".bwameth.c2t.ann", ".bwameth.c2t.bwt", ".bwameth.c2t.pac", ".bwameth.c2t.sa")
    log:
        "logs/bwameth_index.log"
    resources:
        mem_mb=369000
    cache: True
    conda:
        "../envs/bwameth.yaml"
    shell:
        "bwameth.py index {input.ref}"
