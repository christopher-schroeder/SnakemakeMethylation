rule snpeff:
    input:
        calls="results/dmr/metilene/dmrs/{experiment}.bed",
        db="resources/snpeff/hg38"
    output:
        calls="results/dmr/metilene/dmrs/{experiment}.annotated.bed",
        stats="results/dmr/metilene/dmrs/snpeff_stats/{experiment}.html"
    log:
        "logs/snpeff/{experiment}.log"
    resources:
        mem_mb=4096
    params:
        extra="-i bed"
    wrapper:
        "0.77.0/bio/snpeff/annotate"


rule snpeff_download:
    output:
        directory("resources/snpeff/hg38")
    log:
        "logs/snpeff/download/hg38.log"
    params:
        reference="hg38"
    resources:
        mem_mb=1024
    wrapper:
        "0.77.0/bio/snpeff/download"