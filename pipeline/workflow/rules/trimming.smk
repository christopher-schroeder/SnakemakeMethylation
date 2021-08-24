rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq"
    log:
        "logs/get-sra/{accession}.log"
    wrapper:
        "0.56.0/bio/sra-tools/fasterq-dump"


def get_cutadapt_pipe_input(wildcards):
    pattern = units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]
    if "*" in pattern:
        files = sorted(glob.glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]))
        if not files:
            raise ValueError(
                "No raw fastq files found for unit pattern {} (sample {}). "
                "Please check the your sample sheet.".format(wildcards.unit, wildcards.sample)
            )
    else:
        files = [pattern]

    return files


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input
    output:
        pipe('pipe/cutadapt/{sample}/{unit}.{fq}.{ext}')
    log:
        "logs/pipe-fastqs/catadapt/{sample}-{unit}.{fq}.{ext}.log"
    wildcard_constraints:
        ext=r"fastq|fastq\.gz"
    threads: 0 # this does not need CPU
    shell:
        "cat {input} > {output} 2> {log}"


def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    # only paired end supported
    return expand("pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(S=unit.sample_name, U=unit.unit_name, E=ending), read=["fq1","fq2"])


def get_cutadapt_adapters(wildcards):
    adapters = units.loc[wildcards.sample].loc[wildcards.unit, "adapters"]
    if isinstance(adapters, str):
        return adapters
    return ""


rule cutadapt_pe:
    input:
        get_cutadapt_input
    output:
        fastq1="results/trimmed/{sample}/{unit}_R1.fastq.gz",
        fastq2="results/trimmed/{sample}/{unit}_R2.fastq.gz",
        qc="results/trimmed/{sample}/{unit}.paired.qc.txt"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    params:
        others = config["params"]["cutadapt"],
        adapters = get_cutadapt_adapters,
    threads: 8
    wrapper:
        "0.59.2/bio/cutadapt/pe"


def get_fastqs(wc):
    return expand("results/trimmed/{sample}/{unit}_{read}.fastq.gz", unit=units.loc[wc.sample, "unit_name"], sample=wc.sample, read=wc.read)


rule merge_fastqs:
    input:
        get_fastqs
    output:
        "results/merged/{sample}_{read}.fastq.gz"
    log:
        "logs/merge-fastqs/{sample}_{read}.log"
    wildcard_constraints:
        read="single|R1|R2"
    shell:
        "cat {input} > {output} 2> {log}"
