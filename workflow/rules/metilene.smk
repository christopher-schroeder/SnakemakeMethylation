rule methylation_to_bed:
    input:
        "results/methylation/{sample}_CpG.bedGraph"
    output:
        "results/dmr/metilene/bed/{sample}.bg"
    params:
        minimum_coverage=5,
    shell:
        """cat {input} | awk 'BEGIN {{ FS = "\\t"; OFS=FS }} ; {{ if (($5+$6)>={params.minimum_coverage}) {{ print $1,$2,$3,$4/100; }} }}' > {output}"""


rule generate_metilene_input:
    input:
        case=lambda wc: expand("results/dmr/metilene/bed/{sample}.bg", sample=get_case(wc)),
        control=lambda wc: expand("results/dmr/metilene/bed/{sample}.bg", sample=get_control(wc)),
    output:
        "results/dmr/metilene/input/{experiment}.bg"
    params:
        header=metilene_header
    conda:
        "../envs/bedtools.yaml"
    shell:
        """(echo -e {params.header} & bedtools unionbedg -i {input.case} {input.control} -filler "." | cut -f 1,2,4-100) > {output}"""


rule call_dmr_metilene:
    threads:
        8
    input:
        "results/dmr/metilene/input/{experiment}.bg"
    output:
        "results/dmr/metilene/dmrs_unfiltered/{experiment}.tsv"
    log:
        "logs/dmr/metilene/call/{experiment}.log"
    conda:
        "../envs/metilene.yaml"
    params:
        minimum_diff=config["dmrs"]["metilene"]["min_diff"],
        minimum_cgps=config["dmrs"]["metilene"]["min_cpg"],
        min_n_group_A = 2,
        min_n_group_B = 2,
    shell:
        "metilene -t {threads} -X {params.min_n_group_A} -Y {params.min_n_group_B} --mincpgs {params.minimum_cgps} --minMethDiff {params.minimum_diff} -a case_ -b control_ {input} > {output}"


rule filter_pvalue_metilene:
    input:
        "results/dmr/metilene/dmrs_unfiltered/{experiment}.tsv"
    output:
        "results/dmr/metilene/dmrs/{experiment}.tsv"
    shell:
        """cat {input} | awk 'BEGIN {{FS="\\t"; OFS=FS}} {{if ($4<0.05) {{print $0}}}}' > {output}"""


rule dmr_to_bed_metilene:
    input:
        "results/dmr/metilene/dmrs/{experiment}.tsv"
    output:
        "results/dmr/metilene/dmrs/{experiment}.bed"
    shell:
        """cat {input} | awk 'BEGIN {{FS="\\t"; OFS=FS;}} {{j=$4; for(i=5; i <= NF; i++) {{j = j"\;"$(i);}} print $1,$2,$3,j}}' > {output}"""


rule bed_metilene_to_table:
    input:
        "results/dmr/metilene/dmrs/{experiment}.annotated.bed"
    output:
        "results/dmr/metilene/dmrs/{experiment}.table.tsv"
    script:
        "../scripts/bed_to_table.py"