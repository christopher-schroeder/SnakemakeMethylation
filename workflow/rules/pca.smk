rule pca:
    input:
        calls=expand("results/methylation/{sample}_CpG.bedGraph", sample=samples["sample_name"])
    output:
        "results/plots/pca.pdf"
    log:
        "logs/plots/pca.log"
    params:
        pca_params
    shell:
        "python workflow/scripts/pca.py {params} -o {output}"