import glob
import yaml
import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

samples = pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str, "group": str}).set_index("sample_name", drop=False).sort_index()

experiments = config["dmrs"]["experiments"]

groups = samples["group"].unique()

def get_case(wildcards):
    #return experiments[wildcards.experiment]["case"]
    case = config["dmrs"]["experiments"][wildcards.experiment]["case"]
    case_samples = samples[samples.group == case]["sample_name"]
    return case_samples
    

def get_control(wildcards):
    #return experiments[wildcards.experiment]["control"]
    control = config["dmrs"]["experiments"][wildcards.experiment]["control"]
    control_samples = samples[samples.group == control]["sample_name"]
    return control_samples


def metilene_header(wildcards):
    case = get_case(wildcards)
    control = get_control(wildcards)
    return "\\\\t".join(["chr", "pos"] + list(map(lambda x: f"case_{x}", case)) + list(map(lambda x: f"control_{x}", control)))


def get_final_output():
    final_output = expand("results/mapped/{sample}.bam",
        sample=samples["sample_name"])

    if config["qc"]["activate"]:
        final_output.extend(expand(["results/qc/conversion_rate/{sample}.tsv",
            "results/qc/bamqc/{sample}/qualimapReport.html"], 
            sample=samples["sample_name"]))
    if config["meth"]["activate"]:
        final_output.extend(expand("results/methylation/{sample}_CpG.bedGraph", 
            sample=samples["sample_name"]))
    if config["dmrs"]["metilene"]:
        final_output.extend(expand("results/dmr/metilene/dmrs/{experiment}.table.tsv",
            experiment=config["dmrs"]["experiments"]))
    if config["pca"]["activate"]:
        final_output.append("results/plots/pca.pdf")

    return final_output

def _group_or_sample(row):
    group = row.get("group", None)
    if pd.isnull(group):
        return row["sample_name"]
    return group

samples["group"] = [_group_or_sample(row) for _, row in samples.iterrows()]
units = pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str}).set_index(["sample_name", "unit_name"], drop=False).sort_index()


def pca_params(wildcards):
    filenames = " ".join(["--group " + " ".join(expand("results/methylation/{sample}_CpG.bedGraph", sample=samples["sample_name"][samples["group"]==g])) for g in groups])
    symbols = " ".join([" " + " ".join(samples["pca_symbol"][samples["group"]==g]) for g in groups])
    return f"{filenames} --symbols {symbols}"
        


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"--read-group '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=samples.loc[wildcards.sample, "platform"])

def get_read_group_biscuit(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=samples.loc[wildcards.sample, "platform"])

wildcard_constraints:
    sample="|".join(samples["sample_name"]),
    experiment="|".join(experiments),
