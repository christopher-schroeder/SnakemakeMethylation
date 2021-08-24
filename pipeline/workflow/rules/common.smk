import glob
import yaml
import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

samples = pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str, "group": str}).set_index("sample_name", drop=False).sort_index()

with open(config["experiments"]) as file:
    experiments = yaml.load(file, Loader=yaml.FullLoader)


def get_case(wildcards):
    return experiments[wildcards.experiment]["case"]


def get_control(wildcards):
    return experiments[wildcards.experiment]["control"]


def metilene_header(wildcards):
    case = get_case(wildcards)
    control = get_control(wildcards)
    return "\\\\t".join(["chr", "pos"] + list(map(lambda x: f"case_{x}", case)) + list(map(lambda x: f"control_{x}", control)))


def get_final_output():
    final_output = []

    if config["report"]["activate"]:
        final_output.extend(expand("results/vcf-report/all.{event}/",
                            event=config["calling"]["fdr-control"]["events"]))
    else:
        final_output.extend(expand("results/merged-calls/{group}.{event}.fdr-controlled.bcf",
                            group=groups,
                            event=config["calling"]["fdr-control"]["events"]))

    if config["tables"]["activate"]:
        final_output.extend(expand("results/tables/{group}.{event}.fdr-controlled.tsv",
                            group=f,
                            event=config["calling"]["fdr-control"]["events"]))
        if config["tables"].get("generate_excel", False):
            final_output.extend(expand("results/tables/{group}.{event}.fdr-controlled.xlsx",
                            group=groups,
                            event=config["calling"]["fdr-control"]["events"]))
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

wildcard_constraints:
    sample="|".join(samples["sample_name"]),
    experiment="|".join(experiments),
