filename = snakemake.input[0]
outfile = snakemake.output[0]

# filename = "results/dmr/metilene/dmrs/ko_vs_upmm2wt.annotated.bed"
# outfile = "/dev/stdout"

i=0
with open(filename, "r") as f, open(outfile, "w") as o:
    print("i", "chrom", "start", "end", "q", "diff", "p_mwu", "p_2d_ks", "mean_case", "mean_control", "effect", "transcript", "gene", "gene_type", sep="\t", file=o)
    for line in f:
        if line.startswith("#"):
            continue
        i+=1
        chrom, start, end, data = line.strip().split("\t")
        q, diff, n, p_mwu, p_2d_ks, mean_case, mean_control = data.split(";")[0:7]
        annotation = data.split(";")[7:]
        for a in annotation:
            a_split = a.split("|")
            if len(a_split) == 3:
                effect, transcript, gene = a.split("|")
            else:
                effect = ""
                transcript = ""
                gene = ""
                for field in a_split:
                    if field.startswith("Transcript"):
                        transcript = field
                    elif field.startswith("Gene"):
                        gene = field
                    else:
                        effect = field

            gene = gene[len("Gene:"):].split(":")
            transcript = transcript[len("Transcript:"):]

            print(i, f"chr{chrom}", start, end, q, diff, p_mwu, p_2d_ks, mean_case, mean_control, effect, transcript, *gene, sep="\t", file=o)
            break # only first transcript

