rule plot:
    input:
        dmr="results/dmr/metilene/dmrs/{experiment}.bed",
        mvalues="results/dmr/metilene/input/{experiment}.bg"
    output:
        outdir=directory("results/dmr/metilene/plots/{experiment}")
    script:
        "../scripts/dendogram.py"