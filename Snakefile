rule all:
    input:
        "plots/cms_ditau_reconstruction.png",
        "plots/roofit_mvis_fit.png"

rule reconstruct_ditau:
    input:
        config="configs/workflow.yaml"
    output:
        root="plots/ditau_histograms.root",
        plot="plots/cms_ditau_reconstruction.png"
    shell:
        "python3 src/reconstruct_ditau.py"

rule roofit_mass_fit:
    input:
        root="plots/ditau_histograms.root"
    output:
        plot="plots/roofit_mvis_fit.png"
    shell:
        "python3 src/fit_ditau_mass.py"
