rule seurat:
    input:
        sc_data = OUTPUTDIR + "00_cellranger/Mix_MM_lines/outs/filtered_feature_bc_matrix/"
    output:
        seurat_report = report(OUTPUTDIR + "01_seurat/seurat.html", caption = ROOTDIR + REPORT + "seurat.rst", category="01 seurat")
    conda:
        CONTAINER + "seurat.yaml"
    params:
        sample = config["sample"]["sname"],
    message: 
        "Run Seurat for the clustering"
    script:
        SCRIPTDIR + "seurat_reports_compilation.R"