SAMPLE = config["sample"]["sname"]

rule seurat:
    input:
        sc_data = OUTPUTDIR + "00_cellranger/Mix_MM_lines/outs/filtered_feature_bc_matrix/",
        tsne = OUTPUTDIR + "00_cellranger/Mix_MM_lines/outs/analysis/tsne/2_components/projection.csv",
        demuxlet = OUTPUTDIR + "00_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.best"
    output:
        tabdemuxlet = OUTPUTDIR + "00_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.tsv",
        seurat_report = report(OUTPUTDIR + "01_seurat/" + SAMPLE + "_seurat_report.html", caption = ROOTDIR + REPORT + "seurat.rst", category="01 seurat"),
        #seurat_object = OUTPUTDIR + "01_seurat/" + SAMPLE + "_seurat_object.rds",
    conda:
        CONTAINER + "seurat.yaml"
    params:
        sample = config["sample"]["sname"],
        biop = config["sample"]["biop"],
        min_cells = config["seurat"]["min_cells"],
        min_features = config["seurat"]["min_features"],
        minFeature_RNA = config["seurat"]["minFeature_RNA"],
        maxFeature_RNA = config["seurat"]["maxFeature_RNA"],
        percent_mt = config["seurat"]["percent_mt"],
        normwf = config["seurat"]["normwf"],
        normalization_method = config["seurat"]["normalization_method"],
        scale_factor = config["seurat"]["scale_factor"],
        nfeatures = config["seurat"]["nfeatures"],
        selection_method = config["seurat"]["selection_method"],
        nHVG = config["seurat"]["nHVG"],
        dims =  config["seurat"]["dims"],
        resolution = config["seurat"]["resolution"],
    message: 
        "Run Seurat for the clustering"
    script:
        SCRIPTDIR + "seurat_reports_compilation.R"