NPROJ = config["sample"]["nproject"]

rule seurat:
    input:
        sc_data = OUTPUTDIR + "01_cellranger/" + NPROJ + "/outs/filtered_feature_bc_matrix/",
        aggrcsv = ROOTDIR + "/aggregation.csv",
        # tsne = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/analysis/tsne/2_components/projection.csv",
        # demuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.best"
    output:
        #tabdemuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.tsv",
        seurat_report = report(OUTPUTDIR + "02_seurat/" + NPROJ + "_seurat_report.html", caption = ROOTDIR + REPORT + "seurat.rst", category="02 seurat"),
        seurat_object = OUTPUTDIR + "02_seurat/" + NPROJ + "_seurat_object.rds",
        count_matrix = report(OUTPUTDIR + "02_seurat/" + NPROJ + "_count_matrix.csv", caption = ROOTDIR + REPORT + "data_matrix.rst", category="02 seurat"),
        data_matrix = report(OUTPUTDIR + "02_seurat/" + NPROJ + "_data_matrix.csv", caption = ROOTDIR + REPORT + "data_matrix.rst", category="02 seurat"),
        scale_data_matrix = report(OUTPUTDIR + "02_seurat/" + NPROJ + "_scale_data_matrix.csv", caption = ROOTDIR + REPORT + "data_matrix.rst", category="02 seurat"),
    conda:
        CONTAINER + "seurat.yaml"
    params:
        nproj = config["sample"]["nproject"],
        samples = config["cellranger"]["sagrr"],
        wt = config["seurat"]["wt"],
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
        cluster_ids = config["seurat"]["cluster_ids"],
    message: 
        "Run Seurat for the clustering"
    script:
        SCRIPTDIR + "seurat_reports_compilation.R"