#----------------------------------------
# Seurat
#----------------------------------------

rule seurat:
    input:
        cellranger_output = expand(OUTPUTDIR + "01_cellranger/cellranger_output.txt"),

    output:
        seurat_output = expand(OUTPUTDIR + "02_seurat/seurat_output.txt"),
        seurat_report = report(expand(OUTPUTDIR + "02_seurat/{sample_id}/{sample_id}_seurat_report.html", sample_id = SAMPLE_ID), caption = REPORT + "seurat.rst", category = "02 seurat"),

    conda:
        CONTAINER + "seurat.yaml"

    params:
        run_demultiplex = config["run_demultiplex"],
        sc_data = expand(OUTPUTDIR + "01_cellranger/{sample_id}/count/sample_filtered_feature_bc_matrix/",sample_id = SAMPLE_ID),
        cell_ranger_count_path = config["seurat"]["cell_ranger_count_path"],
        sample_id = expand("{sample_id.id}", sample_id = sample_id.itertuples()),
        multiplex = expand("{multiplex.multiplex}", multiplex = multiplex.itertuples()),
        hto = expand("{hto.hto}", hto = hto.itertuples()),
        plot_raster_nbcells_threshold = config["seurat"]["plot_raster_nbcells_threshold"],
        qc_exploration_mode = config["seurat"]["qc_exploration_mode"],
        min_cells = config["seurat"]["min_cells"],
        min_features = config["seurat"]["min_features"],
        filter_umi_min = config["seurat"]["filter_umi_min"].split(','),
        filter_umi_max = config["seurat"]["filter_umi_max"].split(','),
        filter_feature_min = config["seurat"]["filter_feature_min"].split(','),
        filter_feature_max = config["seurat"]["filter_feature_max"].split(','),
        filter_percent_mt = config["seurat"]["filter_percent_mt"].split(','),
        filter_percent_mt_min = config["seurat"]["filter_percent_mt_min"].split(','),
        filter_percent_rb = config["seurat"]["filter_percent_rb"].split(','),
        pattern_mt = config["seurat"]["pattern_mt"],
        pattern_rb = config["seurat"]["pattern_rb"],
        norm_method = config["seurat"]["norm_method"],
        norm_scale_factor = config["seurat"]["norm_scale_factor"],
        feature_select_method = config["seurat"]["feature_select_method"],
        variable_features = config["seurat"]["variable_features"],
        variable_features_showtop = config["seurat"]["variable_features_showtop"],
        dims = config["seurat"]["dims"],
        pca_npc = config["seurat"]["pca_npc"],
        pca_plots_nbdims = config["seurat"]["pca_plots_nbdims"],
        pca_plot_nbfeatures = config["seurat"]["pca_plot_nbfeatures"],
        dimreduc_use_pca_nbdims = config["seurat"]["dimreduc_use_pca_nbdims"],
        findclusters_use_pca_nbdims = config["seurat"]["findclusters_use_pca_nbdims"],
        findneighbors_k = config["seurat"]["findneighbors_k"],
        findclusters_resolution = config["seurat"]["findclusters_resolution"],
        findclusters_algorithm = config["seurat"]["findclusters_algorithm"],

    message:
        "Run Seurat for the clustering"

    script:
        SCRIPTDIR + "seurat_reports_compilation.R"

#----------------------------------------
# Label Transfert with Seurat
#----------------------------------------

rule seurat_labelTransfert:
    input:
        seurat_output = expand(OUTPUTDIR + "02_seurat/seurat_output.txt"),

    output:
        seurat_labelTransfert_output = expand(OUTPUTDIR + "02_seurat/seurat_labelTransfert_output.txt"),

    conda:
        CONTAINER + "seurat.yaml"

    params:
        sample_id = expand("{sample_id.id}", sample_id = sample_id.itertuples()),

    message:
        "Run Seurat for the transfert data"

    script:
        SCRIPTDIR + "labelTransfert_seurat.R"

#----------------------------------------
# DEG analysis with Seurat
#----------------------------------------

rule seurat_DEG:
    input:
        seurat_labelTransfert_output = expand(OUTPUTDIR + "02_seurat/seurat_labelTransfert_output.txt"),

    output:
        seurat_DEG_output = expand(OUTPUTDIR + "04_diffexp/seurat_DEG_output.txt"),

    conda:
        CONTAINER + "seurat.yaml"

    params:
        time_point = config["diffexp"]["time_point"],

    message:
        "Run Seurat for the DE analysis"

    script:
        SCRIPTDIR + "DEG_analysis.R"
