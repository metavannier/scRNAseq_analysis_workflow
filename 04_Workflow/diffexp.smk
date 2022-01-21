NPROJ = config["sample"]["nproject"]
FEATURES = config["diffexp"]["features"].split(',')
CELLMARKER = config["diffexp"]["cell_marker"].split(',')
WT = config["seurat"]["wt"]

rule diffexp:
    input:
        seurat_object = OUTPUTDIR + "02_seurat/" + NPROJ + "_seurat_object.rds",
    output:
        diffexp_report = report(OUTPUTDIR + "03_diffexp/" + NPROJ + "_diffexp_report.html", caption = REPORT + "diffexp.rst", category="03 differential expression", subcategory="Report on the results of the differential expression analysis"),
        violinplot = expand(OUTPUTDIR + "03_diffexp/violin_plot/{features}_violin_plot.pdf", features=FEATURES),
        umapfeature = expand(OUTPUTDIR + "03_diffexp/umap_plot/{features}_umapfeature_plot.pdf", features=FEATURES),
        tsnefeature = expand(OUTPUTDIR + "03_diffexp/tsne_plot/{features}_tsnefeature_plot.pdf", features=FEATURES),
        ridgefeature = expand(OUTPUTDIR + "03_diffexp/ridge_plot/{features}_ridgefeature_plot.pdf", features=FEATURES),
        heatmapfeature = OUTPUTDIR + "03_diffexp/heatmap/heatmapfeature.pdf",
        defile_allcells = expand(OUTPUTDIR + "03_diffexp/differential_expression_features/{wt}_vs_AllCellMarker_DE.csv", wt=WT),
        defile = expand(OUTPUTDIR + "03_diffexp/differential_expression_features/{wt}_vs_{cellmarker}_DE.csv", wt=WT, cellmarker=CELLMARKER),
        volcano_allcluster = expand(OUTPUTDIR + "03_diffexp/volcanoplot/volcano_{wt}_vs_AllCellMarker.pdf", wt=WT),
        volcano = expand(OUTPUTDIR + "03_diffexp/volcanoplot/volcano_{wt}_vs_{cellmarker}.pdf", wt=WT, cellmarker=CELLMARKER),
        sign_up_allcell = expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_AllCellMarker_signif-up-regulated.txt", wt=WT),
        sign_up = expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_{cellmarker}_signif-up-regulated.txt", wt=WT, cellmarker=CELLMARKER),
        sign_down_allcell = expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_AllCellMarker_signif-down-regulated.txt", wt=WT),
        sign_down = expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_{cellmarker}_signif-down-regulated.txt", wt=WT, cellmarker=CELLMARKER),
    conda:
        CONTAINER + "seurat.yaml"
    params:
        nproj = config["sample"]["nproject"],
        sample = config["sample"]["sname"],
        features = config["diffexp"]["features"],
        wt = config["seurat"]["wt"],
        cellmarker = config["diffexp"]["cell_marker"],
        test = config["diffexp"]["test"],
        threshold = config["diffexp"]["threshold"],
        FCcutoff = config["diffexp"]["FCcutoff"],
        pCutoff = config["diffexp"]["pCutoff"],
    message: 
        "Run Seurat for differential genes expression across clusters"
    script:
        SCRIPTDIR + "diffexp_reports_compilation.R"