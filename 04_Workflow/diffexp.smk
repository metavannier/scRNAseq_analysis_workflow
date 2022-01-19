NPROJ = config["sample"]["nproject"]
FEATURES = config["diffexp"]["features"].split(',')
CELLMARKER = config["diffexp"]["cell_marker"].split(',')
WT = config["seurat"]["wt"]

rule diffexp:
    input:
        seurat_object = OUTPUTDIR + "02_seurat/" + NPROJ + "_seurat_object.rds",
    output:
        diffexp_report = report(OUTPUTDIR + "03_diffexp/" + NPROJ + "_diffexp_report.html", caption = ROOTDIR + REPORT + "diffexp.rst", category="03 diffexpression"),
        violinplot = report(expand(OUTPUTDIR + "03_diffexp/violin_plot/{features}_violin_plot.pdf", features=FEATURES), caption = ROOTDIR + REPORT + "violin.rst", category="03 diffexpression"),
        umapfeature = report(expand(OUTPUTDIR + "03_diffexp/umap_plot/{features}_umapfeature_plot.pdf", features=FEATURES), caption = ROOTDIR + REPORT + "umapfeature.rst", category="03 diffexpression"),
        tsnefeature = report(expand(OUTPUTDIR + "03_diffexp/tsne_plot/{features}_tsnefeature_plot.pdf", features=FEATURES), caption = ROOTDIR + REPORT + "tsnefeature.rst", category="03 diffexpression"),
        ridgefeature = report(expand(OUTPUTDIR + "03_diffexp/ridge_plot/{features}_ridgefeature_plot.pdf", features=FEATURES), caption = ROOTDIR + REPORT + "ridgefeature.rst", category="03 diffexpression"),
        heatmapfeature = report(OUTPUTDIR + "03_diffexp/heatmap/heatmapfeature.pdf", caption = ROOTDIR + REPORT + "heatmapfeature.rst", category="03 diffexpression"),
        defile_allcells = report(expand(OUTPUTDIR + "03_diffexp/differential_expression_features/{wt}_vs_AllCellMarker_DE.csv", wt=WT), caption = ROOTDIR + REPORT + "de_tests.rst", category="03 diffexpression"),
        defile = report(expand(OUTPUTDIR + "03_diffexp/differential_expression_features/{wt}_vs_{cellmarker}_DE.csv", wt=WT, cellmarker=CELLMARKER), caption = ROOTDIR + REPORT + "de_tests.rst", category="03 diffexpression"),
        volcano_allcluster = report(expand(OUTPUTDIR + "03_diffexp/volcanoplot/volcano_{wt}_vs_AllCellMarker.pdf", wt=WT), caption = ROOTDIR + REPORT + "volcanoplot.rst", category="03 diffexpression"),
        volcano = report(expand(OUTPUTDIR + "03_diffexp/volcanoplot/volcano_{wt}_vs_{cellmarker}.pdf", wt=WT, cellmarker=CELLMARKER), caption = ROOTDIR + REPORT + "volcanoplot.rst", category="03 diffexpression"),
        sign_up_allcell = report(expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_AllCellMarker_signif-up-regulated.txt", wt=WT), caption = ROOTDIR + REPORT + "de_significatif.rst", category="03 diffexpression"),
        sign_up = report(expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_{cellmarker}_signif-up-regulated.txt", wt=WT, cellmarker=CELLMARKER), caption = ROOTDIR + REPORT + "de_significatif.rst", category="03 diffexpression"),
        sign_down_allcell = report(expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_AllCellMarker_signif-down-regulated.txt", wt=WT), caption = ROOTDIR + REPORT + "de_significatif.rst", category="03 diffexpression"),
        sign_down = report(expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_{cellmarker}_signif-down-regulated.txt", wt=WT, cellmarker=CELLMARKER), caption = ROOTDIR + REPORT + "de_significatif.rst", category="03 diffexpression"),
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