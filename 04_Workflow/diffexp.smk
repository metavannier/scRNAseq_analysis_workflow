NPROJ = config["sample"]["nproject"]
FEATURES = config["diffexp"]["features"].split(',')

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
    conda:
        CONTAINER + "seurat.yaml"
    params:
        nproj = config["sample"]["nproject"],
        sample = config["sample"]["sname"],
        features = config["diffexp"]["features"],
    message: 
        "Run Seurat for violin plot of expression distributions across clusters"
    script:
        SCRIPTDIR + "diffexp_reports_compilation.R"