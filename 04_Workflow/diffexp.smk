NPROJ = config["sample"]["nproject"]

rule diffexp:
    input:
        seurat_object = OUTPUTDIR + "02_seurat/" + NPROJ + "_seurat_object.rds",
    output:
        diffexp_report = report(OUTPUTDIR + "03_diffexp/" + NPROJ + "_diffexp_report.html", caption = ROOTDIR + REPORT + "diffexp.rst", category="03 diffexpression"),
        violinplot = report(OUTPUTDIR + "03_diffexp/violin_plot.pdf", caption = ROOTDIR + REPORT + "violin.rst", category="03 diffexpression"),
        umapfeature = report(OUTPUTDIR + "03_diffexp/umapfeature_plot.pdf", caption = ROOTDIR + REPORT + "umapfeature.rst", category="03 diffexpression"),
        tsnefeature = report(OUTPUTDIR + "03_diffexp/tsnefeature_plot.pdf", caption = ROOTDIR + REPORT + "tsnefeature.rst", category="03 diffexpression"),
        ridgefeature = report(OUTPUTDIR + "03_diffexp/ridgefeature_plot.pdf", caption = ROOTDIR + REPORT + "ridgefeature.rst", category="03 diffexpression"),
        heatmapfeature = report(OUTPUTDIR + "03_diffexp/heatmapfeature.pdf", caption = ROOTDIR + REPORT + "heatmapfeature.rst", category="03 diffexpression"),
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