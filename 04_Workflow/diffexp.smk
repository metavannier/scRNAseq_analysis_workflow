SAMPLE = config["sample"]["sname"]

rule diffexp:
    input:
        seurat_object = OUTPUTDIR + "01_seurat/" + SAMPLE + "_seurat_object.rds",
    output:
        diffexp_report = report(OUTPUTDIR + "02_diffexp/" + SAMPLE + "_diffexp_report.html", caption = ROOTDIR + REPORT + "diffexp.rst", category="02 diffexpression"),
        violinplot = report(OUTPUTDIR + "02_diffexp/violin_plot.pdf", caption = ROOTDIR + REPORT + "violin.rst", category="02 diffexpression"),
        umapfeature = report(OUTPUTDIR + "02_diffexp/umapfeature_plot.pdf", caption = ROOTDIR + REPORT + "umapfeature.rst", category="02 diffexpression"),
        tsnefeature = report(OUTPUTDIR + "02_diffexp/tsnefeature_plot.pdf", caption = ROOTDIR + REPORT + "tsnefeature.rst", category="02 diffexpression"),
        ridgefeature = report(OUTPUTDIR + "02_diffexp/ridgefeature_plot.pdf", caption = ROOTDIR + REPORT + "ridgefeature.rst", category="02 diffexpression"),
        heatmapfeature = report(OUTPUTDIR + "02_diffexp/heatmapfeature.pdf", caption = ROOTDIR + REPORT + "heatmapfeature.rst", category="02 diffexpression"),
    conda:
        CONTAINER + "seurat.yaml"
    params:
        sample = config["sample"]["sname"],
        features = config["diffexp"]["features"],
    message: 
        "Run Seurat for violin plot of expression distributions across clusters"
    script:
        SCRIPTDIR + "diffexp_reports_compilation.R"