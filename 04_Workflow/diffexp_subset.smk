NPROJ = config["sample"]["nproject"]
FEATURES = config["diffexp"]["features"].split(',')
CELLMARKER = config["diffexp"]["cell_marker"].split(',')
WT = config["seurat"]["wt"]
MARKERGENE = config["diffexpsubset"]["marker_gene"]

rule diffexp_subset:
    input:
        seurat_object = OUTPUTDIR + "02_seurat/" + NPROJ + "_seurat_object.rds",
    output:
        diffexp_subset_report = report(OUTPUTDIR + "04_diffexp_subset/" + NPROJ + "_diffexp_subset_report.html", caption = REPORT + "diffexp.rst", category="04 differential expression after subset", subcategory="Report on the results of the differential expression analysis after subset"),
        defile = expand(OUTPUTDIR + "04_diffexp_subset/differential_expression/{cluster}_subcluster_DE.csv", cluster=CLUSTER),
        volcano = expand(OUTPUTDIR + "04_diffexp_subset/volcanoplot/volcano_{cluster}.pdf", cluster=CLUSTER),
        sign_up = expand(OUTPUTDIR + "04_diffexp_subset/de_significant/{cluster}_signif-up-regulated.txt", cluster=CLUSTER),
        sign_down = expand(OUTPUTDIR + "04_diffexp_subset/de_significant/{cluster}_signif-down-regulated.txt", cluster=CLUSTER),
        violinplot = expand(OUTPUTDIR + "04_diffexp_subset/violin_plot/{features}_violin_plot.pdf", features=FEATURES),
    conda:
        CONTAINER + "seurat.yaml"
    params:
        markergene = config["diffexpsubset"]["marker_gene"],
        genetoremove = config["diffexpsubset"]["genetoremove"],
        nproj = config["sample"]["nproject"],
        sample = config["sample"]["sname"],
        features = config["diffexp"]["features"],
        wt = config["seurat"]["wt"],
        cluster = config["diffexpsubset"]["cluster"],
        test = config["diffexp"]["test"],
        threshold = config["diffexp"]["threshold"],
        FCcutoff = config["diffexp"]["FCcutoff"],
        pCutoff = config["diffexp"]["pCutoff"],
    message: 
        "Run Seurat for differential genes expression inside clusters"
    script:
        SCRIPTDIR + "diffexp_subset_reports_compilation.R"