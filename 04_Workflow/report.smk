rule archive_report:
    input:
        count_matrix = OUTPUTDIR + "02_seurat/data_matrix/" + NPROJ + "_count_matrix.csv",
        data_matrix = OUTPUTDIR + "02_seurat/data_matrix/" + NPROJ + "_data_matrix.csv",
        scale_data_matrix = OUTPUTDIR + "02_seurat/data_matrix/" + NPROJ + "_scale_data_matrix.csv",
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
        defile_subset = expand(OUTPUTDIR + "04_diffexp_subset/differential_expression/{cluster}_subcluster_DE.csv", cluster=CLUSTER),
        volcano_subset = expand(OUTPUTDIR + "04_diffexp_subset/volcanoplot/volcano_{cluster}.pdf", cluster=CLUSTER),
        sign_up_subset = expand(OUTPUTDIR + "04_diffexp_subset/de_significant/{cluster}_signif-up-regulated.txt", cluster=CLUSTER),
        sign_down_subset = expand(OUTPUTDIR + "04_diffexp_subset/de_significant/{cluster}_signif-down-regulated.txt", cluster=CLUSTER),
        violinplot_subset = expand(OUTPUTDIR + "04_diffexp_subset/violin_plot/{features}_violin_plot.pdf", features=FEATURES),
    output:
        data_matrix_tar = OUTPUTDIR + "02_seurat/data_matrix.tar.gz",
        violinplot_tar = report(OUTPUTDIR + "03_diffexp/violin_plot.tar.gz", caption = REPORT + "violin.rst", category="03 differential expression", subcategory="violin"),
        umapfeature_tar = report(OUTPUTDIR + "03_diffexp/umapfeature_plot.tar.gz", caption = REPORT + "umapfeature.rst", category="03 differential expression", subcategory="umap feature"),
        tsnefeature_tar = report(OUTPUTDIR + "03_diffexp/tsnefeature_plot.tar.gz", caption = REPORT + "tsnefeature.rst", category="03 differential expression", subcategory="tsne feature"),
        ridgefeature_tar = report(OUTPUTDIR + "03_diffexp/ridgefeature_plot.tar.gz", caption = REPORT + "ridgefeature.rst", category="03 differential expression", subcategory="ridge feature"),
        heatmapfeature_tar = report(OUTPUTDIR + "03_diffexp/heatmapfeature.tar.gz", caption = REPORT + "heatmapfeature.rst", category="03 differential expression", subcategory="heatmap feature"),
        defile_tar = report(OUTPUTDIR + "03_diffexp/differencial_expression_tests.tar.gz", caption = REPORT + "de_tests.rst", category="03 differential expression", subcategory="differencial expression tests"),
        volcano_tar = report(OUTPUTDIR + "03_diffexp/volcano_plot.tar.gz", caption = REPORT + "volcanoplot.rst", category="03 differential expression", subcategory="volcano plot"),
        sign_up_down_tar = report(OUTPUTDIR + "03_diffexp/up_down_regulated_genes_list.tar.gz", caption = REPORT + "de_significatif.rst", category="03 differential expression", subcategory="significant regulated genes list"),
        defile_subset_tar = report(OUTPUTDIR + "04_diffexp_subset/differencial_expression_tests.tar.gz", caption = REPORT + "de_tests.rst", category="04 differential expression after subset", subcategory="differencial expression tests"),
        volcano_subset_tar = report(OUTPUTDIR + "04_diffexp_subset/volcano_plot.tar.gz", caption = REPORT + "volcanoplot.rst", category="04 differential expression after subset", subcategory="volcano plot"),
        sign_up_subset_tar = report(OUTPUTDIR + "04_diffexp_subset/up_down_regulated_genes_list.tar.gz", caption = REPORT + "de_significatif.rst", category="04 differential expression after subset", subcategory="significant regulated genes list"),
        violinplot_subset_tar = report(OUTPUTDIR + "04_diffexp_subset/violin_plot.tar.gz", caption = REPORT + "violin.rst", category="04 differential expression after subset", subcategory="violin"),
    params:
        outputdir = OUTPUTDIR
    message: 
        "Compressing the files for snakemake report"
    shell:
        """
        outputdir=({params.outputdir})
        cd ${{outputdir}}02_seurat/
        tar -czvf data_matrix.tar.gz -C ${{outputdir}}02_seurat/ data_matrix
        cd ${{outputdir}}03_diffexp/
        tar -czvf violin_plot.tar.gz -C ${{outputdir}}03_diffexp/ violin_plot
        tar -czvf umapfeature_plot.tar.gz -C ${{outputdir}}03_diffexp/ umap_plot
        tar -czvf tsnefeature_plot.tar.gz -C ${{outputdir}}03_diffexp/ tsne_plot
        tar -czvf ridgefeature_plot.tar.gz -C ${{outputdir}}03_diffexp/ ridge_plot
        tar -czvf heatmapfeature.tar.gz -C ${{outputdir}}03_diffexp/ heatmap
        tar -czvf differencial_expression_tests.tar.gz -C ${{outputdir}}03_diffexp/ differential_expression_features
        tar -czvf volcano_plot.tar.gz -C ${{outputdir}}03_diffexp/ volcanoplot
        tar -czvf up_down_regulated_genes_list.tar.gz -C ${{outputdir}}03_diffexp/ de_significant
        cd ${{outputdir}}04_diffexp_subset/
        tar -czvf violin_plot.tar.gz -C ${{outputdir}}04_diffexp_subset/ violin_plot
        tar -czvf differencial_expression_tests.tar.gz -C ${{outputdir}}04_diffexp_subset/ differential_expression
        tar -czvf volcano_plot.tar.gz -C ${{outputdir}}04_diffexp_subset/ volcanoplot
        tar -czvf up_down_regulated_genes_list.tar.gz -C ${{outputdir}}04_diffexp_subset/ de_significant
        """

rule cleaning:
    input:
        data_matrix_tar = OUTPUTDIR + "02_seurat/data_matrix.tar.gz",
        violinplot_tar = OUTPUTDIR + "03_diffexp/violin_plot.tar.gz",
        umapfeature_tar = OUTPUTDIR + "03_diffexp/umapfeature_plot.tar.gz",
        tsnefeature_tar = OUTPUTDIR + "03_diffexp/tsnefeature_plot.tar.gz",
        ridgefeature_tar = OUTPUTDIR + "03_diffexp/ridgefeature_plot.tar.gz",
        heatmapfeature_tar = OUTPUTDIR + "03_diffexp/heatmapfeature.tar.gz",
        defile_tar = OUTPUTDIR + "03_diffexp/differencial_expression_tests.tar.gz",
        volcano_tar = OUTPUTDIR + "03_diffexp/volcano_plot.tar.gz",
        sign_up_down_tar = OUTPUTDIR + "03_diffexp/up_down_regulated_genes_list.tar.gz",
        violinplot_tar_sub = OUTPUTDIR + "04_diffexp_subset/violin_plot.tar.gz",
        defile_tar_sub = OUTPUTDIR + "04_diffexp_subset/differencial_expression_tests.tar.gz",
        volcano_tar_sub = OUTPUTDIR + "04_diffexp_subset/volcano_plot.tar.gz",
        sign_up_down_tar_sub = OUTPUTDIR + "04_diffexp_subset/up_down_regulated_genes_list.tar.gz",
    output:
        clean = OUTPUTDIR + "03_diffexp/clean.txt"
    params:
        outputdir = OUTPUTDIR
    message: 
        "Compressing the files for snakemake report"
    shell:
        """
        outputdir=({params.outputdir})
        cd ${{outputdir}}03_diffexp/
        rm -r ${{outputdir}}03_diffexp/violin_plot
        rm -r ${{outputdir}}03_diffexp/umap_plot
        rm -r ${{outputdir}}03_diffexp/tsne_plot
        rm -r ${{outputdir}}03_diffexp/ridge_plot
        rm -r ${{outputdir}}03_diffexp/heatmap
        rm -r ${{outputdir}}03_diffexp/differential_expression_features
        rm -r ${{outputdir}}03_diffexp/volcanoplot
        rm -r ${{outputdir}}03_diffexp/de_significant
        cd ${{outputdir}}04_diffexp_subset/
        rm -r ${{outputdir}}04_diffexp_subset/violin_plot
        rm -r ${{outputdir}}04_diffexp_subset/differential_expression
        rm -r ${{outputdir}}04_diffexp_subset/volcanoplot
        rm -r ${{outputdir}}04_diffexp_subset/de_significant
        echo "cleaning ok" > ${{outputdir}}03_diffexp/clean.txt
        """