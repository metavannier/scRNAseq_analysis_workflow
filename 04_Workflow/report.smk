rule archive_report:
    input:
        count_matrix = OUTPUTDIR + "02_seurat/" + NPROJ + "_count_matrix.csv",
        data_matrix = OUTPUTDIR + "02_seurat/" + NPROJ + "_data_matrix.csv",
        scale_data_matrix = OUTPUTDIR + "02_seurat/" + NPROJ + "_scale_data_matrix.csv",
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
    output:
        data_matrix_tar = report(OUTPUTDIR + "02_seurat/data_matrix.tar.gz", caption = REPORT + "data_matrix.rst", category="02 seurat"),
        violinplot_tar = report(OUTPUTDIR + "03_diffexp/violin_plot.tar.gz", caption = REPORT + "violin.rst", category="03 differential expression", subcategory="violin"),
        umapfeature_tar = report(OUTPUTDIR + "03_diffexp/umapfeature_plot.tar.gz", caption = REPORT + "umapfeature.rst", category="03 differential expression", subcategory="umap feature"),
        tsnefeature_tar = report(OUTPUTDIR + "03_diffexp/tsnefeature_plot.tar.gz", caption = REPORT + "tsnefeature.rst", category="03 differential expression", subcategory="tsne feature"),
        ridgefeature_tar = report(OUTPUTDIR + "03_diffexp/ridgefeature_plot.tar.gz", caption = REPORT + "ridgefeature.rst", category="03 differential expression", subcategory="ridge feature"),
        heatmapfeature_tar = report(OUTPUTDIR + "03_diffexp/heatmapfeature.tar.gz", caption = REPORT + "heatmapfeature.rst", category="03 differential expression", subcategory="heatmap feature"),
        defile_tar = report(OUTPUTDIR + "03_diffexp/differencial_expression_tests.tar.gz", caption = REPORT + "de_tests.rst", category="03 differential expression", subcategory="differencial expression tests"),
        volcano_tar = report(OUTPUTDIR + "03_diffexp/volcano_plot.tar.gz", caption = REPORT + "volcanoplot.rst", category="03 differential expression", subcategory="volcano plot"),
        sign_up_tar = report(OUTPUTDIR + "03_diffexp/up_regulated_genes_list.tar.gz", caption = REPORT + "de_significatif.rst", category="03 differential expression", subcategory="up regulated genes list"),
        sign_down_tar = report(OUTPUTDIR + "03_diffexp/down_regulated_genes_list.tar.gz", caption = REPORT + "de_significatif.rst", category="03 differential expression", subcategory="down regulated genes list"),
    params:
        outputdir = OUTPUTDIR
    message: 
        "Compressing the files for snakemake report"
    shell:
        """
        outputdir=({params.outputdir})
        cd ${{outputdir}}02_seurat/
        tar -czvf data_matrix.tar.gz -C ${{outputdir}}02_seurat/data_matrix
        cd ${{outputdir}}03_diffexp/
        tar -czvf violin_plot.tar.gz -C ${{outputdir}}03_diffexp/ violin_plot
        rm -r ${{outputdir}}03_diffexp/violin_plot
        tar -czvf umapfeature_plot.tar.gz -C ${{outputdir}}03_diffexp/ umap_plot
        rm -r ${{outputdir}}03_diffexp/umap_plot
        tar -czvf tsnefeature_plot.tar.gz -C ${{outputdir}}03_diffexp/ tsne_plot
        rm -r ${{outputdir}}03_diffexp/tsne_plot
        tar -czvf ridgefeature_plot.tar.gz -C ${{outputdir}}03_diffexp/ ridge_plot
        rm -r ${{outputdir}}03_diffexp/ridge_plot
        tar -czvf heatmapfeature.tar.gz -C ${{outputdir}}03_diffexp/ heatmap
        rm -r ${{outputdir}}03_diffexp/heatmap
        tar -czvf differencial_expression_tests.tar.gz -C ${{outputdir}}03_diffexp/ differential_expression_features
        rm -r ${{outputdir}}03_diffexp/differential_expression_features
        tar -czvf volcano_plot.tar.gz -C ${{outputdir}}03_diffexp/ volcanoplot
        rm -r ${{outputdir}}03_diffexp/volcanoplot
        tar -czvf up_regulated_genes_list.tar.gz -C ${{outputdir}}03_diffexp/ de_significant
        rm -r ${{outputdir}}03_diffexp/de_significant
        tar -czvf down_regulated_genes_list.tar.gz -C ${{outputdir}}03_diffexp/ de_significant
        rm -r ${{outputdir}}03_diffexp/de_significant
        """