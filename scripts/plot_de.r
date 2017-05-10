modules::import_package('dplyr', attach = TRUE)
modules::import('klmr/ggplots', attach = TRUE)

plot_disp_ests = function (dds) {
    dispersion_data = S4Vectors::mcols(dds) %>%
        as.data.frame() %>%
        mutate(disp = DESeq2::dispersions(dds)) %>%
        filter(baseMean > 0)

    ymin = with(dispersion_data,
                10 ^ floor(log10(min(dispGeneEst[dispGeneEst > 0],
                                     na.rm = TRUE)) - 0.1))

    ggplot(dispersion_data) +
        aes(baseMean, pmax(dispGeneEst, ymin)) +
        geom_point(aes(y = disp), filter(dispersion_data, dispOutlier),
                   size = 2.5, color = 'dodgerblue') +
        geom_point(aes(color = 'estimate'), size = 0.75) +
        geom_point(aes(y = disp), filter(dispersion_data, ! dispOutlier),
                   color = 'dodgerblue', size = 0.75) +
        geom_point(aes(y = dispFit), color = 'red', size = 0.75) +
        scale_color_manual('Dispersion value',
                           values = c('black', 'dodgerblue', 'red'),
                           limits = c('estimate', 'fitted', 'final')) +
        scale_x_log10() +
        scale_y_log10() +
        labs(x = 'Mean of normalised counts', y = 'Dispersion') +
        guides(color = guide_legend(override.aes = list(size = 1.5)))
}

plot_ma = function (results) {
    data = results %>%
        as.data.frame() %>%
        mutate(DE = ! is.na(padj) & padj < 0.05)

    ylim = with(data, quantile(abs(log2FoldChange[is.finite(log2FoldChange)]),
                               probs = 0.99) * 1.1 * c(-1, 1))

    ggplot(data) +
        aes(baseMean, pmax(ylim[1], pmin(ylim[2], log2FoldChange)),
            color = DE,
            shape = ifelse(log2FoldChange < ylim[1], 'low',
                           ifelse(log2FoldChange > ylim[2], 'high', 'normal'))) +
        geom_point() +
        geom_segment(aes(x, xend = xend, yend = log2FoldChange),
                     data.frame(x = 0, xend = Inf, log2FoldChange = 0), color = 'red') +
        scale_color_manual(values = c('black', 'red'),
                           guide = guide_legend(title = NULL),
                           breaks = TRUE,
                           labels = 'differentially expressed') +
        scale_shape_manual(values = c(low = 6, high = 2, normal = 16), guide = FALSE) +
        scale_x_log10() +
        labs(x = 'Mean of normalised counts', y = expression(log[2]~FC)) +
        theme(legend.position = 'bottom')
}

plot_pca = function (dds, intgroup) {
    vsd = DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
    pca_data = BiocGenerics::plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
    var_exp = attr(pca_data, 'percentVar') * 100

    modules::import_package('ggrepel', attach = TRUE)

    ggplot(pca_data) +
        aes(PC1, PC2) +
        geom_point(aes_string(color = intgroup[1], shape = intgroup[2])) +
        geom_text_repel(aes(label = name), segment.width = 0, segment.color = 'white') +
        coord_fixed() +
        labs(x = sprintf('PC1 (%.2f%% variance explained)', var_exp[1]),
             y = sprintf('PC1 (%.2f%% variance explained)', var_exp[2]))
}
