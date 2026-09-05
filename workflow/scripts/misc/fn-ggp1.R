## desc: Function to make plot displaying validation R-squared for different p1 thresholds
##
##       Takes as input a data.frame object with columns for:
##       - p1 thresholds
##       - R-squared value
##       - optimum p1 value and corresponding R-squared for labelling


library(ggplot2)
library(ggrepel)



ggp1 <- function(df, p1 = 'p1', rsq = 'rsq', p1_opt, rsq_opt, title = '', transparent = FALSE) {
    ## compute range of p1 for optimum point label
    rg_p1  <- max(df$p1) - min(df$p1)
    rg_rsq <- max(df$rsq)

    p1_plot <- ggplot(data = df, aes_string(x = p1, y = rsq)) +
        geom_vline(xintercept = p1_opt, size = 0.3, linetype = 2, col = '#898989') +
        ## geom_hline(yintercept = rsq_opt, size = 0.3, linetype = 2, col = '#898989') +
        geom_line(size = 0.4, col = '#0095AF') + 
        geom_point(aes(x = p1_opt, y = rsq_opt), size = 0.8, col = '#0095AF') +
        geom_label_repel(data = data.frame(x = p1_opt, y = rsq_opt), aes(x = x, y = y),
                         label = paste0('p-value = ', format(round(p1_opt, digits = 10), scientific = FALSE),
                                        '; R-squared = ', round(rsq_opt, digits = 4)),
                         box.padding = 0.5,
                         nudge_x = 0.3 * rg_p1,
                         nudge_y = -0.3 * rg_rsq,
                         size = 2.6,
                         segment.size = 0.3,
                         color = '#898989') +
        scale_x_continuous(trans = 'log10', breaks = c(1e-9, 1e-7, 1e-5, 1e-3, 1e-1),
                           limits = c(NA, 0.1)) +
        ylim(0, NA) +
        ggtitle(title) +
        labs(x = 'Global p-value threshold',
             y = 'R-squared') +
        theme_bw()

    ## add theme options
    if (transparent) {
        p1_plot <- p1_plot + 
            theme(plot.title = element_text(size = 9),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7),
                  panel.background = element_rect(fill = "transparent", colour = NA),
                  plot.background  = element_rect(fill = "transparent", colour = NA))
    } else {
        p1_plot <- p1_plot + 
            theme(plot.title = element_text(size = 9),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7))
    }

    return(p1_plot)
}




## alternative version for iterative step script
ggp1_iter <- function(df, p1 = 'p1', rsq = 'rsq_vali', p1_opt, rsq_opt, title = '', transparent = FALSE) {
    ## compute range of p1 for optimum point label
    rg_p1  <- max(df$p1) - min(df$p1)

    p1_plot <- ggplot(data = df, aes(x = .data[[p1]], y = .data[[rsq]])) +
        geom_vline(xintercept = p1_opt, size = 0.3, linetype = 2, col = '#898989') +
        ## geom_hline(yintercept = rsq_opt, size = 0.3, linetype = 2, col = '#898989') +
        geom_line(size = 0.4, col = '#0095AF') + 
        geom_point(aes(x = p1_opt, y = rsq_opt), size = 0.8, col = '#0095AF') +
        geom_label_repel(data = data.frame(x = p1_opt, y = rsq_opt), aes(x = x, y = y),
                         label = paste0('p-value = ', format(round(p1_opt, digits = 20), scientific = FALSE),
                                        '; R-squared = ', round(rsq_opt, digits = 10)),
                         box.padding = 0.5,
                         nudge_x = 0.3 * rg_p1,
                         nudge_y = -0.3 * rsq_opt,
                         size = 2.6,
                         segment.size = 0.3,
                         color = '#898989') +
        scale_x_continuous(trans = 'log10', breaks = c(1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 1),
                           limits = c(NA, 1)) +
        ylim(0, NA) +
        ggtitle(title) +
        labs(x = 'p-value threshold',
             y = 'R-squared') +
        theme_bw()

    ## add theme options
    if (transparent) {
        p1_plot <- p1_plot + 
            theme(plot.title = element_text(size = 9),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7),
                  panel.background = element_rect(fill = "transparent", colour = NA),
                  plot.background  = element_rect(fill = "transparent", colour = NA))
    } else {
        p1_plot <- p1_plot + 
            theme(plot.title = element_text(size = 9),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7))
    }

    return(p1_plot)
}
