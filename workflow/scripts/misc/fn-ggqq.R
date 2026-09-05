## desc: Function to make Q-Q plot
##
##       Takes as input a data.frame object with columns for:
##       - chromosome (named CHR by default)
##       - position (named POS by default)
##       - -log10-transformed p-values (named LOG10_P by default)


library(ggplot2)



ggqq <- function(df, chr = 'CHR', pos = 'POS', log_pv = 'LOG10_P',
                 title = '', transparent = FALSE) {
    df <- df[!is.na(df[[log_pv]]),]
    n_vars   <- nrow(df)
    unif_vec <- seq(1 / (n_vars + 1), n_vars / (n_vars + 1), length.out = n_vars)
    df_qq    <- data.frame(the = sort(-log10(unif_vec)),
                           obs = sort(df[[log_pv]]))

    qq <- ggplot(df_qq, aes(the, obs)) +
        geom_abline(slope = 1, intercept = 0, linewidth= 0.4) +
        geom_point(size = 0.2, stroke = 0.2, col = '#0095AF') +
        ggtitle(title) +
        labs(x = 'Uniform quantiles',
             y = 'Empirical p-value quantiles') +
        theme_bw() 

    ## add theme options
    if (transparent) {
        qq <- qq +
            theme(plot.title = element_text(size = 7),
                  axis.title = element_text(size = 6),
                  axis.text  = element_text(size = 5),
                  legend.position = 'none',
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.background = element_rect(fill = 'transparent', colour = NA),
                  plot.background  = element_rect(fill = 'transparent', colour = NA))
    } else {
        qq <- qq +
            theme(plot.title = element_text(size = 7),
                  axis.title = element_text(size = 6),
                  axis.text  = element_text(size = 5),
                  legend.position = 'none',
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank())
    }

    return(qq)
}


ggqq_line <- function(df, chr = 'CHR', pos = 'POS', log_pv = 'LOG10_P',
                 title = '', transparent = FALSE) {
    df <- df[!is.na(df[[log_pv]]),]
    n_vars   <- nrow(df)
    unif_vec <- seq(1 / (n_vars + 1), n_vars / (n_vars + 1), length.out = n_vars)
    df_qq    <- data.frame(the = sort(-log10(unif_vec)),
                           obs = sort(df[[log_pv]]))

    qq <- ggplot(df_qq, aes(the, obs)) +
        geom_abline(slope = 1, intercept = 0, linewidth= 0.4) +
        geom_line(linewidth= 0.5, col = '#0095AF') +
        ggtitle(title) +
        labs(x = 'Uniform quantiles',
             y = 'Empirical p-value quantiles') +
        theme_bw() 

    ## add theme options
    if (transparent) {
        qq <- qq +
            theme(plot.title = element_text(size = 7),
                  axis.title = element_text(size = 6),
                  axis.text  = element_text(size = 5),
                  legend.position = 'none',
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.background = element_rect(fill = 'transparent', colour = NA),
                  plot.background  = element_rect(fill = 'transparent', colour = NA))
    } else {
        qq <- qq +
            theme(plot.title = element_text(size = 7),
                  axis.title = element_text(size = 6),
                  axis.text  = element_text(size = 5),
                  legend.position = 'none',
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank())
    }

    return(qq)
}
