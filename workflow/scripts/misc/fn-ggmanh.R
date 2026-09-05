## desc: Function to make Manhattan plot
##
##       Takes as input a data.frame object with columns for:
##       - chromosome (named CHR by default)
##       - position (named POS by default)
##       - minus-log10-transformed p-values (named LOG10_P by default)


library(dplyr)
library(ggplot2)



ggmanh <- function(df, chr = "CHR", pos = "POS", log_pv = "LOG10_P",
                   sig_line_thr = 5e-8, sig_line_col = "black",
                   title = "", gap = 2e7, transparent = FALSE) {
    ## make the first SNP of each chromosome have position 0
    df <- df %>%
        arrange(.data[[chr]], .data[[pos]]) %>%
        group_by(.data[[chr]]) %>%
        mutate(POS_0 = .data[[pos]] - .data[[pos]][1])

    df <- df %>%
        ## compute chromosome size
        group_by(.data[[chr]]) %>%
        dplyr::summarise(chr_len = max(POS_0)) %>%
        ## calculate cumulative position of each chromosome
        mutate(gap_add = seq(0, gap * (n() - 1), gap)) %>%
        mutate(tot = cumsum(as.numeric(chr_len)) - chr_len + gap_add) %>%
        select(-chr_len) %>%  
        ## add this info to the initial dataset
        left_join(df, ., by = chr) %>%
        ## add a cumulative position of each SNP
        arrange(.data[[chr]], POS_0) %>%
        mutate(POScum = POS_0 + tot)

    ## get chromosome center positions for x-axis
    axisdf <- df %>%
        group_by(.data[[chr]]) %>%
        dplyr::summarise(centre = (max(POScum) + min(POScum)) / 2)

    ## plot
    cols = c("#00468B", "#0095AF")  # colours
    manh <- ggplot(df, aes(x = POScum, y = .data[[log_pv]])) +
        ## show all points
        geom_point(aes(color=as.factor(.data[[chr]])), alpha=1, size=0.3, stroke=0.2, shape = 19) +
        scale_color_manual(values = rep(cols, 22 )) +
        ## custom X axis
        scale_x_continuous(label = axisdf[[chr]], breaks= axisdf$centre,
                           expand = expansion(mult = c(0.015, 0.015))) +  # remove space at the edges
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.02))) +
        ## add plot and axis titles
        ggtitle(title) +
        labs(x = "Chromosome", y = expression("-log"[10]*"(p-value)")) +
        ## add genome-wide sig line
        geom_hline(yintercept = -log10(sig_line_thr),
                   linewidth = 0.4, linetype = 'solid', colour = sig_line_col) +
        theme_bw()

    ## add theme options
    if (transparent) {
        manh <- manh + 
            theme(plot.title = element_text(size = 9),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7),
                  legend.position = "none",
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.background = element_rect(fill = "transparent", colour = NA),
                  plot.background  = element_rect(fill = "transparent", colour = NA))
    } else {
        manh <- manh + 
            theme(plot.title = element_text(size = 9),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7),
                  legend.position = "none",
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank())
    }
    
    return(manh)
}




ggmanh_filter <- function(df, chr = "CHR", pos = "POS", log_pv = "LOG10_P",
                          sig_line_thr = 5e-8, sig_line_col = "black",
                          filter_var,
                          title = "", gap = 2e7, transparent = FALSE) {
    ## make the first SNP of each chromosome have position 0
    df <- df %>%
        arrange(.data[[chr]], .data[[pos]]) %>%
        group_by(.data[[chr]]) %>%
        mutate(POS_0 = .data[[pos]] - .data[[pos]][1])

    df <- df %>%
        ## compute chromosome size
        group_by(.data[[chr]]) %>%
        dplyr::summarise(chr_len = max(POS_0)) %>%
        ## calculate cumulative position of each chromosome
        mutate(gap_add = seq(0, gap * (n() - 1), gap)) %>%
        mutate(tot = cumsum(as.numeric(chr_len)) - chr_len + gap_add) %>%
        select(-chr_len) %>%  
        ## add this info to the initial dataset
        left_join(df, ., by = chr) %>%
        ## add a cumulative position of each SNP
        arrange(.data[[chr]], POS_0) %>%
        mutate(POScum = POS_0 + tot)

    ## get chromosome center positions for x-axis
    axisdf <- df %>%
        group_by(.data[[chr]]) %>%
        dplyr::summarise(centre = (max(POScum) + min(POScum)) / 2)

    ## plot
    cols = c("#00468B", "#0095AF")  # colours
    manh <- ggplot(df, aes(x = POScum,
                           y = .data[[log_pv]])) +
        ## show all points
        geom_blank() +  # same range as if all points were plotted
        geom_point(data = . %>% filter(.data[[filter_var]]),  # filter by filter_var indicator
                   aes(color=as.factor(.data[[chr]])),
                   alpha=1, size=0.3, stroke=0.2, shape = 19) +
        scale_color_manual(values = rep(cols, 22 )) +
        ## custom X axis
        scale_x_continuous(label = axisdf[[chr]], breaks= axisdf$centre,
                           expand = expansion(mult = c(0.015, 0.015))) +  # remove space at the edges
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.02))) +
        ## add plot and axis titles
        ggtitle(title) +
        labs(x = "Chromosome", y = expression("-log"[10]*"(p-value)")) +
        ## add genome-wide sig line
        geom_hline(yintercept = -log10(sig_line_thr),
                   linewidth = 0.4, linetype = 'solid', colour = sig_line_col) +
        theme_bw()

    ## add theme options
    if (transparent) {
        manh <- manh + 
            theme(plot.title = element_text(size = 8),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7),
                  legend.position = "none",
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.background = element_rect(fill = "transparent", colour = NA),
                  plot.background  = element_rect(fill = "transparent", colour = NA))
    } else {
        manh <- manh + 
            theme(plot.title = element_text(size = 8),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7),
                  legend.position = "none",
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank())
    }
    
    return(manh)
}




ggmanh_plotly <- function(df, chr = "CHR", pos = "POS", log_pv = "LOG10_P",
                          sig_line_thr = 5e-8, sig_line_col = "black",
                          title = "", gap = 2e7, transparent = FALSE) {
    ## make the first SNP of each chromosome have position 0
    df <- df %>%
        arrange(.data[[chr]], .data[[pos]]) %>%
        group_by(.data[[chr]]) %>%
        mutate(POS_0 = .data[[pos]] - .data[[pos]][1])

    df <- df %>%
        ## compute chromosome size
        group_by(.data[[chr]]) %>%
        dplyr::summarise(chr_len = max(POS_0)) %>%
        ## calculate cumulative position of each chromosome
        mutate(gap_add = seq(0, gap * (n() - 1), gap)) %>%
        mutate(tot = cumsum(as.numeric(chr_len)) - chr_len + gap_add) %>%
        select(-chr_len) %>%  
        ## add this info to the initial dataset
        left_join(df, ., by = chr) %>%
        ## add a cumulative position of each SNP
        arrange(.data[[chr]], POS_0) %>%
        mutate(POScum = POS_0 + tot)

    ## get chromosome center positions for x-axis
    axisdf <- df %>%
        group_by(.data[[chr]]) %>%
        dplyr::summarise(centre = (max(POScum) + min(POScum)) / 2)

    ## plot
    cols = c("#00468B", "#0095AF")  # colours
    manh <- ggplot(df, aes(x = POScum,
                           y = .data[[log_pv]])) +
        geom_blank() +  # same range as if all points were plotted
        geom_point(data = . %>% filter(.data[[log_pv]] >= -log10(sig_line_thr)),
                   aes(label1 = ID,
                       label2 = rsID,
                       label3 = A1,
                       label4 = MAF,
                       label5 = orig_sign_LOG10_P,
                       label6 = int_full_sign_LOG10_P,
                       label7 = int_loco_sign_LOG10_P,
                       label8 = ANNOVAR_FN,
                       label9 = ANNOVAR_NEAR_GENES,
                       label10 = VEP_Consequence,
                       label11 = VEP_SYMBOL,
                       color = as.factor(.data[[chr]])),
                   alpha=1, size=1, stroke=0.5, shape = 19) +
        scale_color_manual(values = rep(cols, 22 )) +
        ## custom X axis
        scale_x_continuous(label = axisdf[[chr]], breaks= axisdf$centre,
                           expand = expansion(mult = c(0.015, 0.015))) +  # remove space at the edges
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.02))) +
        ## add plot and axis titles
        ggtitle(title) +
        labs(x = "Chromosome", y = "-log_10(p-value)") +
        ## add genome-wide sig lines
        geom_hline(yintercept = -log10(sig_line_thr),
                   linewidth = 0.4, linetype = 'solid', colour = sig_line_col) +
        theme_bw()

    ## add theme options
    if (transparent) {
        manh <- manh + 
            theme(plot.title = element_text(size = 9),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7),
                  legend.position = "none",
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.background = element_rect(fill = "transparent", colour = NA),
                  plot.background  = element_rect(fill = "transparent", colour = NA))
    } else {
        manh <- manh + 
            theme(plot.title = element_text(size = 9),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7),
                  legend.position = "none",
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank())
    }
    
    return(manh)
}




ggmanh_plotly_filter <- function(df, chr = "CHR", pos = "POS", log_pv = "LOG10_P",
                                 sig_line_thr = 5e-8, sig_line_col = "black",
                                 filter_var,
                                 title = "", gap = 2e7, transparent = FALSE) {
    ## make the first SNP of each chromosome have position 0
    df <- df %>%
        arrange(.data[[chr]], .data[[pos]]) %>%
        group_by(.data[[chr]]) %>%
        mutate(POS_0 = .data[[pos]] - .data[[pos]][1])

    df <- df %>%
        ## compute chromosome size
        group_by(.data[[chr]]) %>%
        dplyr::summarise(chr_len = max(POS_0)) %>%
        ## calculate cumulative position of each chromosome
        mutate(gap_add = seq(0, gap * (n() - 1), gap)) %>%
        mutate(tot = cumsum(as.numeric(chr_len)) - chr_len + gap_add) %>%
        select(-chr_len) %>%  
        ## add this info to the initial dataset
        left_join(df, ., by = chr) %>%
        ## add a cumulative position of each SNP
        arrange(.data[[chr]], POS_0) %>%
        mutate(POScum = POS_0 + tot)

    ## get chromosome center positions for x-axis
    axisdf <- df %>%
        group_by(.data[[chr]]) %>%
        dplyr::summarise(centre = (max(POScum) + min(POScum)) / 2)

    ## plot
    cols = c("#00468B", "#0095AF")  # colours
    manh <- ggplot(df, aes(x = POScum,
                           y = .data[[log_pv]])) +
        geom_blank() +  # same range as if all points were plotted
        geom_point(data = . %>% filter(.data[[filter_var]]),
                   aes(label1 = ID,
                       label2 = rsID,
                       label3 = A1,
                       label4 = MAF,
                       label5 = orig_sign_LOG10_P,
                       label6 = int_full_sign_LOG10_P,
                       label7 = int_loco_sign_LOG10_P,
                       label8 = ANNOVAR_FN,
                       label9 = ANNOVAR_NEAR_GENES,
                       label10 = VEP_Consequence,
                       label11 = VEP_SYMBOL,
                       color = as.factor(.data[[chr]])),
                   alpha=1, size=1, stroke=0.5, shape = 19) +
        scale_color_manual(values = rep(cols, 22 )) +
        ## custom X axis
        scale_x_continuous(label = axisdf[[chr]], breaks= axisdf$centre,
                           expand = expansion(mult = c(0.015, 0.015))) +  # remove space at the edges
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.02))) +
        ## add plot and axis titles
        ggtitle(title) +
        labs(x = "Chromosome", y = "-log_10(p-value)") +
        ## add genome-wide sig lines
        geom_hline(yintercept = -log10(sig_line_thr),
                   linewidth = 0.4, linetype = 'solid', colour = sig_line_col) +
        theme_bw()

    ## add theme options
    if (transparent) {
        manh <- manh + 
            theme(plot.title = element_text(size = 9),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7),
                  legend.position = "none",
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.background = element_rect(fill = "transparent", colour = NA),
                  plot.background  = element_rect(fill = "transparent", colour = NA))
    } else {
        manh <- manh + 
            theme(plot.title = element_text(size = 9),
                  axis.title = element_text(size = 8),
                  axis.text  = element_text(size = 7),
                  legend.position = "none",
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank())
    }
    
    return(manh)
}
