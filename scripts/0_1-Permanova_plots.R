kretxeu.permanova <- function(mypermanova = NULL,pval_thres = 0.05,
                            legend_proportion_size =2,
                            y_vjust=1,size_axis_text=20,
                            size_axis_title=30,size_legend_text=20,
                            size_title_text = 30,
                            size_ticks_x = 2.5, size_ticks_y =2.5, 
                            font_family = "Helvetica",aspect.ratio =3,
                            size_panel_border = 0.1,
                            terms_exclude_plot = NULL){
  pval_thres <- pval_thres
  df_aov <- mypermanova %>% as.data.frame
  colnames(df_aov)[5] <- "pvalue"
  df_aov$Term <- rownames(df_aov) %>% factor
  df_aov$Description <- rep("Term",nrow(df_aov))
  df_aov <- df_aov %>% subset(Term != "Total") %>%
    droplevels
  df_aov <- df_aov %>% subset(pvalue <= pval_thres)
  
  residual <- 1- sum(df_aov$R2)
  df_plot <- df_aov[,c(3,5:7)] %>% droplevels
  df_plot <- with(df_plot,order(R2)) %>% df_plot[.,]
  ord_ele <- df_plot$Term %>% as.character
  df_plot <- data.frame(R2 = residual,pvalue = 1,
                        Term = "Residual",Description = "Term") %>%
    rbind(df_plot,.)
  df_plot$Term <- factor(df_plot$Term,levels = c("Residual",ord_ele))
  rownames(df_plot) <- NULL
  df_plot$R2 <- df_plot$R2*100
  df_plot_original <- df_plot
  if(! is.null(terms_exclude_plot)){
    df_plot <- which(!(df_plot$Term %in% terms_exclude_plot)) %>%
      df_plot[.,] %>% droplevels
  }
  p_var <- ggplot(data = df_plot, aes(x = Description,y = R2, fill = Term)) +
    geom_bar(stat = "identity") + ylab(label = "Variance Explained (%)") +
    theme(aspect.ratio = aspect.ratio) +
    scale_y_continuous(breaks=seq(0,100,10)) +
    coord_cartesian(ylim = c(0,100), expand = FALSE) +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", size = size_panel_border, fill = NA),
      axis.ticks.y =element_line(colour = "black", size = size_ticks_y),
      axis.ticks.x =element_blank(),
      axis.text.x =element_blank(),
      axis.text.y = element_text(family = font_family,face="bold",size=size_axis_text,colour="black",vjust = y_vjust),
      axis.title.x = element_blank(),
      axis.title.y = element_text(family = font_family,face="bold",size=size_axis_title,colour="black"),
      legend.background = element_blank(),legend.key.size = unit(legend_proportion_size,"line"),
      legend.title=element_text(size=size_title_text,
                                family = font_family,face = "plain",colour = "black"),
      legend.key = element_blank(),
      legend.text = element_text(size=size_legend_text,family = font_family,face = "bold",colour = "black"),
      legend.position ="right"
    )
  return(list(df = df_plot_original,p = p_var))
}
