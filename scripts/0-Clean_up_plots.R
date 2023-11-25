#Clean graph
size=25
alpha=1
stroke=1.5
col_shape_background="white"
alpha_shape_background=0
lines_zero = TRUE
size_axis_line=2
type_axis_line = "longdash"
ratio_size_shape_background=1.3
y_vjust=0.5
x_hjust=0.5
size_axis_text=20
size_axis_title=30
size_legend_text=20
size_title_text = 30
legend_proportion_size=0.5
size_lines_panel = 0
size_panel_border = 2
font_family = "Helvetica"

clean <- theme(axis.line = element_blank(),
               panel.background = element_rect(fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(fill=NA,color =  "black",size = size_panel_border),
               axis.ticks = element_line(colour = "black",size = 2.5),
               axis.text.x = element_text(family = font_family,size =size_axis_text,colour="black",hjust = x_hjust),
               axis.text.y = element_text(family = font_family,size=size_axis_text,colour="black",vjust = y_vjust),
               axis.title.x = element_text(family = font_family,size = size_axis_title,colour = "black"),
               axis.title.y = element_text(family = font_family,size=size_axis_title,colour="black"),
               legend.background = element_blank(),legend.key.size = unit(legend_proportion_size,"line"),
               legend.title=element_text(size=size_title_text,
                                         family = font_family,face = "bold",colour = "black"),            
               #legend.title = element_blank(),
               legend.key = element_blank(),
               legend.text = element_text(size=size_legend_text,
                                          family = font_family,face = "bold",colour = "black"),
               legend.position ="right")

# function that draws key without gap
draw_key_polygon2 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}