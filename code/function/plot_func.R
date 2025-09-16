
## Functions for plot 
theme_for_coe_plot = function(x){
  ggplot2::theme_test() + 
    theme(text = element_text(family = 'Arial', face = "plain", 
                              colour = "black", size = 18, lineheight = 0.9, 
                              hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), 
                              debug = FALSE),
          legend.position = c(0.18, 0.15),
          legend.background = element_blank(),
          rect = element_rect(fill = "white", 
                              colour = "black", linewidth = 2.5, linetype = 1),
          legend.key = element_blank(),
          axis.ticks = element_line(colour = "black", linewidth = 1.5),
          axis.text = element_text(size = rel(0.8), 
                                   colour = "black"),
          axis.title = element_text(size = 15, 
                                    colour = "black"),
          axis.title.x = element_text(margin = margin(t = 0.3, r = 0,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          axis.title.y = element_text(margin = margin(t = 0, r = 0.3,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          plot.margin = margin(t = 0.2, r = 0.2, b = 0.5, l = 0.2, 
                               unit = "cm"))
}

theme_regular = function(x){
  ggplot2::theme_test() + 
    theme(text = element_text(family = 'Arial', face = "plain", 
                              colour = "black", size = 14.4, lineheight = 0.9, 
                              hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), 
                              debug = FALSE),
          legend.position = 'None',
          legend.background = element_blank(),
          rect = element_rect(fill = "white", 
                              colour = "black", linewidth = 2.5, linetype = 1),
          legend.key = element_blank(),
          axis.ticks = element_line(colour = "black", linewidth = 1.5),
          axis.text = element_text(size = rel(1), 
                                   colour = "black"),
          axis.title = element_text(size = 15, 
                                    colour = "black"),
          axis.title.x = element_text(margin = margin(t = 0.3, r = 0,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          axis.title.y = element_text(margin = margin(t = 0, r = 0.3,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          plot.margin = margin(t = 0.2, r = 0.2, b = 0.5, l = 0.2, 
                               unit = "cm"))
}


theme_regular_1 = function(x){
    theme(text = element_text(family = 'Arial', face = "plain", 
                              colour = "black", size = 16, lineheight = 0.9, 
                              hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), 
                              debug = FALSE),
          legend.position = 'None',
          legend.background = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_line(colour = "#93785B", linewidth = 0.6),
          axis.text = element_text(size = 14, 
            colour = "black"),
          axis.title = element_text(size = 14, 
            colour = "black"),
          axis.title.x = element_text(margin = margin(t = 0.3, r = 0,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          axis.title.y = element_text(margin = margin(t = 0, r = 0.3,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          plot.margin = margin(t = 0.2, r = 0.2, b = 0.5, l = 0.2, 
                               unit = "cm"))
}



theme_regular_2 = function(x){
  ggplot2::theme_test() + 
    theme(text = element_text(family = 'Arial', face = "plain", 
                              colour = "black", size = 14, lineheight = 0.9, 
                              hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), 
                              debug = FALSE),
          legend.position = 'None',
          legend.background = element_blank(),
          panel.border = element_rect(color = "#93785B",
                                      fill = NA,
                                      linewidth = 1),
          rect = element_rect(fill = "white", 
                              colour = "black", linewidth = 2.5, linetype = 1),
          legend.key = element_blank(),
          axis.ticks = element_line(colour = "#93785B", linewidth = 0.6),
          axis.text = element_text(#size = rel(1), 
                                   colour = "black"),
          axis.title = element_text(#size = 15, 
                                    colour = "black"),
          axis.title.x = element_text(margin = margin(t = 0.3, r = 0,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          axis.title.y = element_text(margin = margin(t = 0, r = 0.3,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          plot.margin = margin(t = 0.2, r = 0.2, b = 0.5, l = 0.2, 
                               unit = "cm"))
}

theme_for_sem = function(x){
    theme(legend.position = 'None',
          legend.background = element_blank(),
          rect = element_rect(fill = "white", 
                              colour = "black", linewidth = 2.5, linetype = 1),
          legend.key = element_blank(),
          #axis.ticks = element_line(colour = "black", linewidth = 0.8),
          axis.text = element_text(size = rel(0.7), 
                                   colour = "black"),
          axis.title=element_blank(),
          strip.text=element_text(face="plain", 
                                  colour = "black"),
          plot.margin=unit(c(0.4,0.4,0.2,-0.2),units="lines"))
}

CMYKtoRGB = function(C,M,Y,K){
  Rc <- (1-C)*(1-K)
  Gc <- (1-M)*(1-K)
  Bc <- (1-Y)*(1-K)
  R <- as.character(as.hexmode((ceiling(Rc*255))))
  G <- as.character(as.hexmode((ceiling(Gc*255))))
  B <- as.character(as.hexmode((ceiling(Bc*255))))
  if(nchar(R) == 1){
    R <- paste0('0',R)
  }
  if(nchar(G) == 1){
    G <- paste0('0',G)
  }
  if(nchar(B) == 1){
    B <- paste0('0',B)
  }
  RGBColor <- paste0('#',R,G,B)
  return(RGBColor)
}

##==== Function ====##