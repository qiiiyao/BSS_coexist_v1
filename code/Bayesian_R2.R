load("D:/R projects/BSS/results/fit_results/BSS_exculde_trees_raw/bayesian_r2_dat.rdata")
library(ggplot2)
bayesian_r2_dat$name = 'Original'
bayesian_r2 = ggplot(data = bayesian_r2_dat, aes(y = bayesian_r2, x = name, 
                                                 color = name, fill = name)) +
  geom_violin(width = 3) +
  geom_boxplot(width = 1.2, color="black", alpha=1) +
  #geom_hline(yintercept = mean(p_rant_mean_sp_ul_r$r_mean,
  #                             na.rm = T),
  #          color = 'black', linetype = 2, size = 2) +
  scale_y_continuous(limits = c(0.4, 0.9), breaks = seq(0, 1, length.out = 11)) +
  scale_fill_brewer(palette="Accent")+
  scale_color_brewer(palette="Accent")+
  labs(x = '', y = expression(paste('R'^2)), 
       title = '') +
  theme(legend.position = 'none',
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.x = element_text(margin = margin(t = 20, r = 20, b = 20, l = 20),
                                   color = "black"),
        axis.text.y = element_text(margin = margin(t = 20, r = 20, b = 20, l = 20),
                                   color = "black"),
        plot.title = element_text(hjust = 0.5,
                                  vjust = 0.5, size=80, face='bold'),
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 180),
        axis.ticks.y = element_line(linetype=1,color="black",size=1.2),
        axis.ticks.x = element_line(linetype=1,color="black",size=1.2),
        axis.ticks.length.y = unit(-0.3,'cm'),
        axis.ticks.length.x = unit(-0.3,'cm'))

ggsave('results/figures_alltime_raw/Fig.S1.svg',
       plot=bayesian_r2, 
       width=105, height=120, dpi=600, units='cm',
       limitsize=F)
