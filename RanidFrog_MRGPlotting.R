source('RanidFrog_MRGAnalysis.R')

#Bar Plot
ggplot(melt.percent.fam[melt.percent.fam$variable %in% c("GF", "PF", "WF", "Tad"),], aes(fill=Group.1, y=variable, x=value)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic() + ylab('') + xlab('Mean Relative Abundance')+ 
  theme(legend.title = element_blank()) +   scale_fill_viridis(option="viridis", discrete = TRUE, direction = -1) +
  scale_y_discrete(labels=c( "L. clamitans", 
                             "L. palustris", "L. sylvaticus", "Lithobates tadpole")) +
  guides(fill=guide_legend(ncol=1))


#Bray/Jac Plots -- Figure 2
bray.plot <- ggplot() + 
  geom_point(data=frog2, aes(x=bray_x.nmds,y=bray_y.nmds,
                             color=species, shape = as.factor(infect)),size=5) +
  scale_color_manual(breaks = c("GF", "PF", "Tad", "WF"),
                     values=c("Green", "red", "orange", "blue"),
                     labels = c('L. clamitans', 'L. palustris', 'Lithobates tadpole',
                                'L. sylvaticus')) + 
  labs(col = "Amphibian") +
  coord_equal() +
  theme_half_open(12) + 
  scale_shape_manual(values=c(1, 16)) +
  theme(
    axis.title.x = element_text(size=axis.text.size),# remove x-axis labels
    axis.title.y = element_text(size=axis.text.size),# remove y-axis labels
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),#remove major-grid labels
    panel.grid.minor = element_blank(),#remove minor-grid labels
    plot.background = element_blank(),
    plot.margin = unit(c(.5, 0, 0, 0), "cm"),
    legend.title = element_text(size = legend.title.size),
    legend.text = element_text(size = legend.text.size)) +
  xlab('NMDS 1') + ylab('NMDS 2') + xlim(-2,2) + ylim(-1.5,2) + 
  labs(col="Amphibian", shape= NULL) + 
  guides(shape = FALSE, color = guide_legend(override.aes = list(shape = 1))); bray.plot


jac.plot <- ggplot() + 
  geom_point(data=frog2, aes(x=jac_x.nmds,y=jac_y.nmds,
                             color=species, shape = as.factor(infect)),size=5) +
  scale_color_manual(breaks = c("GF", "PF", "Tad", "WF"),
                     values=c("Green", "red", "orange", "blue"),
                     labels = c('L. clamitans', 'L. palustris', 'Lithobates tadpole',
                                'L. sylvaticus')) + 
  labs(col = "Amphibian") +
  coord_equal() +
  theme_half_open(12) + 
  scale_shape_manual(values=c(1, 16)) +
  theme(
    axis.title.x = element_text(size=axis.text.size),# remove x-axis labels
    axis.title.y = element_text(size=axis.text.size),# remove y-axis labels
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),#remove major-grid labels
    panel.grid.minor = element_blank(),#remove minor-grid labels
    plot.background = element_blank(),
    plot.margin = unit(c(.5, 0, 0, 0), "cm"),
    legend.title = element_text(size = legend.title.size),
    legend.text = element_text(size = legend.text.size)) +
  xlab('NMDS 1') + ylab('') + xlim(-2,2) + ylim(-1.5,2) + 
  labs(col="Amphibian", shape= NULL) + 
  guides(shape = FALSE, color = guide_legend(override.aes = list(shape = 16)));jac.plot


figure <- ggarrange(bray.plot, jac.plot,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
figure


legend.jac <- get_legend(
  # create some space to the left of the legend
  jac.plot + theme(legend.box.margin = margin(0, 0, 0, 10))
)

legend.loc <- get_legend(
  # create some space to the left of the legend
  Loc.plot + theme(legend.box.margin = margin(0, 0, 0, 5))
)


bray.jac.plot <- plot_grid(bray.plot + theme(legend.position = 'none'),
                           jac.plot + theme(legend.position = 'none'),
                           legend.jac,
                           rel_widths = c(1, 1, .7),
                           nrow = 1, align = 'vh', hjust = -1)


bray.jac.plot


#Alpha diversity species/sites plots -- Figure 4

axis.text.size <- 16
axis.text.lab <- 16

alpha.amph$species <- as.character(alpha.amph$species)

alpha.amph$species <- factor(alpha.amph$species, levels=c("GF", "PF", "WF", "Tad"))


sp.rich.plot <- ggplot(alpha.amph, aes(x = species, y = observed_OTU)) +
  theme_classic()+
  labs(x = 'Species', y = 'ASV Richness') +
  theme(plot.title = element_text(hjust = 0, face='bold'),
        axis.title.x = element_text(size=axis.text.size),
        axis.title.y = element_text(size=axis.text.size),
        axis.text.x = element_text(size = axis.text.lab, angle = 30, vjust = 1, 
                                   hjust = 1),
        axis.text.y = element_text(size = axis.text.lab)) +
  geom_boxplot(stat = "boxplot") +
  scale_x_discrete(labels= c('L. clamitans', 'L. palustris', 'L. sylvaticus', 'Lithobates tadpole'))


sp.sh.plot <- ggplot(alpha.amph, aes(x = species, y = shannon)) +
  theme_classic()+
  labs(x = 'Species', y = 'Shannon Diversity') +
  theme(plot.title = element_text(hjust = 0.5, face='bold'),
        axis.title.x = element_text(size=axis.text.size),
        axis.title.y = element_text(size=axis.text.size),
        axis.text.x = element_text(size = axis.text.lab, angle = 30, vjust = 1,
                                   hjust = 1),
        axis.text.y = element_text(size = axis.text.lab)) +
  geom_boxplot(stat = "boxplot") +
  scale_x_discrete(labels= c('L. clamitans', 'L. palustris','L. sylvaticus', 'Lithobates tadpole'))

sp.faith.plot <- ggplot(alpha.amph, aes(x = species, y = faith)) +
  theme_classic()+
  labs(x = 'Species', y = 'Faith Phlyogenetic Diversity') +
  theme(plot.title = element_text(hjust = 0.5, face='bold'),
        axis.title.x = element_text(size=axis.text.size),
        axis.title.y = element_text(size=axis.text.size),
        axis.text.x = element_text(size = axis.text.lab, angle = 30, vjust = 1, 
                                   hjust = 1),
        axis.text.y = element_text(size = axis.text.lab)) +
  geom_boxplot(stat = "boxplot") +
  scale_x_discrete(labels= c('L. clamitans', 'L. palustris', 'L. sylvaticus', 'Lithobates tadpole'))

ggarrange(sp.rich.plot, sp.sh.plot,
          labels = c("A", "B"),
          nrow = 1)



site.rich.plot <- ggplot(alpha, aes(x = site, y = observed_OTU)) +
  theme_classic()+
  labs(x = 'Site', y = 'ASV Richness') +
  theme(plot.title = element_text(hjust = 0.5, face='bold'),
        axis.title.x = element_text(size=axis.text.size),
        axis.title.y = element_text(size=axis.text.size),
        axis.text.x = element_text(size = axis.text.lab, angle = 30, 
                                   vjust = 1, hjust = 1),
        axis.text.y = element_text(size = axis.text.lab)) +
  geom_boxplot(stat = "boxplot")

site.sh.plot <- ggplot(alpha, aes(x = site, y = shannon)) +
  theme_classic()+
  labs(x = 'Site', y = 'Shannon Diversity') +
  theme(plot.title = element_text(hjust = 0.5, face='bold'),
        axis.title.x = element_text(size=axis.text.size),
        axis.title.y = element_text(size=axis.text.size),
        axis.text.x = element_text(size = axis.text.lab, angle = 30, 
                                   vjust = 1, hjust = 1),
        axis.text.y = element_text(size = axis.text.lab)) +  
  geom_boxplot(stat = "boxplot") 

site.faith.plot <- ggplot(alpha, aes(x = site, y = faith)) +
  theme_classic()+
  labs(x = 'Site', y = 'Faith Phylogenetic Diversity') +
  theme(plot.title = element_text(hjust = 0.5, face='bold'),
        axis.title.x = element_text(size=axis.text.size),
        axis.title.y = element_text(size=axis.text.size),
        axis.text.x = element_text(size = axis.text.lab, angle = 30, 
                                   vjust = 1, hjust = 1),
        axis.text.y = element_text(size = axis.text.lab)) +
  geom_boxplot(stat = "boxplot") 

#pdf('~/Desktop/MRGAlphaPlots_18Apr22.pdf',
#   width = 11, height = 8)
ggarrange(sp.rich.plot, sp.sh.plot, site.rich.plot, site.sh.plot,
          labels = c("A", "B", "C", "D"),
          nrow = 1)
#dev.off()

#Location Plot -- Figure 5
Loc.plot <- ggplot() + 
  geom_point(data=frog, aes(x=bray_x.nmds,y=bray_y.nmds,
                            color = location),size=5) + # add the point markers
  scale_color_viridis(discrete=TRUE, option = 'viridis') +
  guides(color = guide_legend("Location")) +
  coord_equal() +
  theme_half_open(12) + 
  theme(
    axis.title.x = element_text(size=axis.text.size), # remove x-axis labels
    axis.title.y = element_text(size=axis.text.size), # remove y-axis labels
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),  #remove major-grid labels
    panel.grid.minor = element_blank(),  #remove minor-grid labels
    plot.background = element_blank(),
    plot.margin = unit(c(.5, 0, 0, 0), "cm"),
    legend.title = element_text(size = legend.title.size),
    legend.text = element_text(size = legend.text.size)) +
  scale_shape_manual(breaks = c('Environment'),
                     values=c(1,16),
                     labels = c('Enviroment')) +
  xlab('NMDS 1') + ylab('NMDS 2') + xlim(-2,2) + ylim(-1.5,2) + labs(shape = 'Sample Type'); Loc.plot 




#Extra Unused Plots
axis.text.size <- 16
legend.title.size <- 16
legend.text.size <- 16

bray.plot <- ggplot() + 
  geom_point(data=frog2, aes(x=bray_x.nmds,y=bray_y.nmds,
                             color=species, shape = type),size=5) +
  scale_color_manual(breaks = c("GF", "PF", "Tad", "WF"),
                     values=c("Green", "red", "orange", "blue"),
                     labels = c('L. clamitans', 'L. palustris', 'Lithobates tadpole',
                                'L. sylvaticus')) + 
  labs(col = "Amphibian") +
  new_scale_color() +
  geom_point(data=frog1, aes(x=bray_x.nmds,y=bray_y.nmds,
                             color=species, shape = type),size=5) + 
  scale_color_viridis(discrete=TRUE, option = 'viridis',labels = c("Ann's Meadow",
                                                                   "Bog",
                                                                   "Tunnel","Stream",
                                                                   "Vernal Pool")) +
  coord_equal() +
  theme_half_open(12) + 
  scale_shape_manual(values=c(16, 1)) +
  theme(
    axis.title.x = element_text(size=axis.text.size),# remove x-axis labels
    axis.title.y = element_text(size=axis.text.size),# remove y-axis labels
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),#remove major-grid labels
    panel.grid.minor = element_blank(),#remove minor-grid labels
    plot.background = element_blank(),
    plot.margin = unit(c(.5, 0, 0, 0), "cm"),
    legend.title = element_text(size = legend.title.size),
    legend.text = element_text(size = legend.text.size)) +
  xlab('NMDS 1') + ylab('NMDS 2') + xlim(-2,2) + ylim(-1.5,2) + 
  labs(col="Environment", shape= NULL) + 
  guides(shape = FALSE, color = guide_legend(override.aes = list(shape = 1))); bray.plot


jac.plot <- ggplot() + 
  geom_point(data=frog2, aes(x=jac_x.nmds,y=jac_y.nmds,
                             color=species, shape = type),size=5) +
  scale_color_manual(breaks = c("GF", "PF", "Tad", "WF"),
                     values=c("Green", "red", "orange", "blue"),
                     labels = c('L. clamitans', 'L. palustris', 'Lithobates tadpole',
                                'L. sylvaticus')) + 
  labs(col = "Amphibian") +
  new_scale_color() +
  geom_point(data=frog1, aes(x=bray_x.nmds,y=bray_y.nmds,
                             color=species, shape = type),size=5) + 
  scale_color_viridis(discrete=TRUE, option = 'viridis',labels = c("Ann's Meadow",
                                                                   "Bog",
                                                                   "Tunnel","Stream",
                                                                   "Vernal Pool")) +
  coord_equal() +
  theme_half_open(12) + 
  scale_shape_manual(values=c(16, 1)) +
  theme(
    axis.title.x = element_text(size=axis.text.size),# remove x-axis labels
    axis.title.y = element_text(size=axis.text.size),# remove y-axis labels
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),#remove major-grid labels
    panel.grid.minor = element_blank(),#remove minor-grid labels
    plot.background = element_blank(),
    plot.margin = unit(c(.5, 0, 0, 0), "cm"),
    legend.title = element_text(size = legend.title.size),
    legend.text = element_text(size = legend.text.size)) +
  xlab('NMDS 1') + ylab('') + xlim(-2,2) + ylim(-1.5,2) + 
  labs(col="Environment", shape= NULL) + 
  guides(shape = FALSE, color = guide_legend(override.aes = list(shape = 1)))


frog.pf$infect <- as.factor(frog.pf$infect)
Inf.plot <- ggplot() + 
  geom_point(data=frog.pf, aes(x=bray_x.nmds,y=bray_y.nmds,
                               col = infect),size=5) + 
  coord_equal() +
  theme_half_open(12) + 
  guides(color=guide_legend("Infection Status")) +
  scale_color_manual(breaks = c(0, 1),
                     values=c("grey", "black"),
                     labels = c("Not Infected","Infected")) + 
  theme(
    axis.title.x = element_text(size=axis.text.size), # remove x-axis labels
    axis.title.y = element_text(size=axis.text.size), # remove y-axis labels
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),  #remove major-grid labels
    panel.grid.minor = element_blank(),  #remove minor-grid labels
    plot.background = element_blank(),
    plot.margin = unit(c(.5, 0, 0, 0), "cm"),
    legend.position = c(.8,.8),
    legend.title = element_text(size = legend.title.size),
    legend.text = element_text(size = legend.text.size)) +
  xlab('NMDS 1') + ylab('NMDS 2') + xlim(-2,2) + ylim(-1.5,2) 

Loc.plot <- ggplot() + 
  geom_point(data=frog, aes(x=bray_x.nmds,y=bray_y.nmds,
                            color = location, shape = type),size=5) + # add the point markers
  scale_color_viridis(discrete=TRUE, option = 'viridis') +
  guides(color = guide_legend("Location")) +
  coord_equal() +
  theme_half_open(12) + 
  theme(
    axis.ticks = element_blank(),  # remove axis ticks
    axis.title.x = element_text(size=axis.text.size), # remove x-axis labels
    axis.title.y = element_text(size=axis.text.size), # remove y-axis labels
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),  #remove major-grid labels
    panel.grid.minor = element_blank(),  #remove minor-grid labels
    plot.background = element_blank(),
    plot.margin = unit(c(.5, 0, 0, 0), "cm"),
    legend.title = element_text(size = legend.title.size),
    legend.text = element_text(size = legend.text.size)) +
  scale_shape_manual(breaks = c('Environment', 'Amphibian'),
                     values=c(1,16),
                     labels = c('Enviroment', 'Amphibian')) +
  xlab('NMDS 1') + ylab('NMDS 2') + xlim(-2,2) + ylim(-1.5,2) + labs(shape = 'Sample Type')  
