# This script performs all calculations and reproduces all figures from the manuscript and Supplementary Material of:
# Buchwald SZ, Birgel D, Senger K, Mosociova T, Pei Y, Zuchuat V, Tarhan L, Frank AB, Galasso F,
# Gómez Correa MA, Koşun E, Karapunar B, Wang X, Kustatscher E, Prinoth H, Lahajnar N, Steinkrauss R,
# Peckmann J, and Foster WJ
# "Phytoplankton blooms on the Barents Shelf, Svalbard, associated with the Permian–Triassic mass extinction"

# This script and the script "Extinction_marker_functions.R" need to be in the working directory to be properly loaded.
# Create a sub-folder titled "Raw_data" in the working directory, which should contain all raw and metadata, including
# "Extinction_Marker_Metadata.xlsx", "MPI_raw.xlsx", "Nabbefeld.xlsx", and the folders "area_mol_sieve"
# and "area_n_alk", in which the integrated peak areas of the compounds of interest are saved in an
# individual .txt-files per sample.

# When initiating the project, a folder "Output" is created in the working directory,
# that will contain all data produced.

setwd("") # set your working directory
source('Extinction_marker_functions.R')
init_project()

############################### DATA IMPORT ####################################
# post molecular sieve
a_mol <- read_raw_data("Raw_data/area_mol_sieve")
q_mol_ug_TOC <- quantify_compounds(a_mol, "meta_mol")

################################## PLOTS #######################################
meta_mol <- read.xlsx("Raw_data/Extinction_Marker_Metadata.xlsx", sheet = "meta_mol") # import metadata for post molecular sieve

sample_order <- as.data.frame(select(meta_mol, sample, log_height)) # log height
mol <- merge(q_mol_ug_TOC, sample_order, by = "sample")

mol <- mol[order(mol$log_height), ] # order by log height
write.xlsx(mol, file = "Output/FID_mol_quantified_ug_TOC_Chol.xlsx")

############ C33-n-ACH ############ 
ACH33 <- as.data.frame(select(mol, log_height, "C33-ACH"))

ACH33_plot <- melt(ACH33,  id.vars = 'log_height', variable.name = 'compound') #long format
names(ACH33_plot)
names(ACH33_plot) <- c("log_height", "compound", "mass")

plot_33ACH <- ggplot(ACH33_plot, aes(log_height, mass)) +
  geom_point(aes(colour = compound), pch = 1, size = 2.5, show.legend = FALSE) +
  geom_line(aes(colour = compound), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(y_name = expression(paste(mu, "g/g TOC")),
                         x_name = "height [m]") +
  scale_x_break(c(-25, -3), ticklabels = seq(-5, 30, by = 5)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90))
plot_33ACH
ggsave(plot_33ACH, file = "Output/C33-nACH.pdf", width = 8, height = 15, units = "cm")

############ All n-ACHs ############ 
all_ACHs <- as.data.frame(mol[,c("log_height", colnames(mol)[grep("ACH",colnames(mol))])]) # select all columns with "ACH" in col_name
all_ACHs <- all_ACHs[, !names(all_ACHs) %in% "C33-ACH"] # remove C33-ACH column

all_ACHs[is.na(all_ACHs)] <- 0 # replace NAs by 0

############ Sum all n-ACHs ############ 
log_height <- all_ACHs$log_height
all_ACHs <- subset(all_ACHs, select = -log_height) # remove log_height
all_ACHs_sum <- data.frame(log_height, rowSums(all_ACHs))

ACH_sum_plot <- melt(all_ACHs_sum,  id.vars = 'log_height', variable.name = 'compound') # long format
names(ACH_sum_plot)
names(ACH_sum_plot) <- c("log_height", "compound", "mass")

plot_ACH_sum <- ggplot(ACH_sum_plot, aes(log_height, mass)) +
  geom_point(aes(colour = compound), pch = 1, size = 2.5, show.legend = FALSE) +
  geom_line(aes(colour = compound), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(y_name = expression(paste(mu, "g/g TOC")),
                         x_name = "height [m]") +
  scale_x_break(c(-25, -3), ticklabels = seq(-5, 30, by = 5)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90))
plot_ACH_sum
ggsave(plot_ACH_sum, file = "Output/ACHs_sum.pdf", width = 8, height = 15, units = "cm")

############ Pristane and phytane in one panel ############ 
Pr_Ph <- as.data.frame(select(mol, log_height, Pristane, Phytane))

Pr_Ph_plot <- melt(Pr_Ph,  id.vars = 'log_height', variable.name = 'compound') # long format
names(Pr_Ph_plot)
names(Pr_Ph_plot) <- c("log_height", "compound", "mass")

plot_PrPh <- ggplot(Pr_Ph_plot, aes(log_height, mass)) +
  geom_point(aes(colour = compound, shape = compound), size = 2.5, show.legend = FALSE) +
  geom_line(aes(colour = compound), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(name="",values = c("olivedrab4", "olivedrab2"))+
  plot_common_parameters(y_name = expression(paste(mu, "g/g TOC")),
                         x_name = "height [m]")+
  scale_x_break(c(-25, -3), ticklabels = seq(-5, 30, by = 5)) +
  scale_shape_manual(values = c(4, 1)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90))
plot_PrPh
ggsave(plot_PrPh, file = "Output/Pristane_Phytane.pdf", width = 8, height = 15, units = "cm")

############ Phytanyl toluene ############ 
phyt <- as.data.frame(select(mol, log_height, "phytanyltoluene"))

phyt_plot <- melt(phyt,  id.vars = 'log_height', variable.name = 'compound') # long format
phyt_plot[is.na(phyt_plot)] <- 0 # replace NAs by 0
names(phyt_plot)
names(phyt_plot) <- c("log_height", "compound", "mass")

plot_phyt <- ggplot(phyt_plot, aes(log_height, mass)) +
  geom_point(aes(colour = compound), pch = 1, size = 2.5, show.legend = FALSE) +
  geom_line(aes(colour = compound), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(y_name = expression(paste(mu, "g/g TOC")),
                         x_name = "height [m]") +
  scale_x_break(c(-25, -3), ticklabels = seq(-5, 30, by = 5)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90))
plot_phyt
ggsave(plot_phyt, file = "Output/Phytanyl_toluene.pdf", width = 8, height = 15, units = "cm")

############ Maturity ############
MPI_raw_data <- read.xlsx("Raw_data/MPI_raw.xlsx", sheet = "F2")
class(MPI_raw_data$location)
MPI_raw_data$location <- as.factor(MPI_raw_data$location)
MPI_raw_data$section <- as.factor(MPI_raw_data$section)

# methylphenanthrene index:
# calculated after Cassini et al., 1988
# MPI = (1.89 * (2-MP + 3-MP))/(phenanthrene + 1.26 * (1-MP + 9-MP))
# 1-MP = 1-methylphenanthrene (from m/z = 192)
# 2-MP = 2-methylphenanthrene (from m/z = 192)
# 3-MP = 3-methylphenanthrene (from m/z = 192)
# 9-MP = 9-methylphenanthrene (from m/z = 192)
# phenanthrene from m/z = 178

MPI_raw_data <- na.omit(MPI_raw_data) # delete rows with NA (otherwise MPI cannot be calculated)

MPI <- data.frame(MPI_raw_data$location,
                  MPI_raw_data$section)
names(MPI) <- c("location", "section")
MPI$MPI <- (1.89 * (MPI_raw_data$'2MP' + MPI_raw_data$'3MP'))/(MPI_raw_data$phenanthrene + 1.26 * (MPI_raw_data$'1MP' + MPI_raw_data$'9MP'))

MPI$section <- with(MPI, factor(section, levels = unique(section))) # to have order as in data frame later

col <- paletteer_d("MetBrewer::Egypt") # define color palette

MPI_plot <- ggplot(MPI, aes(x = section, y = MPI, colour = location)) +
  geom_boxplot(lwd = 0.3) +
  scale_y_continuous(name = "methylphenanthrene index (MPI)",
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25),
                     limits = c(0, 1.35)) +
  xlab("") +
  scale_color_manual(values = col) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.85, 0.86),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))
MPI_plot
ggsave(MPI_plot, file = "Output/MPI_boxplot.pdf", width = 15, height = 12, units = "cm")

############################## SUPPLEMENT #####################################
# SUPPLEMENTARY FIG. S1
#TOC
plot_TOC <- ggplot(meta_mol, aes(log_height, TOC)) +
  geom_point(aes(colour = "grey13"), pch = 1, size = 2.5, stroke = 1, show.legend = FALSE) +
  geom_line(aes(colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(y_name = expression(paste("TOC (wt %)")),
                         x_name = "height (m)")+
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90),
        legend.position = "right")

# MPI
MPI_Lusitaniadalen <- MPI_raw_data[MPI_raw_data$section %in% "Lusitaniadalen",]
MPI_Lusitaniadalen$MPI <- (1.89 * (MPI_Lusitaniadalen$'2MP' + MPI_Lusitaniadalen$'3MP'))/(MPI_Lusitaniadalen$phenanthrene + 1.26 * (MPI_Lusitaniadalen$'1MP' + MPI_Lusitaniadalen$'9MP'))
MPI_Lusitaniadalen <- merge(MPI_Lusitaniadalen, sample_order, by = "sample")

plot_MPI_Lusitaniadalen <- ggplot(MPI_Lusitaniadalen, aes(log_height, MPI)) +
  geom_point(aes(colour = "grey13"), pch = 1, size = 2.5, stroke = 1, show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(y_name = expression(paste(mu, "g"["compound"]," / g"["TOC"])),
                         x_name = "height [m]")+
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90),
        legend.position = "right") +
  ylim(c(0.5, 1))

# all n-alkanes (quantified)
a_nalk <- read_raw_data("Raw_data/area_n_alk") # import n-alkane raw data
q_nalk_ug_TOC <- quantify_compounds(a_nalk, "meta_nalk") # quantify n-alkanes

meta_nalk <- read.xlsx("Raw_data/Extinction_Marker_Metadata.xlsx", sheet = "meta_nalk") # import metadata for n-alkanes

sample_order <- as.data.frame(select(meta_nalk, sample, log_height)) # get df with sample ID and log height
nalk <- merge(q_nalk_ug_TOC, sample_order, by = "sample")

nalk <- nalk[order(nalk$log_height), ] # order by log height
write.xlsx(nalk, file = "Output/FID_nalk_quantified_ug_TOC_Chol.xlsx")

# plot
names(nalk)
nalk <- nalk[, !names(nalk) %in% c("sample", "Cholestane")] # removes columns
nalk[is.na(nalk)] <- 0 # replace NAs by 0

nalk_plot <- melt(nalk,  id.vars = 'log_height', variable.name = 'compound') # long format
names(nalk_plot)
names(nalk_plot) <- c("log_height", "compound", "mass")

col <- paletteer_d("ggthemes::manyeys") # define color palette

plot_nalk <- ggplot(nalk_plot, aes(log_height, mass)) +
  geom_point(aes(colour = compound), pch = 1, size = 2.5, show.legend = FALSE) +
  plot_common_parameters(y_name = expression(paste(mu, "g"["compound"]," / g"["TOC"])),
                         x_name = "height [m]")+
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90),
        legend.position = "right")

nC34 <- as.data.frame(cbind(log_height = nalk$log_height,
                            nC34 = nalk$`n-C34`))
nC34 <- merge(nC34, ACH33, by = "log_height")
nC34$ACH33_nC34 <- nC34$`C33-ACH`/nC34$nC34

plot_ACH33_nC34 <- ggplot(nC34, aes(log_height, ACH33_nC34)) +
  geom_point(aes(colour = "grey13"), pch = 1, size = 2.5, stroke = 1, show.legend = FALSE) +
  geom_line(aes(colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(y_name = expression(paste("C"[33], "-", italic("n"), "-ACH/", italic("n"), "-C"[34])),
                         x_name = "height (m)")+
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90),
        legend.position = "right")

combi_nalk <- ggarrange(plot_TOC, plot_MPI_Lusitaniadalen, plot_nalk, plot_ACH33_nC34,
                        align = "hv",
                        nrow = 1, labels = c("a","b", "c"))
combi_nalk
ggsave(combi_nalk, file = "Output/Supplementary_Figure1.pdf", width = 25, height = 15, units = "cm")

# SUPPLEMENTARY FIG. S2
MPI_F1_raw_data <- read.xlsx("Raw_data/MPI_raw.xlsx", sheet = "F1")
class(MPI_F1_raw_data$location)
MPI_F1_raw_data$location <- as.factor(MPI_F1_raw_data$location)
MPI_F1_raw_data$section <- as.factor(MPI_F1_raw_data$section)

# methylphenanthrene index:
# calculated after Cassini et al., 1988
# MPI = (1.89 * (2-MP + 3-MP))/(phenanthrene + 1.26 * (1-MP + 9-MP))

MPI_F1_raw_data <- na.omit(MPI_F1_raw_data) # delete rows with NA (otherwise MPI cannot be calculated)

MPI_F1 <- data.frame(MPI_F1_raw_data$location,
                     MPI_F1_raw_data$section)
names(MPI_F1) <- c("location", "section")
MPI_F1$MPI <- (1.89 * (MPI_F1_raw_data$'2MP' + MPI_F1_raw_data$'3MP'))/(MPI_F1_raw_data$phenanthrene + 1.26 * (MPI_F1_raw_data$'1MP' + MPI_F1_raw_data$'9MP'))

MPI_F1$section <- with(MPI_F1, factor(section, levels = unique(section))) # to have order as in data frame later

col <- paletteer_d("MetBrewer::Egypt") # define color palette

MPI_F1_plot <- ggplot(MPI_F1, aes(x = section, y = MPI, colour = location)) +
  geom_boxplot(lwd = 0.3) +
  scale_y_continuous(name = "methylphenanthrene index (MPI)",
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25),
                     limits = c(0, 1.35)) +
  xlab("") +
  scale_color_manual(values = col) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.85, 0.86),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))
MPI_F1_plot
ggsave(MPI_plot, file = "Output/Supplementary_Figure2.pdf", width = 15, height = 12, units = "cm")

# SUPPLEMENTARY FIG. S3
#correlations of content of all ACHs
all_ACHs <- as.data.frame(mol[,c("log_height", colnames(mol)[grep("ACH",colnames(mol))])]) # select all columns with "ACH" in col_name
names(all_ACHs)
all_ACHs[is.na(all_ACHs)] <- 0 # replace NAs by 0

all_ACHs_pre <- subset(all_ACHs, all_ACHs$log_height < 0) # select pre-extinction samples only
all_ACHs_pre <- subset(all_ACHs_pre, select = -log_height)

all_ACHs_post <- subset(all_ACHs, all_ACHs$log_height > 0) # select post-extinction samples only
all_ACHs_post <- subset(all_ACHs_post, select = -log_height) 

all_ACHs_pre_long <- melt(all_ACHs_pre,  id.vars = 'C33-ACH', variable.name = 'compound') # long format
all_ACHs_post_long <- melt(all_ACHs_post,  id.vars = 'C33-ACH', variable.name = 'compound') # long format

col <- paletteer::paletteer_d("ggsci::category20b_d3") # define color palette
scatterplot_33ACH_ACHs_pre <- ggplot(all_ACHs_pre_long, aes(`C33-ACH`, value, colour = compound)) +
  geom_smooth(formula = y ~ x, method = lm, se = FALSE, show.legend = FALSE) +
  geom_point(aes(colour = compound), pch = 16, size = 1.5, show.legend = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, show.legend = FALSE) +
  xlab(expression(paste("C"[33], "-", italic("n"), "-ACH (", mu, "g/g TOC)"))) +
  ylab(expression(paste(italic("n"), "-ACHs (", mu, "g/g TOC)"))) +
  scale_color_manual(values = col) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

scatterplot_33ACH_ACHs_post <- ggplot(all_ACHs_post_long, aes(`C33-ACH`, value, colour = compound)) +
  geom_smooth(formula = y ~ x, method = lm, se = FALSE, show.legend = FALSE) +
  geom_point(aes(colour = compound), pch = 16, size = 1.5, show.legend = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, show.legend = FALSE) +
  xlab(expression(paste("C"[33], "-", italic("n"), "-ACH (", mu, "g/g TOC)"))) +
  ylab(expression(paste(italic("n"), "-ACHs (", mu, "g/g TOC)"))) +
  scale_color_manual(values = col) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# pristane and phytane vs. C33-n-ACH 
C33ACH_Pr_Ph <- as.data.frame(cbind('C33-ACH' = q_mol_ug_TOC$`C33-ACH`,
                                    Pristane = q_mol_ug_TOC$Pristane,
                                    Phytane = q_mol_ug_TOC$Phytane))
C33ACH_Pr_Ph_long <- melt(C33ACH_Pr_Ph,  id.vars = 'C33-ACH', variable.name = 'compound') #long format

col <- c("yellowgreen", "seagreen")
scatterplot_33ACH_pristane <- ggplot(C33ACH_Pr_Ph_long, aes(`C33-ACH`, value, color = compound)) +
  geom_smooth(formula = y ~ x, method = lm, show.legend = FALSE) +
  geom_point(aes(colour = compound, shape = compound), size = 1.5, stroke = 1.5, show.legend = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, show.legend = FALSE) +
  xlab(expression(paste("C"[33], "-", italic("n"), "-ACH (", mu, "g/g TOC)"))) +
  ylab(expression(paste("Pr or Ph (", mu, "g/g TOC)"))) +
  scale_color_manual(values = col) +
  scale_shape_manual(values = c(4, 16)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# phytanyl toluene vs. C33-n-ACH 
scatterplot_33ACH_phyt <- ggplot(q_mol_ug_TOC, aes(`C33-ACH`, phytanyltoluene)) +
  geom_smooth(formula = y ~ x, method = lm, color = "grey33") +
  geom_point(aes(colour = "black"), pch = 16, size = 1.5, show.legend = FALSE) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) +
  xlab(expression(paste("C"[33], "-", italic("n"), "-ACH (", mu, "g/g TOC)"))) +
  ylab(expression(paste("phytanyl toluene (", mu, "g/g TOC)"))) +
  scale_color_manual(values = "grey13") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

combi_plot <- ggarrange(scatterplot_33ACH_ACHs_pre, scatterplot_33ACH_ACHs_post,
                        scatterplot_33ACH_phyt, scatterplot_33ACH_pristane,
                        nrow = 2, ncol = 2, align = "hv", labels = c("a","b", "c", "d"))
combi_plot
ggsave(combi_plot, file = "Output/Supplementary_Figure3.pdf", width = 20, height = 20, units = "cm")

# SUPPLEMENTARY FIG. S4
Nabbefeld <- read.xlsx("Raw_data/Nabbefeld.xlsx")
Nabbefeld$log_height <- Nabbefeld$`log_height_Nabbefeld`-16.90 # unify log_height
Nabbefeld$`C33-ACH` <- Nabbefeld$`C33-ACH`/1000 # from n/g TOC to ug/g TOC
Nabbefeld$pristane <- Nabbefeld$pristane/1000 # from n/g TOC to ug/g TOC
Nabbefeld$phytane <- Nabbefeld$phytane/1000 # from n/g TOC to ug/g TOC

plot_C33ACH_Nabbefeld <- ggplot(Nabbefeld, aes(log_height, `C33-ACH`)) +
  geom_point(aes(colour = "deepskyblue"), pch = 16, size = 2.5, stroke = 1, show.legend = FALSE) +
  geom_point(data = ACH33, aes(log_height, `C33-ACH`, colour = "grey13"), pch = 1, size = 2.5, stroke = 1, show.legend = FALSE) +
  scale_color_manual(values = c("deepskyblue", "grey13")) +
  plot_common_parameters(y_name = expression(paste(mu, "g/g TOC")),
                         x_name = "height (m)") +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90)) +
  xlim(c(-2, 1.2)) + ylim(c(0, 320))

plot_TOC_Nabbefeld <- ggplot(Nabbefeld, aes(log_height, TOC)) +
  geom_point(aes(colour = "deepskyblue"), pch = 16, size = 2.5, stroke = 1, show.legend = FALSE) +
  geom_point(data = meta_mol, aes(log_height, TOC, colour = "grey13"), pch = 1, size = 2.5, stroke = 1, show.legend = FALSE) +
  scale_color_manual(values = c("deepskyblue", "grey13")) +
  plot_common_parameters(y_name = expression(paste("TOC (wt%)")),
                         x_name = "height (m)") +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90)) +
  xlim(c(-2, 1.2)) + ylim(c(0, 1))

Pr_Ph_Nabbefeld <- as.data.frame(cbind(log_height = Nabbefeld$log_height,
                                       pristane = Nabbefeld$pristane,
                                       phytane = Nabbefeld$phytane))
Pr_Ph_Nabbefeld_long <- melt(Pr_Ph_Nabbefeld,  id.vars = 'log_height', variable.name = 'compound') #long format

plot_PrPh_Nabbefeld <- ggplot(Pr_Ph_Nabbefeld_long, aes(log_height, value)) +
  geom_point(aes(colour = compound, shape = compound), size = 2.5, stroke = 1, show.legend = FALSE) +
  geom_point(data = Pr_Ph_plot, aes(log_height, mass, colour = compound, shape = compound), size = 2.5, stroke = 1, show.legend = FALSE) +
  scale_color_manual(name="",values = c("olivedrab2", "olivedrab2", "olivedrab4", "olivedrab4"))+
  plot_common_parameters(y_name = expression(paste(mu, "g/g TOC")),
                         x_name = "height (m)")+
  scale_shape_manual(values = c(1, 16, 2, 4)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90)) +
  xlim(c(-2, 1.2)) + ylim(c(0, 1250))

combi_Nabbefeld <- ggarrange(plot_TOC_Nabbefeld, plot_C33ACH_Nabbefeld,
                             plot_PrPh_Nabbefeld,
                        nrow = 1, labels = c("a", "b", "c"))
combi_Nabbefeld
ggsave(combi_Nabbefeld, file = "Output/Supplementary_Figure4.pdf", width = 20, height = 15, units = "cm")

# SUPPLEMENTARY FIG. S5
log_height_nalk <- nalk$log_height
nalk <- nalk[, !names(nalk) %in% "log_height"] # removes column
nalk_rel <- (nalk / rowSums(nalk))*100 # relative abundance
nalk_rel$log_height <- log_height_nalk

plot_nC33_rel <- ggplot(nalk_rel, aes(log_height, `n-C33`)) +
  geom_point(aes(colour = "grey13"), pch = 1, size = 2.5, stroke = 1, show.legend = FALSE) +
  geom_line(aes(colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(y_name = expression(paste(mu, "g"["compound"]," / g"["TOC"])),
                         x_name = "height [m]")+
  scale_x_break(c(-25, -3), ticklabels = seq(-5, 30, by = 5)) +
  theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y = element_text( angle = 90),
        legend.position = "right")

ACH_nalk_rel <- merge(ACH33, nalk_rel, by = "log_height")
ACH_nalk_rel <- as.data.frame(cbind(ACH33 = ACH_nalk_rel$`C33-ACH`,
                                    nC33 = ACH_nalk_rel$`n-C33`,
                                    nC34 = ACH_nalk_rel$`n-C34`))

scatterplot_ACH_nC33_rel <- ggplot(ACH_nalk_rel, aes(ACH33, nC33)) +
  geom_smooth(formula = y ~ x, method = lm, color = "grey13",  show.legend = FALSE) +
  geom_point(aes(colour = "black"), pch = 16, size = 1.5, show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, show.legend = FALSE) +
  xlab(expression(paste("C"[33], "-", italic("n"), "-ACH (", mu, "g/g TOC)"))) +
  ylab(expression(paste("relative ", italic("n"), "-C"[33], " (%)"))) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

scatterplot_ACH_nC34_rel <- ggplot(ACH_nalk_rel, aes(ACH33, nC34)) +
  geom_smooth(formula = y ~ x, method = lm, color = "grey13",  show.legend = FALSE) +
  geom_point(aes(colour = "black"), pch = 16, size = 1.5, show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, show.legend = FALSE) +
  xlab(expression(paste("C"[33], "-", italic("n"), "-ACH (", mu, "g/g TOC)"))) +
  ylab(expression(paste("relative ", italic("n"), "-C"[34], " (%)"))) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

combi_nalk_rel_scatterplot <- ggarrange(scatterplot_ACH_nC33_rel, scatterplot_ACH_nC34_rel,
                                        nrow = 1, labels = c("a","b"))
combi_nalk_rel_scatterplot
ggsave(combi_nalk_rel_scatterplot, file = "Output/Supplementary_Figure5.pdf", width = 20, height = 10, units = "cm")
