# This script performs all calculations and reproduces all figures from the manuscript:
# Buchwald SZ, Birgel D, Senger K, Mosočiová T, Pei Y, Frank AB, Galasso F, Gómez Correa MA,
# Koşun E, Karapunar B, Wang X, Kustatscher E, Prinoth H, Steinkrauss R, Peckmann J, and Foster WJ
# "Primary productivity blooms on the Barents Shelf, Svalbard, associated with the Permian–Triassic mass extinction"

# This script and the script "Extinction_marker_functions.R" need to be in the working directory to be properly loaded.
# Create a sub-folder titled "Raw_data" in the working directory, which should contain all raw and metadata, including
# "Extinction_marker_functions.xlsx", "MPI_raw.xlsx" and the folder "area_molar_sieve", in which the integrated
# peak areas of the compounds of interest are saved in an individual .txt-files per sample.

# When initiating the project, a folder "Output" is created in the working directory,
# that will contain all data produced.

setwd("") # set your working directory
source('Extinction_Marker_functions.R')
init_project()

############################### DATA IMPORT ####################################
# post molar sieve
a_mol <- read_raw_data("Raw_data/area_mol_sieve")
q_mol_ug_TOC <- quantify_compounds(a_mol, "meta_mol")

################################## PLOTS #######################################
meta_mol <- read.xlsx("Raw_data/Extinction_Marker_Metadata.xlsx", sheet = "meta_mol") # import metadata for post molar sieve

sample_order <- as.data.frame(select(meta_mol, sample, log_height)) # get df with sample ID and log height
mol <- merge(q_mol_ug_TOC, sample_order, by = "sample")

mol <- mol[order(mol$log_height), ] # order by log height
write.xlsx(mol, file = "Output/FID_mol_quantified_ug_TOC_Chol.xlsx")

############ C33-n-ACH ############ 
ACH33 <- as.data.frame(select(mol, log_height, "C33-ACH"))

ACH33_plot <- melt(ACH33,  id.vars = 'log_height', variable.name = 'compound') #melt data frame into long format
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
all_ACHs <- as.data.frame(mol[,c("log_height", colnames(mol)[grep("ACH",colnames(mol))])]) # select all columns with "ACH" in colname
all_ACHs <- all_ACHs[, !names(all_ACHs) %in% "C33-ACH"] # remove C33-ACH column

all_ACHs[is.na(all_ACHs)] <- 0 # replace NAs by 0

############ Sum all n-ACHs ############ 
log_height <- all_ACHs$log_height
all_ACHs <- subset(all_ACHs, select = -log_height) # remove log_height
all_ACHs_sum <- data.frame(log_height, rowSums(all_ACHs))

ACH_sum_plot <- melt(all_ACHs_sum,  id.vars = 'log_height', variable.name = 'compound') #melt data frame into long format
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

Pr_Ph_plot <- melt(Pr_Ph,  id.vars = 'log_height', variable.name = 'compound') # melt data frame into long format
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

phyt_plot <- melt(phyt,  id.vars = 'log_height', variable.name = 'compound') #melt data frame into long format
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

############################## CORRELATIONS ####################################

############ correlations of content of all ACHs ############
all_ACHs <- as.data.frame(mol[,c("log_height", colnames(mol)[grep("ACH",colnames(mol))])]) # select all columns with "ACH" in colname
names(all_ACHs)
all_ACHs <- subset(all_ACHs, select = -log_height) # remove log_height
all_ACHs[is.na(all_ACHs)] <- 0 # replace NAs by 0

ACH_corr <- rcorr(as.matrix(all_ACHs), type = "pearson")
ACH_corr$r # correlation coefficients
ACH_corr$P #p-values

ACH_corr_mat <- as.matrix(ACH_corr$r)
ACH_sign_mat <- as.matrix(ACH_corr$P)

corrplot(ACH_corr_mat, tl.col = "black", method = "circle", # methods: "circle", color", "number", "shade"
         type = "upper")

############ correlation all n-ACHs pre- vs. post-extinction ############
# pre-extinction
names(all_ACHs)
all_ACHs$log_height <- mol$log_height

all_ACHs_pre <- subset(all_ACHs, all_ACHs$log_height < 0) # select pre-extinction samples only
all_ACHs_pre <- subset(all_ACHs_pre, select = -log_height) 

ACH_pre_corr <- rcorr(as.matrix(all_ACHs_pre), type = "pearson")

ACH_pre_corr_mat <- as.matrix(ACH_pre_corr$r) # correlation coefficients
ACH_pre_sign_mat <- as.matrix(ACH_pre_corr$P) # p-values

# post-extinction
all_ACHs$log_height <- mol$log_height

all_ACHs_post <- subset(all_ACHs, all_ACHs$log_height > 0) # select post-extinction samples only
all_ACHs_post <- subset(all_ACHs_post, select = -log_height) 

ACH_post_corr <- rcorr(as.matrix(all_ACHs_post), type = "pearson")

ACH_post_corr_mat <- as.matrix(ACH_post_corr$r) # correlation coefficients
ACH_post_sign_mat <- as.matrix(ACH_post_corr$P) # p-values

############ correlations including other markers ############
cor.test(all_ACHs$`C33-ACH`, mol$Pristane) # rho = 0.9096739; p = 3.369e-12
cor.test(all_ACHs$`C33-ACH`, mol$Phytane) # rho = 0.9029993; p = 8.761e-12

cor.test(all_ACHs$`C33-ACH`, mol$phytanyltoluene)
cor.test(mol$phytanyltoluene, mol$Pristane)
cor.test(mol$phytanyltoluene, mol$Phytane)

############ Maturity ############
MPI_raw_data <- read.xlsx("Raw_data/MPI_raw.xlsx")
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

MPI_raw_data <- na.omit(MPI_raw_data) # delete rows with NA (MPI cannot be calculated)

MPI <- data.frame(MPI_raw_data$location,
                  MPI_raw_data$section)
names(MPI) <- c("location", "section")
MPI$MPI <- (1.89 * (MPI_raw_data$'2MP' + MPI_raw_data$'3MP'))/(MPI_raw_data$phenanthrene + 1.26 * (MPI_raw_data$'1MP' + MPI_raw_data$'9MP'))

MPI$section <- with(MPI, factor(section, levels = unique(section))) # to have order as in data frame later

col <- paletteer_d("MetBrewer::Egypt") # define color palette

MPI_plot<- ggplot(MPI, aes(x = section, y = MPI, colour = location)) +
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
