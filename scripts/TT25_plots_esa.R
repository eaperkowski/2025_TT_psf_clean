##############################################################################
# Prep
##############################################################################
# Libraries
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(MuMIn)
library(ggpubr)

# Read compiled dataset
df <- read.csv("../data/TT25_tri_compiled.csv")

# Remove sterile plant treatment from analyses
df_noSterile <- df %>%
  filter(ExpFungSource != "Wsterile" & ExpFungSource != "NWsterile") %>%
  mutate(FullExpTrt = gsub("-", "_", FullExpTrt),
         wet_root_shoot = ifelse(wet_root_shoot == "Inf", NA,
                                 ifelse(wet_root_shoot > 20, NA, wet_root_shoot)),
         plot_trt_comb = str_c(PlantGMTrt, "_", ExpFungSource))

gm.colors <- c("#F1B700", "#00B2BE")

##############################################################################
# Net photosynthesis
##############################################################################
# Model
vcmax_tri <- lmer(log(vcmax25) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                   wet_rhizome_mass_g + (1 | machine) + (1 | plot_field), 
                 data = df_noSterile)

# Check model result
Anova(vcmax_tri)

# Create compact letters
vcmax_tri_cld <- cld(
  emmeans(vcmax_tri, ~PlantGMTrt * ExpFungSource), 
  Letters = LETTERS, alpha = 0.1) %>%
  mutate(.group = trimws(.group, which = "both"))

# Create plot
png("../drafts/figs/evan_esa/TT25_vcmax_legacy_AM_plot.png", width = 6, height = 5,
    units = "in", res = 600)
ggplot(data = vcmax_tri_cld, 
       aes(x = PlantGMTrt, y = emmean, fill = ExpFungSource)) +
  geom_rect(aes(fill = PlantGMTrt),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = PlantGMTrt),
            xmin = 1.5, xmax = Inf,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_jitter(data = df_noSterile, 
             aes(x = PlantGMTrt, y = log(vcmax25), fill = ExpFungSource),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), alpha = 0.1) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(label = .group, y = 4.5), position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  scale_y_continuous(limits = c(0, 4.5), breaks = seq(0, 4, 1)) +
  labs(x = "Plant treatment legacy",
       y = expression(bold("ln ")*bolditalic("V")[bold("cmax25")]* bold(" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

##############################################################################
# Leaf d15N
##############################################################################
# Model
leaf15n_tri <- lmer(leaf_d15n ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                    wet_rhizome_mass_g + (1 | plot_field), 
                  data = df_noSterile)

# Check model result
Anova(leaf15n_tri)

# Create compact letters
d15n_tri_cld <- cld(
  emmeans(leaf15n_tri, ~PlantGMTrt * ExpFungSource , type = "response"), 
  Letters = LETTERS) %>%
  mutate(.group = trimws(.group, which = "both"))

# Create plot
png("../drafts/figs/evan_esa/TT25_d15n_legacy_AM_plot.png", width = 6, height = 5,
    units = "in", res = 600)
ggplot(data = d15n_tri_cld, 
       aes(x = PlantGMTrt, y = emmean, fill = ExpFungSource)) +
  geom_rect(aes(fill = PlantGMTrt),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = PlantGMTrt),
            xmin = 1.5, xmax = Inf,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_jitter(data = df_noSterile, 
              aes(x = PlantGMTrt, y = leaf_d15n, fill = ExpFungSource),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), alpha = 0.1) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(label = .group, y = 3), position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  scale_y_continuous(limits = c(0, 3.1), breaks = seq(0, 3, 1)) +
  labs(x = "Plant treatment legacy",
       y = expression(bold("Leaf "*delta^"15"*"N (‰)")),
       fill = "AM fungal source") +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

##############################################################################
# Leaf d15N - AM fungal length
##############################################################################
# Model
leaf15n_tri_amf <- lmer(leaf_d15n ~ ExpFungSource * AMF_asep_plant + (1 | plot_field), 
                        data = df_noSterile)
Anova(leaf15n_tri_amf)
test(emtrends(leaf15n_tri_amf, ~ExpFungSource, "AMF_asep_plant"))

leaf15n_tri_amf_regline <- data.frame(
  emmeans(leaf15n_tri_amf, ~ExpFungSource, "AMF_asep_plant",
          at = list(AMF_asep_plant = seq(0.25, 3.75, 0.01)))) %>%
  mutate(linetype = ifelse(ExpFungSource == "NW", "solid", "dashed"))

png("../drafts/figs/evan_esa/TT25_d15n_AMhyphae_plot.png", width = 8, height = 6,
    units = "in", res = 600)
leaf15n_tri_amf_plot <- ggplot(data = df_noSterile, 
       aes(x = AMF_asep_plant, y = leaf_d15n, fill = ExpFungSource)) +
  geom_point(shape = 21) +
  geom_ribbon(data = leaf15n_tri_amf_regline, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), alpha = 0.3) +
  geom_smooth(data = leaf15n_tri_amf_regline, 
              aes(y = emmean, color = ExpFungSource)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_color_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +  
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  scale_y_continuous(limits = c(-1, 4), breaks = seq(-1, 4, 1)) +
  labs(x = expression(bold("Aseptate AM hyphal length (m g"^"-1"*" dry soil)")),
       y = expression(bold("Leaf "*delta^"15"*"N (‰)")),
       fill = "AM fungal source", color = "AM fungal source")  +
  guides(linetype = "none") +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

##############################################################################
# Vcmax25 - AM fungal length
##############################################################################
# Model
vcmax25_tri_amf <- lmer(log(vcmax25) ~ ExpFungSource * AMF_asep_plant + (1 | machine) + (1 | plot_field), 
                        data = df_noSterile)
Anova(vcmax25_tri_amf)
test(emtrends(vcmax25_tri_amf, ~ExpFungSource, "AMF_asep_plant"))

vcmax25_tri_amf_regline <- data.frame(
  emmeans(vcmax25_tri_amf, ~ExpFungSource, "AMF_asep_plant",
          at = list(AMF_asep_plant = seq(0.25, 3.75, 0.01)))) %>%
  mutate(linetype = ifelse(ExpFungSource == "NW", "solid", "dashed"))

# Plot
png("../drafts/figs/evan_esa/TT25_vcmax_AMhyphae_plot.png", width = 8, height = 6,
    units = "in", res = 600)
vcmax25_tri_amf_plot <- ggplot(data = df_noSterile, 
       aes(x = AMF_asep_plant, y = log(vcmax25), fill = ExpFungSource)) +
  geom_point(shape = 21) +
  geom_ribbon(data = vcmax25_tri_amf_regline, 
              aes(y = emmean, ymin = lower.CL, ymax = upper.CL), alpha = 0.3) +
  geom_smooth(data = vcmax25_tri_amf_regline, 
              aes(y = emmean, linetype = linetype, color = ExpFungSource)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_color_manual(values = gm.colors, 
                     labels = c("ambient", "weeded")) +  
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  scale_y_continuous(limits = c(-1, 4.5), breaks = seq(-1, 4, 1)) +
  labs(x = expression(bold("Aseptate AM hyphal length (m g"^"-1"*" dry soil)")),
       y = expression(bold("ln ")*bolditalic("V")[bold("cmax25")]* bold(" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source", color = "AM fungal source")  +
  guides(linetype = "none") +
  theme_classic(base_size = 20) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()


png("../drafts/figs/evan_esa/TT25_d15n_vcmax_hyphae_plots.png",
    height = 12, width = 7, units = "in", res = 600)
ggarrange(leaf15n_tri_amf_plot, vcmax25_tri_amf_plot, nrow = 2,
          common.legend = TRUE, legend = "bottom")
dev.off()

