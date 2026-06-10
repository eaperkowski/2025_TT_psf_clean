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
df <- read.csv("../../data/2025_2026/TT25_tri_compiled.csv")

# Remove sterile plant treatment from analyses
df_noSterile <- df %>%
  filter(ExpFungSource != "Wsterile" & ExpFungSource != "NWsterile") %>%
  mutate(FullExpTrt = gsub("-", "_", FullExpTrt),
         wet_root_shoot = ifelse(wet_root_shoot == "Inf", NA,
                                 ifelse(wet_root_shoot > 20, NA, wet_root_shoot)))

# Some plot aesthetics
gm.colors <- c("#F1B700", "#00B2BE")
facet_lab <- c("Plant history: ambient", "Plant history: weeded")
names(facet_lab) <- c("NW", "W")

###############
# Anet
###############
df_noSterile$anet[89] <- NA

anet_tri <- lmer(log(anet + 1) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                   wet_rhizome_mass_g + (1 | machine) + (1 | plot_field), 
                 data = df_noSterile)

# Check normality assumptions
plot(anet_tri)
qqnorm(residuals(anet_tri))
qqline(residuals(anet_tri))
hist(residuals(anet_tri))
shapiro.test(residuals(anet_tri))
outlierTest(anet_tri)

# Model output
summary(anet_tri)
Anova(anet_tri)
r.squaredGLMM(anet_tri)

# Pairwise comparisons
cld(emmeans(anet_tri, pairwise~ExpSoilSource*ExpFungSource))
## Weeded plants inoculated with AM fungi from ambient plots exhibited a 
## significant reduction in net photosynthesis compared to weeded plants
## inoculated with AM fungi from weeded plots (home-field advantage?)

cld(emmeans(anet_tri, pairwise~PlantGMTrt*ExpFungSource))
## Weeded fungal community increased Anet by 150% in plants growing in ambient
## plots, while AMF decreased Anet by 53% in plants growing in weeded plots


## Plants with history of growing in weeded plot experience significant
## net photosynthesis boost when inoculated with ambient AMF; no fungal
## 

# Plot prep (full)
anet_tri_results <- cld(emmeans(anet_tri, 
                                ~PlantGMTrt*ExpSoilSource*ExpFungSource, 
                                type = "response"), 
                        Letters = LETTERS, reversed = TRUE, alpha = 0.1) %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant", PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant", PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded"))) %>% 
  data.frame()

# Plot prep (ExpSoilSource * ExpFungSource interaction)
anet_soilFung_int <- cld(emmeans(anet_tri, ~ExpSoilSource*ExpFungSource), 
                                 Letters = LETTERS, reversed = TRUE) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()


## Full plot
png("../../drafts/figs/TT25_tri_anet_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = anet_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource, y = emmean,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, y = emmean,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 3, 
                group = full_trt, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("A"["net"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

## Soil source x AM fungal source interaction plot
tri_anet_soilAMFint_plot <-
  ggplot(data = anet_soilFung_int,
         x = ExpSoilSource, 
         y = anet) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource, y = emmean,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, y = emmean,
                 fill = ExpFungSource),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpSoilSource, y = 2, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("A"["net"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")

###############
# gsw
###############
gsw_tri <- lmer(log(gsw) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                  wet_rhizome_mass_g + (1 | machine) + (1 | plot_field), 
                data = filter(df_noSterile, gsw > 0))

# Check normality assumptions
plot(gsw_tri)
qqnorm(residuals(gsw_tri))
qqline(residuals(gsw_tri))
hist(residuals(gsw_tri))
shapiro.test(residuals(gsw_tri))
outlierTest(gsw_tri)

# Model output
summary(gsw_tri)
Anova(gsw_tri)
r.squaredGLMM(gsw_tri)

# Pairwise comparisons
cld(emmeans(gsw_tri, pairwise~ExpSoilSource*ExpFungSource, type = "response"), alpha = 0.1)
## Plants growing in non-weeded soil have greater stomatal conductance than
## plants growing in weeded soil; however, this response is only observed
## when plants are inoculated with nonweeded AMF community

# Plot prep (full)
gsw_tri_results <- cld(emmeans(gsw_tri, 
                               ~PlantGMTrt*ExpSoilSource*ExpFungSource, 
                               type = "response"), 
                       Letters = LETTERS, reversed = TRUE, alpha = 0.1) %>%
  data.frame() %>% 
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant", PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant", PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded")))

# Plot prep (ExpSoilSource * ExpFungSource interaction)
gsw_soilFung_int <- cld(emmeans(gsw_tri, ~ExpSoilSource*ExpFungSource, 
                                 type = "response"), 
                         Letters = LETTERS, reversed = TRUE) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()


## Full plot
png("../../drafts/figs/TT25_tri_gsw_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = gsw_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3, ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, 
                width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, 
                y = 0.03, 
                group = full_trt, 
                label = .group), 
            size = 6, position = position_dodge(width = 0.75), 
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 0.03), breaks = seq(0, 0.03, 0.01)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("g"["sw"]* " (mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

## Soil source x AM fungal source interaction plot
tri_gsw_soilAMFint_plot <-
  ggplot(data = gsw_soilFung_int,
         x = ExpSoilSource) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource, y = response,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, 
                 y = response,
                 fill = ExpFungSource),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpSoilSource, y = 0.02, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(0, 0.02), breaks = seq(0, 0.02, 0.005)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("g"["sw"]* " (mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
tri_gsw_soilAMFint_plot

###############
# Vcmax25
###############
vcmax_tri <- lmer(log(vcmax25) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                    wet_rhizome_mass_g + (1 | machine) + (1 | plot_field), 
                  data = df_noSterile)

# Check normality assumptions
plot(vcmax_tri)
qqnorm(residuals(vcmax_tri))
qqline(residuals(vcmax_tri))
hist(residuals(vcmax_tri))
shapiro.test(residuals(vcmax_tri))
outlierTest(vcmax_tri)

# Model output
summary(vcmax_tri)
Anova(vcmax_tri)
r.squaredGLMM(vcmax_tri)

# Pairwise comparisons
emmeans(vcmax_tri, pairwise~ExpFungSource, type = "response")
## Weeded fungal source has greater Vcmax25 than ambient fungal source

cld(emmeans(vcmax_tri, pairwise~ExpSoilSource*ExpFungSource))
## Weeded fungal source has greater Vcmax25 than ambient fungal source,
## but only in plants conditioned with weeded soils

cld(emmeans(vcmax_tri, pairwise~PlantGMTrt*ExpSoilSource))
## Plants conditioned with ambient soils have greater Vcmax25 than plants
## conditioned with weeded soils, but only in plants with legacy of growing
## in weeded treatment

cld(emmeans(vcmax_tri, pairwise~PlantGMTrt*ExpFungSource), alpha = 0.1)
## Weeded fungal source has marginally greater Vcmax than ambient fungal
## source, but only in plants with legacy of growing in ambient plots.

# Plot prep (full)
vcmax_tri_results <- cld(emmeans(vcmax_tri, 
                                 ~PlantGMTrt*ExpSoilSource*ExpFungSource, 
                                 type = "response"), 
                         Letters = LETTERS, reversed = TRUE) %>%
  data.frame() %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded")))

# Plot prep (ExpSoilSource * ExpFungSource interaction)
vcmax_soilFung_int <- cld(emmeans(vcmax_tri, ~ExpSoilSource*ExpFungSource, 
                                type = "response"), 
                        Letters = LETTERS) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()

# Plot prep (ExpFungSource)
vcmax_soilFung_ind <- cld(emmeans(vcmax_tri, ~ExpFungSource, 
                                  type = "response"), 
                          Letters = LETTERS, alpha = 0.1) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()

## Full plot
png("../../drafts/figs/TT25_tri_vcmax_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = vcmax_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf,xmax = 1.5,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5,xmax = 3,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 30, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("V"["cmax25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

## Soil source x AM fungal source interaction plot
tri_vcmax_soilAMFint_plot <-
  ggplot(data = vcmax_soilFung_int,
         x = ExpSoilSource) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource, 
                    y = response,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, 
                 y = response,
                 fill = ExpFungSource),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpSoilSource, y = 20, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("V"["cmax25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
tri_vcmax_soilAMFint_plot

## Individual AMF response
vcmax_AMF_ind_plot <- ggplot(data = vcmax_soilFung_ind) +
  geom_errorbar(aes(x = ExpFungSource,
                    y = response,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5) +
  geom_point(aes(x = ExpFungSource,
                 y = response,
                 fill = ExpFungSource),
             size = 6, shape = 21) +
  geom_text(aes(x = ExpFungSource, y = 20, label = .group), 
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Fungal Source", 
       y = expression(bold("V"["cmax25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")"))) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")

###############
# Jmax
###############
jmax_tri <- lmer(log(jmax25) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                   wet_rhizome_mass_g + (1 | machine) + (1 | plot_field), 
                 data = df_noSterile)

# Check normality assumptions
plot(jmax_tri)
qqnorm(residuals(jmax_tri))
qqline(residuals(jmax_tri))
hist(residuals(jmax_tri))
shapiro.test(residuals(jmax_tri))
outlierTest(jmax_tri)

# Model output
summary(jmax_tri)
Anova(jmax_tri)
r.squaredGLMM(jmax_tri)

# Pairwise comparisons
emmeans(jmax_tri, pairwise~ExpFungSource)
## Weeded fungal source has greater Jmax25 than non-weeded fungal source

cld(emmeans(jmax_tri, pairwise~PlantGMTrt*ExpFungSource))
## Weeded fungal source has greater Jmax25 than non-weeded fungal source, but
## this response is only observed in plants with legacy of growing in ambient
## plots

cld(emmeans(jmax_tri, pairwise~ExpSoilSource*ExpFungSource))
## Weeded fungal source has greater Jmax25 than non-weeded fungal source, but
## this response is only observed in plants conditioned with weeded soils

cld(emmeans(jmax_tri, pairwise~PlantGMTrt*ExpSoilSource), alpha = 0.1)
## Plants conditioned with ambient soil experienced a marginal increase in
## Jmax25 compared to plants conditioned with weeded soil, but this response
## was only observed in plants with legacy of growing in weeded plots

# Plot prep (full)
jmax_tri_results <- cld(emmeans(jmax_tri, 
                                ~PlantGMTrt*ExpSoilSource*ExpFungSource, 
                                type = "response"), 
                        Letters = LETTERS, reversed = TRUE) %>%
  data.frame() %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded")))

# Plot prep (ExpSoilSource * ExpFungSource interaction)
jmax_soilFung_int <- cld(emmeans(jmax_tri, ~ExpSoilSource*ExpFungSource, 
                                  type = "response"), 
                          Letters = LETTERS) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()

# Plot prep (ExpFungSource)
jmax_soilFung_ind <- cld(emmeans(jmax_tri, ~ExpFungSource, 
                                  type = "response"), 
                          Letters = LETTERS, alpha = 0.1) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()


# Full plot
png("../../drafts/figs/TT25_tri_jmax_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = jmax_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3, ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 60, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 20)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("J"["max25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

## Soil source x AM fungal source interaction plot
tri_jmax_soilAMFint_plot <-
  ggplot(data = jmax_soilFung_int,
         x = ExpSoilSource) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf, xmax = 1.5,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource, 
                    y = response,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, 
                 y = response,
                 fill = ExpFungSource),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpSoilSource, y = 40, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("J"["max25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")

## AM fungal source (individual) plot
jmax_AMF_ind_plot <- ggplot(data = jmax_soilFung_ind) +
  geom_errorbar(aes(x = ExpFungSource,
                    y = response,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5) +
  geom_point(aes(x = ExpFungSource,
                 y = response,
                 fill = ExpFungSource),
             size = 6, shape = 21) +
  geom_text(aes(x = ExpFungSource, y = 40, label = .group), 
            size = 6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Fungal Source", 
       y = expression(bold("J"["max25"]* " ("*mu*"mol m"^"-2"*" s"^"-1"*")"))) +
  guides(fill = "none") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")

###############
# Jmax:Vcmax
###############
jmaxvcmax_tri <- lmer(jmax25_vcmax25 ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                        wet_rhizome_mass_g + (1 | machine) + (1 | plot_field), 
                   data = df_noSterile)

# Check normality assumptions
plot(jmaxvcmax_tri)
qqnorm(residuals(jmaxvcmax_tri))
qqline(residuals(jmaxvcmax_tri))
hist(residuals(jmaxvcmax_tri))
shapiro.test(residuals(jmaxvcmax_tri))
outlierTest(jmaxvcmax_tri)

# Model output
summary(jmaxvcmax_tri)
Anova(jmaxvcmax_tri)
r.squaredGLMM(jmaxvcmax_tri)

# Pairwise comparisons
cld(emmeans(jmaxvcmax_tri, pairwise~PlantGMTrt*ExpSoilSource, type = "response"))
## Plants with legacy of growing in ambient plots have greater Jmax25:Vcmax25
## compared to plants with legacy of growing in weeded plots, but this response
## is only observed when plants were conditioned with soils from ambient plots

# Plot prep (full)
jmaxvcmax_tri_results <- cld(emmeans(jmaxvcmax_tri, 
                                     ~PlantGMTrt*ExpSoilSource*ExpFungSource, type = "response"), 
                             Letters = LETTERS, reversed = TRUE) %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded"))) %>% 
  data.frame()

## Full plot
png("../../drafts/figs/TT25_tri_jmaxvcmax_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = jmaxvcmax_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf,xmax = 1.5,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5,xmax = 3,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = emmean,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = emmean,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 3.5, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(1.5, 3.5), breaks = seq(1.5, 3.5, 1)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("J"["max25"]* ": V"["cmax25"]*" (unitless)")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

###############
# iWUE
###############
df_noSterile$iwue[c(45, 89)] <- NA

iwue_tri <- lmer(iwue ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                   wet_rhizome_mass_g + (1 | machine) + (1 | plot_field), 
                 data = filter(df_noSterile))

# Check normality assumptions
plot(iwue_tri)
qqnorm(residuals(iwue_tri))
qqline(residuals(iwue_tri))
hist(residuals(iwue_tri))
shapiro.test(residuals(iwue_tri))
outlierTest(iwue_tri)

# Model output
summary(iwue_tri)
Anova(iwue_tri)
r.squaredGLMM(iwue_tri)

# Post-hoc comparisons
cld(emmeans(iwue_tri, pairwise~PlantGMTrt*ExpSoilSource), alpha = 0.2)

# Plot prep (full)
iwue_tri_results <- cld(emmeans(iwue_tri, 
                              ~PlantGMTrt*ExpSoilSource*ExpFungSource, type = "response"), 
                      Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded"))) %>% 
  data.frame()

## Full plot
png("../../drafts/figs/TT25_tri_iwue_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = iwue_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf,xmax = 1.5,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5,xmax = 3,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = emmean,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = emmean,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 150, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("iWUE ("*mu*"mol mol"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

###############
# PNUE
###############
df_noSterile$pnue[c(44)] <- NA

pnue_tri <- lmer(pnue ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                   (1 | machine) + (1 | plot_field), 
                 data = df_noSterile)

# Check normality assumptions
plot(pnue_tri)
qqnorm(residuals(pnue_tri))
qqline(residuals(pnue_tri))
hist(residuals(pnue_tri))
shapiro.test(residuals(pnue_tri))
outlierTest(pnue_tri)

# Model output
summary(pnue_tri)
Anova(pnue_tri)
r.squaredGLMM(pnue_tri)

# Post-hoc comparisons
cld(emmeans(pnue_tri, pairwise~ExpSoilSource*ExpFungSource), alpha = 0.3)

cld(emmeans(pnue_tri, ~ExpSoilSource * ExpFungSource | PlantGMTrt))
## No difference in experimental fungal or soil source in plants growing
## in ambient plots
##
## In plants with history of growing in weeded plots, PNUE was greatest
## in plants inoculated with ambient AMF communities, but only in ambient
## soils
## 
## Additionally, PNUE is greater in ambient soils than weeded soils, but only
## in plants inoculated with ambient fungal communities

###############
# Nmass
###############
df_noSterile$pnue[c(44)] <- NA

nmass_tri <- lmer(nmass ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                   (1 | plot_field), data = df_noSterile)

# Check normality assumptions
plot(nmass_tri)
qqnorm(residuals(nmass_tri))
qqline(residuals(nmass_tri))
hist(residuals(nmass_tri))
shapiro.test(residuals(nmass_tri))
outlierTest(nmass_tri)

# Model output
summary(nmass_tri)
Anova(nmass_tri)

# Post-hoc comparisons
cld(emmeans(nmass_tri, pairwise~ExpFungSource*PlantGMTrt))

###############
# Marea
###############
df_noSterile$marea[c(11, 15, 83, 101)] <- NA

marea_tri <- lmer(log(marea) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                    (1 | plot_field), data = df_noSterile)

# Check normality assumptions
plot(marea_tri)
qqnorm(residuals(marea_tri))
qqline(residuals(marea_tri))
hist(residuals(marea_tri))
shapiro.test(residuals(marea_tri))
outlierTest(marea_tri)

# Model output
summary(marea_tri)
Anova(marea_tri)

# Post-hoc comparisons
cld(emmeans(marea_tri, ~ExpFungSource*PlantGMTrt), alpha = 0.1)

cld(emmeans(marea_tri, ~ExpSoilSource*PlantGMTrt), alpha = 0.1)

cld(emmeans(marea_tri, ~ExpSoilSource*ExpFungSource), alpha = 0.2)

###############
# Narea
###############
df_noSterile$narea[c(11, 15, 83, 101)] <- NA

narea_tri <- lmer(log(narea) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                    (1 | plot_field), data = df_noSterile)

# Check normality assumptions
plot(narea_tri)
qqnorm(residuals(narea_tri))
qqline(residuals(narea_tri))
hist(residuals(narea_tri))
shapiro.test(residuals(narea_tri))
outlierTest(narea_tri)

# Model output
summary(narea_tri)
Anova(narea_tri)

# Post-hoc comparisons
cld(emmeans(narea_tri, pairwise~ExpSoilSource*PlantGMTrt), alpha = 0.1)

cld(emmeans(narea_tri, ~ExpSoilSource * ExpFungSource), alpha = 0.2)

 ###############
# Leaf d15N
###############
d15n_tri <- lmer(leaf_d15n ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                   (1 | plot_field), data = df_noSterile)

# Check normality assumptions
plot(d15n_tri)
qqnorm(residuals(d15n_tri))
qqline(residuals(d15n_tri))
hist(residuals(d15n_tri))
shapiro.test(residuals(d15n_tri))
outlierTest(d15n_tri)

# Model output
summary(d15n_tri)
Anova(d15n_tri)
r.squaredGLMM(d15n_tri)

# Post-hoc comparisons
emmeans(d15n_tri, pairwise~ExpFungSource)
## Increased 15N in plants inoculated with weeded fungal communities.
## Implies increased N acquisition from weeded fungal communities

emmeans(d15n_tri, pairwise~ExpSoilSource)
## Increased 15N in plants conditioned with weeded soils

###############
# TLA
###############
tla_tri <- lmer(sqrt(tla_cm2) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + (1 | plot_field), 
              data = df_noSterile)

# Check normality assumptions
plot(tla_tri)
qqnorm(residuals(tla_tri))
qqline(residuals(tla_tri))
hist(residuals(tla_tri))
shapiro.test(residuals(tla_tri))
outlierTest(tla_tri)

# Model output
summary(tla_tri)
Anova(tla_tri) # Exp fungal source effect

# Post-hoc comparisons
emmeans(tla_tri, pairwise~ExpFungSource, type = "response") 
## Greater total leaf area when plants were inoculated with weeded fungal communities

cld(emmeans(tla_tri, ~ExpFungSource*ExpSoilSource|PlantGMTrt))

# Plot prep (full)
tla_tri_results <- cld(emmeans(tla_tri, 
                                ~PlantGMTrt*ExpSoilSource*ExpFungSource, type = "response"), 
                        Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded"))) %>% 
  data.frame()

## Full plot
png("../../drafts/figs/TT25_tri_tla_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = tla_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf,xmax = 1.5,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5,xmax = 3,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 80, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("Total leaf area (cm"^"2"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

###############
# total biomass
###############
tbio_tri <- lmer(log(total_wet_biomass_g) ~ PlantGMTrt * ExpSoilSource * ExpFungSource  + (1 | plot_field), 
               data = df_noSterile)

# Check normality assumptions
plot(tbio_tri)
qqnorm(residuals(tbio_tri))
qqline(residuals(tbio_tri))
hist(residuals(tbio_tri))
shapiro.test(residuals(tbio_tri))
outlierTest(tbio_tri)

# Model output
summary(tbio_tri)
Anova(tbio_tri)

# Post-hoc comparison
emmeans(tbio_tri, pairwise~PlantGMTrt)

# Plot prep (full)
tbio_tri_results <- cld(emmeans(tbio_tri, 
                               ~PlantGMTrt*ExpSoilSource*ExpFungSource, type = "response"), 
                       Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded"))) %>% 
  data.frame()

## Full plot
png("../../drafts/figs/TT25_tri_tbio_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = tbio_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf,xmax = 1.5,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 30, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = "Total biomass (g)",
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

###############
# Root:shoot
###############
df_noSterile$root_shoot[c(25, 31, 81, 92, 107, 109, 111, 113, 114, 115, 116)] <- NA

rootshoot_tri <- lm(wet_root_shoot ~ PlantGMTrt * ExpSoilSource * ExpFungSource, 
               data = df_noSterile)

# Check normality assumptions
plot(rootshoot_tri)
qqnorm(residuals(rootshoot_tri))
qqline(residuals(rootshoot_tri))
hist(residuals(rootshoot_tri))
shapiro.test(residuals(rootshoot_tri))
outlierTest(rootshoot_tri)

# Model output
summary(rootshoot_tri)
Anova(rootshoot_tri)

# Plot prep (full)
rootshoot_tri_results <- cld(emmeans(rootshoot_tri, 
                                ~PlantGMTrt*ExpSoilSource*ExpFungSource, type = "response"), 
                        Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"),
         full_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource, "_Fung", ExpFungSource),
         plot_trt = str_c("Plant",PlantGMTrt, "_Soil", ExpSoilSource),
         facet_label = 
           factor(PlantGMTrt,
                  levels = c("NW", "W"),
                  labels = c("Plant history: ambient", "Plant history: weeded"))) %>% 
  data.frame()

## Full plot
png("../../drafts/figs/TT25_tri_rootshoot_full.png", height = 5, width = 8, 
    units = "in", res = 600)
ggplot(data = rootshoot_tri_results) +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = -Inf,xmax = 1.5,
            ymin = -Inf,ymax = Inf,
            alpha = 0.05, fill = "#F1B700") +
  geom_rect(aes(fill = ExpSoilSource),
            xmin = 1.5, xmax = 3,
            ymin = -Inf, ymax = Inf,
            alpha = 0.05, fill = "#00B2BE") +
  geom_errorbar(aes(x = ExpSoilSource,
                    y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource,
                 y = response,
                 fill = ExpFungSource,
                 group = full_trt),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpFungSource, y = 8, group = full_trt, label = .group), 
            size = 5, position = position_dodge(width = 0.75),
            fontface = "bold") +
  scale_y_continuous(limits = c(2, 8), breaks = seq(2, 8, 2)) +
  scale_fill_manual(values = gm.colors, labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  facet_grid(~facet_label) +
  labs(x = "Experimental Soil Source", 
       y = "Root:shoot (unitless)",
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
dev.off()

###############
# Rhizome mass
###############
rhizome_tri <- lmer(log(wet_rhizome_mass_g) ~ PlantGMTrt * ExpSoilSource * ExpFungSource +
                     (1 | plot_field), 
                    data = df_noSterile)

# Check normality assumptions
plot(rhizome_tri)
qqnorm(residuals(rhizome_tri))
qqline(residuals(rhizome_tri))
hist(residuals(rhizome_tri))
shapiro.test(residuals(rhizome_tri))
outlierTest(rhizome_tri)

# Model output
summary(rhizome_tri)
Anova(rhizome_tri)


##############################
# Plot arrangements
###############
png("../../drafts/figs/TT25_soilAMint_plot.png", height = 8, 
    width = 8, units = "in", res = 600)
ggarrange(tri_anet_soilAMFint_plot, tri_gsw_soilAMFint_plot,
          tri_vcmax_soilAMFint_plot, tri_jmax_soilAMFint_plot,
          nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", "(d)"), align = "hv")
dev.off()

png("../../drafts/figs/TT25_AMind_plot.png", height = 4, 
    width = 10, units = "in", res = 600)
ggarrange(vcmax_AMF_ind_plot, jmax_AMF_ind_plot,
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)"), align = "hv",
          font.label = list(size = 18))
dev.off()


