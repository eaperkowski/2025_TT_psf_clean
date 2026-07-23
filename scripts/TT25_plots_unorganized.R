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
                                 ifelse(wet_root_shoot > 20, NA, wet_root_shoot)))

# Some plot aesthetics
gm.colors <- c("#F1B700", "#00B2BE")
facet_lab <- c("Plant history: ambient", "Plant history: weeded")
names(facet_lab) <- c("NW", "W")

# Remove outliers
df_noSterile$anet[89] <- NA

# Models
anet_tri <- lmer(log(anet + 1) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                   wet_rhizome_mass_g + (1 | machine) + (1 | plot_field), 
                 data = df_noSterile)
anet_tri_amf <- lmer(log(anet + 1) ~ ExpFungSource * AMF_asep_plant + (1 | machine) + (1 | plot_field), 
                     data = df_noSterile)


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
                         Letters = LETTERS, reversed = TRUE, type = "response") %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()


## Full plot
png("../drafts/figs/TT25_tri_anet_full.png", height = 5, width = 8, 
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
  geom_errorbar(aes(x = ExpSoilSource, y = response,
                    group = full_trt,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, y = response,
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
  geom_errorbar(aes(x = ExpSoilSource, y = response,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, y = response,
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
                        Letters = LETTERS, reversed = TRUE, alpha = 0.1) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()


## Full plot
png("../drafts/figs/TT25_tri_gsw_full.png", height = 5, width = 8, 
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
png("../drafts/figs/TT25_tri_vcmax_full.png", height = 5, width = 8, 
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
png("../drafts/figs/TT25_tri_jmax_full.png", height = 5, width = 8, 
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
png("../drafts/figs/TT25_tri_jmaxvcmax_full.png", height = 5, width = 8, 
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
png("../drafts/figs/TT25_tri_iwue_full.png", height = 5, width = 8, 
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

# Plot Prep (ExpSoilSource * ExpFungSource interaction)
iwue_soilFung_int <- cld(emmeans(iwue_tri, ~ExpFungSource*ExpSoilSource, 
                                 type = "response"), 
                         Letters = LETTERS, alpha = 0.1) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()

## Soil source x AM fungal source interaction plot
tri_iwue_soilAMFint_plot <-
  ggplot(data = iwue_soilFung_int,
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
                    y = emmean,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, 
                 y = emmean,
                 fill = ExpFungSource),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpSoilSource, y = 150, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("iWUE ("*mu*"mol mol"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
tri_iwue_soilAMFint_plot


# Plot Prep (ExpSoilSource * ExpFungSource interaction)
narea_soilFung_int <- cld(emmeans(narea_tri, ~ExpFungSource*ExpSoilSource, 
                                  type = "response"), 
                          Letters = LETTERS, alpha = 0.2, reversed = TRUE) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()

## Soil source x AM fungal source interaction plot
tri_narea_soilAMFint_plot <-
  ggplot(data = narea_soilFung_int,
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
  geom_text(aes(x = ExpSoilSource, y = 4.5, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(2, 4.5), breaks = seq(2, 4, 1)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("N"["area"]*" (g m"^"-2"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
tri_narea_soilAMFint_plot

# Plot Prep (ExpSoilSource * ExpFungSource interaction)
d15n_soilFung_int <- cld(emmeans(d15n_tri, ~ExpFungSource*ExpSoilSource, 
                                 type = "response"), 
                         Letters = LETTERS) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()

## Soil source x AM fungal source interaction plot
tri_d15n_soilAMFint_plot <-
  ggplot(data = d15n_soilFung_int,
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
                    y = emmean,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, 
                 y = emmean,
                 fill = ExpFungSource),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpSoilSource, y = 3, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("Leaf "*delta^"15"*"N (‰)")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
tri_d15n_soilAMFint_plot

# Plot Prep (ExpSoilSource * ExpFungSource interaction)
pnue_soilFung_int <- cld(emmeans(pnue_tri, ~ExpFungSource*ExpSoilSource, 
                                 type = "response"), 
                         Letters = LETTERS, alpha = 0.3, reversed = TRUE) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()

## Soil source x AM fungal source interaction plot
tri_pnue_soilAMFint_plot <-
  ggplot(data = pnue_soilFung_int,
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
                    y = emmean,
                    group = ExpFungSource,
                    ymin = lower.CL,
                    ymax = upper.CL),
                linewidth = 1, width = 0.5, 
                position = position_dodge(width = 0.75)) +
  geom_point(aes(x = ExpSoilSource, 
                 y = emmean,
                 fill = ExpFungSource),
             size = 6, shape = 21, 
             position = position_dodge(width = 0.75)) +
  geom_text(aes(x = ExpSoilSource, y = 1.2, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(-0.2, 1.2), breaks = seq(0, 1.2, 0.4)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("PNUE (μmol g"^"-1"*"N s"^"-1"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
tri_pnue_soilAMFint_plot

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
png("../drafts/figs/TT25_tri_tla_full.png", height = 5, width = 8, 
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

# Plot Prep (ExpSoilSource * ExpFungSource interaction)
tla_soilFung_int <- cld(emmeans(tla_tri, ~ExpFungSource*ExpSoilSource, 
                                type = "response"), 
                        Letters = LETTERS, alpha = 0.1) %>% 
  mutate(.group = trimws(.group, "both")) %>% data.frame()

## Soil source x AM fungal source interaction plot
tri_tla_soilAMFint_plot <-
  ggplot(data = tla_soilFung_int,
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
  geom_text(aes(x = ExpSoilSource, y = 60, 
                group = ExpFungSource, label = .group), 
            size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(20, 60), breaks = seq(20, 60, 20)) +
  scale_fill_manual(values = gm.colors, 
                    labels = c("ambient", "weeded")) +
  scale_x_discrete(labels = c("ambient", "weeded")) +
  labs(x = "Experimental Soil Source", 
       y = expression(bold("Total leaf area (cm"^"2"*")")),
       fill = "AM fungal source") +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")
tri_tla_soilAMFint_plot

