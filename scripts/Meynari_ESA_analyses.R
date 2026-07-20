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
         canopy_A_umols = (anet / 10000) * tla_cm2,
         canopy_A_mmols = canopy_A_umols * 1000)
head(df_noSterile)

###############
# Anet model
###############

# Remove outliers (Bonferroni: p < 0.001)
df_noSterile$anet[89] <- NA

anet_tri <- lmer(anet ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
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
cld(emmeans(anet_tri, pairwise~PlantGMTrt*ExpFungSource), alpha = 0.15)

###############
# Anet plot - soil * AMF
###############
anet_plot_prep <- cld(emmeans(anet_tri, ~ExpSoilSource*ExpFungSource), 
                      Letters = LETTERS, reversed = TRUE) %>%
  mutate(.group = trimws(.group, which = "both"))

anet_plot <- ggplot(data = anet_plot_prep, 
                    aes(x = ExpSoilSource, y = emmean, group = ExpFungSource)) +
  geom_rect(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#FCDE0A") +
  geom_rect(xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#668026") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.75), width = 0.2) +
  geom_point(aes(shape = ExpFungSource, fill = ExpFungSource),
             position = position_dodge(width = 0.75), size = 6) +
  geom_text(aes(y = 2.4, label = .group), size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(-0.8, 2.4), breaks = seq(-0.8, 2.4, 0.8)) +
  scale_x_discrete(labels = c("Ambient", "Weeded")) +
  scale_shape_manual(values = c(21, 22), labels = c("ambient", "weeded")) +
  scale_fill_manual(values = c("#FCDE0A", "#668026"), labels = c("ambient", "weeded")) +
  labs(x = "",
       y = expression(bolditalic("A")[bold("net")]*bold(" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       shape = "AM fungal source",
       fill = "AM fungal source") +
  theme_classic(base_size = 22) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
anet_plot

###############
# Anet plot - legacy * soil environment
###############
anet_plot_prep2 <- cld(emmeans(anet_tri, ~ExpSoilSource*ExpFungSource*PlantGMTrt), 
                      Letters = LETTERS) %>%
  filter(ExpSoilSource == ExpFungSource) %>%
  mutate(soil_trt = ExpSoilSource,
         .group = trimws(.group, which = "both"))
  
anet_plot2 <- ggplot(data = anet_plot_prep2, 
                    aes(x = PlantGMTrt, y = emmean, group = soil_trt)) +
  geom_rect(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#FCDE0A", alpha = 0.3) +
  geom_rect(xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#668026", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.75), width = 0.2) +
  geom_point(aes(shape = soil_trt, fill = soil_trt),
             position = position_dodge(width = 0.75), size = 6) +
  geom_text(aes(y = 3, label = .group), size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  scale_x_discrete(labels = c("+GM", "- GM")) +
  scale_shape_manual(values = c(21, 22), labels = c("+GM", "- GM")) +
  scale_fill_manual(values = c("#FCDE0A", "#668026"), labels = c("+GM", "- GM")) +
  labs(x = "Plant treatment legacy",
       y = expression(bolditalic("A")[bold("net")]*bold(" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       shape = "Soil + AMF source",
       fill = "Soil + AMF source") +
  theme_classic(base_size = 22) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
anet_plot2

###############
# TLA model
###############
tla_tri <- lmer(sqrt(tla_cm2) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                  (1 | plot_field), 
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
Anova(tla_tri)

# Post-hoc comparisons
emmeans(tla_tri, pairwise~ExpFungSource, type = "response") 
cld(emmeans(tla_tri, ~ExpFungSource*ExpSoilSource*PlantGMTrt))

###############
# TLA plot
###############
tla_plot_prep <- cld(emmeans(tla_tri, ~ExpSoilSource*ExpFungSource), 
                      Letters = LETTERS, alpha = 0.15) %>%
  mutate(.group = trimws(.group, which = "both"))

tla_plot <- ggplot(data = tla_plot_prep, 
                   aes(x = ExpSoilSource, y = emmean, group = ExpFungSource)) +
  geom_rect(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#FCDE0A") +
  geom_rect(xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#668026") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.75), width = 0.2) +
  geom_point(aes(shape = ExpSoilSource, fill = ExpFungSource),
             position = position_dodge(width = 0.75), size = 6) +
  geom_text(aes(y = 8, label = .group), size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(4, 8), breaks = seq(4, 8, 1)) +
  scale_x_discrete(labels = c("Ambient", "Weeded")) +
  scale_shape_manual(values = c(21, 22), labels = c("ambient", "weeded")) +
  scale_fill_manual(values = c("#FCDE0A", "#668026"), labels = c("ambient", "weeded")) +
  labs(x = "Soil Source",
       y = expression(bold("Total leaf area (cm"^"2"*")")),
       shape = "AM fungal source",
       fill = "AM fungal source") +
  theme_classic(base_size = 22) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

###############
# TLA plot - legacy * soil environment
###############
tla_plot_prep2 <- cld(emmeans(tla_tri, ~ExpSoilSource*ExpFungSource*PlantGMTrt), 
                       Letters = LETTERS) %>%
  filter(ExpSoilSource == ExpFungSource) %>%
  mutate(soil_trt = ExpSoilSource,
         .group = trimws(.group, which = "both"))

tla_plot2 <- ggplot(data = tla_plot_prep2, 
                     aes(x = PlantGMTrt, y = emmean, group = soil_trt)) +
  geom_rect(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#FCDE0A", alpha = 0.3) +
  geom_rect(xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#668026", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.75), width = 0.2) +
  geom_point(aes(shape = soil_trt, fill = soil_trt),
             position = position_dodge(width = 0.75), size = 6) +
  geom_text(aes(y = 9, label = .group), size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(3, 9), breaks = seq(3, 9, 1.5)) +
  scale_x_discrete(labels = c("+GM", "- GM")) +
  scale_shape_manual(values = c(21, 22), labels = c("+GM", "- GM")) +
  scale_fill_manual(values = c("#FCDE0A", "#668026"), labels = c("+GM", "- GM")) +
  labs(x = "Plant treatment legacy",
       y = expression(bold("Total leaf area (cm"^"2"*")")),
       shape = "Soil + AMF source",
       fill = "Soil + AMF source") +
  theme_classic(base_size = 22) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
tla_plot2

###############
# Canopy photosynthesis model
###############
df_noSterile$canopy_A_mmols[78] <- NA

canopyA_tri <- lmer(canopy_A_mmols ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                 (1 | machine) + (1 | plot_field), data = df_noSterile)

# Check normality assumptions
plot(canopyA_tri)
qqnorm(residuals(canopyA_tri))
qqline(residuals(canopyA_tri))
hist(residuals(canopyA_tri))
shapiro.test(residuals(canopyA_tri))
outlierTest(canopyA_tri)

# Model output
summary(canopyA_tri)
Anova(canopyA_tri)
performance::performance(canopyA_tri)

# Post-hoc comparisons
cld(emmeans(canopyA_tri, ~ExpFungSource*ExpSoilSource), alpha = 0.15, Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"))


###############
# Canopy A plot
###############
canA_plot_prep <- cld(emmeans(canopyA_tri, ~ExpSoilSource*ExpFungSource), 
                     Letters = LETTERS, alpha = 0.15) %>%
  mutate(.group = trimws(.group, which = "both"))

canA_plot <- ggplot(data = canA_plot_prep, 
                    aes(x = ExpSoilSource, y = emmean, group = ExpFungSource)) +
  geom_rect(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#FCDE0A", alpha = 0.2) +
  geom_rect(xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#668026", alpha = 0.2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.75), width = 0.2) +
  geom_point(aes(shape = ExpSoilSource, fill = ExpFungSource),
             position = position_dodge(width = 0.75), size = 6) +
  geom_text(aes(y = 9, label = .group), size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(-3, 9), breaks = seq(-3, 9, 3)) +
  scale_x_discrete(labels = c("Ambient", "Weeded")) +
  scale_shape_manual(values = c(21, 22), labels = c("ambient", "weeded")) +
  scale_fill_manual(values = c("#FCDE0A", "#668026"), labels = c("ambient", "weeded")) +
  labs(x = "",
       y = expression(bold("Canopy ")*bolditalic("A")[bold("net")]*bold(" (mmol s"^"-1"*")")),
       shape = "AM fungal source",
       fill = "AM fungal source") +
  theme_classic(base_size = 22) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
canA_plot

###############
# Canopy A plot - legacy * soil environment
###############
canA_plot_prep2 <- cld(emmeans(canopyA_tri, ~ExpSoilSource*ExpFungSource*PlantGMTrt), 
                      Letters = LETTERS, alpha = 0.15) %>%
  filter(ExpSoilSource == ExpFungSource) %>%
  mutate(soil_trt = ExpSoilSource,
         .group = trimws(.group, which = "both"))

canA_plot2 <- ggplot(data = canA_plot_prep2, 
                    aes(x = PlantGMTrt, y = emmean, group = soil_trt)) +
  geom_rect(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#FCDE0A", alpha = 0.3) +
  geom_rect(xmin = 1.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#668026", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.75), width = 0.2) +
  geom_point(aes(shape = soil_trt, fill = soil_trt),
             position = position_dodge(width = 0.75), size = 6) +
  geom_text(aes(y = 12, label = .group), size = 6, fontface = "bold",
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(-3, 12), breaks = seq(-3, 12, 3)) +
  scale_x_discrete(labels = c("+GM", "- GM")) +
  scale_shape_manual(values = c(21, 22), labels = c("+GM", "- GM")) +
  scale_fill_manual(values = c("#FCDE0A", "#668026"), labels = c("+GM", "- GM")) +
  labs(x = "Plant treatment legacy",
       y = expression(bold("Canopy ")*bolditalic("A")[bold("net")]*bold(" (mmol s"^"-1"*")")),
       shape = "Soil + AMF source",
       fill = "Soil + AMF source") +
  theme_classic(base_size = 22) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
canA_plot2

png("../drafts/figs/niyomi_esa/TT25_NM_soil_AMF_int.png", width = 16, height = 6,
    units = "in", res = 600)
ggarrange(anet_plot, tla_plot, canA_plot, ncol = 3, common.legend = TRUE, legend = "bottom",
          labels = c("(a)", "(b)", "(c)"), font.label = list(size = 22))
dev.off()


png("../drafts/figs/niyomi_esa/TT25_NM_legacy_soil_int.png", width = 12, height = 6,
    units = "in", res = 600)
ggarrange(anet_plot2, tla_plot2, ncol = 2, common.legend = TRUE, 
          legend = "bottom", font.label = list(size = 22))
dev.off()




png("../drafts/figs/niyomi_esa/TT25_NM_legacy_soil_canopyA.png",
    height = 5, width = 9, units = "in", res = 600)
canA_plot2
dev.off()

png("../drafts/figs/niyomi_esa/TT25_NM_legacy_soil_Anet_tla.png",
    height = 5, width = 9, units = "in", res = 600)
canA_plot2
dev.off()
