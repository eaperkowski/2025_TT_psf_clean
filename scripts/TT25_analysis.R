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

###############
# Anet - general
###############
# Remove outliers (Bonferroni: p < 0.001)
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
cld(emmeans(anet_tri, pairwise~PlantGMTrt*ExpFungSource), alpha = 0.15)

###############
# Anet - AMF
###############
# Model
anet_tri_amf <- lmer(log(anet + 1) ~ ExpFungSource * AMF_asep_plant + (1 | machine) + (1 | plot_field), 
                 data = df_noSterile)

# Check normality assumptions
plot(anet_tri_amf)
qqnorm(residuals(anet_tri_amf))
qqline(residuals(anet_tri_amf))
hist(residuals(anet_tri_amf))
shapiro.test(residuals(anet_tri_amf))
outlierTest(anet_tri_amf)

# Model output
summary(anet_tri_amf)
Anova(anet_tri_amf)


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
cld(emmeans(vcmax_tri, pairwise~ExpSoilSource*ExpFungSource))
cld(emmeans(vcmax_tri, pairwise~PlantGMTrt*ExpSoilSource), alpha = 0.1)
cld(emmeans(vcmax_tri, pairwise~PlantGMTrt*ExpFungSource), alpha = 0.1)

###############
# Vcmax25 - AM hyphae
###############
vcmax_tri_amf <- lmer(log(vcmax25) ~ ExpFungSource * AMF_asep_plant + (1 | machine) + (1 | plot_field), 
                      data = df_noSterile)

# Check normality assumptions
plot(vcmax_tri_amf)
qqnorm(residuals(vcmax_tri_amf))
qqline(residuals(vcmax_tri_amf))
hist(residuals(vcmax_tri_amf))
shapiro.test(residuals(vcmax_tri_amf))
outlierTest(vcmax_tri_amf)

# Model results
Anova(vcmax_tri_amf)

# Post-hoc tests
test(emtrends(vcmax_tri_amf, ~ExpFungSource, "AMF_asep_plant"))

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
cld(emmeans(jmax_tri, pairwise~PlantGMTrt*ExpFungSource))
cld(emmeans(jmax_tri, pairwise~ExpSoilSource*ExpFungSource))
cld(emmeans(jmax_tri, pairwise~PlantGMTrt*ExpSoilSource), alpha = 0.2)

###############
# Jmax25 - AM hyphae
###############
jmax_tri_amf <- lmer(log(jmax25) ~ ExpFungSource * AMF_asep_plant + (1 | machine) + (1 | plot_field), 
                      data = df_noSterile)

# Check normality assumptions
plot(jmax_tri_amf)
qqnorm(residuals(jmax_tri_amf))
qqline(residuals(jmax_tri_amf))
hist(residuals(jmax_tri_amf))
shapiro.test(residuals(jmax_tri_amf))
outlierTest(jmax_tri_amf)

# Model results
Anova(jmax_tri_amf)

# Post-hoc tests
test(emtrends(jmax_tri_amf, ~ExpFungSource, "AMF_asep_plant"))

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



###############
# Nmass
###############
df_noSterile$nmass[c(9, 11)] <- NA

nmass_tri <- lmer(log(nmass) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                    wet_rhizome_mass_g + (1 | plot_field), data = df_noSterile)

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
cld(emmeans(nmass_tri, pairwise~ExpFungSource*PlantGMTrt, type = "response"))

###############
# Marea
###############
df_noSterile$marea[c(11, 15, 83, 101)] <- NA

marea_tri <- lmer(log(marea) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                    wet_rhizome_mass_g + (1 | plot_field), data = df_noSterile)

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

cld(emmeans(marea_tri, ~ExpSoilSource*ExpFungSource*PlantGMTrt), alpha = 0.2)

###############
# Narea
###############
df_noSterile$narea[c(11, 15, 83,101)] <- NA

narea_tri <- lmer(log(narea) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                    wet_rhizome_mass_g + (1 | plot_field), data = df_noSterile)

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
                   wet_rhizome_mass_g + (1 | plot_field), data = df_noSterile)

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
emmeans(d15n_tri, pairwise~ExpSoilSource)
cld(emmeans(d15n_tri, pairwise~ExpSoilSource * ExpFungSource))

###############
# d15N - AM hyphae
###############
d15n_tri_amf <- lmer(leaf_d15n ~ PlantGMTrt * ExpFungSource * AMF_asep_plant + (1 | plot_field), 
                      data = df_noSterile)

# Check normality assumptions
plot(d15n_tri_amf)
qqnorm(residuals(d15n_tri_amf))
qqline(residuals(d15n_tri_amf))
hist(residuals(d15n_tri_amf))
shapiro.test(residuals(d15n_tri_amf))
outlierTest(d15n_tri_amf)

# Model results
Anova(d15n_tri_amf)

# Post-hoc comparisons
test(emtrends(d15n_tri_amf, ~PlantGMTrt * ExpFungSource, "AMF_asep_plant"))


###############
# PNUE
###############
pnue_tri <- lmer(sqrt(pnue) ~ PlantGMTrt * ExpSoilSource * ExpFungSource + 
                  wet_rhizome_mass_g + (1 | machine) + (1 | plot_field), 
                 data = subset(df_noSterile, pnue > 0))

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
cld(emmeans(pnue_tri, ~ExpSoilSource * ExpFungSource | PlantGMTrt))

###############
# TLA
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
cld(emmeans(tla_tri, ~ExpFungSource*ExpSoilSource|PlantGMTrt))

###############
# Total leaf area - AM hyphae
###############
tla_tri_amf <- lmer(log(tla_cm2) ~ ExpFungSource * AMF_asep_plant + (1 | machine) + (1 | plot_field), 
                      data = df_noSterile)

# Check normality assumptions
plot(tla_tri_amf)
qqnorm(residuals(tla_tri_amf))
qqline(residuals(tla_tri_amf))
hist(residuals(tla_tri_amf))
shapiro.test(residuals(tla_tri_amf))
outlierTest(tla_tri_amf)

# Model results
Anova(tla_tri_amf)

# Post-hoc tests
test(emtrends(tla_tri_amf, ~ExpFungSource, "AMF_asep_plant"))

###############
#  Total wet biomass
###############
tbio_tri <- lmer(log(total_wet_biomass_g) ~ PlantGMTrt * ExpSoilSource * ExpFungSource +
                  (1 | plot_field), data = df_noSterile)

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

# Post hoc comparisons
emmeans(tbio_tri, pairwise~PlantGMTrt, type = "response")

###############
# Rhizome mass
###############
rhizome_tri <- lmer(log(wet_rhizome_mass_g) ~ PlantGMTrt * ExpSoilSource * ExpFungSource +
                     (1 | plot_field), data = df_noSterile)

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

# Posthoc comparisons
emmeans(rhizome_tri, pairwise~PlantGMTrt, type = "response")

###############
# Leaf area ratio
###############
df_noSterile$lar[c(83)] <- NA

lar_tri <- lmer(sqrt(lar) ~ PlantGMTrt * ExpSoilSource * ExpFungSource +
                  (1 | plot_field), data = df_noSterile)

# Check normality assumptions
plot(lar_tri)
qqnorm(residuals(lar_tri))
qqline(residuals(lar_tri))
hist(residuals(lar_tri))
shapiro.test(residuals(lar_tri))
outlierTest(lar_tri)

# Model output
summary(lar_tri)
Anova(lar_tri)

# Post hoc comparisons
emmeans(lar_tri, pairwise~ExpSoilSource, type = "response")
