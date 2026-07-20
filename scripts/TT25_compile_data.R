
##############################################################################
# Libraries
##############################################################################
library(tidyverse)

##############################################################################
# Read data
##############################################################################

# Plant IDs and their treatments
df_ids <- read.csv("../data/TT25_tri_trt_list.csv")
head(df_ids)

# Read photosynthetic data
df_photo <- read.csv("../data/TT25_tri_photo_traits.csv")
head(df_photo)

# Read harvest data
df_harvest <- read.csv("../data/TT25_tri_harvest.csv") %>%
  dplyr::select(id:gasket, wet_shoot_mass_g = shoot_mass_g, 
                wet_root_mass_g = root_mass_g, wet_rhizome_mass_g = rhizome_mass_g)
head(df_harvest)

# Read leaf area data
df_la <- read.csv("../data/TT25_leaf_area.csv") %>%
  filter(id != "5904a" & id != "7574a") %>%
  mutate(id = ifelse(id == "5904b", "5904", 
                     ifelse(id == "7574b", "7574", id)))
head(df_la)

# Read dry mass data (only focal leaf so far)
df_drymass <- read.csv("../data/TT25_dry_masses.csv")
head(df_drymass)

# Read leaf isotope data, merge with leaf area data, then calculate additional
# leaf nutrient data
df_isotope <- read.csv("../data/TT25_leaf_isotope_data.csv") %>%
  dplyr::select(id = Identifier.1, nmass_perc = X.N, d15N = d15Nair, 
                cmass_perc = X.C, d13C = d13Cvpdb, run.comment) %>%
  filter(id != "UT729" & id != "LGlu4510" & id != "ureaJT" & 
           id != "5904a" & id != "7574a") %>%
  mutate(id = ifelse(id == "5904b", "5904", 
                     ifelse(id == "7574b", "7574", id)))
head(df_isotope)

# Read AMF dataset (for aseptate/septate hyphal lengths on plant side)
df_AMhyphae_plant <- read.csv("../data/TT25_AMhyphae.csv") %>%
  filter(Plant == "Plant") %>%
  dplyr::select(id = SAMPLE_number, 
                AMF_asep_plant = Aseptate, 
                AMF_sep_plant = Septate)

# Read AMF dataset (for aseptate/septate hyphal lengths on non-plant side)
df_AMhyphae_noplant <- read.csv("../data/TT25_AMhyphae.csv") %>%
  filter(Plant == "No plant") %>%
  dplyr::select(id = SAMPLE_number, 
                AMF_asep_noplant = Aseptate, 
                AMF_sep_noplant = Septate)
  

##############################################################################
# Compile data into single dataset
##############################################################################
# Merge photosynthetic and harvest dataset
df <- df_ids %>%
  full_join(df_photo, by = "id") %>%
  full_join(df_harvest, by = "id") %>%
  full_join(df_isotope, by = "id") %>%
  full_join(df_la, by = "id") %>%
  full_join(df_drymass, by = "id") %>%
  full_join(df_AMhyphae_plant, by = "id") %>%
  full_join(df_AMhyphae_noplant, by = "id") %>%
  mutate(jmax25_vcmax25 = Jmax / Vcmax,
         nmass = nmass_perc / 100,
         cmass = cmass_perc / 100,
         leaf_cn = nmass / cmass,
         marea = focal_dryMass_g / isotope_la_cm2 * 10000,
         narea = nmass * marea,
         pnue = anet / narea,
         pnue = ifelse(pnue < 0, NA, pnue),
         total_wet_biomass_g = wet_shoot_mass_g + wet_root_mass_g + wet_rhizome_mass_g,
         wet_root_shoot = (wet_root_mass_g + wet_rhizome_mass_g) / wet_shoot_mass_g,
         rmf = (wet_root_mass_g + wet_rhizome_mass_g) / total_wet_biomass_g,
         lar = tla_cm2 / total_wet_biomass_g,
         canopy_A_umols = (anet / 10000) * tla_cm2,
         canopy_A_mmols = canopy_A_umols * 1000) %>%
  dplyr::select(id:ExpFungSource, gasket, alive_dead, senescing, flowering,
                machine, anet:ci.ca, vcmax25 = Vcmax, jmax25 = Jmax, 
                jmax25_vcmax25, rd25 = Rd, tpu = TPU, focal_area_cm2 = isotope_la_cm2,
                focal_dryMass_g, marea, nmass, narea, leaf_d15n = d15N,
                cmass, leaf_cn, leaf_d13c = d13C, pnue, wet_shoot_mass_g:wet_rhizome_mass_g,
                phosphorus_la_cm2, extra_la_cm2, tla_cm2, total_wet_biomass_g:lar,
                canopy_A_mmols, AMF_asep_plant:AMF_sep_noplant)
head(df)

write.csv(df, "../data/TT25_tri_compiled.csv", row.names = F)
