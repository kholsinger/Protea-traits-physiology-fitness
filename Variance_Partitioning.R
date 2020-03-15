

library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(rstanarm)
library(shinystan)
library(tidyr)
library(lme4)



### Step 1: Get Data in Place

dat_full <- read.csv("Individual_Data_for_Analysis.csv", na.strings = ".")


### Step 2: Remove Obs from LDSA, McGregor, and JK2
dat_sub <- dat_full[dat_full$Species != "LDSA",]
dat_sub <- dat_sub[dat_sub$Site_ID != "McGregor",]
dat_sub <- dat_sub[dat_sub$Site_ID != "Jonaskop_2",]
dat_sub <- droplevels(dat_sub)

## Make Ind_No a Factor (not important for this, but for other analysis)
dat_sub$Ind_No <- as.factor(dat_sub$Ind_No)

## Add WUE_Instantaneous
dat_sub$WUE_Instantaneous <- dat_sub$Photo/dat_sub$Trmmol


### Step 3: Remove Individuals that do not have all observations

## Try this with tidyr

# First, select the columns of interest
dat_sub <- dat_sub %>% dplyr::select(Species, Site, Age, Height_in, Canopy_area_in2,
                                     Total_Fruit_No, HV_m2_cm2, WD, Ks, KL_Units, TotalAssim,
                                     BW_dry_avg, Photo, Trmmol, Cond, WUE_Intrinsic, Thick_cm,
                                     Area_cm2, LMA_g.cm2, LWR, LD_g.cm3, Length_mm_Top,
                                     Density_mm2_Top, SPI_Top, WUE_Instantaneous)

# Remove Observations that Have Missing Data
dat_sub <- dat_sub %>% drop_na()



### Step 4: Run Spp. Random Intercept Model to get Among vs. Within Spp. Variation

options(mc.cores = parallel::detectCores())

# Run following model using different traits
mod.rstan <- stan_glmer(Height_in ~ (1|Species), family = gaussian(link = "identity"),
                        data = dat_sub)
summary(mod.rstan, digits = 5)


### Step 5: Determine the percent of within species variation

# From Kent: With a small number of species, the better estimate of among is
# (variance of the intercepts)/(variance of the intercepts + square of sigma)


mod.rstan.df <- as.data.frame(mod.rstan)

among_var <- apply(mod.rstan.df[,2:6], 1, var)
within_var <- (mod.rstan.df$sigma)^2
within_var_percent <- (within_var/(among_var + within_var))*100


## Bark Thickness
BT_within_var <- cbind(rep("Bark_Thick", times = length(within_var_percent)),
                        rep("Structural", times = length(within_var_percent)),
                        as.data.frame(within_var_percent))
colnames(BT_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")


## WD
WD_within_var <- cbind(rep("WD", times = length(within_var_percent)),
                       rep("Structural", times = length(within_var_percent)),
                       as.data.frame(within_var_percent))
colnames(WD_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")


## Leaf Area
Area_within_var <- cbind(rep("Leaf_Area", times = length(within_var_percent)),
                        rep("Structural", times = length(within_var_percent)),
                        as.data.frame(within_var_percent))
colnames(Area_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")


## LMA
LMA_within_var <- cbind(rep("LMA", times = length(within_var_percent)),
                        rep("Structural", times = length(within_var_percent)),
                        as.data.frame(within_var_percent))
colnames(LMA_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")



## LWR
LWR_within_var <- cbind(rep("LWR", times = length(within_var_percent)),
                        rep("Structural", times = length(within_var_percent)),
                        as.data.frame(within_var_percent))
colnames(LWR_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")


## LD
LD_within_var <- cbind(rep("LD", times = length(within_var_percent)),
                        rep("Structural", times = length(within_var_percent)),
                        as.data.frame(within_var_percent))
colnames(LD_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")


## Stomatal Length
Stom_Length_within_var <- cbind(rep("Stom_Length", times = length(within_var_percent)),
                                 rep("Structural", times = length(within_var_percent)),
                                 as.data.frame(within_var_percent))
colnames(Stom_Length_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")



## Stomatal Density
Stom_Density_within_var <- cbind(rep("Stom_Density", times = length(within_var_percent)),
                       rep("Structural", times = length(within_var_percent)),
                       as.data.frame(within_var_percent))
colnames(Stom_Density_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")


## Stomatal Pore Index
SPI_within_var <- cbind(rep("SPI", times = length(within_var_percent)),
                                 rep("Structural", times = length(within_var_percent)),
                                 as.data.frame(within_var_percent))
colnames(SPI_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")





## Ks
Ks_within_var <- cbind(rep("Ks", times = length(within_var_percent)),
                       rep("Physiological", times = length(within_var_percent)),
                       as.data.frame(within_var_percent))
colnames(Ks_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")


## KL
KL_within_var <- cbind(rep("KL", times = length(within_var_percent)),
                       rep("Physiological", times = length(within_var_percent)),
                       as.data.frame(within_var_percent))
colnames(KL_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")



## WUE_Instantaneous
WUE_within_var <- cbind(rep("WUE", times = length(within_var_percent)),
                        rep("Physiological", times = length(within_var_percent)),
                        as.data.frame(within_var_percent))
colnames(WUE_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")



## A_area
A_area_within_var <- cbind(rep("A_area", times = length(within_var_percent)),
                       rep("Physiological", times = length(within_var_percent)),
                       as.data.frame(within_var_percent))
colnames(A_area_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")


## Stomatal Conductance
Stom_Cond_within_var <- cbind(rep("Stom_Cond", times = length(within_var_percent)),
                           rep("Physiological", times = length(within_var_percent)),
                           as.data.frame(within_var_percent))
colnames(Stom_Cond_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")



## Total_Assimilation
TotalAssim_within_var <- cbind(rep("Total_Assim", times = length(within_var_percent)),
                        rep("Physiological", times = length(within_var_percent)),
                        as.data.frame(within_var_percent))
colnames(TotalAssim_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")




## Height
Height_within_var <- cbind(rep("Height", times = length(within_var_percent)),
                               rep("Size", times = length(within_var_percent)),
                               as.data.frame(within_var_percent))
colnames(Height_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")


## Canopy
Canopy_within_var <- cbind(rep("Canopy Area", times = length(within_var_percent)),
                           rep("Size", times = length(within_var_percent)),
                           as.data.frame(within_var_percent))
colnames(Canopy_within_var) <- c("Trait_Name", "Trait_Class", "Within_Spp_Var")



## Merge Dataframes
var_df <- rbind(BT_within_var, WD_within_var, Area_within_var, LMA_within_var,
                LWR_within_var, LD_within_var, Stom_Length_within_var,
                Stom_Density_within_var, SPI_within_var,
                Ks_within_var, KL_within_var, WUE_within_var, A_area_within_var,
                Stom_Cond_within_var, TotalAssim_within_var,
                Height_within_var, Canopy_within_var)

write.csv(var_df, "within_spp_variance_posterior_distributions.csv")


## Plot the Within Species Percent Variance

# First, Violin Plot of the Whole Posterior

p1 <- ggplot(data = var_df, aes(x = Trait_Name, y = Within_Spp_Var))
p1 + geom_violin() + facet_wrap(~Trait_Class) + theme_bw()

# Second, try boxplot

p2 <- ggplot(data = var_df, aes(x = Trait_Name, y = Within_Spp_Var))
p2 + geom_boxplot() + facet_grid(~Trait_Class) + theme_bw()



### Separate Structural and Physiological, then Merge plots

tiff("within_spp_var_phys_size.tiff", units="in", width=3, height=5, res=300)

p.phys <- ggplot(data = subset(var_df, Trait_Class == "Physiological"), aes(x = Trait_Name, y = Within_Spp_Var))
p.phys <- p.phys + geom_violin() + theme_bw() + ylim(0,105) + ggtitle("Physiological Traits") + coord_flip()

p.size <- ggplot(data = subset(var_df, Trait_Class == "Size"), aes(x = Trait_Name, y = Within_Spp_Var))
p.size <- p.size + geom_violin() + theme_bw() + ylim(0,105) + ggtitle("Size Traits") + coord_flip()

grid.arrange(p.phys, p.size, nrow=2)

dev.off()


tiff("within_spp_var_structural.tiff", units="in", width=3, height=5, res=300)

p.structural <- ggplot(data = subset(var_df, Trait_Class == "Structural"), aes(x = Trait_Name, y = Within_Spp_Var))
p.structural <- p.structural + geom_violin() + theme_bw() + ylim(0,105) + ggtitle("Structural Traits") + coord_flip()
p.structural

dev.off()



#################################################################
## Variance Partitioning Using Jane et al. Full Protea Dataset ##
#################################################################

## Step 1: Read in the data

## Note, this is the full data file from Mitchell et al. (2015)
##
trait_full <- read.csv("Protea_Trait_Data_Full.csv", header = TRUE, na.strings = ".")

names(trait_full)

# Subset the Data to Include the Traits and Columns of Interest

trait_sub <- trait_full %>% dplyr::select(Site_location, Region, Site_number, Species_name,
                                          Plant_no, Plant_age, Height, Age_of_oldest_leaf, area_lam_pet,
                                          lamina_thickness, atlas_name, resprouter, polymorphic,
                                          pollinator, dominant_color, wood_density, LMA,
                                          Canopy_area, LWratio, FWC, LDMC, Succulence,
                                          Stomatal_size_top, Stomatal_density_top, Stomatal_pore_index_top,
                                          Stomatal_size_bottom, Stomatal_density_bottom, Stomatal_pore_index_bottom)



mod.rstan <- stan_glmer(LMA ~ (1|Species_name), family = gaussian(link = "identity"),
                        data = trait_sub)
summary(mod.rstan, digits = 5)

mod.rstan.df.all <- as.data.frame(mod.rstan)

within_var_all <- (mod.rstan.df.all[,61])^2
region_var_all <- (mod.rstan.df.all[,153])^2
species_var_all <- (mod.rstan.df.all[,62])^2

within_var_all_percent <- (within_var_all/(within_var_all + region_var_all + species_var_all))*100
region_var_all_percent <- (region_var_all/(within_var_all + region_var_all + species_var_all))*100
species_var_all_percent <- (species_var_all/(within_var_all + region_var_all + species_var_all))*100


within.tmp <- ((within_var_all)/(within_var_all + species_var_all))*100


p1 <- ggplot(data = subset(trait_sub, Region == "Barberton"), aes(x = Species_name, y = LMA))
p1 <- p1 + geom_point() + ggtitle("Barberton") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p2 <- ggplot(data = subset(trait_sub, Region == "Baviaanskloof"), aes(x = Species_name, y = LMA))
p2 <- p2 + geom_point() + ggtitle("Baviaanskloof") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p3 <- ggplot(data = subset(trait_sub, Region == "Blyde_River"), aes(x = Species_name, y = LMA))
p3 <- p3 + geom_point() + ggtitle("Blyde_River") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p4 <- ggplot(data = subset(trait_sub, Region == "Cape_Peninsula"), aes(x = Species_name, y = LMA))
p4 <- p4 + geom_point() + ggtitle("Cape_Peninsula") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p5 <- ggplot(data = subset(trait_sub, Region == "CapeTown"), aes(x = Species_name, y = LMA))
p5 <- p5 + geom_point() + ggtitle("CapeTown") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p6 <- ggplot(data = subset(trait_sub, Region == "Cedarberg"), aes(x = Species_name, y = LMA))
p6 <- p6 + geom_point() + ggtitle("Cedarberg") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p7 <- ggplot(data = subset(trait_sub, Region == "De_Hoop"), aes(x = Species_name, y = LMA))
p7 <- p7 + geom_point() + ggtitle("De_Hoop") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p8 <- ggplot(data = subset(trait_sub, Region == "Garcia_Pass"), aes(x = Species_name, y = LMA))
p8 <- p8 + geom_point() + ggtitle("Garcia_Pass") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p9 <- ggplot(data = subset(trait_sub, Region == "Hottentots_Hollands"), aes(x = Species_name, y = LMA))
p9 <- p9 + geom_point() + ggtitle("Hottentots_Hollands") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p10 <- ggplot(data = subset(trait_sub, Region == "Kogelberg"), aes(x = Species_name, y = LMA))
p10 <- p10 + geom_point() + ggtitle("Kogelberg") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p11 <- ggplot(data = subset(trait_sub, Region == "Long_Tom_pass"), aes(x = Species_name, y = LMA))
p11 <- p11 + geom_point() + ggtitle("Long_Tom_pass") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p12 <- ggplot(data = subset(trait_sub, Region == "Nieuwoudtville"), aes(x = Species_name, y = LMA))
p12 <- p12 + geom_point() + ggtitle("Nieuwoudtville") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p13 <- ggplot(data = subset(trait_sub, Region == "Royal_Natal"), aes(x = Species_name, y = LMA))
p13 <- p13 + geom_point() + ggtitle("Royal_Natal") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p14 <- ggplot(data = subset(trait_sub, Region == "Swartberg"), aes(x = Species_name, y = LMA))
p14 <- p14 + geom_point() + ggtitle("Swartberg") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())

p15 <- ggplot(data = subset(trait_sub, Region == "Swellendam"), aes(x = Species_name, y = LMA))
p15 <- p15 + geom_point() + ggtitle("Swellendam") + ylim(0,0.085) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title.x = element_blank())


grid.arrange(p1,p2,p3,p6,p7,p8,p9,p10,p11,p12,p13,p14,nrow=4)



##############################
## With lme4, lmer function ##
##############################

Model <- lmer(Stomatal_density_top ~ (1|Species_name/Region), data = trait_sub)
summary(Model)

tmp <- as.data.frame(VarCorr(Model), comp="Variance")

var_region <- tmp[1,4]/sum(tmp[1,4], tmp[2,4], tmp[3,4])
var_spp <- tmp[2,4]/sum(tmp[1,4], tmp[2,4], tmp[3,4])
var_resid <- tmp[3,4]/sum(tmp[1,4], tmp[2,4], tmp[3,4])

var_region
var_spp
var_resid

