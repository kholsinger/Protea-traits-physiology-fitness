library(tidyverse)
library(rstanarm)
library(ggplot2)
library(dplyr)
library(ggcorrplot)

rm(list = ls())

options(mc.cores = parallel::detectCores())


## Get Data and Format

dat <- read.csv("Individual_Data_for_Analysis.csv", na = ".") %>%
  filter(Species != "LDSA") %>%
  filter(Site != "McGregor") %>%
  filter(Site != "Jonaskop_2") %>%
  mutate(WUE_Instantaneous = Photo/Trmmol) %>%
  select(Species, Site, Age, Height_in, Canopy_area_in2, Total_Fruit_No, HV_m2_cm2, WD,
         Ks, KL_Units, TotalAssim, BW_dry_avg, Photo, Trmmol, Cond, WUE_Intrinsic,
         Thick_cm, Area_cm2, LMA_g.cm2, LWR, LD_g.cm3, Length_mm_Top, Density_mm2_Top,
         SPI_Top, WUE_Instantaneous) %>%
  drop_na()

dat <- droplevels(dat)


Morph_Traits <- c("BW_dry_avg", "WD", "Area_cm2", "LMA_g.cm2",
                  "LD_g.cm3", "LWR", "Length_mm_Top", "Density_mm2_Top",
                  "SPI_Top")

Performance_Traits <- c("Ks", "KL_Units", "Photo", "TotalAssim", "Cond", "WUE_Instantaneous")

Size_Traits <- c("Height_in", "Canopy_area_in2")


####################################################################################
## First, let's look at the trait-phys relationships (6 phys x 9 traits = 54 R2s) ##
####################################################################################

phys_trait_df <- data.frame(matrix(nrow = (length(Performance_Traits)*(length(Morph_Traits))), ncol = 3))
colnames(phys_trait_df) <- c("Phys", "Trait", "R2")

ct <- 0
for(i in 1:length(Performance_Traits)){

  phys <- Performance_Traits[i]

  for(j in 1:length(Morph_Traits)){

    ct <- ct + 1
    morph <- Morph_Traits[j]
    mod <- stan_lmer(scale(dat[[phys]]) ~ scale(dat[[morph]]) + (1|Species) + (1|Site),
                     data = dat, adapt_delta = 0.99)
    r_square <- mean(bayes_R2(mod))

    phys_trait_df[ct,1] <- phys
    phys_trait_df[ct,2] <- morph
    phys_trait_df[ct,3] <- r_square

  }

  print(phys_trait_df)

}

phys_trait_df$model <- factor(rep("Phys-Trait", times = nrow(phys_trait_df)))
colnames(phys_trait_df) <- c("Response", "Predictor", "R2", "Model")

#####################################################################################
## Second, let's look at the trait-size relationships (2 size x 9 traits = 18 R2s) ##
#####################################################################################

size_trait_df <- data.frame(matrix(nrow = (length(Size_Traits)*(length(Morph_Traits))), ncol = 3))
colnames(size_trait_df) <- c("Size", "Trait", "R2")

ct <- 0
for(i in 1:length(Size_Traits)){

  size <- Size_Traits[i]

  for(j in 1:length(Morph_Traits)){

    ct <- ct + 1
    morph <- Morph_Traits[j]
    mod <- stan_lmer(scale(dat[[size]]) ~ scale(dat[[morph]]) + (1|Species) + (1|Site),
                     data = dat, adapt_delta = 0.99)
    r_square <- mean(bayes_R2(mod))

    size_trait_df[ct,1] <- size
    size_trait_df[ct,2] <- morph
    size_trait_df[ct,3] <- r_square

  }

  print(size_trait_df)

}

size_trait_df$model <- factor(rep("Size-Trait", times = nrow(size_trait_df)))
colnames(size_trait_df) <- c("Response", "Predictor", "R2", "Model")


#################################################################################
## Third, let's look at the phys-size relationships (2 size x 6 phys = 12 R2s) ##
#################################################################################

size_phys_df <- data.frame(matrix(nrow = (length(Size_Traits)*(length(Performance_Traits))), ncol = 3))
colnames(size_phys_df) <- c("Response", "Predictor", "R2")

ct <- 0
for(i in 1:length(Size_Traits)){

  size <- Size_Traits[i]

  for(j in 1:length(Performance_Traits)){

    ct <- ct + 1
    phys <- Performance_Traits[j]
    mod <- stan_lmer(scale(dat[[size]]) ~ scale(dat[[phys]]) + (1|Species) + (1|Site),
                     data = dat, adapt_delta = 0.99)
    r_square <- mean(bayes_R2(mod))

    size_phys_df[ct,1] <- size
    size_phys_df[ct,2] <- phys
    size_phys_df[ct,3] <- r_square

  }

  print(size_phys_df)

}

size_phys_df$Model <- factor(rep("Size-Phys", times = nrow(size_phys_df)))
colnames(size_phys_df) <- c("Response", "Predictor", "R2", "Model")



##############
## Plotting ##
##############

all <- rbind(phys_trait_df, size_trait_df, size_phys_df)


p1 <- ggplot(all, aes(x = R2, color = Model)) +
  geom_histogram(fill="white", alpha=0.5, position="identity", binwidth = 0.01)
p1 + theme_bw()



##########################################################################
## Check the Single Multivariate Models w/ Coefficients from Path Model ##
##########################################################################

height_str_lmer <- stan_lmer(scale(Height_in) ~ scale(BW_dry_avg) + scale(WD) + scale(Area_cm2) +
                          scale(LMA_g.cm2) + scale(LD_g.cm3) + scale(LWR) +
                          scale(Length_mm_Top) + scale(Density_mm2_Top) +
                          scale(SPI_Top)
                        + (1|Species),
                        data = dat, adapt_delta = 0.99)

summary(height_str_lmer, digits = 5, probs = c(0.05,0.1,0.9,0.95))
mean(bayes_R2(height_str_lmer))


canopy_str_lmer <- stan_lmer(scale(Canopy_area_in2) ~ scale(BW_dry_avg) + scale(WD) + scale(Area_cm2) +
                               scale(LMA_g.cm2) + scale(LD_g.cm3) + scale(LWR) +
                               scale(Length_mm_Top) + scale(Density_mm2_Top) +
                               scale(SPI_Top)
                             + (1|Species),
                             data = dat, adapt_delta = 0.99)

summary(canopy_str_lmer, digits = 5, probs = c(0.05,0.1,0.9,0.95))
mean(bayes_R2(canopy_str_lmer))

#################
## Size ~ Phys ##
#################

height_phys_lmer <- stan_lmer(scale(Height_in) ~ scale(Ks) + scale(KL_Units) + scale(TotalAssim) +
                                scale(Photo) + scale(Cond) + scale(WUE_Instantaneous) +
                                (1|Species), data = dat, adapt_delta = 0.99)

summary(height_phys_lmer, digits = 5, probs = c(0.05,0.1,0.9,0.95))
mean(bayes_R2(height_phys_lmer))


canopy_phys_lmer <- stan_lmer(scale(Canopy_area_in2) ~ scale(Ks) + scale(KL_Units) + scale(TotalAssim) +
                                scale(Photo) + scale(Cond) + scale(WUE_Instantaneous) +
                                (1|Species), data = dat, adapt_delta = 0.99)

summary(canopy_phys_lmer, digits = 5, probs = c(0.05,0.1,0.9,0.95))
mean(bayes_R2(canopy_phys_lmer))





tmp_mod <- stan_lmer(scale(Photo) ~ scale(LD_g.cm3) + (1|Site) + (1|Species), data = dat,
                     adapt_delta = 0.99)
mean(bayes_R2((tmp_mod)))
summary(tmp_mod, digits = 5, probs = c(0.05,0.1,0.9,0.95))

