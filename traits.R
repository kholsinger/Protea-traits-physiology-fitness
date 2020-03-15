## traits defined here for consistent sharing across scripts
##

if (STOMATES_STRUCT) {
  Performance_Traits <- c("Ks", "KL", "Photo", "TotAssim", "Cond", "wue")
  Size_Traits <- c("Height_in", "Canopy_area_in2")
  Morph_Traits <- c("BW_dry_avg", "WD", "Area_cm2", "LMA_g.cm2",
                    "LD_g.cm3", "LWR", "Length_mm_Top", "Density_mm2_Top",
                    "SPI_Top")
} else {
  Performance_Traits <- c("Ks", "KL", "Photo", "TotAssim", "Cond", "wue",
                          "Length_mm_Top", "Density_mm2_Top", "SPI_Top")
  Size_Traits <- c("Height_in", "Canopy_area_in2")
  Morph_Traits <- c("BW_dry_avg", "WD", "Area_cm2", "LMA_g.cm2",
                    "LD_g.cm3", "LWR")
}

