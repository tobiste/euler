## code to prepare `DATASET` dataset goes here
library(dplyr)
### Load data
plates <- readRDS("data/plates.rds")


frame <- rnaturalearth::ne_download(category = "physical", type = "wgs84_bounding_box", scale = "small", returnclass = "sf")
graticules <- rnaturalearth::ne_download(category = "physical", type = "graticules_30", scale = "small", returnclass = "sf")
world <- rnaturalearth::ne_countries(returnclass = "sf", scale = 50)

gsrm <- readxl::read_xlsx("D:/Global_data/Recent_euler_poles.xlsx", sheet = "GSRM v.2.1") %>%
  mutate(plate.fix = "NNR") %>%
  select(-c(lat.IGS08, lon.IGS08, rate.IGS08, xx.IGS08, xy.IGS08, xz.IGS08, yy.IGS08, yz.IGS08, zz.IGS08, Source)) %>%
  rename(lat = lat.NNR, lon = lon.NNR, angle = rate.NNR)

morvel <- readxl::read_xlsx("D:/Global_data/Recent_euler_poles.xlsx", sheet = "NNR-MORVEL56") %>%
  mutate(
    plate.rot = ifelse(plate.rot == "nb", "af", plate.rot),
    plate.rot = ifelse(plate.rot == "sm", "so", plate.rot),
    plate.rot = toupper(plate.rot)
  ) %>%
  select(-c(Plate.model, rate.std, RMS, area)) %>%
  rename(angle = rate)

### Save data
usethis::use_data(plates, overwrite = TRUE)
usethis::use_data(frame, overwrite = TRUE)
usethis::use_data(graticules, overwrite = TRUE)
usethis::use_data(world, overwrite = TRUE)
usethis::use_data(gsrm, overwrite = TRUE)
usethis::use_data(morvel, overwrite = TRUE)
