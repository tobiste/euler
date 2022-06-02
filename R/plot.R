library(ggplot2)
library(sf)

frame <- rnaturalearth::ne_download(category = 'physical', type = 'wgs84_bounding_box', scale = 'small', returnclass = 'sf')
graticules <- rnaturalearth::ne_download(category = 'physical', type = 'graticules_30', scale = 'small', returnclass = 'sf')
world <- rnaturalearth::ne_countries(returnclass = "sf", scale = 50)
data("PB2002", package = "tectonicr")

plates  <- readRDS('data/plates.rds')
in.plate <- plates %>%
  subset(Code == "IN")
som.plate <- plates %>%
  subset(Code == "SO")

in.eu <- c(27.12746847, 17.32482497, 0.402388191)
so.eu <- c(22.2078593, -92.40545103, 0.085835298)
in.eu.cart <- to_euler(in.eu)
so.eu.cart <- to_euler(so.eu)
in.so <- pole_migration(in.eu, so.eu)
so.in <- pole_migration(so.eu, in.eu)
in.so.pole <- relative_euler2(in.eu.cart, so.eu.cart)

gc <- common_greatcircle(in.eu, so.eu)

in.eu.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(in.eu[1], in.eu[2]), n = 90) %>% subset(n%%30 == 0)
so.eu.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(so.eu[1], so.eu[2]), n = 90) %>% subset(n%%30 == 0)
in.so.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(in.so.pole$axis.fin[1], in.so.pole$axis.fin[2]), n = 90) %>% subset(n%%30 == 0)

cgc <- common_greatcircle(in.eu, so.eu)
cgc.gc <- tectonicr::eulerpole_smallcircles(data.frame(lat = cgc[1], lon = cgc[2])) %>% subset(n == 90)

ggplot() +
  geom_sf(data = world, color = NA) +
  geom_sf(data = in.plate, aes(fill = 'India'), alpha=.5) +
  geom_sf(data = som.plate, aes(fill = 'Somalia'), alpha=.5) +
  geom_sf(data = graticules, color = 'grey80', lwd = .1) +
  geom_sf(data = PB2002, lwd = .25) +
  geom_sf(data = in.eu.sm, aes(color = 'India'), lty = 3) +
  geom_sf(data = so.eu.sm, aes(color = 'Somalia'), lty = 2) +
  geom_sf(data = in.so.sm, aes(color = 'India-Somalia'), lty = 4) +
  geom_sf(data = cgc.gc) +
  #borders(fill = 'grey90') +
  geom_point(aes(c(in.eu[2], in.eu[2]+180), c(in.eu[1], -in.eu[1]), color = 'India'), size = 2) +
  geom_point(aes(x = c(so.eu[2], so.eu[2]+180), y = c(so.eu[1], -so.eu[1]), color = 'Somalia'), size = 2) +
  geom_point(data = in.so, aes(axis.fin.lon, axis.fin.lat, color = 'India-Somalia'), shape = 1, size = 2) +
  geom_point(data = in.so, aes(axis.fin.lon+180, -axis.fin.lat, color = 'India-Somalia'), shape = 1, size = 2) +
  scale_color_discrete(name = 'Plate motion') +
  scale_fill_discrete(name = 'Plate') +


  ggnewscale::new_scale_color() +
  geom_point(data = in.so, aes(axis.inf.lon, axis.inf.lat, color = time), size = 2) +
  scale_color_viridis_c(name = 'Model time (Myr)') +

  geom_sf(data = frame, fill = NA, color = 'black', lwd = .5) +
  coord_sf(default_crs = st_crs(4326), crs = st_crs("ESRI:54030")) +


  ggthemes::theme_map() +
  labs(
    title = "Migration path of relative Euler pole",
    subtitle = '(Absolute motion is motion relative to fixed Eurasia)',
    caption = bquote(
    omega[In] == .(round(in.eu[3], 3)) |
    omega[So] == .(round(so.eu[3], 3)))
    )

