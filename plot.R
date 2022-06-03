library(ggplot2)
library(sf)
library(gt)

frame <- rnaturalearth::ne_download(category = "physical", type = "wgs84_bounding_box", scale = "small", returnclass = "sf")
graticules <- rnaturalearth::ne_download(category = "physical", type = "graticules_30", scale = "small", returnclass = "sf")
world <- rnaturalearth::ne_countries(returnclass = "sf", scale = 50)
data("PB2002", package = "tectonicr")

plates <- readRDS("data/plates.rds")
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

in.so <- cbind(in.so, pole_migration_stats(in.so, in.eu, so.eu))

# Small circles
in.eu.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(in.eu[1], in.eu[2]), n = 90) %>% subset(n %% 30 == 0)
so.eu.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(so.eu[1], so.eu[2]), n = 90) %>% subset(n %% 30 == 0)
in.so.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(in.so.pole$axis.fin[1], in.so.pole$axis.fin[2]), n = 90) %>% subset(n %% 30 == 0)

# Common great circle
cgc <- common_greatcircle(in.eu, so.eu)
cgc.gc <- tectonicr::eulerpole_smallcircles(data.frame(lat = cgc[1], lon = cgc[2])) %>% subset(n == 90)

# Rotated plates

india.rotated.ineu <- plate_rotation(
  x = data.frame(time = in.so$time,
                 lat = in.eu[1],
                 lon = in.eu[2],
                 angle = in.so$time*in.eu[3]),
  p = in.plate
) %>% arrange(desc(time))

india.rotated.inso <- plate_rotation(
  x = data.frame(time = in.so$time,
                 lat = in.so$axis.inf.lat,
                 lon = in.so$axis.inf.lon,
                 angle = in.so$angle.inf),
  p = in.plate
) %>% arrange(desc(time))

somalia.rotated.soeu <- plate_rotation(
  x = data.frame(time = in.so$time,
                 lat = so.eu[1],
                 lon = so.eu[2],
                 angle = in.so$time*so.eu[3]),
  p = som.plate
) %>% arrange(desc(time))



# Plot --------------------
ggplot() +
  geom_sf(data = world, color = NA) +
  geom_sf(data = PB2002, color = "grey40", lwd = .25) +

  geom_sf(data = somalia.rotated.soeu, aes(fill = "Somalia"), color = NA, alpha = .1) +
  geom_sf(data = india.rotated.ineu, aes(fill = "India"),  color = NA, alpha = .1) +
  geom_sf(data = india.rotated.inso, aes(fill = "India-Somalia"), color = NA, alpha = .1) +
  geom_sf(data = in.plate, aes(fill = "India"), alpha = .5) +
  geom_sf(data = som.plate, aes(fill = "Somalia"), alpha = .5) +

  geom_sf(data = graticules, color = "grey80", lwd = .1) +

  geom_sf(data = in.eu.sm, aes(color = "India"), lty = 3) +
  geom_sf(data = so.eu.sm, aes(color = "Somalia"), lty = 2) +
  geom_sf(data = in.so.sm, aes(color = "India-Somalia"), lty = 4) +
  geom_sf(data = cgc.gc) +
  # borders(fill = 'grey90') +
  geom_point(aes(c(in.eu[2], in.eu[2] + 180), c(in.eu[1], -in.eu[1]), color = "India"), size = 3) +
  geom_point(aes(x = c(so.eu[2], so.eu[2] + 180), y = c(so.eu[1], -so.eu[1]), color = "Somalia"), size = 3) +
  geom_point(data = in.so, aes(axis.fin.lon, axis.fin.lat, color = "India-Somalia"), shape = 1, size = 3) +
  geom_point(data = in.so, aes(axis.fin.lon + 180, -axis.fin.lat, color = "India-Somalia"), shape = 1, size = 3) +
  ggthemes::scale_color_colorblind(name = "Plate motion") +
  ggthemes::scale_fill_colorblind(name = "Plate") +

  ggnewscale::new_scale_color() +
  geom_point(data = in.so, aes(axis.inf.lon, axis.inf.lat, color = time), size = 2) +
  geom_point(data = in.so, aes(axis.inf.lon +180, -axis.inf.lat, color = time), size = 2) +
  scale_color_viridis_c(name = expression(Time ~ tau ~ (l) ~ (Myr)), limits = c(0, 300)) +
  geom_sf(data = frame, fill = NA, color = "black", lwd = .5) +
  coord_sf(default_crs = st_crs(4326), crs = st_crs("ESRI:54030")) +
  #ggthemes::theme_map() +
  theme_void() +
  labs(
    title = "Migration path of relative Euler pole",
    subtitle = "(Absolute motion is motion relative to fixed Eurasia)",
    caption = bquote("Present-day rotation rates:"~
      omega["In"] == .(round(in.eu[3], 3)) ~ degree ~ "Myr"^-1 ~ "|" ~
      omega["So"] == .(round(so.eu[3], 3))~degree ~ "Myr"^-1
    )
  )
ggsave("output/map1.png", bg = "white", width = 8.5, height = 5, scale = 1.5)

# Table --------------------
in.so %>%
  mutate(
    angle.in = time * in.eu[3],
    angle.so = time * so.eu[3]
    ) %>%
  select(time, angle.in, angle.so, axis.inf.lat, axis.inf.lon, angle.inf, axis.fin.lat, axis.fin.lon, angle.fin, d, eta) %>%
  gt() %>%
  tab_header(
    subtitle = "Comparison of 'as-if-infinitesimal' composition of rotations and proper
    concatenation of finite rotations of motion between India (In) and Somalia (So)",
    title = "Euler pole migration"
  ) %>%
  tab_spanner(
    label = "Proper concatenation",
    columns = c(axis.inf.lat, axis.inf.lon, angle.inf)
  ) %>%
  tab_spanner(
    label = "'as-if-infinitesimal' composition",
    columns = c(axis.fin.lat, axis.fin.lon, angle.fin)
  ) %>%
  cols_label(
    time = html("&tau;<sub>l</sub>"),
    angle.in = html("&omega;<sub>In</sub>(&tau;<sub>l</sub>)"),
    angle.so = html("&omega;<sub>So</sub>(&tau;<sub>l</sub>)"),
    axis.inf.lat = html("e<sub>&phi;</sub>(&tau;<sub>l</sub>)"),
    axis.inf.lon = html("e<sub>&lambda;</sub>(&tau;<sub>l</sub>)"),
    angle.inf = html("&omega;<sub>In-So</sub>(&tau;<sub>l</sub>)"),
    axis.fin.lat = html("e<sub>&phi;</sub>"),
    axis.fin.lon = html("e<sub>&lambda;</sub>"),
    angle.fin = html("&omega;<sub>So</sub>(&tau;<sub>l</sub>) - &omega;<sub>In</sub>(&tau;<sub>l</sub>)"),
    d = html("||d(&tau;<sub>l</sub>)||"),
    eta = html("&eta;(&tau;<sub>l</sub>)")
  ) %>%
  fmt_number(
    columns = c(axis.inf.lat, axis.inf.lon, axis.fin.lat, axis.fin.lon),
    decimals = 2
  ) %>%
  fmt_number(
    columns = c(angle.in, angle.so, angle.inf, angle.fin, d, eta),
    decimals = 3
  ) %>%
  tab_source_note(
    source_note = html("Acronyms: &tau; - Time [Myr]; &phi; - Latitude [&deg;N]; &lambda; - Longitude [&deg;E]; &omega; - Angle of rotation [&deg;] (rotation angles are positive CCW); ||d(&tau;<sub>l</sub>)|| - difference of rotation vectors; &eta;(&tau;<sub>l</sub>) - angle of the Euler poles [&deg;].")
  ) %>%
  gtsave("output/tab1.html")
