reticulate::py_run_file(system.file("python", "quaternions.py", package = "euler"), convert = FALSE)
# reticulate::source_python("inst/python/quaternions.py", convert = FALSE)

#' Vector norm
#'
#' Vector norm or magnitude
#'
#' @param v Numeric vector
vector_norm <- function(v) sqrt(sum(v^2))

#' Normalization of a vector
#'
#' normalizes a vector to unit length (length = 1)
#'
#' @inheritParams vector_norm
normalize_vector <- function(v) v / vector_norm(v)


#' Euler class
#'
#' Converts Euler pole in from geographic to Cartesian coordinates and Euler
#' angle from degrees to radians
#'
#' @param x three-column vector or 3*n matrix of the
#' geographic coordinates latitude and longitude, and the amount of rotation in
#' degrees
#' @importFrom tectonicr geographical_to_cartesian
#' @importFrom dplyr %>%
#' @export
#' @examples
#' euler1 <- c(90, 0, 10)
#' to_euler(euler1)
to_euler <- function(x) {
  stopifnot(is.numeric(x) & length(x) == 3)
  cart <- tectonicr::geographical_to_cartesian(c(x[1], x[2])) %>%
    normalize_vector()
  angle <- x[3] * pi / 180
  c(x = cart[1], y = cart[2], z = cart[3], angle = angle)
}


#' Relative Rotation associated to two absolute rotations
#'
#' Calculates Euler rotation axis and angle for two given absolute rotations using the
#' infinitesimal (Schaeben et al. 2021) and the finite approach
#'
#' @param x,y three-column vectors  giving the geographic coordinates latitude
#' and longitude, and the amount of rotation in degrees for first rotation (\code{x})
#' and subsequent second rotation (\code{y})
#' @param infinitesimal,finite logical. Should the rotation be calculated using
#' the infinitesimal and/or the finite rotation approach?
#' @note \code{x} will be considered as the fixed plate for the relative plate
#' motion between \code{x} and \code{y}
#' @references Schaeben, H., Kroner, U., Stephan, T. (2021). Euler Poles
#' of Tectonic Plates. In B. S. Daza Sagar, Q. Cheng, J. McKinley,; F. Agterberg
#' (Eds.), Encyclopedia of Mathematical Geosciences. Encyclopedia of Earth
#' Sciences Series (pp. 1--7). Springer Nature Switzerland AG 2021.
#' \doi{10.1007/978-3-030-26050-7_435-1}
#' @returns \code{list}. Infinitesimal and the finite approach Euler axes
#' (geographical coordinates) and Euler angles (in degrees)
#' @export
#' @examples
#' x <- c(90, 0, 0.7)
#' y <- c(45, 30, 0.15)
#' relative_rotation(x, y)
relative_rotation <- function(x, y, infinitesimal = TRUE, finite = TRUE) {
  xe <- to_euler(x)
  ye <- to_euler(y)

  # "as-if-infinitesimal" Rotation axis:
  if (finite) {
    res.fin <- as_if_infinitesimal_euler(xe, ye)
  }

  # transform to py
  if (infinitesimal) {
    res.inf <- infinitesimal_quaternion(xe, ye)
  }

  if (infinitesimal & finite) {
    res <- append(res.inf, res.fin)
  } else if (infinitesimal & !finite) {
    res <- res.inf
  } else if (!infinitesimal & finite) {
    res <- res.fin
  } else {
    res <- NULL
  }

  return(res)
}


#' Euler pole migration
#'
#' Migration of the relative Euler pole associated with two absolute plate motions
#'
#' @inheritParams relative_rotation
#' @param steps numeric vector of time increments. The default is a sequence
#' from 1 to 10 by an incremental step of 1 (e.g. Myr)
#' @note \code{x} is considered to be the "fixed" for the relative motion
#' between \code{x} and \code{y}
#' @export
#' @examples
#' in.eu <- c(27.12746847, 17.32482497, 0.402388191)
#' som.eu <- c(22.2078593, -92.40545103, 0.085835298)
#' pole_migration(som.eu, in.eu)
#' pole_migration(in.eu, som.eu)
pole_migration <- function(x, y, steps = c(1, seq(25, 300, 25)), infinitesimal = TRUE, finite = TRUE) {
  res <- data.frame(time = NULL, axis.inf.lat = NULL, axis.inf.lon = NULL, angle.inf = NULL, axis.fin.lat = NULL, axis.fin.lon = NULL, angle.fin = NULL)
  rate.x <- x[3]
  rate.y <- y[3]
  for (i in steps) {
    x[3] <- rate.x * i
    y[3] <- rate.y * i

    # y is fix, x is rotating
    rel.i <- relative_rotation(y, x, infinitesimal, finite)

    if (infinitesimal & finite) {
      res <- rbind(
        res,
        data.frame(
          time = i,
          axis.inf.lat = rel.i$axis.inf[1],
          axis.inf.lon = rel.i$axis.inf[2],
          angle.inf = rel.i$angle.inf,
          axis.fin.lat = rel.i$axis.fin[1],
          axis.fin.lon = rel.i$axis.fin[2],
          angle.fin = rel.i$angle.fin
        )
      )
    } else if (infinitesimal & !finite) {
      res <- rbind(
        res,
        data.frame(
          time = i,
          axis.inf.lat = rel.i$axis.inf[1],
          axis.inf.lon = rel.i$axis.inf[2],
          angle.inf = rel.i$angle.inf
        )
      )
    } else if (!infinitesimal & finite) {
      res <- rbind(
        res,
        data.frame(
          time = i,
          axis.fin.lat = rel.i$axis.fin[1],
          axis.fin.lon = rel.i$axis.fin[2],
          angle.fin = rel.i$angle.fin
        )
      )
    } else {
      res <- NULL
    }
  }
  return(res)
}

#' Great circle of Euler poles
#'
#' Calculates the pole to the great circle of Euler poles
#'
#' @param x,y two-column vectors giving the geographic coordinates latitude
#' and longitude in degrees
common_greatcircle <- function(x, y) {
  x.cart <- to_euler(x)
  y.cart <- to_euler(y)
  pracma::cross(
    c(x.cart[[1]], x.cart[[2]], x.cart[[3]]), c(y.cart[[1]], y.cart[[2]], y.cart[[3]])
  ) %>% tectonicr::cartesian_to_geographical()
}

#' Rotate a set a vector
#'
#' Rotate a vector, set of points, or polygon by a given rotation
#'
#' @inheritParams to_euler
#' @param p \code{sf} object
#' @importFrom tectonicr geographical_to_cartesian cartesian_to_geographical
#' @importFrom reticulate r_to_py py_to_r
#' @importFrom dplyr %>%
#' @importFrom sf st_wrap_dateline
#' @return \code{sf} object
#' @export
#' @examples
#' plates <- readRDS("data/plates.rds")
#' india <- subset(plates, Code == "IN")
#' euler <- c(90, 0, 90)
#' rotate_vector(euler, india) %>% plot()
rotate_vector <- function(x, p) {
  stopifnot(("sf" %in% class(p)) & is.numeric(x))

  crs <- sf::st_crs(p)

  p <- sf_to_vector(p)
  lons <- p[, 1]
  lats <- p[, 2]
  p[, 1] <- lats
  p[, 2] <- lons

  p.geo <- c()
  for (i in seq_along(p[, 1])) {
    p.geo <- rbind(
      p.geo,
      tectonicr::geographical_to_cartesian(c(p[i, 1], p[i, 2]))
    )
  }

  u <- p.geo %>%
    reticulate::r_to_py()

  euler <- to_euler(x)
  q <- reticulate::r_to_py(euler) %>% euler2quat()

  w.cart <- rotate_vector_quat(u, q) %>%
    reticulate::py_to_r()


  w <- c()
  for (i in seq_along(w.cart[, 1])) {
    w <- rbind(
      w,
      tectonicr::cartesian_to_geographical(c(w.cart[i, 1], w.cart[i, 2], w.cart[i, 3]))
    )
  }

  lats <- w[, 1]
  lons <- w[, 2]

  p.rot <- cbind(X = lons, Y = lats, L1 = p[, 3], L2 = p[, 4]) %>%
    vector_to_sf() %>%
    sf::st_wrap_dateline()
  sf::st_crs(p.rot) <- sf::st_crs(crs)
  return(p.rot)
}

#' sf object to vector
#'
#' Converts a sf object into a two-column vector
#' @param x \code{sf} object
#' @importFrom sf st_coordinates st_as_sf st_sfc st_polygon
#' @importFrom dplyr %>%
#' @examples
#' readRDS("data/plates.rds")
#' in.plate <- subset(plates, Code == "IN")
#' sf_to_vector(in.plate)
sf_to_vector <- function(x) {
  sf::st_as_sf(x) %>% sf::st_coordinates()
}

#' vector to sf object
#'
#' Converts a four-column vector into a sf object into
#' @param x vector
#' @param multi logical. Is x a MULTIPOLYGON or MULTILINESTRING.
#' @importFrom sf st_cast st_polygon st_sfc st_sf st_as_sf
#' @importFrom dplyr %>% group_by summarise
#' @examples
#' readRDS("data/plates.rds")
#' in.plate <- subset(plates, Code == "IN")
#' sf_to_vector(in.plate) %>% vector_to_sf()
vector_to_sf <- function(x, multi = FALSE) {
  if (multi) {
    x %>%
      st_as_sf(coords = c("X", "Y")) %>%
      group_by(L1, L2) %>%
      summarise(do_union = FALSE) %>%
      st_cast("POLYGON")
  } else {
    st_polygon(x = list(x[, 1:2])) %>%
      st_sfc() %>%
      st_sf()
  }
}

#' Rotate plate using a set of migrated poles
#'
#' @param x \code{data.frame} with time, lat, lon, and angle
#' @param p \code{sf} object. Plate with location at 0 Ma
#' @importFrom dplyr %>% group_by
#' @export
plate_rotation <- function(x, p) {
  p.rot <- c()
  for (i in seq_along(x$time)) {
    euler.i <- c(x$lat[i], x$lon[i], x$angle[i])

    p.i <- rotate_vector(euler.i, p)
    p.i$time <- x$time[i]
    p.rot <- rbind(
      p.rot,
      p.i
    )
  }
  return(p.rot %>% dplyr::group_by(time))
}

twe <- function(time, e) {
  time * e[4] * c(e[1], e[2], e[3])
}


#' Euler pole migration statistics
#'
#' Migration of the relative Euler pole associated with two absolute plate motions
#'
#' @param x \code{data.frame}. Output from \code{pole_migration}
#' @param euler1,euler2 three-column vectors  giving the geographic coordinates latitude
#' and longitude, and the amount of rotation in degrees for first rotation (\code{x})
#' and subsequent second rotation (\code{y})
#' @importFrom tectonicr geographical_to_cartesian angle_vectors deviation_norm
#' @export
#' @examples
#' in.eu <- c(27.12746847, 17.32482497, 0.402388191)
#' som.eu <- c(22.2078593, -92.40545103, 0.085835298)
#' pm <- pole_migration(som.eu, in.eu)
#' pole_migration_stats(pm, som.eu, in.eu)
pole_migration_stats <- function(x, euler1, euler2) {
  euler1.cart <- to_euler(euler1)
  euler2.cart <- to_euler(euler2)

  d <- c()
  eta <- c()
  for (i in seq_along(x$time)) {
    e.mep.i <- tectonicr::geographical_to_cartesian(c(x$axis.inf.lat[i], x$axis.inf.lon[i])) %>%
      normalize_vector()

    twe1 <- twe(x$time[i], euler1.cart)
    twe2 <- twe(x$time[i], euler2.cart)


    d[i] <- vector_norm((twe2 - twe1) - x$angle.inf[i] * e.mep.i)

    eta[i] <- tectonicr::angle_vectors(
      normalize_vector(twe2 - twe1),
      e.mep.i
    ) %>%
      tectonicr::deviation_norm()
  }
  return(cbind(d = d, eta = eta))
}

#' Plate motion model
#'
#' Load global plate motions
#'
#' @param model Model to choose from. Either \code{"GSRM"} for the "Global
#' Strain Rate Model" v2.1 by Kreemer et all. 2014, or \code{"MORVEL"} for the
#' "NNR-MORVEL56" model by Argus et al. 2011
#' @param plateA,plateB plates to extract
#' @param fix Reference system that is considered to be fixed.
#' @importFrom dplyr %>% filter select
#' @importFrom tectonicr equivalent_rotation
#' @export
#' @examples
#' load_plate_motions(model = "GSRM", plateA = "IN", plateB = "SO", fix = "EU")
load_plate_motions <- function(model = c("GSRM", "MORVEL"), plateA, plateB, fix) {
  stopifnot(is.character(plateA) & is.character(plateB) & is.character(fix))

  model <- match.arg(model)
  if (model == "GSRM") {
    load("data/gsrm.rda")
    rots <- gsrm
  } else {
    load("data/morvel.rda")
    rots <- morvel
  }

  tectonicr::equivalent_rotation(rots, fixed = fix) %>%
    dplyr::filter(plate.rot %in% c(plateA, plateB)) %>%
    dplyr::select(lat, lon, angle, plate.rot, plate.fix)
}



#' Quick analysis
#'
#' Analysis the data and returns a ggplot and a gt table
#'
#' @inheritParams load_plate_motions
#' @returns list containing a plot, a gt table and data.frame with the results
#' @import gt
#' @import ggplot2
#' @import sf
#' @importFrom dplyr %>% filter select
#' @importFrom tectonicr eulerpole_smallcircles euler_pole
#' @importFrom ggthemes scale_color_colorblind scale_fill_colorblind
#' @importFrom ggnewscale new_scale_color
#' @export
#' @examples
#' quick_analysis("GSRM", "IN", "SO", "EU")
quick_analysis <- function(model = c("GSRM", "MORVEL"), plateA, plateB, fix) {
  A.plate <- plates %>%
    subset(Code == plateA)
  B.plate <- plates %>%
    subset(Code == plateB)

  rots <- load_plate_motions(model, plateA, plateB, fix)
  a.fix <- rots %>%
    filter(plate.rot == plateA) %>%
    select(lat, lon, angle) %>%
    as.numeric()
  b.fix <- rots %>%
    filter(plate.rot == plateB) %>%
    select(lat, lon, angle) %>%
    as.numeric()

  a.fix.cart <- to_euler(a.fix)
  b.fix.cart <- to_euler(b.fix)
  a.b <- pole_migration(a.fix, b.fix)
  b.a <- pole_migration(b.fix, a.fix)
  a.b.pole.fin <- finite_euler(a.fix.cart, b.fix.cart)
  a.b.asisinf <- as_if_infinitesimal_euler(a.fix.cart, b.fix.cart)

  a.b <- cbind(a.b, pole_migration_stats(a.b, a.fix, b.fix))

  # Small circles
  a.fix.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(a.fix[1], a.fix[2]), n = 90) %>% subset(n %% 30 == 0)
  b.fix.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(b.fix[1], b.fix[2]), n = 90) %>% subset(n %% 30 == 0)
  a.b.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(a.b.asisinf$axis.fin[1], a.b.asisinf$axis.fin[2]), n = 90) %>% subset(n %% 30 == 0)

  # Common great circle
  cgc <- common_greatcircle(a.fix, b.fix)
  cgc.gc <- tectonicr::eulerpole_smallcircles(data.frame(lat = cgc[1], lon = cgc[2])) %>% subset(n == 90)

  # Rotated plates

  A.rotated.a.fix <- plate_rotation(
    x = data.frame(
      time = a.b$time,
      lat = a.fix[1],
      lon = a.fix[2],
      angle = a.b$time * a.fix[3]
    ),
    p = A.plate
  ) %>% arrange(desc(time))

  A.rotated.a.b <- plate_rotation(
    x = data.frame(
      time = a.b$time,
      lat = a.b$axis.inf.lat,
      lon = a.b$axis.inf.lon,
      angle = a.b$angle.inf
    ),
    p = A.plate
  ) %>% arrange(desc(time))

  B.rotated.b.fix <- plate_rotation(
    x = data.frame(
      time = a.b$time,
      lat = b.fix[1],
      lon = b.fix[2],
      angle = a.b$time * b.fix[3]
    ),
    p = B.plate
  ) %>% arrange(desc(time))



  plot <- ggplot() +
    geom_sf(data = world, color = NA) +
    geom_sf(data = PB2002, color = "grey40", lwd = .25) +
    geom_sf(data = B.rotated.b.fix, aes(fill = plateB), color = NA, alpha = .1) +
    geom_sf(data = A.rotated.a.fix, aes(fill = plateA), color = NA, alpha = .1) +
    geom_sf(data = A.rotated.a.b, aes(fill = paste0(plateA, "-", plateB)), color = NA, alpha = .1) +
    geom_sf(data = A.plate, aes(fill = plateA), alpha = .5) +
    geom_sf(data = B.plate, aes(fill = plateB), alpha = .5) +
    geom_sf(data = graticules, color = "grey80", lwd = .1) +
    geom_sf(data = a.fix.sm, aes(color = plateA), lty = 3) +
    geom_sf(data = b.fix.sm, aes(color = plateB), lty = 2) +
    geom_sf(data = a.b.sm, aes(color = paste0(plateA, "-", plateB)), lty = 4) +
    geom_sf(data = cgc.gc) +
    # borders(fill = 'grey90') +
    geom_point(aes(c(a.fix[2], a.fix[2] + 180), c(a.fix[1], -a.fix[1]), color = plateA), size = 3) +
    geom_point(aes(x = c(b.fix[2], b.fix[2] + 180), y = c(b.fix[1], -b.fix[1]), color = plateB), size = 3) +
    geom_point(data = a.b, aes(axis.fin.lon, axis.fin.lat, color = paste0(plateA, "-", plateB)), shape = 1, size = 3) +
    geom_point(data = a.b, aes(axis.fin.lon + 180, -axis.fin.lat, color = paste0(plateA, "-", plateB)), shape = 1, size = 3) +
    ggthemes::scale_color_colorblind(name = "Plate motion") +
    ggthemes::scale_fill_colorblind(name = "Plate") +
    ggnewscale::new_scale_color() +
    geom_point(data = a.b, aes(axis.inf.lon, axis.inf.lat, color = time), size = 2) +
    geom_point(data = a.b, aes(axis.inf.lon + 180, -axis.inf.lat, color = time), size = 2) +
    scale_color_viridis_c(name = expression(Time ~ tau ~ (l) ~ (Myr)), limits = c(0, 300)) +
    geom_sf(data = frame, fill = NA, color = "black", lwd = .5) +
    coord_sf(default_crs = st_crs(4326), crs = st_crs("ESRI:54030")) +
    # ggthemes::theme_map() +
    theme_void() +
    labs(
      title = "Migration path of relative Euler pole",
      subtitle = paste0("(Absolute motion is motion relative to ", fix, ")"),
      caption = bquote("Present-day rotation rates:" ~
      omega[.(plateA)] == .(round(a.fix[3], 3)) ~ degree ~ "Myr"^-1 ~ "|" ~
      omega[.(plateB)] == .(round(b.fix[3], 3)) ~ degree ~ "Myr"^-1)
    )

  table <- a.b %>%
    mutate(
      angle.a = time * a.fix[3],
      angle.b = time * b.fix[3]
    ) %>%
    select(time, angle.a, angle.b, axis.inf.lat, axis.inf.lon, angle.inf, axis.fin.lat, axis.fin.lon, angle.fin, d, eta) %>%
    gt() %>%
    tab_header(
      subtitle = paste0("Comparison of 'as-if-infinitesimal' composition of rotations and proper
    concatenation of finite rotations of motion between ", plateA, " and ", plateB),
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
      angle.a = html(paste0("&omega;<sub>", plateA, "</sub>(&tau;<sub>l</sub>)")),
      angle.b = html(paste0("&omega;<sub>", plateB, "</sub>(&tau;<sub>l</sub>)")),
      axis.inf.lat = html("e<sub>&phi;</sub>(&tau;<sub>l</sub>)"),
      axis.inf.lon = html("e<sub>&lambda;</sub>(&tau;<sub>l</sub>)"),
      angle.inf = html(paste0("&omega;<sub>", plateA, "-", plateB, "</sub>(&tau;<sub>l</sub>)")),
      axis.fin.lat = html("e<sub>&phi;</sub>"),
      axis.fin.lon = html("e<sub>&lambda;</sub>"),
      angle.fin = html(paste0("&omega;<sub>", plateB, "</sub>(&tau;<sub>l</sub>) - &omega;<sub>", plateA, "</sub>(&tau;<sub>l</sub>)")),
      d = html("||d(&tau;<sub>l</sub>)||"),
      eta = html("&eta;(&tau;<sub>l</sub>)")
    ) %>%
    fmt_number(
      columns = c(axis.inf.lat, axis.inf.lon, axis.fin.lat, axis.fin.lon),
      decimals = 2
    ) %>%
    fmt_number(
      columns = c(angle.a, angle.b, angle.inf, angle.fin, d, eta),
      decimals = 3
    ) %>%
    tab_source_note(
      source_note = html("Acronyms: &tau; - Time [Myr]; &phi; - Latitude [&deg;N]; &lambda; - Longitude [&deg;E]; &omega; - Angle of rotation [&deg;] (rotation angles are positive CCW); ||d(&tau;<sub>l</sub>)|| - difference of rotation vectors; &eta;(&tau;<sub>l</sub>) - angle of the Euler poles [&deg;].")
    )


  list(
    plot = plot,
    table = table,
    data = a.b
  )
}
