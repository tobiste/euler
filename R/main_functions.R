#' Relative Rotation associated to two absolute rotations
#'
#' Calculates Euler rotation axis and angle for two given absolute rotations using the
#' infinitesimal (Schaeben et al., 2021) and the finite approach (Cox \& Hart, 1989)
#'
#' @param x,y three-column vectors  giving the geographic coordinates latitude
#' and longitude, and the amount of rotation in degrees for first rotation (\code{x})
#' and subsequent second rotation (\code{y})
#' @param infinitesimal,finite logical. Should the rotation be calculated using
#' the infinitesimal and/or the finite rotation approach?
#' @note \code{x} will be considered as the fixed plate for the relative plate
#' motion between \code{x} and \code{y}
#' @references Schaeben, H., Kroner, U. and Stephan, T. (2021). Euler Poles of Tectonic
#' Plates. In B. S. Daza Sagar, Q. Cheng, J. McKinley and F. Agterberg (Eds.),
#' *Encyclopedia of Mathematical Geosciences. Encyclopedia of Earth Sciences Series*
#' (pp. 1--7). Springer Nature Switzerland AG 2021.
#' @returns \code{list}. Infinitesimal and the finite approach Euler axes
#' (geographical coordinates) and Euler angles (in degrees)
#' @export
#' @seealso [approximation_euler()] and [relative_euler_schaeben()] for quasi-infinitesimal and infinitesimal rotation, respectively.
#' @examples
#' x <- c(90, 0, 0.7)
#' y <- c(45, 30, 0.15)
#' relative_rotation(x, y)
relative_rotation <- function(x, y, infinitesimal = TRUE, finite = TRUE) {
  xe <- to_euler(x)
  ye <- to_euler(y)

  # "as-if-infinitesimal" rotation:\
  if (finite) {
    res.fin <- approximation_euler(xe, ye)
  }

  # infinitesimal rotation
  if (infinitesimal) {
    res.inf <- relative_euler_schaeben2(xe, ye)
  }

  if (infinitesimal & finite) {
    res <- list("infinitesimal" = res.inf, "finite" = res.fin)
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
#' @seealso [pole_migration_stats] for some additional statistics on the pole migration,
#' [relative_rotation()] for calculating the relative rotation,
#' [approximation_euler()] and [relative_euler_schaeben()] for quasi-infinitesimal and infinitesimal rotation, respectively.
#' @examples
#' in.eu <- c(27.12746847, 17.32482497, 0.402388191)
#' som.eu <- c(22.2078593, -92.40545103, 0.085835298)
#' pole_migration(som.eu, in.eu)
#' pole_migration(in.eu, som.eu)
#'
#' # example from Cox and Hart (1989):
#' naeu <- c(65.9, 132.4, 0.231)
#' na <- c(-58.3, 319.3, 0.247)
#' eu <- c(0.7, 336.8, 0.038)
#' pole_migration(eu, na, steps = c(1, seq(20, 200, 20))) # dplyr::%>% pole_migration_stats(eu, na)
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
          axis.inf.lat = rel.i$infinitesimal$axis[1],
          axis.inf.lon = rel.i$infinitesimal$axis[2],
          angle.inf = rel.i$infinitesimal$angle,
          axis.fin.lat = rel.i$finite$axis[1],
          axis.fin.lon = rel.i$finite$axis[2],
          angle.fin = rel.i$finite$angle
        )
      )
    } else if (infinitesimal & !finite) {
      res <- rbind(
        res,
        data.frame(
          time = i,
          axis.inf.lat = rel.i$axis[1],
          axis.inf.lon = rel.i$axis[2],
          angle.inf = rel.i$angle
        )
      )
    } else if (!infinitesimal & finite) {
      res <- rbind(
        res,
        data.frame(
          time = i,
          axis.fin.lat = rel.i$axis[1],
          axis.fin.lon = rel.i$axis[2],
          angle.fin = rel.i$angle
        )
      )
    } else {
      res <- NULL
    }
  }
  return(res)
}

#' Euler pole migration statistics
#'
#' Rates of the migration of the relative Euler pole associated with two absolute plate motions
#'
#' @param x \code{data.frame}. Output from \code{pole_migration}
#' @param euler1,euler2 three-column vectors  giving the geographic coordinates latitude
#' and longitude, and the amount of rotation in degrees for first rotation (\code{x})
#' and subsequent second rotation (\code{y})
#' @importFrom tectonicr geographical_to_cartesian angle_vectors deviation_norm
#' @importFrom magrittr %>%
#' @return data.frame containing magnitude of the vector between the
#' Euler poles \eqn{||d||} (\code{d}),
#' ... in degree \eqn{\eta} (\code{eta}),
#' and the change of their great circle distance in degree \eqn{\Delta} (\code{delta}).
#' @export
#' @seealso [pole_migration()] for calculating the migration of the relative Euler pole.
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
  gc <- c()
  for (i in seq_along(x$time)) {
    e.mep.i.geo <- c(x$axis.inf.lat[i], x$axis.inf.lon[i])
    e.mep.i <- tectonicr::geographical_to_cartesian(e.mep.i.geo) %>%
      normalize_vector()

    twe1 <- twe(x$time[i], euler1.cart)
    twe2 <- twe(x$time[i], euler2.cart)


    d[i] <- vector_norm((twe2 - twe1) - deg2rad(x$angle.inf[i]) * e.mep.i)

    eta[i] <- tectonicr::angle_vectors(
      normalize_vector(twe2 - twe1),
      e.mep.i
    ) %>%
      tectonicr::deviation_norm()


    gc[i] <- gc_distance(euler1, e.mep.i.geo)
    if (i == 1) gc0 <- gc[i]
  }
  return(data.frame(d = d, eta = eta, delta = gc - gc0))
}

#' Rotate a set a vector
#'
#' Rotate a vector, set of points, or polygon by a given rotation
#'
#' @inheritParams to_euler
#' @param p \code{sf} object
#' @param ... additional arguments to daughter functions
#' @importFrom tectonicr geographical_to_cartesian cartesian_to_geographical
#' @importFrom reticulate r_to_py py_to_r source_python
#' @importFrom magrittr %>%
#' @importFrom sf st_wrap_dateline st_geometry_type st_crs
#' @return \code{sf} object
#' @export
#' @examples
#' plates <- readRDS("data/plates.rds")
#' india <- subset(plates, Code == "IN")
#' euler <- c(90, 0, 90)
#' rotate_vector(euler, india) %>% plot()
rotate_vector <- function(x, p, ...) {
  stopifnot(inherits(p, "sf") & is.numeric(x))
  reticulate::source_python(system.file("python", "quaternions.py", package = "euler"), convert = FALSE)

  crs_p <- sf::st_crs(p)
  sf_class_p <- sf::st_geometry_type(p)
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
  q <- reticulate::r_to_py(euler) %>% py_euler2quat()

  w.cart <- py_rotate_vector_quat(u, q) %>%
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
  if (ncol(p) == 5) {
    p.rot.vec <- cbind(X = lons, Y = lats, L1 = p[, 3], L2 = p[, 4], L3 = p[, 5])
  } else if (ncol(p) == 4) {
    p.rot.vec <- cbind(X = lons, Y = lats, L1 = p[, 3], L2 = p[, 4])
  } else if (ncol(p) == 3) {
    p.rot.vec <- cbind(X = lons, Y = lats, L1 = p[, 3])
  } else if (ncol(p) == 2) {
    p.rot.vec <- cbind(X = lons, Y = lats, L1 = rep(1, length(lats)))
  }

  vector_to_sf(p.rot.vec, class = sf_class_p, ...) %>%
    sf::st_set_crs(crs_p) %>%
    sf::st_wrap_dateline(quiet = TRUE)
}


#' Plate rotation
#'
#' Rotate sf object using a sequence of rotations
#'
#' @param x Rotation sequence. \code{data.frame} with time, lat, lon, and angle
#' @param p \code{sf} object. Plate with location at 0 Ma
#' @param ... additional arguments to daughter functions
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @export
plate_rotation <- function(x, p, ...) {
  p.rot <- c()
  for (i in seq_along(x$time)) {
    euler.i <- c(x$lat[i], x$lon[i], x$angle[i])

    p.i <- rotate_vector(euler.i, p, ...)
    p.i$time <- x$time[i]
    p.rot <- rbind(
      p.rot,
      p.i
    )
  }
  return(p.rot %>% dplyr::group_by(time))
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
#' @importFrom dplyr filter select
#' @importFrom magrittr %>%
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
#' @returns list containing a plot, a gt table and \code{data.frame} with the results
#' @import gt
#' @import ggplot2
#' @import sf
#' @importFrom dplyr filter select mutate
#' @importFrom magrittr %>%
#' @importFrom tectonicr eulerpole_smallcircles euler_pole
#' @importFrom ggthemes scale_color_colorblind scale_fill_colorblind
#' @importFrom ggnewscale new_scale_color
#' @export
#' @examples
#' \dontrun{
#' quick_analysis("GSRM", "IN", "SO", "EU")
#' }
quick_analysis <- function(model = c("GSRM", "MORVEL"), plateA, plateB, fix) {
  data(plates)
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
  # a.b.pole.fin <- finite_euler(a.fix.cart, b.fix.cart)
  a.b.asisinf <- approximation_euler(a.fix.cart, b.fix.cart)

  a.b <- cbind(a.b, pole_migration_stats(a.b, a.fix, b.fix))

  # Small circles
  a.fix.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(a.fix[1], a.fix[2], angle = a.fix[3]), n = 6)
  b.fix.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(b.fix[1], b.fix[2], angle = b.fix[3]), n = 6)
  a.b.sm <- tectonicr::eulerpole_smallcircles(tectonicr::euler_pole(a.b.asisinf$axis[1], a.b.asisinf$axis[2], angle = a.b.asisinf$angle), n = 6)

  # Common great circle
  cgc <- common_greatcircle(a.fix, b.fix)
  cgc.gc <- tectonicr::eulerpole_smallcircles(data.frame(lat = cgc[1], lon = cgc[2])) %>% subset(d == 90)

  csc <- common_smallcircle(a.fix.cart, b.fix.cart)

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


  data(plates, package = "tectonicr")

  plot <- ggplot() +
    geom_sf(data = world, color = NA) +
    geom_sf(data = plates, color = "grey40", lwd = .25) +
    geom_sf(data = B.rotated.b.fix, aes(fill = plateB), color = NA, alpha = .1) +
    geom_sf(data = A.rotated.a.fix, aes(fill = plateA), color = NA, alpha = .1) +
    geom_sf(data = A.rotated.a.b, aes(fill = paste0(plateA, "-", plateB)), color = NA, alpha = .1) +
    geom_sf(data = A.plate, aes(fill = plateA), alpha = .5) +
    geom_sf(data = B.plate, aes(fill = plateB), alpha = .5) +
    geom_sf(data = graticules, color = "grey80", lwd = .1) +
    geom_sf(data = a.fix.sm, aes(color = plateA), lty = 3) +
    geom_sf(data = b.fix.sm, aes(color = plateB), lty = 2) +
    geom_sf(data = a.b.sm, aes(color = paste0(plateA, "-", plateB)), lty = 4) +
    geom_sf(data = csc, aes(color = plateA)) +
    geom_sf(data = cgc.gc, color = "grey") +

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
    select(time, angle.a, angle.b, axis.inf.lat, axis.inf.lon, angle.inf, axis.fin.lat, axis.fin.lon, angle.fin, d, eta, delta) %>%
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
      eta = html("&eta;(&tau;<sub>l</sub>)"),
      delta = html("&Delta;(&tau;<sub>l</sub>)")
    ) %>%
    fmt_number(
      columns = c(axis.inf.lat, axis.inf.lon, axis.fin.lat, axis.fin.lon),
      decimals = 2
    ) %>%
    fmt_number(
      columns = c(angle.a, angle.b, angle.inf, angle.fin, d, eta, delta),
      decimals = 3
    ) %>%
    tab_source_note(
      source_note = html("Acronyms: &tau; - Time [Myr]; &phi; - Latitude [&deg;N]; &lambda; - Longitude [&deg;E]; &omega; - Angle of rotation [&deg;] (rotation angles are positive CCW); ||d(&tau;<sub>l</sub>)|| - difference of rotation vectors; &eta;(&tau;<sub>l</sub>) - angle of the Euler poles [&deg;]: &Delta; - change of the angle between the Euler poles [&deg;].")
    )

  list(
    plot = plot,
    table = table,
    data = a.b
  )
}
