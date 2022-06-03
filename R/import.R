# global reference to scipy (will be initialized in .onLoad)
quaternion <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  quaternion <<- reticulate::import("quaternion", delay_load = TRUE)
}
reticulate::source_python("src/quaternions.py", convert = FALSE)

# Global
deg2rad <- pi / 180

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
#' @param x Vector for Euler pole position, three-column vector of the
#' geographic coordinates latitude and longitude, and the amount of rotation in
#' degrees
#' @importFrom tectonicr geographical_to_cartesian
#' @importFrom dplyr %>%
#' @export
#' @examples
#' euler1 <- to_euler(c(90, 0, 10))
#' euler2 <- to_euler(c(45, 30, 20))
to_euler <- function(x) {
  stopifnot(is.numeric(x) & length(x) == 3)
  cart <- tectonicr::geographical_to_cartesian(c(x[1], x[2])) %>%
    normalize_vector()
  angle <- x[3] * deg2rad
  data.frame(x = cart[1], y = cart[2], z = cart[3], angle)
}


#' As-if-infinitesimal rotation
#'
#' @inheritParams relative_quat
#' @importFrom dplyr %>%
#' @importFrom tectonicr cartesian_to_geographical
relative_euler <- function(r1, r2) {
  e <- (r1$angle * c(r1$x, r1$y, r1$z) + r2$angle * c(r2$x, r2$y, r2$z)) %>%
    tectonicr::cartesian_to_geographical()

  angle <- (r1$angle - r2$angle) / deg2rad

  list(
    axis.fin = e,
    angle.fin = angle
  )
}

relative_euler2 <- function(r1, r2) {
  r1.pole <- tectonicr::euler_pole(r1$x, r1$y, r1$z, geo = FALSE)
  r2.pole <- tectonicr::euler_pole(r2$x, r2$y, r2$z, geo = FALSE)

  r1.rot <- tectonicr::euler_rot(r1.pole, r1$angle / deg2rad)
  r2.rot <- tectonicr::euler_rot(r2.pole, r2$angle / deg2rad)

  e <- r1.rot %*% r2.rot %>% tectonicr::euler_from_rot()

  list(
    axis.fin = c(e$pole$lat, e$pole$lon),
    angle.fin = e$psi
  )
}



#' Infinitesimal rotation
#'
#' Infinitesimal rotation using Quaternions
#'
#' @param r1,r2 three-column vectors  giving the geographic coordinates latitude
#' and longitude, and the amount of rotation in degrees for first rotation (\code{r1})
#' and subsequent second rotation (\code{r2})
#' @importFrom reticulate r_to_py py_to_r
#' @importFrom dplyr %>%
#' @importFrom tectonicr cartesian_to_geographical
relative_quat <- function(r1, r2) {
  R1 <- reticulate::r_to_py(r1) %>% euler2quat()
  R2 <- reticulate::r_to_py(r2) %>% euler2quat()

  angle <- euler_angle(R1, R2) %>%
    reticulate::py_to_r()
  axis <- euler_axis(R1, R2, angle) %>%
    reticulate::py_to_r() %>%
    tectonicr::cartesian_to_geographical()

  list(
    axis.inf = axis,
    angle.inf = angle / deg2rad
  )
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
  xe <- to_euler(c(x[1], x[2], x[3]))
  ye <- to_euler(c(y[1], y[2], y[3]))

  # "as-if-infinitesimal" Rotation axis:
  if (finite) {
    res.fin <- relative_euler(xe, ye)
  }

  # transform to py
  if (infinitesimal) {
    res.inf <- relative_quat(xe, ye)
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
rotate_vector <- function(x, p){
  stopifnot(("sf" %in% class(p)) & is.numeric(x))

  crs <- sf::st_crs(p)

  p <- sf_to_vector(p)
  lons <- p[, 1]
  lats <- p[, 2]
  p[, 1] <- lats
  p[, 2] <- lons

  p.geo <- c()
  for(i in seq_along(p[, 1])){
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
  for(i in seq_along(w.cart[, 1])){
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
#' @importFrom sf st_coordinates st_as_sf
#' @importFrom dplyr %>%
#' @examples
#' readRDS("data/plates.rds")
#' in.plate <- subset(plates, Code == "IN")
#' sf_to_vector(in.plate)
sf_to_vector <- function(x){
  sf::st_as_sf(x) %>% sf::st_coordinates()
}

#' vector to sf object
#'
#' Converts a four-column vector into a sf object into
#' @param x vector
#' @param multi logical. Is x a MULTIPOLYGON or MULTILINESTRING.
#' @importFrom sf st_cast st_polygon st_sfc st_sf
#' @importFrom dplyr %>% group_by summarise
#' @examples
#' readRDS("data/plates.rds")
#' in.plate <- subset(plates, Code == "IN")
#' sf_to_vector(in.plate) %>% vector_to_sf()
vector_to_sf <- function(x, multi = FALSE){
  if(multi){
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
plate_rotation <- function(x, p){
  p.rot <- c()
  for(i in seq_along(x$time)){
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
  time * e[1, 4] * c(e[1, 1], e[1, 2], e[1, 3])
  }


#' Euler pole migration statistics
#'
#' Migration of the relative Euler pole associated with two absolute plate motions
#'
#' @param x \code{data.frame}. Output from \code{pole_migration}
#' @param euler1,euler2 three-column vectors  giving the geographic coordinates latitude
#' and longitude, and the amount of rotation in degrees for first rotation (\code{x})
#' and subsequent second rotation (\code{y})
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
  for(i in seq_along(x$time)){
    e.mep.i <- tectonicr::geographical_to_cartesian(c(x$axis.inf.lat[i], x$axis.inf.lon[i]))

    twe1 <- twe(x$time[i], euler1.cart)
    twe2 <- twe(x$time[i], euler2.cart)


    d[i] <- vector_norm((twe1 + twe1) - x$angle.inf[i] * e.mep.i)

    eta[i] <- tectonicr::angle_vectors(
      normalize_vector(twe1 + twe2),
      e.mep.i
    )
  }
  return(cbind(d = d, eta = eta))
}
