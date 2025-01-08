rad2deg <- function(x) {
  x * 180 / pi
}

deg2rad <- function(x) {
  x * pi / 180
}

#' Euclidean normalization of a vector
#'
#' magnitude or length of a vector
#'
#' @param v Numeric vector
#' @export
vector_norm <- function(v) sqrt(sum(v^2))

#' Normalization of a vector
#'
#' normalizes a vector to unit length (length = 1)
#'
#' @inheritParams vector_norm
#' @export
normalize_vector <- function(v) v / vector_norm(v)

#' Euler class
#'
#' Converts Euler pole in from geographic to Cartesian coordinates and Euler
#' angle from degrees to radians and back
#'
#' @param g three-column vector or 3*n matrix of the
#' geographic coordinates latitude and longitude, and the amount of rotation in
#' degrees
#' @param x Object of class \code{"euler"}, i.e. 4-column vector or 3*n matrix of the
#' Cartesian coordinates and the amount of rotation in radians
#' @importFrom tectonicr geographical_to_cartesian cartesian_to_geographical
#' @importFrom magrittr %>%
#' @name euler
#' @examples
#' euler1 <- c(90, 0, 10)
#' to_euler(euler1)
#' to_euler(euler1) %>% from_euler()
NULL

#' @rdname euler
#' @export
to_euler <- function(g) {
  stopifnot(is.numeric(g) & length(g) == 3)
  cart <- tectonicr::geographical_to_cartesian(c(g[1], g[2])) %>%
    normalize_vector()
  e <- c(x = cart[1], y = cart[2], z = cart[3], angle = deg2rad(g[3]))
  class(e) <- append(class(e), "euler")
  return(e)
}

#' @rdname euler
#' @export
from_euler <- function(x) {
  stopifnot(inherits(x, "euler"))
  names(x) <- NULL
  geo <- tectonicr::cartesian_to_geographical(c(x[1], x[2], x[3]))
  c(lat = geo[1], lon = geo[2], angle = rad2deg(x[4]))
}


#' Convert between sf objects and point vectors
#'
#' Converts a sf object into a two-column vector or vice versa
#' @param sf \code{sf} object
#' @param x numeric vector
#' @param multi logical. Is x MULTIPOLYGON or MULTILINESTRING.
#' @importFrom sf st_coordinates st_as_sf st_sfc st_polygon st_cast st_sf
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#' @name sf_conversion
#' @examples
#' readRDS("data/plates.rds")
#' in.plate <- subset(plates, Code == "IN")
#' sf_to_vector(in.plate) %>% vector_to_sf()
NULL

#' @rdname sf_conversion
#' @export
sf_to_vector <- function(sf) {
  sf::st_as_sf(sf) %>% sf::st_coordinates()
}

#' @rdname sf_conversion
#' @export
vector_to_sf <- function(x, class, multi = FALSE) {
  if (multi) {
    if (class == "POLYGON" & ncol(x) == 5) {
      x %>%
        as.data.frame() %>%
        st_as_sf(coords = c("X", "Y")) %>%
        group_by(L1, L2, L3) %>%
        summarise(do_union = FALSE) %>%
        st_cast("POLYGON")
    } else if (class == "POLYGON" & ncol(x) == 4) {
      x %>%
        as.data.frame() %>%
        st_as_sf(coords = c("X", "Y")) %>%
        group_by(L1, L2) %>%
        summarise(do_union = FALSE) %>%
        st_cast("POLYGON")
    } else if (class == "MULTILINESTRING" & ncol(x) == 3) {
      x %>%
        as.data.frame() %>%
        st_as_sf(coords = c("X", "Y")) %>%
        group_by(L1) %>%
        summarise(do_union = FALSE) %>%
        st_cast("MULTILINESTRING")
    } else if (class == "MULTIPOINT") {
      x %>%
        as.data.frame() %>%
        st_as_sf(coords = c("X", "Y")) %>%
        group_by(L1) %>%
        summarise(do_union = FALSE) %>%
        st_cast("POINT")
    }
  } else {
    if (class == "POLYGON") {
      st_polygon(x = list(x[, 1:2])) %>%
        st_sfc() %>%
        st_sf()
    } else if (class == "LINESTRING") {
      st_linestring(st_point(c(x[, 1], x[, 2]))) %>%
        st_sfc() %>%
        st_sf()
    } else if (class == "POINT") {
      st_point(c(x[, 1], x[, 2])) %>%
        st_sfc() %>%
        st_sf()
    }
  }
}

#' Inverse rotation
#'
#' @param x Object of class \code{"euler"}, i.e. 4-column vector or 3*n matrix of the
#' Cartesian coordinates and the amount of rotation in radians
#' @export
inverse_euler <- function(x) {
  stopifnot(inherits(x, "euler"))
  x[4] <- -x[4]
  return(x)
}

#' Helper function for stats
twe <- function(time, e) {
  time * e[4] * c(e[1], e[2], e[3])
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
  tectonicr::vcross(
    c(x.cart[[1]], x.cart[[2]], x.cart[[3]]), c(y.cart[[1]], y.cart[[2]], y.cart[[3]])
  ) %>% tectonicr::cartesian_to_geographical()
}

#' Common small circle of absolute and relative pole
#'
#' @inheritParams approximation_euler
common_smallcircle <- function(r1, r2) {
  r1.r2 <- relative_euler_schaeben(r1, r2)

  angle <- tectonicr::angle_vectors(
    tectonicr::geographical_to_cartesian(r1.r2$axis),
    r1[1:3]
  )

  r1.ep <- tectonicr::euler_pole(r1[1], r1[2], r1[3], geo = FALSE)

  sm_np <- data.frame(
    lon = c(seq(-180, 180, 2), seq(-180, 180, 2)),
    lat = c(rep(90 - angle, 181), rep(-90 + angle, 181))
  ) %>%
    st_as_sf(coords = c("lon", "lat")) %>%
    summarise(do_union = FALSE) %>%
    st_cast("MULTILINESTRING") %>%
    sf::st_wrap_dateline(
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
      quiet = TRUE
    )

  tectonicr::PoR_to_geographical_sf(x = sf::st_as_sf(sm_np), euler = r1.ep) %>%
    sf::st_wrap_dateline(
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
      quiet = TRUE
    )
}

#' Great circle distance
#'
#' Angle between two points on a sphere, measured along great circle passing
#' through both point vectors
#' @param a,b two-column vectors giving the geographic coordinates latitude
#' and longitude in degrees
gc_distance <- function(a, b) {
  a <- deg2rad(a)
  b <- deg2rad(b)

  acos(
    sin(a[1]) * sin(b[1]) + cos(a[1]) * cos(b[1]) * cos(abs(a[2] - b[2]))
  ) %>% rad2deg()
}


rotation_rate <- function(wx, wy, wz) {
  sqrt(wx^2 + wy^2 + wz^2)
}

#' Calculate the position and magnitude of the rotation pole from angular
#' velocity components.
#'
#' @param wx,wy,wz numeric. Cartesian velocity components (in any consistent
#' units).
#' @param lat,lon numeric.  Latitude and longitude of the Euler pole in degrees
#' @param mag numeric. Magnitude of rotation in any consistent units
#'
#' @return [euler_cart2geo()] returns three element vector containing the
#' geographic position of the Euler
#' pole and the rotation magnitudes in the units given by the Cartesian
#' components.
#'
#' [euler_geo2cart()] returns the Cartesian components of
#' the rotation around the unit vectors.
#'
#' @name cartesian_comp
#'
#' @examples
#' euler_cart2geo(wx = -0.001583, wy = 0.005064, wz = -0.010430) # should be: c(-63.036, 107.359, deg2rad(0.6705))
#' euler_geo2cart(-63.036, 107.359, deg2rad(0.6705))
NULL

#' @rdname cartesian_comp
#' @export
euler_cart2geo <- function(wx, wy, wz) {
  # if (is.null(dim(w)) & length(w) == 3) w <- t(w)
  #
  # wx <- w[, 1]
  # wy <- w[, 2]
  # wz <- w[, 3]

  rate_rad <- rotation_rate(wx, wy, wz)

  # Calculate the direction cosines of the rotation axis
  l <- wx / rate_rad # x-component
  m <- wy / rate_rad # y-component
  n <- wz / rate_rad # z-component

  # Calculate latitude (φ) of the pole
  # φ = arcsin(n) where n is the z-component of the unit vector
  latitude <- asin(n)

  # Calculate longitude (λ) of the pole
  # λ = arctan(m/l) with appropriate quadrant correction
  longitude <- atan2(m, l)

  # Convert to degrees
  latitude_deg <- rad2deg(latitude)
  longitude_deg <- rad2deg(longitude) |> tectonicr::longitude_modulo()
  mag <- rate_rad

  euler <- cbind(latitude_deg, longitude_deg, mag)
  colnames(euler) <- c("lat", "lon", "mag")
  return(euler)
}

#' @rdname cartesian_comp
#' @export
euler_geo2cart <- function(lat, lon, mag) {
  lat_rad <- deg2rad(lat)
  lon_rad <- deg2rad(lon)

  # Calculate direction cosines of the rotation axis
  # These are the components of a unit vector pointing to the pole
  l <- cos(lat_rad) * cos(lon_rad) # x-component
  m <- cos(lat_rad) * sin(lon_rad) # y-component
  n <- sin(lat_rad) # z-component

  wx <- mag * l
  wy <- mag * m
  wz <- mag * n

  return(cbind(wx, wx, wy = wy, wz = wz))
}
