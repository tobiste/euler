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

#' Le Pichon method
#' @import onion
#' @importFrom tectonicr geographical_to_cartesian cartesian_to_geographical
#' @importFrom dplyr %>%
#' @name lepichon
NULL

#' @rdname lepichon
as_quaternion2 <- function(x) {
  stopifnot(is.numeric(x))

  x <- x * (pi / 180) # to radians
  lat_c <- (pi / 2) - x[1] # colatitude

  omega <- cos(x[3] / 2)
  chi <- sin(x[3] / 2) * sin(lat_c) * cos(x[2])
  eta <- sin(x[3] / 2) * sin(lat_c) * sin(x[2])
  zeta <- sin(x[3] / 2) * cos(lat_c)

  onion::quaternion(Re = omega, i = chi, j = eta, k = zeta)
}

#' @rdname lepichon
quat_composition <- function(q1, q2) {
  stopifnot(onion::is.quaternion(q1) & onion::is.quaternion(q2))
  omega_t <- Re(q1) * Re(q2) - onion::i(q1) * onion::i(q2) - onion::j(q1) * onion::j(q2) - onion::k(q1) * onion::k(q2)
  chi_t <- Re(q1) * onion::i(q2) + onion::i(q1) * Re(q2) - onion::j(q1) * onion::k(q2) + onion::k(q1) * onion::j(q2)
  eta_t <- Re(q1) * onion::j(q2) + onion::i(q1) * onion::k(q2) + onion::j(q1) * Re(q2) - onion::k(q1) * onion::i(q2)
  zeta_t <- Re(q1) * onion::k(q2) - onion::i(q1) * onion::j(q2) + onion::j(q1) * onion::i(q2) + onion::k(q1) * Re(q2)

  qt <- c(omega_t, chi_t, eta_t, zeta_t)
  names(qt) <- NULL
  if (omega_t < 0) {
    qt <- -1 * qt
  }
  onion::quaternion(Re = omega_t, i = chi_t, j = eta_t, k = zeta_t)
}


#' @rdname lepichon
quat_2_angles <- function(q) {
  stopifnot(onion::is.quaternion(q))
  theta <- 2 * acos(Re(q))
  names(theta) <- NULL

  lat <- (pi / 2) - acos(onion::i(q) / (sin(theta / 2)))
  lon <- atan(onion::j(q) / onion::k(q))

  axis <- (c(lat, lon) / (pi / 180)) # %>%
  # tectonicr::geographical_to_cartesian() %>%
  # tectonicr::cartesian_to_geographical()
  names(axis) <- NULL

  list(
    axis.lep = axis,
    angle.lep = (theta / (pi / 180))
  )
}

#' Quaternion class
#'
#' Build an unit real quaternion
#'
#' @param x Object of class \code{"euler"}
#' @seealso [to_euler()] for class \code{"euler"}
#' @export
#' @examples
#' euler1 <- c(90, 0, 10)
#' euler1 dplyr::%>%
#'   to_euler() dplyr::%>%
#'   as_quaternion()
as_quaternion <- function(x) {
  stopifnot("euler" %in% class(x))
  x <- as.numeric(x)

  scalar <- cos(x[4] / 2)
  vector <- c(x[1], x[2], x[3]) * sin(x[4] / 2)

  q <- c(Re = scalar, i = vector[1], j = vector[2], k = vector[3])
  class(q)[2] <- "quaternion"
  return(q)
}

quaternion_as_vector_part <- function(x) {
  stopifnot("quaternion" %in% class(x))
  qvec <- c(x[2], x[3], x[4])
  names(qvec) <- NULL
  return(qvec)
}


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
#' @importFrom dplyr %>%
#' @name euler
#' @examples
#' euler1 <- c(90, 0, 10)
#' to_euler(euler1)
#' to_euler(euler1) dplyr::%>% from_euler()
NULL

#' @rdname euler
#' @export
to_euler <- function(g) {
  stopifnot(is.numeric(g) & length(g) == 3)
  cart <- tectonicr::geographical_to_cartesian(c(g[1], g[2])) %>%
    normalize_vector()
  angle <- g[3] * pi / 180
  e <- c(x = cart[1], y = cart[2], z = cart[3], angle = angle)
  class(e)[2] <- "euler"
  return(e)
}

#' @rdname euler
#' @export
from_euler <- function(x) {
  stopifnot("euler" %in% class(x))
  geo <- tectonicr::cartesian_to_geographical(c(x[1], x[2], x[3]))
  angle <- x[4] / (pi / 180)
  c(lat = geo[1], lon = geo[2], angle = angle)
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
#' sf_to_vector(in.plate) dplyr::%>% vector_to_sf()
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
  pracma::cross(
    c(x.cart[[1]], x.cart[[2]], x.cart[[3]]), c(y.cart[[1]], y.cart[[2]], y.cart[[3]])
  ) %>% tectonicr::cartesian_to_geographical()
}

#' Common small circle of absolute and relative pole
#'
#' @inheritParams quasi_infinitesimal_euler
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

  tectonicr::PoR_to_geographical(x = sf::st_as_sf(sm_np), ep = r1.ep) %>%
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
  a <- a * (pi / 180)
  b <- b * (pi / 180)

  acos(
    sin(a[1]) * sin(b[1]) + cos(a[1]) * cos(b[1]) * cos(abs(a[2] - b[2]))
  ) / (pi / 180)
}
