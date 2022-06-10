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
#' to_euler(euler1) %>% from_euler()
NULL

#' @rdname euler
#' @export
to_euler <- function(g) {
  stopifnot(is.numeric(g) & length(g) == 3)
  cart <- tectonicr::geographical_to_cartesian(c(g[1], g[2])) %>%
    normalize_vector()
  angle <- g[3] * pi / 180
  e <- c(x = cart[1], y = cart[2], z = cart[3], angle = angle)
  class(e) <- append(class(e), "euler")
  return(e)
}

#' @rdname euler
#' @export
from_euler <- function(x) {
  stopifnot(inherits(x, "euler"))
  geo <- tectonicr::cartesian_to_geographical(c(x[1], x[2], x[3]))
  angle <- x[4] / (pi / 180)
  c(lat = geo[1], lon = geo[2], angle = angle)
}


#' Convert between sf objects and point vectors
#'
#' Converts a sf object into a two-column vector or vice versa
#' @param sf \code{sf} object
#' @param x numeric vector
#' @param multi logical. Is x MULTIPOLYGON or MULTILINESTRING.
#' @importFrom sf st_coordinates st_as_sf st_sfc st_polygon st_cast st_sf
#' @importFrom dplyr %>% group_by summarise
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
