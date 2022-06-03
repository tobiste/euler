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
#' @inheritParams as_if_infinitesimal_euler
#' @examples
#' r1 <- c(90, 0, 0.7) %>% to_euler()
#' r2 <- c(45, 30, 0.15) %>% to_euler()
#' common_smallcircle(r1, r2)
common_smallcircle <- function(r1, r2){
  r1.r2 <- infinitesimal_quaternion(r1, r2)

  angle <- tectonicr::angle_vectors(
    tectonicr::geographical_to_cartesian(r1.r2$axis.inf),
    r1[1:3]
  )

r1.ep <- tectonicr::euler_pole(r1[1], r1[2], r1[3], geo = FALSE)

  sm_np <- data.frame(
    lon = c(seq(-180, 180, 2), seq(-180, 180, 2)),
    lat = c(rep(90-angle, 181), rep(-90+angle, 181))
  ) %>%
    st_as_sf(coords = c("lon", "lat")) %>%
    summarise(do_union = FALSE) %>%
    st_cast("MULTILINESTRING") %>%
    sf::st_wrap_dateline(
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
      quiet = TRUE)

  tectonicr::PoR_to_geographical(x = sf::st_as_sf(sm_np), ep = r1.ep) %>%
    sf::st_wrap_dateline(
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
      quiet = TRUE)
}



