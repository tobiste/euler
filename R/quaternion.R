#' Le Pichon method
#' @import onion
#' @importFrom tectonicr geographical_to_cartesian cartesian_to_geographical
#' @importFrom magrittr %>%
#' @name lepichon
NULL

#' @rdname lepichon
as_quaternion2 <- function(x) {
  stopifnot(is.numeric(x))

  x <- deg2rad(x) # to radians
  colat <- (pi / 2) - x[1] # colatitude (= 90 - lat)

  omega <- cos(x[3] / 2)
  chi <- sin(x[3] / 2) * sin(colat) * cos(x[2])
  eta <- sin(x[3] / 2) * sin(colat) * sin(x[2])
  zeta <- sin(x[3] / 2) * cos(colat)

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
  onion::quaternion(Re = qt[1], i = qt[2], j = qt[3], k = qt[4])
}


#' @rdname lepichon
quat_2_angles <- function(q) {
  stopifnot(onion::is.quaternion(q))
  theta <- acos(Re(q))
  names(theta) <- NULL

  lat <- (pi / 2) - acos(onion::k(q) / (sin(theta)))
  lon <- atan(onion::j(q) / onion::i(q))

  axis <- c(lat, lon) # %>%
  # tectonicr::geographical_to_cartesian() %>%
  # tectonicr::cartesian_to_geographical()
  names(axis) <- NULL

  list(
    axis.lep = rad2deg(axis),
    angle.lep = rad2deg(2 * theta)
  )
}

#' Quaternion class
#'
#' Build an unit real quaternion
#'
#' @param x Object of class \code{"euler"}
#' @seealso [to_euler()] for class \code{"euler"}
#' @importFrom magrittr %>%
#' @export
#' @examples
#' euler1 <- c(90, 0, 10)
#' euler1 %>%
#'   to_euler() %>%
#'   as_quaternion()
as_quaternion <- function(x) {
  stopifnot(inherits(x, "euler"))
  x <- as.numeric(x)

  scalar <- cos(x[4] / 2)
  vector <- c(x[1], x[2], x[3]) * sin(x[4] / 2)

  q <- c(Re = scalar, i = vector[1], j = vector[2], k = vector[3])
  class(q) <- append(class(q), "quaternion")

  return(q)
}

quaternion_as_vector_part <- function(x) {
  stopifnot(inherits(x, "quaternion"))
  qvec <- c(x[2], x[3], x[4])
  names(qvec) <- NULL
  return(qvec)
}

#' Concatenation of rotations
#'
#' Accumulation of the two rotations r1 followed by r2
#'
#' @param r1,r2 Objects of class \code{"euler"}, i.e. four-column vectors
#' giving the Cartesian coordinates of the
#' Euler vector and the amount of rotation in radians for first rotation
#' (\code{r1}) and subsequent second rotation (\code{r2})
#' @details \deqn{\text{euler_concatenation(r1, r2)} = R_2 R_1}
#' @returns \code{list}. Euler axes (geographical coordinates) and Euler
#' angles (in degrees)
#' @export
#' @examples
#' x <- c(27.1275, 17.3248, 0.4024) %>% to_euler()
#' y <- c(22.2079, -92.4055, 0.0858) %>% to_euler()
#' euler_concatenation(x, y)
euler_concatenation <- function(r1, r2) {
  names(r1) <- names(r2) <- NULL

  w1 <- r1[4]
  w2 <- r2[4]

  e1 <- c(r1[1], r1[2], r1[3])
  e2 <- c(r2[1], r2[2], r2[3])

  angle <- 2 * acos(
    cos(w2 / 2) * cos(w1 / 2) - sin(w2 / 2) * sin(w1 / 2) * e2 %*% e1
  ) %>% as.numeric()

  a <- 1 / sin(angle / 2)
  b <- cos(w1 / 2) * sin(w2 / 2) * e2 + cos(w2 / 2) * sin(w1 / 2) * e1 + sin(w2 / 2) * sin(w1 / 2) * tectonicr::vcross(e2, e1)

  axis <- a * b

  list(
    axis = axis %>% tectonicr::cartesian_to_geographical(),
    angle = rad2deg(angle)
  )
}
