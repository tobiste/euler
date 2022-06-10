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