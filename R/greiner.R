#' @title Rotation Matrix
#'
#' @description Calculates the rotation matrix using the rotation axis and the
#' angle of rotation
#' @param x Object of class \code{"euler"}
#' @return \code{matrix}
#' @seealso [rotation_axis()] and [rotation_angle()] to extract the axis and the
#' angle from the rotation matrix.
#' @export
#' @examples
#' w <- c(0, 90, 90) %>% to_euler()
#' rotation_matrix(w)
rotation_matrix <- function(x) {
  stopifnot("euler" %in% class(x))
  alpha <- x[4]

  n <- normalize_vector(c(x[1], x[2], x[3]))
  R <- matrix(nrow = 3, ncol = 3)
  R[1, 1] <- n[1]^2 * (1 - cos(alpha)) + cos(alpha)
  R[1, 2] <- n[1] * n[2] * (1 - cos(alpha)) - n[3] * sin(alpha)
  R[1, 3] <- n[1] * n[3] * (1 - cos(alpha)) + n[2] * sin(alpha)
  R[2, 1] <- n[2] * n[1] * (1 - cos(alpha)) + n[3] * sin(alpha)
  R[2, 2] <- n[2]^2 * (1 - cos(alpha)) + cos(alpha)
  R[2, 3] <- n[2] * n[3] * (1 - cos(alpha)) - n[1] * sin(alpha)
  R[3, 1] <- n[3] * n[1] * (1 - cos(alpha)) - n[2] * sin(alpha)
  R[3, 2] <- n[3] * n[2] * (1 - cos(alpha)) + n[1] * sin(alpha)
  R[3, 3] <- n[3]^2 * (1 - cos(alpha)) + cos(alpha)
  R
}


#' @title Rotation angle from rotation matrix
#' @description Extracts the rotation angle from rotation matrix
#' @param A 3x3 matrix
#' @note Infinitesimal small rotation (i.e. small angles) will cause an Error
#' due to round-off errors. In order to avoid the error, these rotations will
#' be treated as equal rotations. The function will print a warning message when
#' this is the case.
#' @return numeric angle in radians
#' @seealso [rotation_axis()] to extract the axis of the rotation matrix.
#' @export
#' @examples
#' w <- c(0, 90, 90) %>% to_euler()
#' rot <- rotation_matrix(w)
#'
#' rotation_angle(rot)
rotation_angle <- function(A) {
  stopifnot(is.matrix(A))

  a <- (sum(diag(A)) - 1) / 2
  if (a >= 1) {
    warning("introduces round-off")
    a <- 1
  }
  acos(a)
}

#' @title Rotation axis from rotation matrix
#' @description Extracts the rotation axis from rotation matrix
#' @param A 3x3 matrix
#' @param psi angle in radians
#' @return Euler vector in Cartesian coordinates
#' @seealso [rotation_angle()] to extract the angle of the rotation matrix.
#' @export
#' @examples
#' w <- c(0, 90, 90) %>% to_euler()
#' rot <- rotation_matrix(w)
#'
#' rotation_axis(rot, rotation_angle(rot))
rotation_axis <- function(A, psi) {
  stopifnot(is.matrix(A))

  e1 <- (A[3, 2] - A[2, 3]) / 2 * sin(psi)
  e2 <- (A[1, 3] - A[3, 1]) / 2 * sin(psi)
  e3 <- (A[2, 1] - A[1, 2]) / 2 * sin(psi)
  c(e1, e2, e3)
}

#' @title Euler rotation matrix
#' @description Creates a matrix from the given set of values.
#' @param x Object of class \code{"euler"}
#' @return \code{matrix}
#' @references Greiner, B. (1999). Euler rotations in plate-tectonic
#' reconstructions. *Computers and Geosciences*, 25(3), 209--216.
#' \doi{10.1016/S0098-3004(98)00160-5}
#' @export
#' @examples
#' ep <- c(0, 90, 90) %>% to_euler()
#' euler_matrix(ep)
euler_matrix <- function(x) {
  stopifnot("euler" %in% class(x))
  x <- from_euler(x) %>% deg2rad()

  mat <- matrix(nrow = 3, ncol = 3)
  mat[1, 1] <- sin(x[1]) * cos(x[2])
  mat[1, 2] <- sin(x[1]) * sin(x[2])
  mat[1, 3] <- -cos(x[1])
  mat[2, 1] <- -sin(x[2])
  mat[2, 2] <- cos(x[2])
  mat[2, 3] <- 0
  mat[3, 1] <- cos(x[1]) * cos(x[2])
  mat[3, 2] <- cos(x[1]) * sin(x[2])
  mat[3, 3] <- sin(x[1])

  R <- matrix(nrow = 3, ncol = 3)
  R[1, 1] <- cos(x[3])
  R[1, 2] <- -sin(x[3])
  R[1, 3] <- 0
  R[2, 1] <- sin(x[3])
  R[2, 2] <- cos(x[3])
  R[2, 3] <- 0
  R[3, 1] <- 0
  R[3, 2] <- 0
  R[3, 3] <- 1

  solve(mat) %*% R %*% mat
}

#' @title Euler axis and angle from Euler matrix
#' @description Extracts the coordinates of an Euler pole and the angle from a
#' Euler matrix
#' @param A 3x3 matrix
#' @details If there is no rotation (i.e. angle = 0), the coordinates of the
#' axis are equal to Earth's spin axis (according to the GPLATES convention).
#' @param as_euler logical. Whether the output should be an object of
#' \code{"euler.pole"} or a list.
#' @return axis and angle in degrees
#' @export
#' @importFrom tectonicr cartesian_to_geographical euler_pole
#' @examples
#' c(90, 0, 90) %>%
#'   to_euler() %>%
#'   euler_matrix() %>%
#'   matrix_2_angles()
matrix_2_angles <- function(A, as_euler = FALSE) {
  stopifnot(is.matrix(A))

  psi <- rotation_angle(A) # radians
  if (psi != 0) {
    ra <- rotation_axis(A, psi)
  } else {
    ra <- c()
    ra[1] <- 0
    ra[2] <- 0
    ra[3] <- 1
  }
  if (!as_euler) {
    list(
      axis = tectonicr::cartesian_to_geographical(ra),
      angle = rad2deg(psi)
    )
  } else {
    tectonicr::euler_pole(ra[1], ra[2], ra[3], geo = FALSE, angle = rad2deg(psi))
  }
}
