#' Relative rotation of between two rotations
#'
#' Calculates the relative rotation between two (absolute) rotations, i.e. the
#' difference from rotation 1 to rotation 2.
#'
#' @param r1,r2 Objects of class \code{"euler"}, i.e. four-column vectors  giving the Cartesian coordinates of the
#' Euler vector and the amount of rotation in radians for first rotation
#' (\code{r1}) and subsequent second rotation (\code{r2})
#' @importFrom reticulate r_to_py py_to_r source_python
#' @importFrom dplyr %>%
#' @importFrom tectonicr cartesian_to_geographical
#' @importFrom pracma cross
#' @details
#' Giving two "absolute"rotations \eqn{R_i = R(\omega_i, \mathbf{e_1}),\; i = 1, 2} and their unit quaternions
#' \eqn{q_1 = q(R_1),\; i = 1, 2}.
#'
#' Their relative rotation is the difference from the absolute rotation
#' \eqn{R_1} to the absolute rotation \eqn{R_2}:
#'
#' \deqn{R_{21} = R_2 - R_1}
#'
#'
#' Because of \eqn{R_1^{-1} = -R(e_1, \omega_1) = R(e_1, -\omega_1)}, the rotation of 2 relative
#' to 1 (1 is "fixed") can be expressed by
#'
#'  \deqn{R(e_{21}, \omega_{21}) = R(e_1, -\omega_1) + R(e_2, \omega_2)}
#'
#' or
#'
#' \deqn{R(e_{21}, \omega_{21}) = R_2R_1^{-1}}
#'
#' In terms of quaternions the angle and axis of the rotation associated with \eqn{q_2 q_1^{-1} = q_2 q_1^*}
#' with \eqn{\textrm{Vec}(q_1^*) = -\textrm{Vec}(q_1)} is expressed as the concatenation of both quaternions
#' (\code{relative_euler_schaeben()}):
#' \deqn{
#'    \omega(R_2R_1^{-1}) = 2 \arccos (q_{20}q_{10} + \mathbf{q_2} \cdot \mathbf{q_1})
#'    }
#'  \deqn{
#'    \mathbf{e}(R_2R_1^{-1}) = \frac{1}{\sin \frac{\omega}{2}} \left(-q_{20}\mathbf{q_1} + q_{10}q_2 - \mathbf{q_2} \times \mathbf{q_1} \right)
#'  }
#'  Using angles and axes instead (\code{relative_euler_schaeben2()}), the same rotation is given by:
#'  \deqn{
#'  \omega(R_2R^{-1}) = 2 \arccos \left( \cos \frac{\omega_2}{2} \cos \frac{\omega_1}{2} + \sin \frac{\omega_2}{2} \mathbf{e_2} \cdot \sin \frac{\omega_1}{2} \mathbf{e_1}   \right)
#'  }
#'  \deqn{
#'  \mathbf{e}(R_2R^{-1}) = \frac{1}{\sin \frac{\omega}{2}} \left(  - \cos \frac{\omega_2}{2} \sin \frac{\omega_1}{2} \mathbf{e_1} + \cos \frac{\omega_1}{2} \sin \frac{\omega_2}{2} \mathbf{e_2} - \sin \frac{\omega_2}{2} \mathbf{e_2} \times \sin\frac{\omega_1}{2} \mathbf{e_1}    \right)
#'  }
#'
#'  "As-if-infinitesimal" angle and axis (\code{quasi_infinitesimal_euler()}) for \eqn{\omega_i << 1, \; i = 1, 2} are:
#'  \deqn{
#'  \omega(R_2R^{-1}) \approx  \omega_2 - \omega_1
#'  }
#'  \deqn{
#'  \mathbf{e}(R_2R^{-1}) \approx \frac{\omega_2 \mathbf{e_2} - \omega_1 \mathbf{e_1}}{ \lVert \omega_2 \mathbf{e_2} - \omega_1 \mathbf{e_1} \rVert }
#'  }
#' @references Cox, A. and Hart, R. B. (1986). *Finite rotations. In Plate Tectonics. How it works*.
#' Wiley-Blackwell.
#'
#' Greiner, B. (1999). Euler rotations in plate-tectonic reconstructions.
#' *Computers and Geosciences*, 25(3), 209--216.
#' \doi{10.1016/S0098-3004(98)00160-5}
#'
#' Le Pichon, X., Francheteau, J. and Bonnin, J. (1973). Kinematics of relative
#' movements. In *Plate Tectonics* (1st ed., pp. 19--125). Elsevier.
#' \doi{10.1016/B978-0-444-41094-8.50010-8}
#'
#' Schaeben, H., Kroner, U. and Stephan, T. (2021). Euler Poles of Tectonic
#' Plates. In B. S. Daza Sagar, Q. Cheng, J. McKinley and F. Agterberg (Eds.),
#' *Encyclopedia of Mathematical Geosciences. Encyclopedia of Earth Sciences Series*
#' (pp. 1--7). Springer Nature Switzerland AG 2021.
#' \doi{10.1007/978-3-030-26050-7_435-1}
#' @returns \code{list}. Infinitesimal and the finite approach Euler axes
#' (geographical coordinates) and Euler angles (in degrees)
#' @name rotation
#' @aliases rotation quaternion
#' @seealso [to_euler()] for class \code{"euler"}
#' @examples
#' x <- c(33, 15, -10) %>% to_euler()
#' y <- c(65, -100, 2) %>% to_euler()
#' # expected results: c(42.35960, 10.60413, 10.85098)
#' quasi_infinitesimal_euler(x, y) # Schaeben
#' relative_euler_greiner(x, y) # Greiner in terms of rotation matrix
#' relative_euler_lepichon(x, y) # LePichon using quaternions
#' relative_euler_schaeben2(x, y) # Schaeben in terms of angle and axis
#' relative_euler_py_schaeben(x, y) # Schaeben using quaternions
#' relative_euler_schaeben(x, y) # Schaeben using quaternions
NULL

#' @rdname rotation
#' @export
relative_euler_schaeben <- function(r1, r2) {
  # names(r1) <- names(r2) <- NULL
  q <- as_quaternion(r1)
  p <- as_quaternion(r2)
  q.vec <- quaternion_as_vector_part(q)
  p.vec <- quaternion_as_vector_part(p)


  angle <- 2 * acos(
    (p[1] * q[1]) + (p.vec %*% q.vec)
  ) %>% as.vector()

  a <- -p[1] * q.vec + q[1] * p.vec - pracma::cross(p.vec, q.vec)

  axis <- 1 / sin(angle / 2) * (a)
  names(axis) <- NULL

  list(
    axis = axis %>% tectonicr::cartesian_to_geographical(),
    angle = angle / (pi / 180)
  )
}

#' @rdname rotation
#' @export
relative_euler_schaeben2 <- function(r1, r2) {
  names(r1) <- names(r2) <- NULL

  w1 <- r1[4]
  w2 <- r2[4]

  e1 <- c(r1[1], r1[2], r1[3])
  e2 <- c(r2[1], r2[2], r2[3])

  angle <- 2 * acos(
    cos(w2 / 2) * cos(w1 / 2) + (sin(w2 / 2) * e2) %*% (sin(w1 / 2) * e1)
  ) %>% as.numeric()

  a <- 1 / sin(angle / 2)
  b <- -cos(w2 / 2) * sin(w1 / 2) * e1 + cos(w1 / 2) * sin(w2 / 2) * e2 - pracma::cross(sin(w2 / 2) * e2, sin(w1 / 2) * e1)

  axis <- a * b

  list(
    axis = axis %>% tectonicr::cartesian_to_geographical(),
    angle = angle / (pi / 180)
  )
}

#' @rdname rotation
#' @export
relative_euler_py_schaeben <- function(r1, r2) {
  # reticulate::py_run_file(system.file("python", "quaternions.py", package = "euler"), convert = FALSE)
  reticulate::source_python(system.file("python", "quaternions.py", package = "euler"), convert = FALSE)
  # reticulate::py_run_file(system.file("python", "quaternions.py", package = "euler"), convert = FALSE)
  # reticulate::source_python(system.file("python", "quaternions.py", package = "euler"), convert = FALSE)
  # reticulate::source_python("inst/python/quaternions.py", convert = FALSE)

  R1 <- reticulate::r_to_py(r1) %>% euler2quat()
  R2 <- reticulate::r_to_py(r2) %>% euler2quat()

  angle <- euler_angle(R1, R2) %>%
    reticulate::py_to_r()
  axis <- euler_axis(R1, R2, angle) %>%
    reticulate::py_to_r() %>%
    tectonicr::cartesian_to_geographical()

  list(
    axis = axis,
    angle = angle / (pi / 180)
  )
}

#' @rdname rotation
#' @export
relative_euler_greiner <- function(r1, r2) {
  r1[4] <- -1*r1[4] # reverse rot1

  r1.rot <- euler_matrix(r1)
  r2.rot <- euler_matrix(r2)

  e <- (r2.rot %*% r1.rot) %>%
    matrix_2_angles()

  list(
    axis = e$axis,
    angle = e$angle
  )
}

#' @rdname rotation
#' @export
relative_euler_lepichon <- function(r1, r2) {
  r1 <- from_euler(r1)
  r2 <- from_euler(r2)
  qt <- as_quaternion2(r2) * as_quaternion2(r1)
  quat_2_angles(qt)
}

#' @rdname rotation
#' @export
quasi_infinitesimal_euler <- function(r1, r2) {
  axis <- (r2[4] * c(r2[1], r2[2], r2[3]) - r1[4] * c(r1[1], r1[2], r1[3])) %>%
    normalize_vector() %>%
    tectonicr::cartesian_to_geographical()

  angle <- abs((r2[4] - r1[4]) / (pi / 180))

  # if(angle < 0){
  #   angle <- abs(angle)
  #   axis[1] <- -1 * axis[1]
  #   axis[2] <- axis[2]+180
  # }

  list(
    axis = as.numeric(axis),
    angle = as.numeric(angle)
  )
}

check_rotation <- function(w1, w2, w3) {
  abs(w3) <= abs(w1) + abs(w2)
}
