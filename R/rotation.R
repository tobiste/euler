#' Rotations
#'
#' @param r1,r2 four-column vectors  giving the Cartesian coordinates of the
#' Euler vector and the amount of rotation in radians for first rotation
#' (\code{r1}) and subsequent second rotation (\code{r2})
#' @importFrom reticulate r_to_py py_to_r source_python
#' @importFrom dplyr %>%
#' @importFrom tectonicr cartesian_to_geographical euler_pole euler_from_rot
#' @importFrom pracma cross
#' @details
#' \eqn{R_i = R(\omega_i, \mathbf{e_1}),\; i = 1, 2} and their unit quaternions
#' \eqn{q_1 = q(R_1),\; i = 1, 2}.
#'
#' Angle and axis of the rotation associated with \eqn{q_2 q_1^{-1} = q_2 q_1^*}
#' with \eqn{\textrm{Vec}(q_1^*) = -\textrm{Vec}(q_1)}
#' (\code{infinitesimal_quaternion()}):
#' \deqn{
#'    \omega(R_2R^{-1}) = 2 \arccos (q_{20}q_{10} + \mathbf{q_2} \cdot \mathbf{q_1})
#'    }
#'  \deqn{
#'    \mathbf{e}(R_2R^{-1}) = \frac{1}{\sin \frac{\omega}{2}} \left(-q_{20}\mathbf{q_1} + q_{10}q_2 - \mathbf{q_2} \times \mathbf{q_1} \right)
#'  }
#'  In terms of angles and axes (\code{infinitesimal_euler()}):
#'  \deqn{
#'  \omega(R_2R^{-1}) = 2 \arccos \left( \cos \frac{\omega_2}{2} \cos \frac{\omega_1}{2} + \sin \frac{\omega_2}{2} \mathbf{e_2} \cdot \sin \frac{\omega_1}{2} \mathbf{e_1}   \right)
#'  }
#'  \deqn{
#'  \mathbf{e}(R_2R^{-1}) = \frac{1}{\sin \frac{\omega}{2}} \left(  - \cos \frac{\omega_2}{2} \sin \frac{\omega_1}{2} \mathbf{e_1} + \cos \frac{\omega_1}{2} \sin \frac{\omega_2}{2} \mathbf{e_2} - \sin \frac{\omega_2}{2} \mathbf{e_2} \times \sin\frac{\omega_1}{2} \mathbf{e_1}    \right)
#'  }
#'  "As-if-infinitesimal" angle and axis (\code{as_if_infinitesimal_euler()}) for \eqn{\omega_i << 1, \; i = 1, 2}:
#'  \deqn{
#'  \omega(R_2R^{-1}) \approx  \omega_2 - \omega_1
#'  }
#'  \deqn{
#'  \mathbf{e}(R_2R^{-1}) \approx \frac{\omega_2 \mathbf{e_2} - \omega_1 \mathbf{e_1}}{ \lVert \omega_2 \mathbf{e_2} - \omega_1 \mathbf{e_1} \rVert }
#'  }
#' @name rotation
#' @examples
#' x <- c(90, 0, 0.7) %>% to_euler()
#' y <- c(45, 30, 0.15) %>% to_euler()
#' as_if_infinitesimal_euler(x, y)
#' finite_euler(x, y)
#' infinitesimal_euler2(x, y)
#' infinitesimal_euler(x, y)
NULL

#' @rdname rotation
#' @export
as_if_infinitesimal_euler <- function(r1, r2) {
  axis <- (r2[4] * c(r2[1], r2[2], r2[3]) - r1[4] * c(r1[1], r1[2], r1[3])) %>%
    normalize_vector() %>%
    tectonicr::cartesian_to_geographical()

  angle <- (r2[4] - r1[4]) / (pi / 180)

  list(
    axis.fin = as.numeric(axis),
    angle.fin = as.numeric(angle)
  )
}

#' @rdname rotation
#' @export
finite_euler <- function(r1, r2) {
  r1.pole <- tectonicr::euler_pole(r1[1], r1[2], r1[3], geo = FALSE)
  r2.pole <- tectonicr::euler_pole(r2[1], r2[2], r2[3], geo = FALSE)

  r1.rot <- tectonicr::euler_rot(r1.pole, r1[4] / (pi / 180))
  r2.rot <- tectonicr::euler_rot(r2.pole, r2[4] / (pi / 180))

  e <- (r1.rot %*% r2.rot) %>% tectonicr::euler_from_rot()

  list(
    axis.fin = c(e$pole$lat, e$pole$lon),
    angle.fin = e$psi
  )
}

#' @rdname rotation
#' @export
infinitesimal_euler <- function(r1, r2) {
  names(r1) <- names(r2) <- NULL

  w1 <- r1[4]
  w2 <- r2[4]

  e1 <- c(r1[1], r1[2], r1[3])
  e2 <- c(r2[1], r2[2], r2[3])

  angle <- 2 * acos(
    cos(w2 / 2) * cos(w1 / 2) + (sin(w2 / 2) * e2) %*% (sin(w1 / 2) * e1)
  ) %>% as.numeric()

  temp_1 <- 1 / sin(angle / 2)
  temp_2 <- -cos(w2 / 2) * sin(w1 / 2) * e1 + cos(w1 / 2) * sin(w2 / 2) * e2 - pracma::cross(sin(w2 / 2) * e2, sin(w1 / 2) * e1)

  axis <- temp_1 * temp_2

  list(
    axis.inf = axis %>% tectonicr::cartesian_to_geographical(),
    angle.inf = angle / (pi / 180)
  )
}

#' @rdname rotation
#' @export
infinitesimal_quaternion <- function(r1, r2) {
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
    axis.inf = axis,
    angle.inf = angle / (pi / 180)
  )
}
