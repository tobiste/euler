#' Relative rotation of between two rotations
#'
#' Calculates the relative rotation between two (absolute) rotations, i.e. the
#' difference from rotation 1 to rotation 2.
#'
#' @param r1,r2 four-column vectors  giving the Cartesian coordinates of the
#' Euler vector and the amount of rotation in radians for first rotation
#' (\code{r1}) and subsequent second rotation (\code{r2})
#' @importFrom reticulate r_to_py py_to_r source_python
#' @importFrom dplyr %>%
#' @importFrom tectonicr cartesian_to_geographical euler_pole euler_from_rot
#' @importFrom pracma cross
#' @importFrom onion quaternion
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
#' \deqn{R(e_{21}, \omega_{21}) = R_2R_1^{-1}}.
#'
#' In terms of quaternions the angle and axis of the rotation associated with \eqn{q_2 q_1^{-1} = q_2 q_1^*}
#' with \eqn{\textrm{Vec}(q_1^*) = -\textrm{Vec}(q_1)}
#' (\code{infinitesimal_quaternion()}) is:
#' \deqn{
#'    \omega(R_2R_1^{-1}) = 2 \arccos (q_{20}q_{10} + \mathbf{q_2} \cdot \mathbf{q_1})
#'    }
#'  \deqn{
#'    \mathbf{e}(R_2R_1^{-1}) = \frac{1}{\sin \frac{\omega}{2}} \left(-q_{20}\mathbf{q_1} + q_{10}q_2 - \mathbf{q_2} \times \mathbf{q_1} \right)
#'  }
#'  Using angles and axes instead (\code{infinitesimal_euler()}), the same rotation is given by:
#'  \deqn{
#'  \omega(R_2R^{-1}) = 2 \arccos \left( \cos \frac{\omega_2}{2} \cos \frac{\omega_1}{2} + \sin \frac{\omega_2}{2} \mathbf{e_2} \cdot \sin \frac{\omega_1}{2} \mathbf{e_1}   \right)
#'  }
#'  \deqn{
#'  \mathbf{e}(R_2R^{-1}) = \frac{1}{\sin \frac{\omega}{2}} \left(  - \cos \frac{\omega_2}{2} \sin \frac{\omega_1}{2} \mathbf{e_1} + \cos \frac{\omega_1}{2} \sin \frac{\omega_2}{2} \mathbf{e_2} - \sin \frac{\omega_2}{2} \mathbf{e_2} \times \sin\frac{\omega_1}{2} \mathbf{e_1}    \right)
#'  }
#'
#'  "As-if-infinitesimal" angle and axis (\code{as_if_infinitesimal_euler()}) for \eqn{\omega_i << 1, \; i = 1, 2} are:
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
#' @examples
#' x <- c(90, 0, 0.7) %>% to_euler()
#' y <- c(45, 30, 0.15) %>% to_euler()
#' as_if_infinitesimal_euler(x, y) # Schaeben
#' finite_euler_greiner(x, y) # Greiner in terms of rotation matrix
#' finite_euler_lepichon(x, y) # LePichon using quaternions
#' infinitesimal_euler(x, y) # Schaeben in terms of angle and axis
#' infinitesimal_quaternion(x, y) # Schaeben using quaternions
NULL

#' @rdname rotation
#' @export
as_if_infinitesimal_euler <- function(r1, r2) {
  axis <- (r2[4] * c(r2[1], r2[2], r2[3]) - r1[4] * c(r1[1], r1[2], r1[3])) %>%
    normalize_vector() %>%
    tectonicr::cartesian_to_geographical()

  angle <- (r2[4] - r1[4]) / (pi / 180)

  # if(angle < 0){
  #   angle <- abs(angle)
  #   axis[1] <- -1 * axis[1]
  #   axis[2] <- axis[2]+180
  # }

  list(
    axis.fin = as.numeric(axis),
    angle.fin = as.numeric(angle)
  )
}

#' @rdname rotation
#' @export
finite_euler_greiner <- function(r1, r2) {
  r1.pole <- tectonicr::euler_pole(r1[1], r1[2], r1[3], geo = FALSE)
  r2.pole <- tectonicr::euler_pole(r2[1], r2[2], r2[3], geo = FALSE)

  r1.rot <- tectonicr::euler_rot(r1.pole, -r1[4] / (pi / 180)) # reverse rot1
  r2.rot <- tectonicr::euler_rot(r2.pole, r2[4] / (pi / 180))

  e <- (r2.rot %*% r1.rot) %>% tectonicr::euler_from_rot()

  # if(e$psi < 0){
  #   e$psi <- abs(e$psi)
  #   e$pole$lat <- -1 * e$pole$lat
  #   e$pole$lon <- 180+e$pole$lon
  #   }

  list(
    axis.fin2 = c(e$pole$lat, e$pole$lon),
    angle.fin2 = e$psi
  )
}

#' @rdname rotation
#' @export
finite_euler_lepichon <- function(r1, r2) {
  r1 <- from_euler(r1)
  r2 <- from_euler(r2)
  to_quaternion(r2) * to_quaternion(r1) %>%
   quat_2_angles()
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
    axis.infeul = axis %>% tectonicr::cartesian_to_geographical(),
    angle.infeul = angle / (pi / 180)
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






check_rotation <- function(w1, w2, w3) {
  abs(w3) <= abs(w1) + abs(w2)
}

