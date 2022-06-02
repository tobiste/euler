# global reference to scipy (will be initialized in .onLoad)
quaternion <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  quaternion <<- reticulate::import("quaternion", delay_load = TRUE)
}
reticulate::source_python("src/quaternions.py", convert = FALSE)

# Global
deg2rad <- pi/180



normalize_vector <- function(v) v / sqrt(sum(v^2))



#' Euler class
#'
#' Converts Euler pole in from geographic to Cartesian coordinates and Euler
#' angle from degrees to radians
#'
#' @param x Vector for Euler pole position, three-column vector of the
#' geographic coordinates latitude and longitude, and the amount of rotation in
#' degrees
#' @importFrom tectonicr geographical_to_cartesian
#' @importFrom dplyr %>%
#' @export
#' @examples
#' euler1 <- to_euler(c(90, 0, 10))
#' euler2 <- to_euler(c(45, 30, 20))
to_euler <- function(x) {
  stopifnot(is.numeric(x) & length(x) == 3)
  cart <- tectonicr::geographical_to_cartesian(c(x[1], x[2]))
  angle <- x[3] * deg2rad
  data.frame(x = cart[1], y = cart[2], z = cart[3], angle)
}


#' As-if-infinitesimal rotation
#'
#' @inheritParams relative_quat
#' @importFrom dplyr %>%
#' @importFrom tectonicr cartesian_to_geographical
relative_euler <- function(p, q) {
  e <- (p$angle * c(p$x, p$y, p$z) - q$angle * c(q$x, q$y, q$z)) %>%
    normalize_vector() %>%
    tectonicr::cartesian_to_geographical()

  angle <- (p$angle - q$angle) / deg2rad

  list(
    axis.fin = e,
    angle.fin = angle
  )
}

relative_euler2 <- function(p, q) {
  p.pole <- tectonicr::euler_pole(p$x, p$y, p$z, geo = FALSE)
  q.pole <- tectonicr::euler_pole(q$x, q$y, q$z, geo = FALSE)

  p.rot <- tectonicr::euler_rot(p.pole, p$angle / deg2rad)
  q.rot <- tectonicr::euler_rot(q.pole, q$angle / deg2rad)

  e <- p.rot %*% q.rot %>% tectonicr::euler_from_rot()

  list(
    axis.fin = c(e$pole$lat, e$pole$lon),
    angle.fin = e$psi
  )
}



#' Infinitesimal rotation
#'
#' Infinitesimal rotation using quaternions
#'
#' @param p,q three-column vectors  giving the geographic coordinates latitude
#' and longitude, and the amount of rotation in degrees for first rotation (\code{p})
#' and subsequent second rotation (\code{q})
#' @importFrom reticulate r_to_py py_to_r
#' @importFrom dplyr %>%
#' @importFrom tectonicr cartesian_to_geographical
relative_quat <- function(p, q) {
  Q <- reticulate::r_to_py(q) %>% euler2quat()
  P <- reticulate::r_to_py(p) %>% euler2quat()

  angle <- euler_angle(P, Q) %>%
              reticulate::py_to_r()
  axis <- euler_axis(P, Q, angle) %>%
    reticulate::py_to_r() %>%
    tectonicr::cartesian_to_geographical()

  list(
    axis.inf = axis,
    angle.inf = angle / deg2rad
  )
}


#' Relative Rotation associated to two absolute rotations
#'
#' Calculates Euler rotation axis and angle for two given absolute rotations using the
#' infinitesimal (Schaeben et al. 2021) and the finite approach
#'
#' @param x,y three-column vectors  giving the geographic coordinates latitude
#' and longitude, and the amount of rotation in degrees for first rotation (\code{x})
#' and subsequent second rotation (\code{y})
#' @param infinitesimal,finite logical. Should the rotation be calculated using
#' the infinitesimal and/or the finite rotation approach?
#' @note \code{y} will be considered as the fixed plate for the relative plate
#' motion between \code{x} and \code{y}
#' @references Schaeben, H., Kroner, U., Stephan, T. (2021). Euler Poles
#' of Tectonic Plates. In B. S. Daza Sagar, Q. Cheng, J. McKinley,; F. Agterberg
#' (Eds.), Encyclopedia of Mathematical Geosciences. Encyclopedia of Earth
#' Sciences Series (pp. 1--7). Springer Nature Switzerland AG 2021.
#' \doi{10.1007/978-3-030-26050-7_435-1}
#' @returns \code{list}. Infinitesimal and the finite approach Euler axes
#' (geographical coordinates) and Euler angles (in degrees)
#' @export
#' @examples
#' x <- c(90, 0, 0.7)
#' y <- c(45, 30, 0.15)
#' relative_rotation(x, y)
relative_rotation <- function(x, y, infinitesimal = TRUE, finite = TRUE) {
  xe <- to_euler(c(x[1], x[2], x[3]))
  ye <- to_euler(c(y[1], y[2], y[3]))

  # "as-if-infinitesimal" Rotation axis:
  if (finite) {
    res.fin <- relative_euler(xe, ye)
  }

  # transform to py
  if (infinitesimal) {
    res.inf <- relative_quat(xe, ye)
  }

  if (infinitesimal & finite) {
    res <- append(res.inf, res.fin)
  } else if (infinitesimal & !finite) {
    res <- res.inf
  } else if (!infinitesimal & finite) {
    res <- res.fin
  } else {
    res <- NULL
  }

  return(res)
}


#' Euler pole migration
#'
#' Migration of the relative Euler pole associated with two absolute plate motions
#'
#' @inheritParams relative_rotation
#' @param steps numeric vector of time increments. The default is a sequence
#' from 1 to 10 by an incremental step of 1 (e.g. Myr)
#' @export
#' @examples
#' in.eu <- c(27.12746847, 17.32482497, 0.402388191)
#' som.eu <- c(22.2078593, -92.40545103, 0.085835298)
#' pole_migration(som.eu, in.eu)
#' pole_migration(in.eu, som.eu)
pole_migration <- function(x, y, steps = c(1, seq(25, 300, 25)), infinitesimal = TRUE, finite = TRUE) {
  res <- data.frame(time = NULL, axis.inf.lat = NULL, axis.inf.lon = NULL, angle.inf = NULL, axis.fin.lat = NULL, axis.fin.lon = NULL, angle.fin = NULL)
  rate.x <- x[3]
  rate.y <- y[3]
  for (i in steps) {
    x[3] <- rate.x * i
    y[3] <- rate.y * i

    rel.i <- relative_rotation(x, y, infinitesimal, finite)

    if (infinitesimal & finite) {
      res <- rbind(
        res,
        data.frame(
          time = i,
          axis.inf.lat = rel.i$axis.inf[1],
          axis.inf.lon = rel.i$axis.inf[2],
          angle.inf = rel.i$angle.inf,
          axis.fin.lat = rel.i$axis.fin[1],
          axis.fin.lon = rel.i$axis.fin[2],
          angle.fin = rel.i$angle.fin
        )
      )
    } else if (infinitesimal & !finite) {
      res <- rbind(
        res,
        data.frame(
          time = i,
          axis.inf.lat = rel.i$axis.inf[1],
          axis.inf.lon = rel.i$axis.inf[2],
          angle.inf = rel.i$angle.inf
        )
      )
    } else if (!infinitesimal & finite) {
      res <- rbind(
        res,
        data.frame(
          time = i,
          axis.fin.lat = rel.i$axis.fin[1],
          axis.fin.lon = rel.i$axis.fin[2],
          angle.fin = rel.i$angle.fin
        )
      )
    } else {
      res <- NULL
    }
  }
  return(res)
}

#' Great circle of Euler poles
#'
#' Calculates the pole to the great circle of Euler poles
#'
#' @param x,y two-column vectors giving the geographic coordinates latitude
#' and longitude in degrees
common_greatcircle <- function(x, y){
  x.cart <- to_euler(x)
  y.cart <- to_euler(y)
  pracma::cross(
    c(x.cart[[1]], x.cart[[2]], x.cart[[3]]), c(y.cart[[1]], y.cart[[2]],y.cart[[3]])
  ) %>% tectonicr::cartesian_to_geographical()
}

