# global reference to scipy (will be initialized in .onLoad)
quaternion <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  quaternion <<- reticulate::import("quaternion", delay_load = TRUE)

}
reticulate::source_python("src/quaternions.py", convert = FALSE)


#' Euler class
#'
#' Converts Euler pole in from geographic to Cartesian coordinates and Euler
#' angle from degrees to radians
#'
#' @param x Vector for Euler pole position, three-column vector of the
#' geographic coordinates latitude and longitude, and the amount of rotation in
#'  degrees
#' @importFrom tectonicr geographical_to_cartesian deg2rad
#' @export
#' @examples
#' euler1 <- euler(c(90, 0, 10))
#' euler2 <- euler(c(45, 30, 20))
euler <- function(x){
  stopifnot(is.numeric(x) & length(x)==3)
  cart <- tectonicr::geographical_to_cartesian(c(x[1], x[2]))
  angle <-  tectonicr::deg2rad(x[3])
  data.frame(x = cart[1], y = cart[2], z = cart[3], angle)
}


#' Relative Rotation associated to two absolute rotations
#'
#' Calculates Euler rotation axis and angle for two given absolute rotations using the
#' infinitesimal (Schaeben et al. 2021) and the finite approach
#'
#' @param x,y three-column vectors  giving the geographic coordinates latitude
#' and longitude, and the amount of rotation in degrees for first rotation (\code{x})
#' and subsequent second rotation (\code{y})
#' @importFrom reticulate source_python r_to_py py_to_r
#' @importFrom tectonicr cartesian_to_geographical rad2deg cartesian_to_geographical
#' @references Schaeben, H., Kroner, U., &#38; Stephan, T. (2021). Euler Poles
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
#' relative_euler(x, y)
relative_euler <- function(x, y){
  eulerx <- euler(c(x[1], x[2], x[3]))
  eulery <- euler(c(y[1], y[2], y[3]))

  #"as-if-infinitesimal" Rotation axis:
  t <- eulerx$angle * c(eulerx$x, eulerx$y, eulerx$z) - eulery$angle * c(eulery$x, eulery$y, eulery$z)
  axis.fin.cart <- t / abs(t)
  axis.fin <- tectonicr::cartesian_to_geographical(proxy.cart)
  angle.fin <- tectonicr::rad2deg(eulerx$angle - eulery$angle)

  # transform to py
  py.eulerx <-  reticulate::r_to_py(eulerx)
  py.eulery <-  reticulate::r_to_py(eulery)

    R1 <- euler2quat(py.eulerx)
    R2 <- euler2quat(py.eulery)

    py.angle <- euler_angle(R1, R2)
    py.axis <- euler_axis(R1, R2, py.angle)

    # transform back to R
    angle <- reticulate::py_to_r(py.angle)
    axis <- reticulate::py_to_r(py.axis)
    #quit # end python code

  return(list(
    axis.inf =  tectonicr::cartesian_to_geographical(axis),
    angle.inf = tectonicr::rad2deg(angle),
    axis.fin,
    angle.fin
    )
  )
}


#' Euler pole migration
#'
#' Migration of the relative Euler pole associated with two absolute plate motions
#'
#' @param x,y three-column vectors  giving the geographic coordinates latitude
#' and longitude, and the amount of rotation in degrees for first rotation
#' (\code{x}) and subsequent second rotation (\code{y})
#' @param steps numeric vector of time increments. The default is a sequence
#' from 1 to 10 by an incremental step of 1 (e.g. Myr)
#' @export
#' @examples
#' in.eu <- c(27.12746847, 17.32482497, 0.402388191)
#' som.eu <- c(22.2078593, -92.40545103, 0.085835298)
#' euler_migration(in.eu, som.eu)
euler_migration <- function(x, y, steps = c(1, seq(25, 300, 25))){
  res <- data.frame(time = NULL, axis.inf.lat = NULL, axis.inf.lon = NULL, angle.inf = NULL, axis.fin.lat = NULL, axis.fin.lon = NULL, angle.fin = NULL)
  rate.x <-  x[3]
  rate.y <- y[3]
  y[3]
  for(i in steps){
    x[3] <- rate.x * i
    y[3] <- rate.y * i

    rel.i <- relative_euler(x, y)
    res <- rbind(res,
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

  }
  return(res)
}
