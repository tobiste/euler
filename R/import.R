#' Euler class
#'
#' Converts Euler pole in from geographic to Cartesian coordinates and Euler
#' angle from degrees to radians
#'
#' @param x Vector for Euler pole position, three-column vector of the
#' geographic coordinates latitude and longitude, and the amount of rotation in
#'  degrees
#'  @importFrom tectonicr geographical_to_cartesian deg2rad
#'  @importFrom ptrotR normalize_vector
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


#' Infinitesimal and finite Rotations
#'
#' get both infinitesimal and finite Euler rotation axis and angles
#'
#' @param x,y three-column vectors  giving the geographic coordinates latitude
#' and longitude, and the amount of rotation in degrees for rotation 1 \code{x}
#' and subsequent rotation 2
#'
#' @export
#' @examples
#' x <- c(90, 0, 0.7)
#' y <- c(45, 30, 0.15)
#' rotate_euler(x, y)
rotate_euler <- function(x, y){
  eulerx <- euler(c(x[1], x[2], x[3]))
  eulery <- euler(c(y[1], y[2], y[3]))

  #"as-if-infinitesimal" Rotation axis:
  t <- eulerx$angle * c(eulerx$x, eulerx$y, eulerx$z) - eulery$angle * c(eulery$x, eulery$y, eulery$z)
  proxy.cart <- t / abs(t)
  proxy.axis <- tectonicr::cartesian_to_geographical(proxy.cart)
  proxy.angle <- tectonicr::rad2deg(eulerx$angle - eulery$angle)

  reticulate::source_python("src/quaternions.py", convert = FALSE)

  # transform to py
  py.eulerx <-  reticulate::r_to_py(eulerx)
  py.eulery <-  reticulate::r_to_py(eulery)

    R1 <- R2quat(py.eulerx)
    R2 <- R2quat(py.eulery)

    py.angle <- euler_angle(R1, R2)
    py.axis <- euler_axis(R1, R2, angle)

    # transform back to R
    angle <- reticulate::py_to_r(py.angle)
    axis <- reticulate::py_to_r(py.axis)
    #quit # end python code

  return(list(
    axis =  tectonicr::cartesian_to_geographical(axis),
    angle = tectonicr::rad2deg(angle),
    proxy.axis = proxy.axis,
    proxy.angle = proxy.angle
    )
  )
}

#"as-if-infinitesimal" Rotation:
#  tectonicr::rotation_matrix(py$angle.inf, py$axis.inf)


