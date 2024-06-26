% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rotation.R
\name{rotation}
\alias{rotation}
\alias{quaternion}
\alias{relative_euler_schaeben}
\alias{relative_euler_schaeben2}
\alias{relative_euler_schaeben3}
\alias{relative_euler_py_schaeben}
\alias{relative_euler_py2_schaeben}
\alias{relative_euler_greiner}
\alias{relative_euler_lepichon}
\alias{approximation_euler}
\title{Relative rotation between two rotations}
\usage{
relative_euler_schaeben(r1, r2)

relative_euler_schaeben2(r1, r2)

relative_euler_schaeben3(r1, r2)

relative_euler_py_schaeben(r1, r2)

relative_euler_py2_schaeben(r1, r2)

relative_euler_greiner(r1, r2)

relative_euler_lepichon(r1, r2)

approximation_euler(r1, r2)
}
\arguments{
\item{r1, r2}{Objects of class \code{"euler"}, i.e. four-column vectors  giving the Cartesian coordinates of the
Euler vector and the amount of rotation in radians for first rotation
(\code{r1}) and subsequent second rotation (\code{r2})}
}
\value{
\code{list}. Infinitesimal and the finite approach Euler axes
(geographical coordinates) and Euler angles (in degrees)
}
\description{
Calculates the relative rotation between two (absolute) rotations, i.e. the
difference from rotation 1 to rotation 2.
}
\details{
Given two "absolute"rotations \eqn{R_i = R(\omega_i, \mathbf{e_1}),\; i = 1, 2} and their unit quaternions
\eqn{q_1 = q(R_1),\; i = 1, 2}.

Their relative rotation is the difference from the absolute rotation
\eqn{R_1} to the absolute rotation \eqn{R_2}:

\deqn{R_{21} = R_2 - R_1}

Because of \eqn{R_1^{-1} = -R(e_1, \omega_1) = R(e_1, -\omega_1)}, the rotation of 2 relative
to 1 (1 is "fixed") can be expressed by

\deqn{R(e_{21}, \omega_{21}) = R(e_1, -\omega_1) + R(e_2, \omega_2)}

or

\deqn{R(e_{21}, \omega_{21}) = R_2R_1^{-1}}

In terms of quaternions the angle and axis of the rotation associated with \eqn{q_2 q_1^{-1} = q_2 q_1^*}
with \eqn{\textrm{Vec}(q_1^*) = -\textrm{Vec}(q_1)} is expressed as the concatenation of both quaternions
(\code{relative_euler_schaeben()}):
\deqn{
   \omega(R_2R_1^{-1}) = 2 \arccos (q_{20}q_{10} + \mathbf{q_2} \cdot \mathbf{q_1})
   }
\deqn{
   \mathbf{e}(R_2R_1^{-1}) = \frac{1}{\sin \frac{\omega}{2}} \left(-q_{20}\mathbf{q_1} + q_{10}q_2 - \mathbf{q_2} \times \mathbf{q_1} \right)
 }
Using angles and axes instead (\code{relative_euler_schaeben2()}), the same rotation is given by:
\deqn{
 \omega(R_2R^{-1}) = 2 \arccos \left( \cos \frac{\omega_2}{2} \cos \frac{\omega_1}{2} + \sin \frac{\omega_2}{2} \mathbf{e_2} \cdot \sin \frac{\omega_1}{2} \mathbf{e_1}   \right)
 }
\deqn{
 \mathbf{e}(R_2R^{-1}) = \frac{1}{\sin \frac{\omega}{2}} \left(  - \cos \frac{\omega_2}{2} \sin \frac{\omega_1}{2} \mathbf{e_1} + \cos \frac{\omega_1}{2} \sin \frac{\omega_2}{2} \mathbf{e_2} - \sin \frac{\omega_2}{2} \mathbf{e_2} \times \sin\frac{\omega_1}{2} \mathbf{e_1}    \right)
 }

Using Finite approximation for very small rotation angles ("as-if-infinitesimal"), the angle and axis (\code{quasi_infinitesimal_euler()}) for \eqn{\omega_i << 1, \; i = 1, 2} are:
\deqn{
 \omega(R_2R^{-1}) \approx  \omega_2 - \omega_1
 }
\deqn{
 \mathbf{e}(R_2R^{-1}) \approx \frac{\omega_2 \mathbf{e_2} - \omega_1 \mathbf{e_1}}{ \lVert \omega_2 \mathbf{e_2} - \omega_1 \mathbf{e_1} \rVert }
 }
}
\examples{
x <- c(27.1275, 17.3248, 0.4024) \%>\% to_euler()
y <- c(22.2079, -92.4055, 0.0858) \%>\% to_euler()
# expected results: c(-21.1008, -151.6591, 0.4202)
approximation_euler(x, y) # Schaeben
relative_euler_greiner(x, y) # Greiner in terms of rotation matrix
relative_euler_lepichon(x, y) # LePichon using quaternions
relative_euler_schaeben2(x, y) # Schaeben in terms of angle and axis
relative_euler_py_schaeben(x, y) # Schaeben using quaternions
relative_euler_py2_schaeben(x, y) # Schaeben using quaternions
relative_euler_schaeben(x, y) # Schaeben using quaternions
}
\references{
Cox, A. and Hart, R. B. (1986). \emph{Finite rotations. In Plate Tectonics. How it works}.
Wiley-Blackwell.

Greiner, B. (1999). Euler rotations in plate-tectonic reconstructions.
\emph{Computers and Geosciences}, 25(3), 209--216.
\doi{10.1016/S0098-3004(98)00160-5}

Le Pichon, X., Francheteau, J. and Bonnin, J. (1973). Kinematics of relative
movements. In \emph{Plate Tectonics} (1st ed., pp. 19--125). Elsevier.
\doi{10.1016/B978-0-444-41094-8.50010-8}

Schaeben, H., Kroner, U. and Stephan, T. (2021). Euler Poles of Tectonic
Plates. In B. S. Daza Sagar, Q. Cheng, J. McKinley and F. Agterberg (Eds.),
\emph{Encyclopedia of Mathematical Geosciences. Encyclopedia of Earth Sciences Series}
(pp. 1--7). Springer Nature Switzerland AG 2021.
\doi{10.1007/978-3-030-26050-7_435-1}
}
\seealso{
\code{\link[=to_euler]{to_euler()}} for class \code{"euler"}
}
