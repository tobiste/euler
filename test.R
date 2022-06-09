library(euler)
library(dplyr)

r1 <- c(33, 15, -10) %>% to_euler()
r2 <- c(65, -100, 2) %>% to_euler()
# expected results: c(42.35960, 10.60413, 10.85098)

microbenchmark::microbenchmark(
  relative_euler_greiner(r1, r2),
  relative_euler_schaeben2(r1, r2),
  relative_euler_py_schaeben(r1, r2),
  relative_euler_schaeben(r1, r2)
)


greiner <- relative_euler_greiner(r1, r2)
schaeben2 <- relative_euler_schaeben2(r1, r2)
python <- relative_euler_py_schaeben(r1, r2)
schaeben <- relative_euler_schaeben(r1, r2)

schaeben$angle == schaeben2$angle
schaeben$angle == python$angle
schaeben$angle == greiner$angle


sprintf("%.14f", schaeben$angle)
sprintf("%.14f", greiner$angle)

sprintf("%.3f", 1000 * schaeben$angle * 2 * pi * 63781e6 / 360) # mm
sprintf("%.3f", 1000 * greiner$angle * 2 * pi * 63781e6 / 360) # mm
# 0.01 mm  difference!
