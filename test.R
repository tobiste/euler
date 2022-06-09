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
