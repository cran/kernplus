kp.env <- new.env(parent = environment())
local(.seed <- 0, envir = kp.env)
local(.n.reg <- 2, envir = kp.env)

inc.seed <- function() {
  local(.seed <- .seed + 1, envir = kp.env)
}
