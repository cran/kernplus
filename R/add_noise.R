decimalplaces <- function(num) {
  if ((num%%1) != 0)
    nchar(strsplit(sub("0+$", "", as.character(num)), ".", fixed = TRUE)[[1]][[2]]) else return(0)
}

rand.add <- function(df) {
  n.obs <- nrow(df)
  n.sp <- 1000
  if (nrow(df) < n.sp)
    n.sp <- n.obs
  num.decimal <- sapply(1:ncol(df), function(x) {
    max(sapply(df[1:n.sp, x], decimalplaces))
  })
  id.add <- which(num.decimal < 3)
  if (length(id.add) > 0) {
    df.new <- df
    num.add <- sapply(id.add, function(id.col) {
      inc.seed()
      set.seed(local(.seed, envir = kp.env))
      rng <- 10^(-num.decimal[id.col])/2
      return(stats::runif(n.obs, -rng, rng))
    })
    local(.seed <- 0, envir = kp.env)
    df.new[, id.add] <- df.new[, id.add] + num.add
    return(df.new)
  } else return(df)
}
