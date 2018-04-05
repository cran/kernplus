get.diff <- function(X.tr, x.TS) {
  X.tr <- as.matrix(X.tr)
  n.TR <- dim(X.tr)[1]
  q.TR <- dim(X.tr)[2]
  x.TS <- matrix(x.TS, 1, q.TR)
  oneV <- matrix(1, n.TR, 1)
  diff <- X.tr - (oneV %*% x.TS)
  return(diff)
}

AMK <- function(y.tr, X.tr, X.ts, bw, id.spd, id.dir) {
  X.tr <- as.matrix(X.tr)
  X.ts <- as.matrix(X.ts)
  y.tr <- as.matrix(y.tr)
  p <- ncol(X.tr)
  n.tr <- nrow(X.tr)
  n.ts <- nrow(X.ts)

  id.adp <- bw$id.adp
  id.fix <- 1:p
  if (!is.na(id.adp))
    id.fix <- id.fix[-which(id.fix == id.adp)]
  h.adp <- bw$bw.adp
  h.fix <- bw$bw.fix
  cutpt <- bw$cutpt

  bins <- c()
  if (!is.na(id.adp)) {
    bins <- as.data.frame(matrix(NA, nrow = n.ts, ncol = length(id.adp)))
    for (k in 1:length(id.adp)) {
      bins[, k] <- .bincode(X.ts[, id.adp[k]], cutpt[[k]], include.lowest = T)
      if (id.adp[k] == id.dir & (cutpt[[k]][length(cutpt[[k]])] - cutpt[[k]][1] ==
        360))
        bins[which(bins[, k] == (length(cutpt[[k]]) - 1)), k] <- 1
    }
  }

  cat("Estimating (%)\n")
  cat("0")
  est <- rep(NA, n.ts)
  cnt <- 0
  for (i in 1:n.ts) {
    cnt <- cnt + 1
    if (cnt/n.ts >= 0.1) {
      cat(".")
      cnt <- 0
    }

    diff <- get.diff(X.tr, X.ts[i, ])

    h <- rep(NA, p)
    h[id.fix] <- h.fix
    if (!is.na(id.adp)) {
      for (q in 1:length(id.adp)) h[id.adp[q]] <- h.adp[[q]][bins[i, q]]
    }

    if (p == 2) {
      # Bivariate kernel regression (Jeon&Taylor)
      dir <- circular::circular(diff[, id.dir], units = "degrees")

      kappa.j <- matrix(NA, n.tr, 2)
      kappa.j[, 1] <- stats::dnorm(diff[, id.spd]/h[id.spd])/h[id.spd]
      kappa.j[, 2] <- circular::dvonmises(dir, circular::circular(0), 1/(h[id.dir] *
        pi/180)^2)
      kappa <- kappa.j[, 1] * kappa.j[, 2]
      yhat <- sum(y.tr * kappa/sum(kappa))
    } else if (p > 2) {
      # Additive multivariate kernel method
      id.gau <- 1:p
      id.gau <- id.gau[-c(id.spd, id.dir)]
      dir <- circular::circular(diff[, id.dir], units = "degrees")
      kappa <- rep(0, n.tr)
      yhat <- rep(NA, (p - 2))

      kappa.v <- stats::dnorm(diff[, id.spd]/h[id.spd])/h[id.spd]
      kappa.d <- circular::dvonmises(dir, circular::circular(0), 1/((h[id.dir] *
        pi/180)^2))
      for (j in 3:p) {
        kappa.j <- stats::dnorm(diff[, id.gau[j - 2]]/h[id.gau[j - 2]])/h[id.gau[j -
          2]]
        kappa <- (kappa.v * kappa.d * kappa.j)
        yhat[j - 2] <- sum(y.tr * kappa/sum(kappa))
      }
    }
    est[i] <- mean(yhat)
  }
  cat("100\n")

  return(est)
}
