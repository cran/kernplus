get.dpill <- function(cov, y) {
  bw <- KernSmooth::dpill(cov, y)
  if (is.nan(bw)) {
    par <- 0.06
    while (is.nan(bw)) {
      bw <- KernSmooth::dpill(cov, y, proptrun = par)
      par <- par + 0.01
    }
  }
  return(bw)
}

bw.adp <- function(y, x, id.dir = NA, id.adp = id.dir) {
  x.rnd <- rand.add(x)

  if (is.na(id.adp)) {
    bw.fix <- sapply(1:ncol(x.rnd), function(p) get.dpill(x.rnd[, p], y))
    list(bw.fix = bw.fix, bw.adp = NA, id.adp = NA, cutpt = NA)
  } else {
    cutpt <- lapply(id.adp, function(p) cut.pts(cov = x[, p], circ = (p %in%
      id.dir), .n.reg = 2))
    if (all(sapply(cutpt, function(pts) all(is.na(pts))))) {
      bw.fix <- sapply(c(1:ncol(x.rnd)), function(p) get.dpill(x.rnd[, p],
        y))
      list(bw.fix = bw.fix, bw.adp = NA, id.adp = NA, cutpt = NA)
    } else {
      bins <- sapply(1:length(id.adp), function(p) .bincode(x[, id.adp[p]],
        breaks = cutpt[[p]], include.lowest = TRUE))
      bw.fix <- sapply(c(1:ncol(x.rnd))[-id.adp], function(p) get.dpill(x.rnd[,
        p], y))
      bw.adp <- lapply(1:length(id.adp), function(p) {
        n.bin <- length(unique(bins[, p]))
        if (id.adp[p] %in% id.dir & cutpt[[p]][1] <= 0 & cutpt[[p]][length(cutpt[[p]])] >=
          360) {
          id.adj <- which(bins[, p] == n.bin)
          bins[id.adj, p] <- 1
          n.bin <- n.bin - 1
        }
        bw <- sapply(1:n.bin, function(b) {
          id.bin <- which(bins[, p] == b)
          cov.sub <- x.rnd[id.bin, id.adp[p]]
          if (id.adp[p] %in% id.dir & cutpt[[p]][1] <= 0 & cutpt[[p]][length(cutpt[[p]])] >=
          360)
          cov.sub[which(cov.sub < cutpt[[p]][2])] <- cov.sub[which(cov.sub <
            cutpt[[p]][2])] + 360
          KernSmooth::dpill(cov.sub, y[id.bin])
        })
      })
      list(bw.fix = bw.fix, bw.adp = bw.adp, id.adp = id.adp, cutpt = cutpt)
    }
  }
}
