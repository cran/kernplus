cut.pts <- function(cov, circ = TRUE, .n.reg) {
  hobj <- graphics::hist(cov, breaks = "FD", plot = FALSE)
  brks <- hobj$breaks
  cnts <- hobj$counts
  id.cand <- which(cnts == 0)  #Bins without any data point
  if (length(id.cand) > 0) {
    loc.cut <- c()

    # For a circular variable, take a left-part of the support and attach it to the
    # right of the support
    if (circ & brks[1] <= 0 & brks[length(brks)] >= 360) {
      if (length(id.cand) == 1)
        loc.cut <- (brks[id.cand] + brks[id.cand + 1])/2
      else {
        id.diff <- diff(id.cand)
        n.cons <- rle(id.diff)
        vals <- n.cons$values
        lens <- n.cons$lengths
        if (any(vals == 1)) { # If there is any consecutive missing bins
          id.pos <- which(lens == max(lens[which(vals == 1)]))
          if (id.pos == 1)
            id.rng <- c(id.cand[1], id.cand[sum(lens[1:(id.pos)]) + 1])
          else
            id.rng <- id.cand[c(sum(lens[1:(id.pos - 1)]), sum(lens[1:(id.pos)])) + 1]
        } else id.rng <- rep(id.cand[1], 2)
        loc.cut <- (brks[id.rng[1]] + brks[(id.rng[2] + 1)])/2
      }
      cov[which(cov < loc.cut)] <- cov[which(cov < loc.cut)] + 360
    }

    # Clustering based on Gaussian Mixture model
    gmm <- mixtools::normalmixEM(cov, k = .n.reg, verb = FALSE)
    if(length(which(gmm$posterior[, 1] > gmm$posterior[, 2])) > 0 & length(which(gmm$posterior[, 1] < gmm$posterior[, 2])) > 0) {
      if (gmm$mu[1] < gmm$mu[2])
        id.clust <- which(gmm$posterior[, 1] > gmm$posterior[, 2])
      else
        id.clust <- which(gmm$posterior[, 1] < gmm$posterior[, 2])
      bnd1 <- max(cov[id.clust])
      bnd2 <- min(cov[-id.clust])
      loc.cut <- c(loc.cut, ((bnd1 + bnd2)/2))
    }

    if (circ & any(loc.cut > 360))
      loc.cut[which(loc.cut > 360)] <- loc.cut[which(loc.cut > 360)] - 360
    loc.cut <- sort(loc.cut)
    loc.cut <- c(brks[1], loc.cut, brks[length(brks)])
    return(loc.cut)
  } else {
    # If all bins have some data points
    return(NA)
  }
}
