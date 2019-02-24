#' Predict Wind Power Output by Using a Multivariate Power Curve
#'
#' Takes multiple environmental variable inputs measured on an operating wind
#' farm and predicts the wind power output under the given environmental
#' condition.
#'
#' @param y An \eqn{n}-dimensional vector or a matrix of size \eqn{n} by 1
#'   containing wind power output data. This along with x trains the
#'   multidimensional power curve model.
#' @param x An \eqn{n} by \eqn{p} matrix or a data frame containing the input
#'   data for \eqn{p} predictor variables (wind and weather variables). This
#'   \code{x} must have the same number of rows as \code{y}, i.e., \eqn{n}.
#' @param x.new A matrix or a data frame containing new input conditions of the
#'   \eqn{p} predictor variables for which a prediction of wind power output
#'   will be made. This is an optional parameter and will be set to \code{x} by
#'   default, if it is not supplied.
#' @param id.spd The column number of \code{x} (and of \code{x.new}, if
#'   supplied) indicating wind speed data . Default to \code{1}.
#' @param id.dir The column number of \code{x} (and of \code{x.new}, if
#'   supplied) indicating wind direction data. Default to \code{NA}, but this
#'   parameter needs to be set if \code{x} includes wind direction data.
#' @note \itemize{ \item This function is developed for wind power prediction.
#'   As such, the response \code{y} represents wind power output and the
#'   covariates \code{x} include multiple wind and weather variables that
#'   potentially affect the power output. \item The data matrix \code{x} is
#'   expected to include at least wind speed and wind direction data. As
#'   measurements of other environmental variables become available, they can be
#'   added to the \code{x}. Typically, the first column of \code{x} corresponds
#'   to wind speed data and the second column to wind direction data and, as
#'   such, \code{id.spd = 1} and \code{id.dir = 2}. \item If \code{x} has a
#'   single variable of wind speed, i.e., \eqn{p = 1} and \code{id.spd = 1},
#'   this function returns an estimate (or prediction) of the Nadaraya-Watson
#'   estimator with a Gaussian kernel by using the \code{ksmooth} function in
#'   the \pkg{stats} package. }
#' @return A vector representing the predicted power output for the new
#'   wind/weather condition specified in \code{x.new}. If \code{x.new} is not
#'   supplied, this function returns the fitted power output for the given
#'   \code{x}.
#' @seealso \code{\link{windpw}}, \code{\link[stats]{ksmooth}}
#' @references Lee, G., Ding, Y., Genton, M.G., and Xie, L. (2015) Power Curve
#'   Estimation with Multivariate Environmental Factors for Inland and Offshore
#'   Wind Farms, \emph{Journal of the American Statistical Association}
#'   110(509):56-67.
#' @examples
#' head(windpw)
#'
#'
#' ### Power curve estimation.
#'
#' # By using a single input of wind speed.
#' pwcurv.est <- kp.pwcurv(windpw$y, windpw$V)
#'
#' # By using wind speed and direction: id.dir needs to be set.
#' pwcurv.est <- kp.pwcurv(windpw$y, windpw[, c('V', 'D')], id.dir = 2)
#'
#' # By using full covariates: confirm whether id.spd and id.dir are correctly specified.
#' pwcurv.est <- kp.pwcurv(windpw$y, windpw[, c('V', 'D', 'rho', 'I', 'Sb')], id.spd = 1, id.dir = 2)
#'
#'
#' ### Wind power prediction.
#'
#' # Suppose only 90% of data are available and use the rest 10% for prediction.
#' df.tr <- windpw[1:900, ]
#' df.ts <- windpw[901:1000, ]
#' id.cov <- c('V', 'D', 'rho', 'I', 'Sb')
#' pred <- kp.pwcurv(df.tr$y, df.tr[, id.cov], df.ts[, id.cov], id.dir = 2)
#'
#'
#' ### Evaluation of wind power prediction based on 10-fold cross validation.
#'
#' # Partition the given dataset into 10 folds.
#' index <- sample(1:nrow(windpw), nrow(windpw))
#' n.fold <- round(nrow(windpw) / 10)
#' ls.fold <- rep(list(c()), 10)
#' for(fold in 1:9) {
#'   ls.fold[[fold]] <- index[((fold-1)*n.fold+1):(fold*n.fold)]
#' }
#' ls.fold[[10]] <- index[(9*n.fold+1):nrow(windpw)]
#'
#' # Predict wind power output.
#' pred.res <- rep(list(c()), 10)
#' id.cov <- c('V', 'D', 'rho', 'I', 'Sb')
#' for(k in 1:10) {
#'   id.fold <- ls.fold[[k]]
#'   df.tr <- windpw[-id.fold, ]
#'   df.ts <- windpw[id.fold, ]
#'   pred <- kp.pwcurv(df.tr$y, df.tr[, id.cov], df.ts[, id.cov], id.dir = 2)
#'   pred.res[[k]] <- list(obs = df.ts$y, pred)
#' }
#'
#' # Calculate rmse and its mean and standard deviation.
#' rmse <- sapply(pred.res, function(res) with(res, sqrt(mean((obs - pred)^2))))
#' mean(rmse)
#' sd(rmse)
#' @export
kp.pwcurv <- function(y, x, x.new = x, id.spd = 1, id.dir = NA) {
  if (missing(y) | is.null(y) | !is.numeric(y))
    stop("Numeric y must be supplied.")
  if (missing(x) | is.null(x))
    stop("Covariate(s) must be supplied.")
  if (!is.null(dim(y)))
    if (dim(y)[2] > 1)
      stop("y must be univariate.")
  if (is.null(dim(x))) {
    if (length(y) != length(x))
      stop("y and x must have the same number of data points.")
  } else if (length(y) != dim(x)[1])
    stop("y and x must have the same number of data points.")
  if (!is.null(dim(x))) {
    if (dim(x)[2] < id.spd)
      stop("id.spd cannot be greater than the number of columns in x.")
    if (!is.na(id.dir) & dim(x)[2] < id.dir)
      stop("id.dir cannot be greater than the number of columns in x.")
    if (dim(x)[2] >= 2 & (is.na(id.spd) | is.na(id.dir)))
      stop("For multivariate analysis, wind speed and direction must be included.")
    if (dim(x)[2] != dim(x.new)[2])
      stop("x and x.new must have the same number of columns.")
  } else {
    if (is.na(id.spd) | (id.spd != 1) | (!is.na(id.dir) & id.dir == 1))
      stop("The single predictor must be wind speed.") else if (!is.na(id.dir) & id.dir > 1)
      stop("id.dir cannot be greater than the number of columns in x.")
  }

  if (is.null(dim(x)) | is.vector(x)) {
    bw <- KernSmooth::dpill(x, y)
    id.seq <- order(x.new)
    est <- rep(NA, length(x.new))
    est[id.seq] <- stats::ksmooth(x, y, kernel = "normal", bandwidth = bw, n.points = length(x.new), x.points = x.new)$y
  } else if (is.matrix(x)) {
    if (ncol(x) == 1) {
      bw <- KernSmooth::dpill(x, y)
      id.seq <- order(x.new)
      est <- rep(NA, length(x.new))
      est[id.seq] <- stats::ksmooth(x, y, kernel = "normal", bandwidth = bw, n.points = length(x.new), x.points = x.new)$y
    } else {
      bw <- bw.adp(y, x, id.dir)
      est <- AMK(y, x, x.new, bw, id.spd, id.dir)
    }
  } else {
    bw <- bw.adp(y, x, id.dir)
    est <- AMK(y, x, x.new, bw, id.spd, id.dir)
  }
  return(est)
}
