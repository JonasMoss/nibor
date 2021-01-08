function (x, H, h, gridsize, gridtype, xmin, xmax, supp = 3.7, 
    eval.points, binned = FALSE, bgridsize, positive = FALSE, 
    adj.positive, w, compute.cont = FALSE, approx.cont = TRUE, 
    unit.interval = FALSE, verbose = FALSE) 
{
    if (is.vector(x)) {
        if (missing(H)) {
            d <- 1
            n <- length(x)
        }
        else {
            if (is.vector(H)) {
                d <- 1
                n <- length(x)
            }
            else {
                x <- matrix(x, nrow = 1)
                d <- ncol(x)
                n <- nrow(x)
            }
        }
    }
    else {
        d <- ncol(x)
        n <- nrow(x)
    }
    if (!missing(w)) 
        if (!(identical(all.equal(sum(w), n), TRUE))) {
            warning("Weights don't sum to sample size - they have been scaled accordingly\n")
            w <- w * n/sum(w)
        }
    if (missing(w)) 
        w <- rep(1, n)
    if (d == 1) {
        if (missing(adj.positive)) 
            adj.positive <- abs(min(x))
        if (positive) 
            y <- log(x + adj.positive)
        else y <- x
        if (missing(h)) 
            h <- hpi(x = y, binned = default.bflag(d = d, n = n), 
                bgridsize = bgridsize)
    }
    if (missing(H) & d > 1) {
        if (binned) 
            H <- Hpi.diag(x = x, binned = binned, bgridsize = bgridsize, 
                verbose = verbose)
        else H <- Hpi(x = x, binned = default.bflag(d = d, n = n), 
            bgridsize = bgridsize, verbose = verbose)
    }
    if (binned) {
        if (d > 1) {
            if (!identical(diag(diag(H)), H)) 
                warning("Binned estimation for non-diagonal bandwidth matrix H can be inaccurate.")
        }
        if (missing(bgridsize)) 
            bgridsize <- default.gridsize(d)
        if (positive) {
            fhat <- kdde.binned(x = y, H = H, h = h, bgridsize = bgridsize, 
                xmin = xmin, xmax = xmax, w = w, deriv.order = 0)
            fhat$estimate <- fhat$estimate/(exp(fhat$eval.points))
            fhat$eval.points <- exp(fhat$eval.points) - adj.positive
            fhat$x <- x
        }
        else if (unit.interval) {
            fhat <- kde.unit.interval.1d(x = x, binned = binned)
        }
        else {
            if (missing(h) & d == 1) 
                h <- hpi(x = x, binned = default.bflag(d = d, 
                  n = n), bgridsize = bgridsize)
            fhat <- kdde.binned(x = x, H = H, h = h, bgridsize = bgridsize, 
                xmin = xmin, xmax = xmax, w = w, deriv.order = 0)
        }
        if (!missing(eval.points)) {
            fhat$estimate <- predict(fhat, x = eval.points)
            fhat$eval.points <- eval.points
        }
    }
    else {
        if (missing(gridsize)) 
            gridsize <- default.gridsize(d)
        if (d == 1) {
            if (missing(eval.points)) {
                if (unit.interval) {
                  fhat <- kde.unit.interval.1d(x = x, binned = FALSE)
                }
                else fhat <- kde.grid.1d(x = x, h = h, gridsize = gridsize, 
                  supp = supp, positive = positive, xmin = xmin, 
                  xmax = xmax, adj.positive = adj.positive, gridtype = gridtype, 
                  w = w)
            }
            else fhat <- kde.points.1d(x = x, h = h, eval.points = eval.points, 
                positive = positive, adj.positive = adj.positive, 
                w = w)
        }
        else {
            if (is.data.frame(x)) 
                x <- as.matrix(x)
            if (missing(eval.points)) {
                if (d == 2) 
                  fhat <- kde.grid.2d(x = x, H = H, gridsize = gridsize, 
                    supp = supp, xmin = xmin, xmax = xmax, gridtype = gridtype, 
                    w = w, verbose = verbose)
                else if (d == 3) 
                  fhat <- kde.grid.3d(x = x, H = H, gridsize = gridsize, 
                    supp = supp, xmin = xmin, xmax = xmax, gridtype = gridtype, 
                    w = w, verbose = verbose)
                else stop("need to specify eval.points for more than 3 dimensions")
            }
            else fhat <- kde.points(x = x, H = H, eval.points = eval.points, 
                w = w)
        }
    }
    fhat$binned <- binned
    fhat$names <- parse.name(x)
    fhat$w <- w
    class(fhat) <- "kde"
    if (compute.cont & missing(eval.points)) 
        fhat$cont <- contourLevels(fhat, cont = 1:99, approx = approx.cont)
    return(fhat)
}