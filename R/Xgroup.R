setupG <- function(group, m, bilevel) {
  gf <- factor(group)
  if (any(levels(gf)=='0')) {
    g <- as.integer(gf) - 1
    lev <- levels(gf)[levels(gf)!='0']
  } else {
    g <- as.integer(gf)
    lev <- levels(gf)
  }
  if (is.numeric(group) | is.integer(group)) {
    lev <- paste0("G", lev)
  }
  if (missing(m)) {
    m <- rep(NA, length(lev))
    names(m) <- lev
  } else {
    #if (all.equal(sort(names(m)), sort(group)))
    TRY <- try(as.integer(group)==g)
    if (inherits(TRY, 'try-error') || any(!TRY)) stop('Attempting to set group.multiplier is ambiguous if group is not a factor', call.=FALSE)
    if (length(m) != length(lev)) stop("Length of group.multiplier must equal number of penalized groups", call.=FALSE)
    if (storage.mode(m) != "double") storage.mode(m) <- "double"
    if (any(m < 0)) stop('group.multiplier cannot be negative', call.=FALSE)
  }
  structure(g, levels=lev, m=m)
}

subsetG <- function(g, nz) {
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  new <- g[nz]
  dropped <- setdiff(g, new)
  if (length(dropped)) {
    lev <- lev[-dropped]
    m <- m[-dropped]
    gf <- factor(new)
    new <- as.integer(gf) - 1*any(levels(gf)=='0')
  }
  structure(new, levels=lev, m=m)
}

reorderG <- function(g, m, bilevel) {
  og <- g
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  if (any(g==0)) {
    g <- as.integer(relevel(factor(g), "0"))-1
  }
  if (any(order(g) != 1:length(g))) {
    reorder <- TRUE
    gf <- factor(g)
    if (any(levels(gf)=="0")) {
      gf <- relevel(gf, "0")
      g <- as.integer(gf) - 1
    } else {
      g <- as.integer(gf)
    }
    ord <- order(g)
    ord.inv <- match(1:length(g), ord)
    g <- g[ord]
  } else {
    reorder <- FALSE
    ord <- ord.inv <- NULL
  }
  structure(g, levels=lev, m=m, ord=ord, ord.inv=ord.inv, reorder=reorder)
}

newXG <- function(X, g, m, ncolY, bilevel) {
  # Coerce X to matrix
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (any(is.na(X))) stop("Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg", call.=FALSE)
  if (length(g) != ncol(X)) stop ("Dimensions of group is not compatible with X", call.=FALSE)
  xnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)

  # Setup group
  G <- setupG(g, m, bilevel)

  # Reconfigure for multiple outcomes, if necessary
  if (ncolY > 1) {
    X <- multiX(X, ncolY)
    G <- multiG(G, ncolY)
  }

  # Feature-level standardization
  std <- .Call("standardize", X)
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  nz <- which(scale > 1e-6)                # non-constant columns
  if (length(nz) != ncol(X)) {
    XX <- XX[, nz, drop=FALSE]
    G <- subsetG(G, nz)
  }

  # Reorder groups, if necessary
  G <- reorderG(G, attr(G, 'm'), bilevel)
  if (attr(G, 'reorder')) XX <- XX[, attr(G, 'ord')]

  # Group-level standardization
  if (!bilevel) {
    XX <- orthogonalize(XX, G)
    g <- attr(XX, "group")
  } else {
    g <- as.integer(G)
  }

  # Set group multiplier if missing
  m <- attr(G, 'm')
  if (all(is.na(m))) {
    m <- if (bilevel) rep(1, max(g)) else sqrt(table(g[g!=0]))
  }

  # Return
  return(list(X=XX, g=g, m=m, reorder=attr(G, 'reorder'), ord.inv=attr(G, 'ord.inv'), names=xnames,
              center=center, scale=scale, nz=nz))
}
