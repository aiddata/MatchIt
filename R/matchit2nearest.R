matchit2nearest <-  function(treat, X, data, distance, discarded,
                             ratio=1, replace=FALSE, m.order="largest",
                             caliper=0, calclosest=FALSE, calrandom=FALSE,
                             mahvars=NULL, exact=NULL, static.selection=FALSE,
                             subclass=NULL, verbose=FALSE, sub.by=NULL,
                             is.full.mahalanobis, ...){

  # ---------------------------------------------------------------------------
  # Exceptions for when spatial information is passed to matching functions.
  if (class(distance) == "list") {
    is.spatial <- TRUE
    spatial.threshold <- distance[[4]]
    spatial.decay.function <- distance[[3]]
    spatial.data <- distance[[2]]
    distance <- distance[[1]]
  } else {
    is.spatial <- FALSE
  }
  # ---------------------------------------------------------------------------

  if (verbose) {
    cat("Nearest neighbor matching... \n")
  }
  # replace
  if (!(typeof(replace) == "logical")) {
    warning("replace=", replace, " is invalid; used replace=FALSE instead",
            call.=FALSE)
    replace <- FALSE
  }
  # m.order
  if (!(m.order %in% c("largest", "smallest", "random"))) {
  # if (!(identical(m.order,"largest") | identical(m.order,"smallest") |
  #      identical(m.order,"random"))) {
    warning("m.order=", m.order, " is invalid; used m.order='largest' instead",
            call.=FALSE)
    m.order <- "largest"
  }
  # ratio
  ratio <- round(ratio)
  if (!is.numeric(ratio) | ratio[1] < 1 |
      !identical(round(length(ratio)), 1)) {
    warning("ratio=", ratio, " is invalid; used ratio=1 instead", call.=FALSE)
    ratio <- 1
  }
  # caliper
  if (!is.vector(caliper) | !identical(round(length(caliper)), 1)) {
    warning("caliper=", caliper, " is invalid; Caliper matching not done",
            call.=FALSE)
    caliper <- 0
  }
  if (caliper < 0) {
    warning("caliper=", caliper, " is less than 0; Caliper matching not done",
            call.=FALSE)
    caliper <- 0
  }
  # calclosest
  if (!(typeof(calclosest) == "logical")) {
    warning("calclosest=", calclosest, " is invalid; used calclosest=FALSE
            instead", call.=FALSE)
    calclosest <- FALSE
  }
  # calrandom
  if (!(typeof(calrandom) == "logical")) {
    warning("calrandom=", calrandom,
            " is invalid; used calrandom=FALSE instead", call.=FALSE)
    calrandom <- FALSE
  }
  # mahvars & caliper
  if (!is.null(mahvars) & caliper[1] == 0) {
    warning("No caliper size specified for Mahalanobis matching.
            Caliper = 0.25 used.", call.=FALSE)
    caliper <- 0.25
  }
  # when mahalanobis distance is used for all covars
  if (is.full.mahalanobis) {
    mahvars <- X
    Sigma <- var(X)
    # Note: caliper irrelevant, but triggers mahalanobis matching
    caliper <- 0.25
    # no subclass with full mahalanobis
    if (!is.null(subclass)) {
      warning("No subclassification with pure Mahalanobis distance.",
              call.=FALSE)
      subclass <- NULL
    }
  }

  # count of all treated/control
  n <- length(treat)
  # count of all control
  c.size <- length(treat[treat == 0])
  # count of all treated
  t.size <- length(treat[treat == 1])
  # psm dists for control
  c.pscores <- distance[treat == 0]
  # psm dists for treated
  t.pscores <- distance[treat == 1]

  # assign range id to names if names are missing
  if (is.null(names(treat))) {
    names(treat) <- 1:n
  }

  # all labels, control labels, treated labels
  labels <- names(treat)
  c.labels <- names(treat[treat == 0])
  t.labels <- names(treat[treat == 1])

  # get rid of discarded
  in.sample <- !discarded
  names(in.sample) <- labels

  # 10/1/07: Warning for if fewer control than ratio*treated and matching
  #          without replacement
  if (c.size < ratio * t.size && replace == FALSE) {
  	if (ratio > 1) {
      warning(paste("Not enough control units for ", ratio, " matches for each
                    treated unit when matching without replacement.  Not all
                    treated units will receive", ratio, "matches"))
  	} else {
      warning(paste("Fewer control than treated units and matching without
                    replacement.  Not all treated units will receive a match.
                    Treated units will be matched in the order specified by
                    m.order:", m.order))
    }
  }

  # generating match matrix
  # row exists for all treated
  # column exists for number of units to match based on ratio
  match.matrix <- matrix(0, nrow=t.size, ncol=ratio,
                         dimnames=list(t.labels, 1:ratio))

  # Vectors of whether unit has been matched:
  #    [0]  if not matched (unit # of match if matched)
  #   [-1]  if cannot be matched due to discarded (if in.sample == 0)
  c.matched <- rep(0, length(c.pscores))
  names(c.matched) <- c.labels

  # These are the units that are ineligible because of discard
  # (in.sample == 0)
  c.matched[in.sample[c.labels] == 0] <- -1
  match.matrix[in.sample[t.labels] == 0, ] <- -1
  t.matched <- match.matrix[, 1]
  names(t.matched) <- t.labels

  # total number of matches (including ratios) = ratio * t.size
  tr <- length(match.matrix[match.matrix != -1])
  # starting match ratio
  r <- 1


  # ---------------------------------------------------------------------------
  # get matrix for caliper
  if (is.spatial == TRUE && !is.null(distance) && caliper != 0) {
    caliper.vector <- spatial.effects.pscore.caliper(
                        spatial.threshold, spatial.decay.function,
                        spatial.data, distance[in.sample == 1],
                        caliper, treat)
  } else {
    caliper.vector <- distance[in.sample == 1]
  }

  # caliper for matching (is equal to 0 if caliper matching not done)
  sd.cal <- caliper * sqrt(var(caliper.vector, na.rm=TRUE))

  # used when caliper == 0 to track # of treated units which had no eligible
  # control pool due to caliper
  empty.c.pool <- 0

  if (is.spatial == TRUE) {
    treated.units <- spatial.data[t.labels, ]
    tmp.controls <- spatial.data[c.labels, ]

    for (i in 1:nrow(tmp.controls)){

      tmp.dists <- spDistsN1(treated.units, tmp.controls[i, ],
                             longlat=TRUE)

      tmp.weights <- run.distance.decay(thresh=spatial.threshold,
                                        dist=tmp.dists,
                                        func=spatial.decay.function)

      tmp.controls$sp.dist.weights[i] <- mean(tmp.weights)
    }

    spatial.data$sp.dist.weights <- NA

    tmp.rows <- rownames(spatial.data@data) %in% rownames(tmp.controls@data)
    spatial.data[tmp.rows, 'sp.dist.weights'] <- tmp.controls$sp.dist.weights

  }

  # ---------------------------------------------------------------------------


  # Var-covar matrix for Mahalanobis (currently set for full sample)
  if (!is.null(mahvars) & !is.full.mahalanobis) {
    if (!sum(mahvars %in% names(data)) == length(mahvars)) {
	    warning("Mahvars not contained in data.  Mahalanobis matching not done.",
              call.=FALSE)

	    mahvars <- NULL

    } else {
      ww <- mahvars %in% dimnames(X)[[2]]
      nw <- length(mahvars)
      mahvars <- data[, mahvars, drop=F]
      Sigma <- var(mahvars)
      if (sum(ww) != nw) {
        X <- cbind(X, mahvars[!ww])
      }
      mahvars <- as.matrix(mahvars)

    }
  }

  # Now for exact matching within nearest neighbor.
  # exact should not equal T for this type of matching--that would get
  # sent to matchit2exact
  if (!is.null(exact)) {
    if (!sum(exact %in% names(data)) == length(exact)) {
	    warning("Exact variables not contained in data. Exact matching not
              done.", call.=FALSE)
	    exact=NULL
  	} else {
      ww <- exact %in% dimnames(X)[[2]]
      nw <- length(exact)
      exact <- data[, exact, drop=F]
      if (sum(ww) != nw) {
        X <- cbind(X, exact[!ww])
      }
    }
  }

  # Looping through nearest neighbour matching for all treatment units
  # Only do matching for units with in.sample == 1 (matched != -1)
  if (verbose) {
    trseq <- floor(seq(tr/10, tr, tr/10))
    cat("Matching Treated: ")
  }


  for (i in 1:tr) {

    # Make new matched column to be used for exact matching
    # Will only be 0 (eligible for matching) if it's an exact match
    if (verbose) {
      if (i %in% trseq) {
        cat(10*which(trseq == i), "%...", sep="")
      }
    }

    # tmp copy of c.matched used for iteration
    c.matched2 <- c.matched

    # in cases there is no replacement and all controls have been used up
    if (!(0 %in% c.matched2)) {
      warnings("controls used up during r=",r," without replacements enabled")
      match.matrix[match.matrix[, r] == 0 & !is.na(match.matrix[, r]), r] <- NA
      if (r < ratio) {
        match.matrix[, (r+1):ratio] <- NA
      }
      break
    }

    # in case there is replacement, but all units have been used in
    # previous ratios
    if (sum(!is.na(match.matrix[, r])) == 0) {
      warnings("all replacement units have been used in previous ratios (r=",
               r ,")")
      if (r < ratio) {
        match.matrix[, (r+1):ratio] <- NA
      }
      break
    }

    # check/update which ratio we are on
    if (r != ceiling(i/(tr/ratio))) {
      r <- r + 1
      t.matched <- match.matrix[, r]
    }

    # get treatment pscore for current iteration
    if (m.order == "largest") {
      t.iter.pscore <- max(t.pscores[t.matched == 0], na.rm=T)
    }
    if (m.order == "smallest") {
      t.iter.pscore <- min(t.pscores[t.matched == 0], na.rm=T)
    }
    if (m.order == "random") {
      t.iter.pscore <- sample(
        t.pscores[t.matched == 0][!is.na(t.pscores[t.matched == 0])], 1)
    }

    # get the treatment unit label for this iteration
    t.iter.label <- as.vector(na.omit(t.labels[t.iter.pscore == t.pscores &
                                               t.matched == 0]))

    # resolve treatment selection ties randomly unless static.selection
    # option is enabled
    if (length(t.iter.label) > 1) {
      if (static.selection) {
        t.iter.label <- t.iter.label[1]
      } else {
        t.iter.label <- sample(t.iter.label, 1)
      }
    }

    # calculating all the absolute deviations in propensity scores
    #
    # - calculate only for those eligible to be matched (c.matched == 0)
    #
    # - this first if statement only applies to replacement ratio
    #   matching, so that each treatment unit is matched to a different
    #   control unit than from the previous round
    #
    # - match number = NA if no units within caliper


    # set things up for exact matching
    # - make c.matched2 == -2 if it is not an exact match
    # - there might be a more efficient way to do this, but I could not
    #   figure out another way to compare a vector with the matrix
    if (!is.null(exact)) {
      for (k in 1:dim(exact)[2]) {
        c.matched2[exact[t.iter.label, k] != exact[c.labels, k]] <- -2
      }
    }

    # get available control labels
    c.labels2 <- c.labels[c.matched2 == 0]

    # prevent selecting same control for a treatment multiple times
    if (replace && r != 1) {
      c.labels2 <-
        c.labels2[!(c.labels2 %in% match.matrix[t.iter.label, (1:r-1)])]
    }


    # check if there are any eligible matches left
    if (length(c.labels2) == 0) {
      c.deviations <- NULL
      c.match.good <- NA
    } else {
      c.deviations <- abs(c.pscores[c.labels2] - t.iter.pscore)
    }


    # -------------------------------------------------------------------------
    if (is.spatial == TRUE && !is.null(c.deviations)) {

      # spatial penalties are applied to c.deviations
      # c.deviations <- spatial.effects.pscore.deviation(spatial.threshold,
      #                                                  spatial.decay.function,
      #                                                  spatial.data,
      #                                                  c.deviations,
      #                                                  t.iter.label)


      c.deviations <- 0.5 * (c.deviations +
                       spatial.data@data[names(c.deviations), ]$sp.dist.weights)

      # update c.labels2 to remove units with NA deviations and update
      # deviations
      c.labels2 <- c.labels2[!is.na(c.deviations)]
      c.deviations <- c.deviations[c.labels2]

      if (length(c.deviations) == 0) {
        c.deviations <- NULL
        c.match.good <- NA
      }

    }
    # -------------------------------------------------------------------------


    if (!is.null(c.deviations)) {

      if (caliper != 0) {

        c.match.pool <- c.labels2[c.deviations <= sd.cal]

        if (length(c.match.pool) == 0) {
          empty.c.pool <- empty.c.pool + 1

          if (calclosest == FALSE) {
            c.match.good <- NA
          } else {
            c.match.good <- c.labels2[c.deviations == min(c.deviations)]
          }

        } else if (length(c.match.pool) == 1) {
          c.match.good <- c.match.pool[1]

        } else if (is.null(mahvars)) {
          if (calrandom) {
            c.match.good <- sample(c.match.pool, 1)
          } else {
            min.eligible.deviation <- min(c.deviations[c.match.pool])
            c.match.good <- c.labels2[c.deviations == min.eligible.deviation]
          }

        } else {
          # This has the important vars for the C's within the caliper
          poolvarsC <- mahvars[c.match.pool, , drop=F]
          # Sigma is the full group var/covar matrix of Mahalvars
          mahal <- mahalanobis(poolvarsC, mahvars[t.iter.label, ], Sigma)
          c.match.good <- c.match.pool[mahal == min(mahal)]
        }

      } else {
        c.match.good <- c.labels2[c.deviations == min(c.deviations)]

      }

    }

    # resolving ties in minimum deviation by random draw unless
    # static.selection option is enabled
    if (length(c.match.good) > 1) {
      if (static.selection) {
        c.match.final <- c.match.good[1]
      } else {
        c.match.final <- sample(c.match.good, 1)
      }
    } else {
      c.match.final <- c.match.good
    }

    # Storing which treatment unit has been matched to control, and
    # vice versa
    t.matched[t.labels == t.iter.label] <- c.match.final
    c.matched[c.labels == c.match.final] <- t.iter.label

    # If matching with replacement, set c.matched back to 0 so it can be reused
    if (replace) {
      c.matched[c.labels == c.match.final] <- 0
    }

    # instead of the in.sample, we now have an index with dimensions t.size
    # by number of matches (ratio)
    match.matrix[which(t.labels == t.iter.label), r] <- c.match.final

  }
  if (verbose){
    cat("Done\n")
  }

  if (empty.c.pool > 0) {
    warning(empty.c.pool,
            " treated units had an empty control pool due to caliper.")
  }

  x <- as.matrix(match.matrix)
  x[x == -1] <- NA

  # Calculate weights and return the results
  res <- list(match.matrix = match.matrix,
              weights = weights.matrix(match.matrix, treat, discarded),
              X = X)


  # Subclassifying
  if (!is.null(subclass)) {
    if (is.null(sub.by)) {
      sub.by = "treat"
    }
    psres <- matchit2subclass(treat, X, data, distance, discarded,
                              match.matrix=match.matrix, subclass=subclass,
                              verbose=verbose, sub.by=sub.by, ...)
    res$subclass <- psres$subclass
    res$q.cut <- psres$q.cut
    class(res) <- c("matchit.subclass", "matchit")
  } else {
    class(res) <- "matchit"
  }

  return(res)
}
