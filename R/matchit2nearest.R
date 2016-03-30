matchit2nearest <-  function(treat, X, data, distance, discarded,
                             ratio=1, replace=FALSE, m.order="largest",
                             caliper=0, calclosest=FALSE,
                             mahvars=NULL, exact=NULL,
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
    caliper.matrix <- spatial.effects.pscore.caliper(
                        spatial.threshold, spatial.decay.function,
                        spatial.data, distance[in.sample == 1], caliper, treat)
  } else {
    caliper.matrix <- distance[in.sample == 1]
  }

  # caliper for matching (is equal to 0 if caliper matching not done)
  sd.cal <- caliper * sqrt(var(caliper.matrix, na.rm=TRUE))

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
  # exact should not equal T for this type of matching--that would get sent to
  # matchit2exact
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
    # a counter
    c.matched2 <- c.matched
    # in cases there is no replacement and all controls have been used up
    if (!(0 %in% c.matched2)) {
      match.matrix[match.matrix[, r] == 0 & !is.na(match.matrix[, r]), r] <- NA
      if (r < ratio) {
        match.matrix[, (r+1):ratio] <- NA
      }
      break
    }

    # in case there is replacement, but all units have been used in
    # previous ratios
    if (sum(!is.na(match.matrix[, r])) == 0) {
      if (r < ratio) {
        match.matrix[, (r+1):ratio] <- NA
      }
      break
    }

    # check/update which ratio we are on
    if (r != ceiling(i/(tr/ratio))) {
      r <- r+1
      t.matched <- match.matrix[, r]
    }

    if (m.order == "largest") {
      t.iter.pscore <- max(t.pscores[t.matched == 0], na.rm=T)
    }
    if (m.order == "smallest") {
      t.iter.pscore <- min(t.pscores[t.matched == 0], na.rm=T)
    }
    if (m.order == "random") {
      t.iter.pscore <- sample(t.pscores[t.matched == 0][!is.na(t.pscores[t.matched == 0])], 1)
    }

    # The treatment unit for this iteration, again resolving ties randomly
    t.iter.label <- as.vector(na.omit(t.labels[t.iter.pscore == t.pscores & t.matched == 0]))
    if (length(t.iter.label) > 1) {
      t.iter.label <- sample(t.iter.label, 1)
    }


    # Calculating all the absolute deviations in propensity scores
    # Calculate only for those eligible to be matched (c.matched == 0)
    # this first if statement only applies to replacement ratio
    # matching, so that each treatment unit is matched to a different
    # control unit than from the previous round

    # match number = NA if no units within caliper


    # Set things up for exact matching
    # Make c.matched2 == -2 if it is not an exact match
    # There might be a more efficient way to do this, but I could not figure
    # out another way to compare a vector with the matrix
    if (!is.null(exact)) {
      for (k in 1:dim(exact)[2]) {
        c.matched2[exact[t.iter.label, k] != exact[c.labels, k]] <- -2
      }
    }

    # Need to add a check in case there are not any eligible matches left...
    if (replace & r != 1) {
      if (sum(!c.labels %in% match.matrix[t.iter.label, (1:r-1)] &
              c.matched2 == 0) == 0) {
        deviation <- NULL
        min.dev <- NA
      } else {
        deviation <- abs(c.pscores[!c.labels %in% match.matrix[t.iter.label, (1:r-1)] &
                            c.matched2 == 0] - t.iter.pscore)
      }
    }
    else {
      if (sum(c.matched2 == 0) == 0) {
        deviation <- NULL
        min.dev <- NA
      } else {
        deviation <- abs(c.pscores[c.matched2 == 0] - t.iter.pscore)
      }
    }


    # ???
    # DONT REDO SPATIAL DEV EACH TIME
    # RUN IT ONCE AND FILTER OUT NA VALS
    # THEN USE IT AS LOOKUP EACH TIME FOR SPATIAL WEIGHTS

    # -------------------------------------------------------------------------
    # Spatial penalties are applied to the deviations calculated here.
    # This should be through a function which draws in the spatial data.
    # The current trick is identifying the relevant data for the current
    # iteration, but that should be doable based on the c.pscores and
    # t.pscores row titles.
    if (is.spatial == TRUE && !is.null(deviation)) {

      deviationx <- spatial.effects.pscore.t.iter.label(spatial.threshold,
                                                spatial.decay.function,
                                                spatial.data,
                                                deviation, t.iter.label)
      print(deviationx)
      print("!")
    }
    # -------------------------------------------------------------------------
    if (!is.null(deviation)) {

      if (caliper != 0) {
        if (replace & r != 1) {
          pool <- c.labels[!c.labels %in% match.matrix[t.iter.label, (1:r-1)]
                          & c.matched2 == 0][deviation <= sd.cal]
        } else {
          pool <- c.labels[c.matched2 == 0][deviation <= sd.cal]
        }

        # added for spatial?
        pool <- pool[!is.na(pool)]

        if (length(pool) == 0) {
          if (calclosest == FALSE) {
            min.dev <- NA
          } else {
            if (replace && r != 1) {
              min.dev.check <- (!c.labels %in% match.matrix[t.iter.label, (1:r-1)])
            } else {
              min.dev.check <- (c.matched2 == 0)
            }
            min.dev <- c.labels[min.dev.check][min(deviation) == deviation]

          }
        } else if (length(pool) == 1) {
          min.dev <- pool[1]
        } else if (is.null(mahvars)) {
          min.dev <- sample(pool, 1)
        } else {
          # This has the important vars for the C's within the caliper
          poolvarsC <- mahvars[pool,,drop=F]
          # Sigma is the full group var/covar matrix of Mahalvars
          mahal <- mahalanobis(poolvarsC, mahvars[t.iter.label, ], Sigma)
          min.dev <- pool[mahal == min(mahal)]
        }
      } else {
        if (replace && r != 1) {
          min.dev.check <- (!c.labels %in% match.matrix[t.iter.label, (1:r-1)] &
                           c.matched2 == 0)
        } else {
          min.dev.check <- (c.matched2 == 0)
        }
        min.dev <- c.labels[min.dev.check][min(deviation) == deviation]

      }

    }

    # Resolving ties in minimum deviation by random draw
    if (length(min.dev) > 1) {
      goodmatch <- sample(min.dev, 1)
    } else {
      goodmatch <- min.dev
    }

    # Storing which treatment unit has been matched to control, and
    # vice versa
    t.matched[t.iter.label == t.labels] <- goodmatch
    c.matched[goodmatch == c.labels] <- t.iter.label

    # instead of the in.sample, we now have an index with dimensions t.size
    # by number of matches (ratio)
    match.matrix[which(t.iter.label == t.labels), r] <- goodmatch

    # If matching with replacement, set c.matched back to 0 so it can be reused
    if (replace) {
      c.matched[goodmatch == c.labels] <- 0
    }

  }
  if (verbose){
    cat("Done\n")
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
