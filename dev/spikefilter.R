## spikefilter2d <- function (spcmatrix) {
##   ## expand matrix by mean values
## #  m <- apply (spcmatrix, 1, mean)
## #  spcmatrix <- cbind (m, spcmatrix, m)
##   spcmatrix <- spcmatrix [c(1, seq_len (nrow (spcmatrix), nrow (spcmatrix))), ]
## #  m <- apply (spcmatrix, 2, mean)
## #  spcmatrix <- rbind (m, spcmatrix, m)
##   spcmatrix <- spcmatrix [, c(1, seq_len (ncol (spcmatrix), ncol (spcmatrix))), ]

##   d <- t(apply (spcmatrix[-c(1, nrow (spcmatrix)),], 1, convolve, c(-1, 2, -1), type="f"))
##   d <- d + apply (spcmatrix[,-c(1, ncol (spcmatrix))], 2, convolve, c(-1, 2, -1), type="f") 
  
##   return (d) 
## }

spikefilter2d <- function (spcmatrix) {
  ## expand matrix by one row and col at each side
  spcmatrix <- spcmatrix [c (1, seq_len (nrow (spcmatrix)), nrow (spcmatrix)), ]
  spcmatrix <- spcmatrix [, c (1, seq_len (ncol (spcmatrix)), ncol (spcmatrix))]

  # filter 
  d <- t  (apply (spcmatrix, 1, filter, c(-1, 2, -1)))
  d <- d + apply (spcmatrix, 2, filter, c(-1, 2, -1)) 

  # the extra row and col are now NA, so don't return them
  d [-c (1, nrow (d)), -c (1, ncol (d))]
}

spikefilter <- function (spcmatrix) {
  ## expand matrix 
  spcmatrix <- spcmatrix [c(1, seq_len (nrow (spcmatrix)), nrow (spcmatrix)), ]

  d <- t (apply (spcmatrix, 1, filter, c(-1, 2, -1)))

  d [, -c (1, ncol (d))]
}

spikes.interactive <- function (spc, spikiness, npts = 10, nspc = 1,
                                save.tmp = 20, use.tmp = FALSE, ispikes = NULL, iispikes = NULL) {

  ## TODO: better move the first part to spike suspicio
  wavelength <- spc@wavelength
  spc <- spc@data$spc
  gc ()
  
  dim <- dim (spikiness) 
  dim(spikiness)  <- NULL # make vector
  
  if (is.null (ispikes))
    ispikes <- order (spikiness, na.last = TRUE, decreasing = TRUE)

  if (is.null (iispikes))
    iispikes <- order (ispikes)
  
  start.i = 1

  if (use.tmp & file.exists ("spikefilter.tmp.RData")){
    cat ("load temporary data\n")
    load ("spikefilter.tmp.RData")
    start.i <- i
  }

##  close.screen(all = TRUE) 
##  split.screen (figs = c(1, 2))
  layout (matrix (c (1, 0, 3, 2), nrow = 2))
 # X11 (width = 7, height = 4); wspc <- dev.cur ()
 # X11 (width = 7, height = 4); wdetail <- dev.cur () 

  save.i <- 1
  if (save.tmp > 0)
    save.tmp = save.tmp + 1
  for (i in start.i : length (ispikes)){
    cat ("Spike suspicion: ", i, "\n" )

    if ((save.tmp > 0) & (save.i %% save.tmp == 0)) {
      save (spc, i, ispikes, iispikes, spikiness, file = "spikefilter.tmp.RData", compress = FALSE)
      save.i <- 1
    } 
    
    if (is.na (ispikes[i]))
      next

    save.i <- save.i + 1

    ind <- vec2array (ispikes [i], dim = dim)
    cat ("   Spectrum: ", ind[1], "\n")
 ##   screen (1)
##    erase.screen()    
#    dev.set (wspc)
    k <- ind[1] + (-nspc : nspc) # suspicious spectrum plus the spectra around
    k <- k [k > 0]
    k <- k [k <= nrow (spc)]
    isna <- apply (spc[k,,drop = FALSE], 1, function (x) all (is.na (x)))
    k <- k[! isna]
##     plot (spc[k,                  ], "spc", 
##           col = c("black", "blue", "black"))
##     plot (spc[ind[1]           , , ind[2]           , index = TRUE], "spc", 
##           type = "p", col = "red", pch = 20, add = TRUE)

    plot (wavelength, spc[ind[1],], ylim = range (spc[k,], na.rm = TRUE), type = "n")
    for (l in k)
      lines (wavelength, spc[l,], col = if (l == ind[1]) "blue" else "black")
      
    points (wavelength[ind[2]], spc[ind[1], ind[2]], col = "red", pch = 20)
   
##    dev.set (wdetail)
##    screen (2)
##    erase.screen()    
    j <- ind[2] + (-npts : npts) # suspicious data points plus points around
    j <- j[j > 0]
    j <- j[j <= ncol (spc)]

    x <- range (wavelength[j], na.rm = TRUE)
    x <- c(x[1], (x[2] - x[1]) * 1.1 + x[1])

    yl <- min (spc[k,j,drop = FALSE], na.rm = TRUE)
    print (yl)
    yu <- median (spc[k,j,drop = FALSE], na.rm = TRUE)
    yu <- (yu - yl) * 4 + yl
    print (yu)

##     plot (spc[k, , j, index = TRUE], "spc", 
##           col = c("black", "blue", "black"), pch = 20, type = "p",
##           xlim = x,
##           ylim = y,
##           cex = 0.5)
##     plot (spc[ind[1], , j, index = TRUE], "spc", 
##           col = "blue", pch = 20, type = "p",
##           add = TRUE)
   plot (wavelength[j], spc[ind[1],j], xlim = x, ylim = c (yl, yu), type = "n")
    for (l in k)
      lines (wavelength[j], spc[l,j], pch = 20,
             type = "p",
             cex = 0.5,
             col = if (l == ind[1]) "blue" else "black")
    lines (wavelength[j], spc[ind [1],j], col = "blue")        
    points (wavelength[j], spc[ind[1], j], 
          col = "blue", pch = 20, type = "p"
          )

   y <- range (spc[k,j,drop = FALSE], na.rm = TRUE)
  
   plot (wavelength[j], spc[ind[1],j], xlim = x, ylim = y, type = "n")

    for (l in k)
      lines (wavelength[j], spc[l,j], pch = 20, type = "p", cex = 0.5,
          col = if (l == ind[1]) "blue" else "black")
          
    points (wavelength[j], spc[ind[1], j], 
          col = "blue", pch = 20, type = "p"
          )

    x <- rep (x[2], 3) 
    y <- (10:8)/10 * (y[2] - y [1]) + y[1]

    points (x, y, pch = 20, cex = 1, col = c("black", "#008000", "red"))
    text (x, y, labels = c("end", "spc OK", "bad spc"), pos = 2,
          col = c("black", "#008000", "red"))

    pts <- identify (c(wavelength[        j], x),
                     c(spc  [ind[1], j], y))

    if (length (pts) == 0)
      next

    if (max (pts) > length (j)){
      pts <- max (pts) - length (j)

      if (pts == 1) break

      # pts == 2 || pts == 3
      j <- array2vec (matrix (c(rep (ind[1], dim [2]),
                                1 : dim[2]),
                              ncol = 2),
                      dim)
      ispikes [iispikes[j]] <- NA

      if (pts == 3) # bad spectrum
        spc[ind[1], ] <- NA

      next
    }
      
    pts <- j [pts]
    #pts <- pts
    
    cat ("   Spike points: ", pts, "\n")
  
    spc[ind[1], pts] <- NA 
    
    pts <- array2vec (matrix (c (rep (ind[1], length (pts)), pts), ncol = 2), 
                      dim = dim)
    
    cat ("   Spikiness: ", spikiness [pts], "\n")
    ispikes [iispikes[pts]] <- NA     # do not look at this spike again

  }

##  close.screen(all = TRUE) 
  ##  dev.off (wdetail)
  ##  dev.off (wspc)  


#  i <- apply (spc$spc, 1, function (x) all (is.na (x)))
#  spc <- spc[!i, , short = "spikefilter: bad spectrum"]
  
  return (spc) 
  
}

spikes.NA.linapprox <- function (spc, neighbours = 1, ...){
  ispc <- which (is.na (spc@data$spc), arr.ind = TRUE)
  ispc <- unique (ispc[,"row"]) 

  for (i in ispc){
    nas <- which (is.na (spc@data$spc[i,]))
    start <- c (0, which (diff (nas) > 1)) + 1
    end <- c (start[-1] - 1, length (nas)) 
    for (j in seq (along = start)) {
      pts <- nas[start[j]] : nas[end[j]]
      xneighbours <- c(-(1:neighbours) + nas[start[j]],
                        (1:neighbours) + nas[end  [j]]) 
      xneighbours <- xneighbours[xneighbours > 0]
      xneighbours <- xneighbours[xneighbours < nwl(spc)]

      if (length (xneighbours) == 0)
        stop ("No data to interpolate from.")
      else if (length (xneighbours) == 1)
        spc@data$spc[i, pts] <- spc@data$spc [i, xneighbours]
      else
        spc@data$spc[i, pts] <- approx (x = spc@wavelength  [xneighbours],
                                        y = spc@data$spc [i, xneighbours],
                                        xout = spc@wavelength  [pts],
                                        method = "linear",
                                        rule = 2)$y
      
    }
  } 
  spc@log <- logentry (spc,
                       long = list (
                         "neighbours" = neighbours
                         ),              
                       ...            # name and user?
                       );
  
  return (spc) 
} 
