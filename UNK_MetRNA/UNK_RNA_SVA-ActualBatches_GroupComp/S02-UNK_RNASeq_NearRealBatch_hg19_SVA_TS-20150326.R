############################################################
# TO-DO LIST.
############################################################
DONE. # 1. Time Series Analyses for SVA Stitched UNK RNA-Seq data (SVA).





##################################################################################
##################################################################################
##################################################################################
### QN U-TEST SUMMARY:                                                         ###
###                                                                            ###
### HRV1 vs HEALTHY: FDR_P < 0.10; 171 Hits.                                   ###
### HRV2 vs HEALTHY: P < 0.005; 114 Hits.                                      ###
### HRV3 vs HEALTHY: FDR_P < 0.05; 153 Hits.                                   ###
### RSV  vs HEALTHY: FDR_P < 0.05; 2388 Hits.                                  ###
### ADV1 vs HEALTHY: FDR_P < 0.05; 698 Hits.                                   ###
### ADV2 vs HEALTHY: FDR_P < 0.05; 536 Hits.                                   ###
### HRVALL vs HEALTHY: FDR_P < 0.05; 606 Hits.                                 ###
### ADVALL vs HEALTHY: FDR_P < 0.05; 1056 Hits.                                ###
### INFALL vs HEALTHY: FDR_P < 0.05; 1223 Hits.                                ###
### LSP (Cluster 1): FDR_P < 0.05; 123 Hits.                                   ###
### LSP (Cluster 2): FDR_P < 0.05; 165 Hits.                                   ###
### AUTOCORRELATED (Cluster 1): 1.645 sigma (One-Tailed); 1422 Hits.           ###
### AUTOCORRELATED (Cluster 2): 1.645 sigma (One-Tailed); 633 Hits.            ###
##################################################################################
##################################################################################
##################################################################################





############################################################
# A. Library Used.
############################################################
############################################################
# A.1. Full Libraries.
############################################################
library("gplots")
library("ggplot2")
library("reshape2")
library("preprocessCore")
library("coin")
library("magic")
library("RColorBrewer")
library("limma")
library("LDheatmap")
library("MASS")
library("lomb")
library("cts")
library("Sushi")
library("heatmap3")
library("geneplotter")


############################################################
# A.2. Modified Functions.
############################################################
############################################################
# A.2.1. LSP-MOD.
############################################################
# LSP-MOD (Version 3.2).
lsp_mod <- function (x, ofac = 4, alpha = 0.01, plot = TRUE, ...) 
{
	# x is a matrix with 2 columns: 1st column is time and 2nd column is reading.
	alpha <- alpha
	names <- colnames(x)
    times <- x[, 1]
    x <- x[, 2]
    times <- times[!is.na(x)]
    x <- x[!is.na(x)]
    times <- as.numeric(times)
    o <- order(times)
    times <- times[o]
    x <- x[o]
    y <- cbind(times, x)
    colnames(y) <- names
    datanames <- colnames(y)
    t <- y[, 1]
    tm <- t - t[1]
    y <- y[, 2]
    n <- length(y)
    tspan <- tm[n] - tm[1]
    fr.d <- 1/tspan
    step <- 1/(tspan * ofac)
    f.max <- floor(0.5 * n * ofac) * step
    # CORRECTED FORMULA.
    freq <- seq(fr.d/ofac, f.max, by = step)
    n.out <- length(freq)
    samplemean <- mean(x)
    y <- y - samplemean
    norm <- 1/(2 * var(y))
    ####################################
    # Additional Information
    # Reference: Hocke K and Kampfer N. Gap filling and noise reduction of unevenly sampled data by means of the Lomb-Scargle periodogram. Atmos. Chem. Phys. 2009(9):4197-4206.
    tave <- (tm[1] + tm[n]) / 2
    tdif <- tspan
    pnow <- step
    arg <- 2 * pi * ((tm - tave) * step)
    wpr <- -2 * sin(0.5 * arg)^2
    wpi <- sin(arg)
    wr <- cos(arg)
    wi <- wpi
    ####################################
    w <- 2 * pi * freq
    PN <- rep(0, n.out)
    A <- rep(0, n.out)
    B <- rep(0, n.out)
    C <- rep(0, n.out)
    D <- rep(0, n.out)
    tauall <- rep(0, n.out)
    ph <- c()
    ph1 <- c()
    teta <- c()
    px <- c()
    for (i in 1:n.out) {
        wii <- w[i]
        # CORRECTED FORMULA.
        tauall[i] <- 0.5 * atan2(sum(sin(2 * wii * tm)), sum(cos(2 * wii * tm)))/wii
        tau <- tauall[i]
        arg <- wii * (tm - tau)
        cs <- cos(arg)
        sn <- sin(arg)
        A <- (sum(y * cs))^2
        B <- sum(cs * cs)
        C <- (sum(y * sn))^2
        D <- sum(sn * sn)
#       PN[i] <- A/B + C/D
	    ####################################
	    # Additional Information
	    # Reference: Hocke K and Kampfer N. Gap filling and noise reduction of unevenly sampled data by means of the Lomb-Scargle periodogram. Atmos. Chem. Phys. 2009(9):4197-4206.
	    px[i] <- pnow
	    sumsh <- sum(wr * wi)
	    sumc <- sum((wr - wi) * (wr + wi))
	    wtau <- 0.5 * atan2(2 * sumsh , sumc)
	    swtau <- sin(wtau)
	    cwtau <- cos(wtau)
	    ss <- wi * cwtau - wr * swtau
	    cc <- wr * cwtau + wi * swtau
	    sums <- sum(ss^2)
	    sumc <- sum(cc^2)
	    sumsy <- sum(y * ss)
	    sumcy <- sum(y * cc)
	    wtemp <- wr
	    wr <- wr * wpr - wi * wpi + wr
	    wi <- wi * wpr + wtemp * wpi + wi
	    iy <- sumsy / sqrt(sums) # imaginary part of Lomb-Scargle spectral component
	    ry <- sumcy / sqrt(sumc) # real part 
	    PN[i] <- ry^2 + iy^2 # power
	    # here, the FFT phase is computed from the Lomb-Scargle Phase 
		# at each new frequency 'pnow' by adding the phase shift 'arg0'
		phLS <- atan2(iy, ry)            # phase of Lomb-Scargle spectrum 
		arg0 <- 2 * pi * (tave+t[1]) * pnow + wtau  # phase shift with respect to 0
		arg1 <- 2 * pi * tave * pnow + wtau   # phase shift for FFT reconstruction 
		ph[i] <- (phLS+arg0) %% (2 * pi)  # phase with respect to 0
		ph1[i] <- (phLS+arg1) %% (2 * pi) # phase for complex FFT spectrum	
		pnow <- pnow + 1 / (ofac * tdif)    # next frequency
		teta[i] <- arg1 %% (2 * pi)
	    ####################################
    }
    PN <- norm * PN
    PN.max <- max(PN)
    peak.freq <- freq[PN == PN.max]
    peak.at <- c(peak.freq, 1/peak.freq)
    scanned <- freq
    effm <- 2 * n.out/ofac
    level <- -log(1 - (1 - alpha)^(1/effm))
    exPN <- exp(-PN.max)
    p <- effm * exPN
    if (p > 0.01) 
        p <- 1 - (1 - exPN)^effm
    ####################################
    # Additional Information
    # Reference: Hocke K and Kampfer N. Gap filling and noise reduction of unevenly sampled data by means of the Lomb-Scargle periodogram. Atmos. Chem. Phys. 2009(9):4197-4206.
    wmax <- 2 * pi * max(peak.freq)
    taumax <- 0.5 * atan2(sum(sin(2 * wmax * tm)), sum(cos(2 * wmax * tm)))/wmax
    AFT <- ((2 * n.out + 1) * var(y) * PN / 2)^0.5 # Dimention of the complex Fourier spectrum F is 2 * n.out + 1.
    argmax <- wmax * (tm - taumax)
    csmax <- cos(argmax)
    snmax <- sin(argmax)
    Amax <- (sum(y * csmax))^2
    Cmax <- (sum(y * snmax))^2
    PhiFT = atan2(C^0.5, A^0.5) + w * tave + w * tauall
    RFT <- AFT * cos(ph1)
    IFT <- AFT * sin(ph1)
    ph <- (ph + 10 * pi) %% (2 * pi)
    RFT2 <- rev(RFT)
    IFT2 <- rev(-IFT)
    RFT2p <- RFT2[2:length(RFT2)]
    IFT2p <- IFT2[2:length(IFT2)]
    FT <- complex(real = c(samplemean, RFT, RFT2), imaginary = c(0, IFT, IFT2))
    InvFT <- fft(FT, inverse = TRUE) / length(FT)
    FTp <- complex(real = c(samplemean, RFT, RFT2p), imaginary = c(0, IFT, IFT2p))
    InvFTp <- fft(FTp, inverse = TRUE) / length(FTp)
    xresample <- Re(InvFTp)[(length(InvFTp) - length(t) + 1):length(InvFTp)] # Take the last cycle.
    xresample <- rev(xresample) # It seems that the order has reversed.
    xresample.one.more <- Re(InvFT)[(length(InvFT) - length(t) + 1):length(InvFT)] # Take the last cycle.
    xresample.one.more <- rev(xresample.one.more) # It seems that the order has reversed.
    ####################################
    sp.out <- list(scanned = scanned, power = PN, data = datanames, 
        n = n, type = "frequency", ofac = ofac, n.out = n.out, alpha = alpha, 
        sig.level = level, peak = PN.max, peak.at = peak.at, 
        p.value = p, omega = w, tau = tauall, omegamax = wmax, taumax = taumax, t_average = tave, y = x, y.minus.mean = y, x_resample = xresample, real_part_resample = Re(InvFTp), imaginary_part_resample = Im(InvFTp), FT.Same.Length = FTp, InvFT.Same.Length = InvFTp, x_resample.one.more = xresample.one.more, real_part_resample.one.more = Re(InvFT), imaginary_part_resample.one.more = Im(InvFT), FT.One.More = FT, InvFT.One.More = InvFT)
    class(sp.out) <- "lsp"
    if (plot) {
        plot(sp.out, ...)
        return(invisible(sp.out))
    }
    else return(sp.out)
}


############################################################
# A.2.2. HEATMAP.2-MOD.
############################################################
# APPENDIX.1. Hijacking heatmap.2 to Calculate Row / Column Medians instead of Row / Column Means in Row / Column Scaling
heatmap.median <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- apply(x, 1, median, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- apply(x, 2, median, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!missing(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    x.scaled <- x.scaled - x.scaled[,7]
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}


# APPENDIX.2. Hijacking heatmap.2 to Use Day 255 as Baseline in Row Scaling.
heatmap.d255 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- x[,7]
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- apply(x, 2, median, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!missing(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    x.scaled <- x.scaled - x.scaled[,7]
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}


# APPENDIX.3. Further Play with heatmap.2 to Display Different Colors
heatmap.colorplay <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, perl = TRUE,
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- apply(x, 1, mean, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- apply(x, 2, mean, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 501
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    
    # Gray12-White-Gold
	my_palette1 <- colorRampPalette(c("gray12", "white", "gold"))(n = 500)
    # Gray12-White-Red
	my_palette2 <- colorRampPalette(c("gray12", "white", "red"))(n = 500)
    # Gray12-White-Blue
	my_palette3 <- colorRampPalette(c("gray12", "white", "blue"))(n = 500)
    # Gray12-White-Purple
	my_palette4 <- colorRampPalette(c("gray12", "white", "turquoise"))(n = 500)

    # # Turquoise-White-Gold
	# my_palette1 <- colorRampPalette(c("turquoise", "white", "gold"))(n = 500)
    # # Blue-White-Red
	# my_palette2 <- colorRampPalette(c("blue", "white", "red"))(n = 500)
    # # Green-White-Orange
	# my_palette3 <- colorRampPalette(c("green", "white", "orange"))(n = 500)
    # # Gray12-White-Purple
	# my_palette4 <- colorRampPalette(c("gray12", "white", "purple"))(n = 500)

    # # White-Turquoise
	# my_palette1 <- colorRampPalette(c("white", "gold"))(n = 500)
    # # White-Red
	# my_palette2 <- colorRampPalette(c("white", "red"))(n = 500)
    # # White-Gold
	# my_palette3 <- colorRampPalette(c("white", "blue"))(n = 500)
    # # White-Purple
	# my_palette4 <- colorRampPalette(c("white", "turquoise"))(n = 500)

	xtr <- t(x)
	rownamex <- row.names(xtr)
	color <-  matrix(nrow = nr, ncol = 500)
	par(mfrow = c(nr,1), mar=c(0,5,0,5))
	for (i in nr:1) {
		if (length(name <- grep("HRV", rownamex[i], ignore.case = TRUE))) {
			color[i,] <- my_palette1
		} else if (length(name <- grep("RSV", rownamex[i], ignore.case = TRUE))) {
			color[i,] <- my_palette2
		} else if (length(name <- grep("ADV", rownamex[i], ignore.case = TRUE))) {
			color[i,] <- my_palette3
		} else {
			color[i,] <- my_palette4
		}
		mat <- matrix(rnorm(nr), nrow=nr, ncol=1)
		mat <- as.matrix(x[,i])
		image(mat, axes = FALSE, col = color[i,], breaks = breaks, ...)
	}

	# image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        # c(0, nr), axes = FALSE, xlab = "", ylab = "", col = color[1:nr,], 
        # ...)
    
    
    
    # image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        # c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        # breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}


############################################################
# A.3. OTHERS.
############################################################
# Package A2R:
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
# colored dendrogram op = par(bg = "#EFEFEF") A2Rplot(hc, k = 3, boxes = FALSE, col.up = "gray50", col.down = c("#FF6B6B","#4ECDC4", "#556270"))
# another colored dendrogram
op = par(bg = "gray15") cols = hsv(c(0.2, 0.57, 0.95), 1, 1, 0.8) A2Rplot(hc, k = 3, boxes = FALSE, col.up = "gray50", col.down = cols)


############################################################
# B. QUANTILE NORMALIZED LOG2 DATA LOADING AND RESCALING.
############################################################
# Load data.

dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_NearRealBatch_ComBat_Rescaled.txt",  header =TRUE, sep = "\t")



# Colume Names.
# Column 1: GeneID;
# Columns 2-58: UNK (1-60);
# Columns 59-60: UNK-07 & 08 RiboZero;
# GeneID	U01RNA	U02RNA	U03RNA	U05RNA	U06RNA	U07RNA	U08RNA	U09RNA	U10RNA	U11RNA	U12RNA	U13RNA	U14RNA	U15RNA	U16RNA	U17RNA	U18RNA	U19RNA	U20RNA	U21RNA	U23RNA	U24RNA	U25RNA	U26RNA	U27RNA	U28RNA	U29RNA	U30RNA	U31RNA	U32RNA	U33RNA	U34RNA	U35RNA	U36RNA	U37RNA	U38RNA	U39RNA	U40RNA	U41RNA	U42RNA	U43RNA	U45RNA	U46RNA	U47RNA	U48RNA	U49RNA	U50RNA	U51RNA	U52RNA	U53RNA	U54RNA	U55RNA	U56RNA	U58RNA	U58_5RNA	U59RNA	U60RNA	U07RNARZ	U08RNARZ

# GeneID	D_0	D_4	D_21	D_116	D_185	D_186	D_255	D_289	D_290	D_292	D_294	D_297	D_301	D_307	D_311	D_322	D_329	D_369	D_380	D_400	D_476	D_532	D_546	D_602	D_615	D_616	D_618	D_620	D_625	D_630	D_647	D_679	D_680	D_683	D_688	D_694	D_700	D_711	D_735	D_796	D_840	D_912	D_944	D_945	D_948	D_959	D_966	D_984	D_1029	D_1030	D_1032	D_1038	D_1045	D_1051	D_1060	D_1109	D_1124	D_186RZ	D_255RZ

dat_data <- as.matrix(dat[,c(2:58)])
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
row.names(dat_data) <- dat$GeneID
# Alternative Column Names (by Day):
colID <- c("D_0", "D_4", "D_21", "D_116", "D_185", "D_186", "D_255", "D_289", "D_290", "D_292", "D_294", "D_297", "D_301", "D_307", "D_311", "D_322", "D_329", "D_369", "D_380", "D_400", "D_476", "D_532", "D_546", "D_602", "D_615", "D_616", "D_618", "D_620", "D_625", "D_630", "D_647", "D_679", "D_680", "D_683", "D_688", "D_694", "D_700", "D_711", "D_735", "D_796", "D_840", "D_912", "D_944", "D_945", "D_948", "D_959", "D_966", "D_984", "D_1029", "D_1030", "D_1032", "D_1038", "D_1045", "D_1051", "D_1060", "D_1109", "D_1124")
colnames(dat_data) <- colID

datrownames <- row.names(dat_data)
datcolnames <- colnames(dat)[2:length(dat[1,])]
total_matrix_unk_qn <- data.matrix(dat_data[,1:57])
row.names(total_matrix_unk_qn) <- datrownames
colnames(total_matrix_unk_qn) <- datcolnames[1:57]

max(total_matrix_unk_qn)
[1] 1
min(total_matrix_unk_qn)
[1] -1


############################################################
# C. LOMB-SCARGLE TRANSFORMATION.
############################################################
############################################################
# SAME BATCH:
# 	D: 0 - 602;
# 	D: 615 - 1124;

# GROUP DETERMINATION:
# 	HRV1 (Days 0-21; 3 Columns): Columns 1-3;
# 	HRV2 (Days 615-630; 6 Columns): Columns 25-30;
# 	HRV3 (Days 1029-1060; 7 Columns): Columns 49-55;
# 	RSV (Days 289-311; 8 Columns): Columns 8-15;
# 	ADV1 (Days 679-700; 6 Columns): Columns 32-37;
# 	ADV2 (Days 944-984; 6 Columns): Columns 43-48;
# 	HLT (Days 116-255, 322-602, 647, 711-912, 1109-1124; 21 Columns): Columns 4-7, 16-24, 31, 38-42, 56-57.
############################################################
############################################################
# C.1. LOMB-SCARGLE TRANSFORMATION.
############################################################
# Day Series.
dayseries <- c(0, 4, 21, 116, 185, 186, 255, 289, 290, 292, 294, 297, 301, 307, 311, 322, 329, 369, 380, 400, 476, 532, 546, 602, 615, 616, 618, 620, 625, 630, 647, 679, 680, 683, 688, 694, 700, 711, 735, 796, 840, 912, 944, 945, 948, 959, 966, 984, 1029, 1030, 1032, 1038, 1045, 1051, 1060, 1109, 1124)


############################################################
# C.1.1. SIGNIFICANT PERIODIC HITS WITH LSP_MOD (FDR-ADJUSTED P-VALUE < 0.05).
############################################################
# Calculate for Significant Periodic Hits.
# Will also Calculate All Lines at the Same Time.
periodhits <- c() # Hits that have Lome-Scargle transformed p-value < 0.05.
rownameperiodhits <- c() # Row name of the periodhits.
pv <- c() # LSP p-value.
scannedv <- c() # Value for $scanned.
powerv <- c() # Value for $power.
omegav <- c() # Value for $omega.
tauv <- c() # Value for $tau.
x_resamplev <- c() # Value for $x_resample.
siglevelv <- c() # Value for $sig.level.
peakv <- c() # Value for $peak.
peakatv <- c() # Value for $peak.at.
inverseffthits <- c() # Hits that have Lome-Scargle transformed p-value < 0.05.
new_p <- c() # FDR-Adjusted p-Value.
dat_data <- total_matrix_unk_qn
# 9335 Analytes Total.
for (i in 1:length(total_matrix_unk_qn[,1])) {
	if (sum(abs(dat_data[i,]-mean(dat_data[i,]))) == 0) {next}
	tsmatrix <- cbind(dayseries, dat_data[i,])
	tmpa <- lsp_mod(tsmatrix, ofac = 4, plot=FALSE)
#		plot.lsp(tmpa, main = row.names(dat_data)[i])
		pv <- rbind(pv, tmpa$p.value)
		scannedv <- rbind(scannedv, tmpa$scanned)
		powerv <- rbind(powerv, tmpa$power)
		omegav <- rbind(omegav, tmpa$omega)
		tauv <- rbind(tauv, tmpa$tau)
		siglevelv <- rbind(siglevelv, tmpa$sig.level)
		peakv <- rbind(peakv, tmpa$peak)
		peakatv <- rbind(peakatv, tmpa$peak.at)
		periodhits <- rbind(periodhits, dat_data[i,])
		x_resamplev <- rbind(x_resamplev, tmpa$x_resample)
		rownameperiodhits <- rbind(rownameperiodhits, row.names(dat_data)[i])
}
new_p <- p.adjust(pv, method = "fdr", n = length(total_matrix_unk_qn[,1]))

write(c("GeneID", "P-Value", "FDR-P", "Sig_Level", "Peak", "Peak_At_Frequency", "Peak_At_Period", colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL.txt", sep = "\t", ncolumn = 64)
write.table(cbind(rownameperiodhits, pv, new_p, siglevelv, peakv, peakatv, periodhits), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL.txt", sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
write(c("GeneID", "P-Value", "FDR-P", "Sig_Level", "Peak", "Peak_At_Frequency", "Peak_At_Period", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39", "V40", "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49", "V50", "V51", "V52", "V53", "V54", "V55", "V56", "V57"), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled.txt", sep = "\t", ncolumn = 64)
write.table(cbind(rownameperiodhits, pv, new_p, siglevelv, peakv, peakatv, x_resamplev), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled.txt", sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
write(c("Omega", "GeneID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39", "V40", "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49", "V50", "V51", "V52", "V53", "V54", "V55", "V56", "V57", "V58", "V59", "V60", "V61", "V62", "V63", "V64", "V65", "V66", "V67", "V68", "V69", "V70", "V71", "V72", "V73", "V74", "V75", "V76", "V77", "V78", "V79", "V80", "V81", "V82", "V83", "V84", "V85", "V86", "V87", "V88", "V89", "V90", "V91", "V92", "V93", "V94", "V95", "V96", "V97", "V98", "V99", "V100", "V101", "V102", "V103", "V104", "V105", "V106", "V107", "V108", "V109", "V110", "V111", "V112", "V113", "V114"), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_OmegaTau.txt", sep = "\t", ncolumn = 116)
write.table(cbind(rownameperiodhits, omegav), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_OmegaTau.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
write(c("Tau", "GeneID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39", "V40", "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49", "V50", "V51", "V52", "V53", "V54", "V55", "V56", "V57", "V58", "V59", "V60", "V61", "V62", "V63", "V64", "V65", "V66", "V67", "V68", "V69", "V70", "V71", "V72", "V73", "V74", "V75", "V76", "V77", "V78", "V79", "V80", "V81", "V82", "V83", "V84", "V85", "V86", "V87", "V88", "V89", "V90", "V91", "V92", "V93", "V94", "V95", "V96", "V97", "V98", "V99", "V100", "V101", "V102", "V103", "V104", "V105", "V106", "V107", "V108", "V109", "V110", "V111", "V112", "V113", "V114"), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_OmegaTau.txt", sep = "\t", ncolumn = 116, append = TRUE)
write.table(cbind(rownameperiodhits, tauv), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_OmegaTau.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
write(c("Scanned", "GeneID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39", "V40", "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49", "V50", "V51", "V52", "V53", "V54", "V55", "V56", "V57", "V58", "V59", "V60", "V61", "V62", "V63", "V64", "V65", "V66", "V67", "V68", "V69", "V70", "V71", "V72", "V73", "V74", "V75", "V76", "V77", "V78", "V79", "V80", "V81", "V82", "V83", "V84", "V85", "V86", "V87", "V88", "V89", "V90", "V91", "V92", "V93", "V94", "V95", "V96", "V97", "V98", "V99", "V100", "V101", "V102", "V103", "V104", "V105", "V106", "V107", "V108", "V109", "V110", "V111", "V112", "V113", "V114"), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ScannedPower.txt", sep = "\t", ncolumn = 116)
write.table(cbind(rownameperiodhits, scannedv), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ScannedPower.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
write(c("Power", "GeneID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39", "V40", "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49", "V50", "V51", "V52", "V53", "V54", "V55", "V56", "V57", "V58", "V59", "V60", "V61", "V62", "V63", "V64", "V65", "V66", "V67", "V68", "V69", "V70", "V71", "V72", "V73", "V74", "V75", "V76", "V77", "V78", "V79", "V80", "V81", "V82", "V83", "V84", "V85", "V86", "V87", "V88", "V89", "V90", "V91", "V92", "V93", "V94", "V95", "V96", "V97", "V98", "V99", "V100", "V101", "V102", "V103", "V104", "V105", "V106", "V107", "V108", "V109", "V110", "V111", "V112", "V113", "V114"), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ScannedPower.txt", sep = "\t", ncolumn = 116, append = TRUE)
write.table(cbind(rownameperiodhits, powerv), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ScannedPower.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)

# Output Significant Hits with FDR-Adjusted P Value < 0.05.
n <- length(new_p)
writecontent1 <- cbind(rownameperiodhits[new_p < 0.05], pv[(1:n)[new_p < 0.05],], new_p[new_p < 0.05], siglevelv[new_p < 0.05], peakv[new_p < 0.05], peakatv[(1:n)[new_p < 0.05],], periodhits[(1:n)[new_p < 0.05],])
writecontent2 <- cbind(rownameperiodhits[new_p < 0.05], scannedv[(1:n)[new_p < 0.05],])
writecontent3 <- cbind(rownameperiodhits[new_p < 0.05], powerv[(1:n)[new_p < 0.05],])
writecontent4 <- cbind(rownameperiodhits[new_p < 0.05], pv[(1:n)[new_p < 0.05],], new_p[new_p < 0.05], siglevelv[new_p < 0.05], peakv[new_p < 0.05], peakatv[(1:n)[new_p < 0.05],], x_resamplev[(1:n)[new_p < 0.05],])
writecontent5 <- cbind(rownameperiodhits[new_p < 0.05], omegav[(1:n)[new_p < 0.05],])
writecontent6 <- cbind(rownameperiodhits[new_p < 0.05], tauv[(1:n)[new_p < 0.05],])
write(c("GeneID", "P-Value", "FDR-P", "Sig_Level", "Peak", "Peak_At_Frequency", "Peak_At_Period", colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_Hits.txt", sep = "\t", ncolumn = 64)
write.table(writecontent1, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_Hits.txt", sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
write(c("Scanned", "GeneID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39", "V40", "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49", "V50", "V51", "V52", "V53", "V54", "V55", "V56", "V57", "V58", "V59", "V60", "V61", "V62", "V63", "V64", "V65", "V66", "V67", "V68", "V69", "V70", "V71", "V72", "V73", "V74", "V75", "V76", "V77", "V78", "V79", "V80", "V81", "V82", "V83", "V84", "V85", "V86", "V87", "V88", "V89", "V90", "V91", "V92", "V93", "V94", "V95", "V96", "V97", "V98", "V99", "V100", "V101", "V102", "V103", "V104", "V105", "V106", "V107", "V108", "V109", "V110", "V111", "V112", "V113", "V114"), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_ResampledHits_ScannedPower.txt", sep = "\t", ncolumn = 116)
write.table(writecontent2, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_ResampledHits_ScannedPower.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
write(c("Power", "GeneID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39", "V40", "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49", "V50", "V51", "V52", "V53", "V54", "V55", "V56", "V57", "V58", "V59", "V60", "V61", "V62", "V63", "V64", "V65", "V66", "V67", "V68", "V69", "V70", "V71", "V72", "V73", "V74", "V75", "V76", "V77", "V78", "V79", "V80", "V81", "V82", "V83", "V84", "V85", "V86", "V87", "V88", "V89", "V90", "V91", "V92", "V93", "V94", "V95", "V96", "V97", "V98", "V99", "V100", "V101", "V102", "V103", "V104", "V105", "V106", "V107", "V108", "V109", "V110", "V111", "V112", "V113", "V114"), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_ResampledHits_ScannedPower.txt", sep = "\t", ncolumn = 116, append = TRUE)
write.table(writecontent3, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_ResampledHits_ScannedPower.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
write(c("GeneID", "P-Value", "FDR-P", "Sig_Level", "Peak", "Peak_At_Frequency", "Peak_At_Period", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39", "V40", "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49", "V50", "V51", "V52", "V53", "V54", "V55", "V56", "V57"), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_ResampledHits.txt", sep = "\t", ncolumn = 64)
write.table(writecontent4, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_ResampledHits.txt", sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
write(c("Omega", "GeneID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39", "V40", "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49", "V50", "V51", "V52", "V53", "V54", "V55", "V56", "V57", "V58", "V59", "V60", "V61", "V62", "V63", "V64", "V65", "V66", "V67", "V68", "V69", "V70", "V71", "V72", "V73", "V74", "V75", "V76", "V77", "V78", "V79", "V80", "V81", "V82", "V83", "V84", "V85", "V86", "V87", "V88", "V89", "V90", "V91", "V92", "V93", "V94", "V95", "V96", "V97", "V98", "V99", "V100", "V101", "V102", "V103", "V104", "V105", "V106", "V107", "V108", "V109", "V110", "V111", "V112", "V113", "V114"), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_ResampledHits_OmegaTau.txt", sep = "\t", ncolumn = 116)
write.table(writecontent5, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_ResampledHits_OmegaTau.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
write(c("Tau", "GeneID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30", "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38", "V39", "V40", "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48", "V49", "V50", "V51", "V52", "V53", "V54", "V55", "V56", "V57", "V58", "V59", "V60", "V61", "V62", "V63", "V64", "V65", "V66", "V67", "V68", "V69", "V70", "V71", "V72", "V73", "V74", "V75", "V76", "V77", "V78", "V79", "V80", "V81", "V82", "V83", "V84", "V85", "V86", "V87", "V88", "V89", "V90", "V91", "V92", "V93", "V94", "V95", "V96", "V97", "V98", "V99", "V100", "V101", "V102", "V103", "V104", "V105", "V106", "V107", "V108", "V109", "V110", "V111", "V112", "V113", "V114"), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_ResampledHits_OmegaTau.txt", sep = "\t", ncolumn = 116, append = TRUE)
write.table(writecontent6, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_ResampledHits_OmegaTau.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)


############################################################
# C.1.2.1. HEATMAP OF ORIGINAL SIGNIFICANT HITS (FDR-ADJUSTED P-VALUE < 0.05).
############################################################
# Heatmap of periodic hits with lsp_mod that has FDR-Adjusted p < 0.05.
periodhits <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_Hits.txt", header = TRUE, sep = "\t")
periodhits <- periodhits[,c(1,8:length(periodhits[1,]))]
periodhitsnames <- periodhits[,1]
periodhits <- as.matrix(periodhits[,2:length(periodhits[1,])])
row.names(periodhits) <- periodhitsnames
colnames(periodhits) <- c("D_0", "D_4", "D_21", "D_116", "D_185", "D_186", "D_255", "D_289", "D_290", "D_292", "D_294", "D_297", "D_301", "D_307", "D_311", "D_322", "D_329", "D_369", "D_380", "D_400", "D_476", "D_532", "D_546", "D_602", "D_615", "D_616", "D_618", "D_620", "D_625", "D_630", "D_647", "D_679", "D_680", "D_683", "D_688", "D_694", "D_700", "D_711", "D_735", "D_796", "D_840", "D_912", "D_944", "D_945", "D_948", "D_959", "D_966", "D_984", "D_1029", "D_1030", "D_1032", "D_1038", "D_1045", "D_1051", "D_1060", "D_1109", "D_1124")
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Use hclust for cluster and replace dendrogram. 
hc <- hclust(dist(periodhits, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(periodhits, Rowv=as.dendrogram(hc), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# hc$height
plot(hc$height)
quantile(hc$height, probs = seq(0,1,0.01))
         0%          1%          2%          3%          4%          5%          6%          7%          8%          9%         10%         11%         12% 
  0.0000000   0.8890005   1.0738963   1.1077720   1.1729788   1.2351045   1.2797613   1.2996658   1.3437737   1.4302297   1.4882566   1.5291895   1.5813209 
        13%         14%         15%         16%         17%         18%         19%         20%         21%         22%         23%         24%         25% 
  1.6017506   1.6464111   1.6931171   1.7210090   1.7812732   1.8074904   1.8264726   1.8519941   1.9229393   1.9413383   1.9690534   1.9885255   2.0120460 
        26%         27%         28%         29%         30%         31%         32%         33%         34%         35%         36%         37%         38% 
  2.0551732   2.0828706   2.0980229   2.1361910   2.1617509   2.1807166   2.2218912   2.2528314   2.2599846   2.2979158   2.3343077   2.3871822   2.4219917 
        39%         40%         41%         42%         43%         44%         45%         46%         47%         48%         49%         50%         51% 
  2.4530526   2.4780213   2.5022230   2.5223491   2.5528212   2.5742521   2.5944136   2.6509546   2.6861003   2.7234310   2.7501806   2.7621804   2.7932137 
        52%         53%         54%         55%         56%         57%         58%         59%         60%         61%         62%         63%         64% 
  2.8178375   2.8320941   2.8929812   2.9149022   2.9231514   2.9810759   3.0409399   3.1214957   3.1370153   3.1722773   3.2060617   3.2207203   3.2645856 
        65%         66%         67%         68%         69%         70%         71%         72%         73%         74%         75%         76%         77% 
  3.3137778   3.3829306   3.4485141   3.4969007   3.5454946   3.5826995   3.6214608   3.7010051   3.7685930   3.8848887   3.9259221   4.0209329   4.2052399 
        78%         79%         80%         81%         82%         83%         84%         85%         86%         87%         88%         89%         90% 
  4.2500077   4.3688547   4.5223715   4.5943883   4.6457772   4.7230788   4.8942191   5.0122250   5.2336736   5.4852752   5.6804950   5.8612003   6.2401256 
        91%         92%         93%         94%         95%         96%         97%         98%         99%        100% 
  6.5390095   7.6616149   8.0942498   8.5670650  10.0737438  11.5897825  14.6669928  21.5102503  35.9873038 277.9685109 
clusters <- cutree(hc, h=200)
write(c(names(clusters), clusters), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_Hits_2clusters.txt", sep = "\t", ncolumn = 288)
# Include Cluster Info in Heatmap.
# The input file combines the cluster info with the "UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_Hits.txt" file.
periodhitscluster <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_1stCol2Cluster.txt",  header = TRUE, sep = "\t", row.names = 1)
total_matrix_cluster <- data.matrix(periodhitscluster[,c(1,8:length(periodhitscluster[1,]))])
row.names(total_matrix_cluster) <- row.names(periodhitscluster)
colnames(total_matrix_cluster) <- c("CLUSTER", "D_0", "D_4", "D_21", "D_116", "D_185", "D_186", "D_255", "D_289", "D_290", "D_292", "D_294", "D_297", "D_301", "D_307", "D_311", "D_322", "D_329", "D_369", "D_380", "D_400", "D_476", "D_532", "D_546", "D_602", "D_615", "D_616", "D_618", "D_620", "D_625", "D_630", "D_647", "D_679", "D_680", "D_683", "D_688", "D_694", "D_700", "D_711", "D_735", "D_796", "D_840", "D_912", "D_944", "D_945", "D_948", "D_959", "D_966", "D_984", "D_1029", "D_1030", "D_1032", "D_1038", "D_1045", "D_1051", "D_1060", "D_1109", "D_1124")
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
total_heatmaphc <- heatmap.2(total_matrix_cluster, Rowv=as.dendrogram(hc), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
write.matrix(total_heatmaphc$carpet, "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_1stCol2Cluster_HeatmapCarpet.txt", sep = "\t")


############################################################
# C.1.2.2. LINE PLOTS OF ORIGINAL SIGNIFICANT HITS (FDR-ADJUSTED P-VALUE < 0.05).
############################################################
# Line plots of periodic hits with lsp_mod that has p < 0.01.
# All.
plot(periodhits[1,]~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", ylim = c(-1.2,1.2), xlim = c(-50, 1200))
for (i in 2:length(periodhits[,1])) {
	lines(periodhits[i,]~dayseries, type = "o", pch = as.integer(25*i/86), col = i)
}
# Clusters. Plot both Median Path and Smoothed Median Path by Moving Average.
# Cluster 1.
cluster1 <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_1stCol2Cluster_Cl1.txt",  header = TRUE, sep = "\t", row.names = 1)
clusterforlineplot <- data.matrix(cluster1[8:length(cluster1[1,])])
# Plot scatterSmooth Plots.
dayseriesforsmoothScatter <- rep(dayseries, length(clusterforlineplot[,1]))
clusterforlineplotforsmoothScatter <- as.vector(t(clusterforlineplot))
# Turquoise.
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "turquoise")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
points(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, pch = 20, cex = 0.1)
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "turquoise")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
for (i in 1:length(clusterforlineplot[,1])) {
	lines(clusterforlineplot[i,]~dayseries, type = "o", pch = ".", col = "black", lwd = 0.1, add = TRUE)
}
# Green.
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "green")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
points(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, pch = 20, cex = 0.1)
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "green")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
for (i in 1:length(clusterforlineplot[,1])) {
	lines(clusterforlineplot[i,]~dayseries, type = "o", pch = ".", col = "black", lwd = 0.1, add = TRUE)
}
# Plot Line Plots.
medianpath <- apply(clusterforlineplot, 2, median) ### Median Path for plotting.
# plot(medianpath~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path")
mincurve <- apply(clusterforlineplot, 2, min)
maxcurve <- apply(clusterforlineplot, 2, max)
plot(medianpath~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path", ylim = c(-1.2,1.2), xlim = c(-50, 1200))
polygon(c(dayseries, rev(dayseries)), c(maxcurve, rev(mincurve)), col = "lightyellow", border = NA) # Candidate Colors: lightyellow, khaki1, lightskyblue1
lines(mincurve~dayseries, col = "turquoise")
lines(maxcurve~dayseries, col = "gold")
lines(medianpath~dayseries,col = "gray12")
mv <- function(x,n=5){filter(x,rep(1/n,n), sides=2)} ### Moving Average Function with every 5 points.
medianpathmv <- as.numeric(mv(medianpath, n=5))
medianpathmv[1] <- medianpath[1]
medianpathmv[2] <- medianpath[2]
medianpathmv[(length(medianpathmv)-1):length(medianpathmv)] <- medianpath[(length(medianpath)-1):length(medianpath)]
mincurvemv <- as.numeric(mv(mincurve, n=5)) # Smoothed Median Path for Plotting.
mincurvemv[1] <- mincurve[1]
mincurvemv[2] <- mincurve[2]
mincurvemv[(length(mincurvemv)-1):length(mincurvemv)] <- mincurve[(length(mincurve)-1):length(mincurve)]
maxcurvemv <- as.numeric(mv(maxcurve, n=5)) # Smoothed Median Path for Plotting.
maxcurvemv[1] <- maxcurve[1]
maxcurvemv[2] <- maxcurve[2]
maxcurvemv[(length(maxcurvemv)-1):length(maxcurvemv)] <- maxcurve[(length(maxcurve)-1):length(maxcurve)]
plot(medianpathmv~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path -- Moving Average", ylim = c(-1.2,1.2), xlim = c(-50, 1200))
polygon(c(dayseries, rev(dayseries)), c(maxcurvemv, rev(mincurvemv)), col = "lightyellow", border = NA) # Candidate Colors: lightyellow, khaki1, lightskyblue1
lines(mincurvemv~dayseries, col = "turquoise")
lines(maxcurvemv~dayseries, col = "gold")
lines(medianpathmv~dayseries,col = "gray12")
# Cluster 2.
cluster2 <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_1stCol2Cluster_Cl2.txt",  header = TRUE, sep = "\t", row.names = 1)
clusterforlineplot <- data.matrix(cluster2[8:length(cluster2[1,])])
# Plot scatterSmooth Plots.
dayseriesforsmoothScatter <- rep(dayseries, length(clusterforlineplot[,1]))
clusterforlineplotforsmoothScatter <- as.vector(t(clusterforlineplot))
# Gold.
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "gold")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
points(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, pch = 20, cex = 0.1)
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "gold")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
for (i in 1:length(clusterforlineplot[,1])) {
	lines(clusterforlineplot[i,]~dayseries, type = "o", pch = ".", col = "black", lwd = 0.1, add = TRUE)
}
# Red.
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "red")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
points(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, pch = 20, cex = 0.1)
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "red")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
for (i in 1:length(clusterforlineplot[,1])) {
	lines(clusterforlineplot[i,]~dayseries, type = "o", pch = ".", col = "black", lwd = 0.1, add = TRUE)
}
# Plot Line Plots.
medianpath <- apply(clusterforlineplot, 2, median) ### Median Path for plotting.
# plot(medianpath~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path")
mincurve <- apply(clusterforlineplot, 2, min)
maxcurve <- apply(clusterforlineplot, 2, max)
plot(medianpath~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path", ylim = c(-1.2,1.2), xlim = c(-50, 1200))
polygon(c(dayseries, rev(dayseries)), c(maxcurve, rev(mincurve)), col = "lightyellow", border = NA) # Candidate Colors: lightyellow, khaki1, lightskyblue1
lines(mincurve~dayseries, col = "turquoise")
lines(maxcurve~dayseries, col = "gold")
lines(medianpath~dayseries,col = "gray12")
mv <- function(x,n=5){filter(x,rep(1/n,n), sides=2)} ### Moving Average Function with every 5 points.
medianpathmv <- as.numeric(mv(medianpath, n=5))
medianpathmv[1] <- medianpath[1]
medianpathmv[2] <- medianpath[2]
medianpathmv[(length(medianpathmv)-1):length(medianpathmv)] <- medianpath[(length(medianpath)-1):length(medianpath)]
mincurvemv <- as.numeric(mv(mincurve, n=5)) # Smoothed Median Path for Plotting.
mincurvemv[1] <- mincurve[1]
mincurvemv[2] <- mincurve[2]
mincurvemv[(length(mincurvemv)-1):length(mincurvemv)] <- mincurve[(length(mincurve)-1):length(mincurve)]
maxcurvemv <- as.numeric(mv(maxcurve, n=5)) # Smoothed Median Path for Plotting.
maxcurvemv[1] <- maxcurve[1]
maxcurvemv[2] <- maxcurve[2]
maxcurvemv[(length(maxcurvemv)-1):length(maxcurvemv)] <- maxcurve[(length(maxcurve)-1):length(maxcurve)]
plot(medianpathmv~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path -- Moving Average", ylim = c(-1.2,1.2), xlim = c(-50, 1200))
polygon(c(dayseries, rev(dayseries)), c(maxcurvemv, rev(mincurvemv)), col = "lightyellow", border = NA) # Candidate Colors: lightyellow, khaki1, lightskyblue1
lines(mincurvemv~dayseries, col = "turquoise")
lines(maxcurvemv~dayseries, col = "gold")
lines(medianpathmv~dayseries,col = "gray12")


############################################################
# C.1.2.3. HEATMAP OF RESAMPLED SIGNIFICANT HITS (FDR-ADJUSTED P-VALUE < 0.05).
############################################################
# Heatmap of Resampled hits that has p < 0.01.
periodhits <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_FDRp_0_05_ResampledHits.txt", header = TRUE, sep = "\t")
periodhits <- periodhits[,c(1,8:length(periodhits[1,]))]
periodhitsnames <- periodhits[,1]
periodhits <- as.matrix(periodhits[,2:length(periodhits[1,])])
row.names(periodhits) <- periodhitsnames
for (i in 1:length(colnames(periodhits))) {
	colnames(periodhits)[i] <- paste("RS", colnames(periodhits)[i], sep = "_")
}
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Range of periodhits: -8.390401086 ~ 7.572316183
# Adjusted range to make difference obvious.
total_heatmap <- heatmap.2(periodhits, Rowv=TRUE, Colv=FALSE, col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-9,-3,length=167), seq(-3,3,length=167), seq(3,9,length=167)), keysize = 0.5)


############################################################
# C.2. ALL ANALYTES (USE LSP_MOD).
############################################################
############################################################
# C.2.1. LOAD ALL LOMB-SCARGLE PARAMETERS AND RESAMPLED SAMPLES (USE LSP_MOD).
############################################################
# Load All Analytes (Resampled and Original).
periodhits <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled.txt", header = TRUE, sep = "\t")
periodhits <- periodhits[,c(1,8:length(periodhits[1,]))]
periodhitsnames <- periodhits[,1]
periodhits <- as.matrix(periodhits[,2:length(periodhits[1,])])
row.names(periodhits) <- periodhitsnames
for (i in 1:length(colnames(periodhits))) {
	colnames(periodhits)[i] <- paste("RS", colnames(periodhits)[i], sep = "_")
}
originalhits <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL.txt", header = TRUE, sep = "\t")
originalhits <- originalhits[,c(1,8:length(originalhits[1,]))]
originalhitsnames <- originalhits[,1]
originalhits <- as.matrix(originalhits[,2:length(originalhits[1,])])
row.names(originalhits) <- originalhitsnames
# for (i in 1:length(colnames(originalhits))) {
	# colnames(originalhits)[i] <- paste("D", colnames(originalhits)[i], sep = "_")
# }


############################################################
# C.2.2. AUTOCORRELATION CALCULATION (USE GEORGE'S FORMULA).
############################################################
############################################################
# C.2.2.1. CALCULATE AUTOCORRELATION COEFFICIENT FOR ALL ANALYTES.
############################################################
# ACF of Resampled Analytes:
autocorrelation <- c()
# Lag.
k <- 1
for (i in 1:length(periodhits[,1])) {
	rowmean <- mean(periodhits[i,])
	numerator <- 0
	denominator <- 0
	for (j in 1:(length(periodhits[1,])-k)) {
		numerator <- numerator + (periodhits[i,j] - rowmean) * (periodhits[i, j + k] - rowmean)
	}
	for (j in 1:length(periodhits[1,])) {
		denominator <- denominator + (periodhits[i,j] - rowmean)^2
	}
	rho <- numerator / denominator
	autocorrelation <- c(autocorrelation, rho)
}
# Obtain list of analytes with acf > mean + 2sigma or acf < mean - 2sigma.
meanautocorrelation <- mean(autocorrelation)
meanautocorrelation
[1] 0.06519466
sigma <- sd(autocorrelation)
sigma
[1] 0.141892
hist(autocorrelation, nclass = 200)
# mean + 1.645 sigma for p < 0.05 one-tailed (right tail).
meanautocorrelation + 1.645 * sigma
[1] 0.2986069


############################################################
# C.2.2.2. BOOTSTRAP OF AUTOCORRELATION COEFFICIENT.
############################################################
# Bootstrap through Randomizing.
# sample(x, size, replace = TRUE, prob = NULL)
autocorrelationrand <- c()
# Lag.
k <- 1
set.seed(7)
for (i in 1:1000000) {
	periodhitsrand <- sample(periodhits, 57, replace = TRUE, prob = NULL)
		rowmeanrand <- mean(periodhitsrand)
		numeratorrand <- 0
		denominatorrand <- 0
		for (j in 1:(length(periodhitsrand)-k)) {
			numeratorrand <- numeratorrand + (periodhitsrand[j] - rowmeanrand) * (periodhitsrand[j + k] - rowmeanrand)
		}
		for (j in 1:length(periodhitsrand)) {
			denominatorrand <- denominatorrand + (periodhitsrand[j] - rowmeanrand)^2
		}
		rhorand <- numeratorrand / denominatorrand
		autocorrelationrand <- c(autocorrelationrand, rhorand)
}
# Save Bootstrap ACF Distribution.
write(autocorrelationrand, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ACF_Bootstrap_Distribution_wRplc.txt", ncol = 1)
# Obtain list of analytes with acf > mean + 2sigma or acf < mean - 2sigma.
meanautocorrelationrand <- mean(autocorrelationrand)
meanautocorrelationrand
[1] -0.01758803
sigmarand <- sd(autocorrelationrand)
sigmarand
[1] 0.1212896
hist(autocorrelationrand, nclass = 200)
# mean + 1.645 sigma for p < 0.05 one-tailed (right tail).
meanautocorrelationrand + 1.645 * sigmarand
[1] 0.1819333

# Output ACF Resampled Dataset.
periodallacf <- cbind(row.names(periodhits), autocorrelation, periodhits)
write.table(periodallacf, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ACF_Dataset.txt", sep = "\t", row.names = FALSE)
# Output ACF Original Dataset.
originalallacf <- cbind(row.names(originalhits), autocorrelation, originalhits)
write.table(originalallacf, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Dataset.txt", sep = "\t", row.names = FALSE)


############################################################
# C.2.2.3. SELECTING SIGNIFICANT ACF WITH REAL DISTRIBUTION MEAN AND SIGMA.
############################################################
# Slow step in R, used perl for selection instead.
# Selecting significant ACF.
# mean + 1.645 sigma for p < 0.05 one-tailed (right tail).
periodallacf <- c()
originalallacf <- c()
periodhitsacf <- c()
originalhitsacf <- c()
for (i in 1:length(periodhits[,1])) {
	tmprow0a <- c(row.names(periodhits)[i], autocorrelation[i], periodhits[i,])
	periodallacf <- rbind(periodallacf, tmprow0a)		
	tmprow0b <- c(row.names(originalhits)[i], autocorrelation[i], originalhits[i,])
	originalallacf <- rbind(originalallacf, tmprow0b)		
	if (autocorrelation[i] > meanautocorrelation + 1.645 * sigma) {
		tmprow1 <- c(row.names(periodhits)[i], autocorrelation[i], periodhits[i,])
		periodhitsacf <- rbind(periodhitsacf, tmprow1)		
		tmprow2 <- c(row.names(originalhits)[i], autocorrelation[i], originalhits[i,])
		originalhitsacf <- rbind(originalhitsacf, tmprow2)		
	}
}
colnames(periodallacf) <- c("GeneID", "ACF", colnames(periodhits))
colnames(originalallacf) <- c("GeneID", "ACF", colnames(originalhits))
colnames(periodhitsacf) <- c("GeneID", "ACF", colnames(periodhits))
colnames(originalhitsacf) <- c("GeneID", "ACF", colnames(originalhits))
# Output ACF Resampled Dataset.
write.table(periodallacf, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ACF_Dataset.txt", sep = "\t", row.names = FALSE)
# Output ACF Original Dataset.
write.table(originalallacf, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Dataset.txt", sep = "\t", row.names = FALSE)
# Output ACF Resampled Hits.
write.table(periodhitsacf, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ACF_Dataset_Hits.txt", sep = "\t", row.names = FALSE)
# Output ACF Original Hits.
write.table(originalhitsacf, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Dataset_Hits.txt", sep = "\t", row.names = FALSE)

# Heatmap of ACF Hits.
# Resampled Hits.
periodhitshm <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ACF_Dataset_Hits.txt", header = TRUE, sep = "\t")
periodhitshm <- periodhitshm[,c(1,3:length(periodhitshm[1,]))]
periodhitshmnames <- periodhitshm[,1]
periodhitshm <- as.matrix(periodhitshm[,2:length(periodhitshm[1,])])
row.names(periodhitshm) <- periodhitshmnames
# Range of periodhitshm: -4.802538 - 5.335815
# Adjusted range to make difference obvious.
hc_ts <- hclust(dist(periodhitshm, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(periodhitshm, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-6,-1,length=167), seq(-1,2,length=167), seq(2,6,length=167)), keysize = 0.5)
# Original Hits.
originalhitshm <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Dataset_Hits.txt", header = TRUE, sep = "\t")
originalhitshm <- originalhitshm[,c(1,3:length(originalhitshm[1,]))]
originalhitshmnames <- originalhitshm[,1]
originalhitshm <- as.matrix(originalhitshm[,2:length(originalhitshm[1,])])
row.names(originalhitshm) <- originalhitshmnames
hc_ts <- hclust(dist(originalhitshm, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(originalhitshm, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)


############################################################
# C.2.2.4. SELECTING SIGNIFICANT ACF WITH BOOTSTRAP DISTRIBUTION MEAN AND SIGMA.
############################################################
# Slow step in R, used perl for selection instead.
# Selecting significant ACF.
# mean + 1.645 sigma for p < 0.05 one-tailed (right tail).
periodallacfbs <- c()
originalallacfbs <- c()
periodhitsacfbs <- c()
originalhitsacfbs <- c()
for (i in 1:length(periodhits[,1])) {
	tmprow0abs <- c(row.names(periodhits)[i], autocorrelation[i], periodhits[i,])
	periodallacfbs <- rbind(periodallacfbs, tmprow0abs)		
	tmprow0bbs <- c(row.names(originalhits)[i], autocorrelation[i], originalhits[i,])
	originalallacfbs <- rbind(originalallacfbs, tmprow0bbs)		
	if (autocorrelation[i] > meanautocorrelationrand + 1.645 * sigmarand) {
		tmprow1bs <- c(row.names(periodhits)[i], autocorrelation[i], periodhits[i,])
		periodhitsacfbs <- rbind(periodhitsacfbs, tmprow1bs)		
		tmprow2bs <- c(row.names(originalhits)[i], autocorrelation[i], originalhits[i,])
		originalhitsacfbs <- rbind(originalhitsacfbs, tmprow2bs)		
	}
}
colnames(periodallacfbs) <- c("GeneID", "ACF", colnames(periodhits))
colnames(originalallacfbs) <- c("GeneID", "ACF", colnames(originalhits))
colnames(periodhitsacfbs) <- c("GeneID", "ACF", colnames(periodhits))
colnames(originalhitsacfbs) <- c("GeneID", "ACF", colnames(originalhits))
# Output ACF Resampled Dataset.
write.table(periodallacfbs, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ACF_Bootstrap.txt", sep = "\t", row.names = FALSE)
# Output ACF Original Dataset.
write.table(originalallacfbs, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Bootstrap.txt", sep = "\t", row.names = FALSE)
# Output ACF Resampled Hits.
write.table(periodhitsacfbs, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ACF_Bootstrap_Hits.txt", sep = "\t", row.names = FALSE)
# Output ACF Original Hits.
write.table(originalhitsacfbs, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Bootstrap_Hits.txt", sep = "\t", row.names = FALSE)

# Heatmap of ACF Hits.
# Resampled Hits.
periodhitshmbs <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Resampled_ACF_Bootstrap_Hits.txt", header = TRUE, sep = "\t")
periodhitshmbs <- periodhitshmbs[,c(1,3:length(periodhitshmbs[1,]))]
periodhitshmnamesbs <- periodhitshmbs[,1]
periodhitshmbs <- as.matrix(periodhitshmbs[,2:length(periodhitshmbs[1,])])
row.names(periodhitshmbs) <- periodhitshmnamesbs
# Range of periodhitshmbs: -6.70043 ~ 6.30206
# Adjusted range to make difference obvious.
hc_ts <- hclust(dist(periodhitshmbs, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(periodhitshmbs, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-7,-1,length=167), seq(-1,1,length=167), seq(1,7,length=167)), keysize = 0.5)
# Original Hits.
originalhitshmbs <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Bootstrap_Hits.txt", header = TRUE, sep = "\t")
originalhitshmbs <- originalhitshmbs[,c(1,3:length(originalhitshmbs[1,]))]
originalhitshmnamesbs <- originalhitshmbs[,1]
originalhitshmbs <- as.matrix(originalhitshmbs[,2:length(originalhitshmbs[1,])])
row.names(originalhitshmbs) <- originalhitshmnamesbs
hc_ts <- hclust(dist(originalhitshmbs, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(originalhitshmbs, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
plot(hc_ts$height)
quantile(hc_ts$height, probs = seq(0,1,0.01))
         0%          1%          2%          3%          4%          5%          6%          7%          8%          9%         10%         11%         12% 
   0.000000    1.066448    1.170718    1.276293    1.310528    1.380548    1.442516    1.496939    1.542332    1.588791    1.633763    1.702591    1.731065 
        13%         14%         15%         16%         17%         18%         19%         20%         21%         22%         23%         24%         25% 
   1.784788    1.810592    1.855035    1.890873    1.923753    1.957450    1.986096    2.010870    2.054602    2.089751    2.121387    2.147289    2.178883 
        26%         27%         28%         29%         30%         31%         32%         33%         34%         35%         36%         37%         38% 
   2.208078    2.231077    2.254416    2.283032    2.318481    2.342289    2.380583    2.400323    2.421691    2.448549    2.489595    2.522482    2.545839 
        39%         40%         41%         42%         43%         44%         45%         46%         47%         48%         49%         50%         51% 
   2.570794    2.598927    2.633287    2.657081    2.690804    2.716033    2.745997    2.777244    2.792390    2.820231    2.851089    2.870641    2.910998 
        52%         53%         54%         55%         56%         57%         58%         59%         60%         61%         62%         63%         64% 
   2.942981    2.962207    2.985609    3.017857    3.050680    3.084296    3.119575    3.160905    3.187973    3.223861    3.263068    3.302372    3.351338 
        65%         66%         67%         68%         69%         70%         71%         72%         73%         74%         75%         76%         77% 
   3.395579    3.442661    3.490396    3.557379    3.596331    3.665627    3.719248    3.768054    3.812272    3.911265    3.984049    4.041799    4.118170 
        78%         79%         80%         81%         82%         83%         84%         85%         86%         87%         88%         89%         90% 
   4.244643    4.312765    4.436952    4.575969    4.647425    4.753664    4.857794    5.011537    5.143194    5.326639    5.515692    5.873330    6.104582 
        91%         92%         93%         94%         95%         96%         97%         98%         99%        100% 
   6.512698    7.090220    7.656851    8.166138    9.551777   10.864165   13.023530   18.681675   31.796283 1098.254074 
clusters <- cutree(hc_ts, h=800)
write(c(names(clusters), clusters), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Bootstrap_Hits_2clusters.txt", sep = "\t", ncolumn = 2055)
# Include Cluster Info in Heatmap.
# The input file combines the cluster info with the "UNK_RNASeq_hg19_QNRest_LSPMOD_ALL_Original_ACF_Bootstrap_Hits.txt" file.
periodhitscluster <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Bootstrap_Hits_1stCol2Cluster.txt",  header = TRUE, sep = "\t", row.names = 1)
total_matrix_cluster <- data.matrix(periodhitscluster[,c(1,3:length(periodhitscluster[1,]))])
row.names(total_matrix_cluster) <- row.names(periodhitscluster)
colnames(total_matrix_cluster) <- c("CLUSTER", "D_0", "D_4", "D_21", "D_116", "D_185", "D_186", "D_255", "D_289", "D_290", "D_292", "D_294", "D_297", "D_301", "D_307", "D_311", "D_322", "D_329", "D_369", "D_380", "D_400", "D_476", "D_532", "D_546", "D_602", "D_615", "D_616", "D_618", "D_620", "D_625", "D_630", "D_647", "D_679", "D_680", "D_683", "D_688", "D_694", "D_700", "D_711", "D_735", "D_796", "D_840", "D_912", "D_944", "D_945", "D_948", "D_959", "D_966", "D_984", "D_1029", "D_1030", "D_1032", "D_1038", "D_1045", "D_1051", "D_1060", "D_1109", "D_1124")
hc_ts <- hclust(dist(total_matrix_cluster[,2:length(total_matrix_cluster[1,])], method = "euclidean"), method = "ward.D")
total_heatmaphc <- heatmap.2(total_matrix_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
write.matrix(total_heatmaphc$carpet, "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Bootstrap_Hits_1stCol2Cluster_HeatmapCarpet.txt", sep = "\t")


############################################################
# C.2.2.5. LINE PLOTS OF ORIGINAL HITS WITH SIGNIFICANT ACF WITH BOOTSTRAP DISTRIBUTION MEAN AND SIGMA.
############################################################
# Line plots of ACF hits with ACF bootstrap mean and sigma.
# All.
plot(originalhitshmbs[1,]~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", ylim = c(-1.2,1.2), xlim = c(-50, 1200))
for (i in 2:length(originalhitshmbs[,1])) {
	lines(originalhitshmbs[i,]~dayseries, type = "o", pch = as.integer(25*i/86), col = i)
}
# Clusters. Plot both Median Path and Smoothed Median Path by Moving Average.
# Cluster 1.
cluster1 <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Bootstrap_Hits_1stCol2Cluster_Cl1.txt",  header = TRUE, sep = "\t", row.names = 1)
clusterforlineplot <- data.matrix(cluster1[3:length(cluster1[1,])])
# Plot scatterSmooth Plots.
dayseriesforsmoothScatter <- rep(dayseries, length(clusterforlineplot[,1]))
clusterforlineplotforsmoothScatter <- as.vector(t(clusterforlineplot))
# Turquoise.
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "turquoise")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
points(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, pch = 20, cex = 0.1)
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "turquoise")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
for (i in 1:length(clusterforlineplot[,1])) {
	lines(clusterforlineplot[i,]~dayseries, type = "o", pch = ".", col = "black", lwd = 0.1, add = TRUE)
}
# Green.
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "green")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
points(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, pch = 20, cex = 0.1)
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "green")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
for (i in 1:length(clusterforlineplot[,1])) {
	lines(clusterforlineplot[i,]~dayseries, type = "o", pch = ".", col = "black", lwd = 0.1, add = TRUE)
}
# Plot Line Plots.
medianpath <- apply(clusterforlineplot, 2, median) ### Median Path for plotting.
# plot(medianpath~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path")
mincurve <- apply(clusterforlineplot, 2, min)
maxcurve <- apply(clusterforlineplot, 2, max)
plot(medianpath~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path", ylim = c(-1.2,1.2), xlim = c(-50, 1200))
polygon(c(dayseries, rev(dayseries)), c(maxcurve, rev(mincurve)), col = "lightyellow", border = NA) # Candidate Colors: lightyellow, khaki1, lightskyblue1
lines(mincurve~dayseries, col = "turquoise")
lines(maxcurve~dayseries, col = "gold")
lines(medianpath~dayseries,col = "gray12")
mv <- function(x,n=5){filter(x,rep(1/n,n), sides=2)} ### Moving Average Function with every 5 points.
medianpathmv <- as.numeric(mv(medianpath, n=5))
medianpathmv[1] <- medianpath[1]
medianpathmv[2] <- medianpath[2]
medianpathmv[(length(medianpathmv)-1):length(medianpathmv)] <- medianpath[(length(medianpath)-1):length(medianpath)]
mincurvemv <- as.numeric(mv(mincurve, n=5)) # Smoothed Median Path for Plotting.
mincurvemv[1] <- mincurve[1]
mincurvemv[2] <- mincurve[2]
mincurvemv[(length(mincurvemv)-1):length(mincurvemv)] <- mincurve[(length(mincurve)-1):length(mincurve)]
maxcurvemv <- as.numeric(mv(maxcurve, n=5)) # Smoothed Median Path for Plotting.
maxcurvemv[1] <- maxcurve[1]
maxcurvemv[2] <- maxcurve[2]
maxcurvemv[(length(maxcurvemv)-1):length(maxcurvemv)] <- maxcurve[(length(maxcurve)-1):length(maxcurve)]
plot(medianpathmv~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path -- Moving Average", ylim = c(-1.2,1.2), xlim = c(-50, 1200))
polygon(c(dayseries, rev(dayseries)), c(maxcurvemv, rev(mincurvemv)), col = "lightyellow", border = NA) # Candidate Colors: lightyellow, khaki1, lightskyblue1
lines(mincurvemv~dayseries, col = "turquoise")
lines(maxcurvemv~dayseries, col = "gold")
lines(medianpathmv~dayseries,col = "gray12")
# Cluster 2.
cluster2 <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Bootstrap_Hits_1stCol2Cluster_Cl2.txt",  header = TRUE, sep = "\t", row.names = 1)
clusterforlineplot <- data.matrix(cluster2[3:length(cluster2[1,])])
# Plot scatterSmooth Plots.
dayseriesforsmoothScatter <- rep(dayseries, length(clusterforlineplot[,1]))
clusterforlineplotforsmoothScatter <- as.vector(t(clusterforlineplot))
# Gold.
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "gold")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
points(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, pch = 20, cex = 0.1)
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "gold")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
for (i in 1:length(clusterforlineplot[,1])) {
	lines(clusterforlineplot[i,]~dayseries, type = "o", pch = ".", col = "black", lwd = 0.1, add = TRUE)
}
# Red.
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "red")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
points(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, pch = 20, cex = 0.1)
smoothScatter(dayseriesforsmoothScatter, clusterforlineplotforsmoothScatter, nbin = 512, nrpoints = 0, colramp=colorRampPalette(c("white", "white", "red")), xlab="Day Number", ylab="RU", ylim=c(-1.2,1.2), xlim = c(-50, 1200))
for (i in 1:length(clusterforlineplot[,1])) {
	lines(clusterforlineplot[i,]~dayseries, type = "o", pch = ".", col = "black", lwd = 0.1, add = TRUE)
}
# Plot Line Plots.
medianpath <- apply(clusterforlineplot, 2, median) ### Median Path for plotting.
# plot(medianpath~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path")
mincurve <- apply(clusterforlineplot, 2, min)
maxcurve <- apply(clusterforlineplot, 2, max)
plot(medianpath~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path", ylim = c(-1.2,1.2), xlim = c(-50, 1200))
polygon(c(dayseries, rev(dayseries)), c(maxcurve, rev(mincurve)), col = "lightyellow", border = NA) # Candidate Colors: lightyellow, khaki1, lightskyblue1
lines(mincurve~dayseries, col = "turquoise")
lines(maxcurve~dayseries, col = "gold")
lines(medianpath~dayseries,col = "gray12")
mv <- function(x,n=5){filter(x,rep(1/n,n), sides=2)} ### Moving Average Function with every 5 points.
medianpathmv <- as.numeric(mv(medianpath, n=5))
medianpathmv[1] <- medianpath[1]
medianpathmv[2] <- medianpath[2]
medianpathmv[(length(medianpathmv)-1):length(medianpathmv)] <- medianpath[(length(medianpath)-1):length(medianpath)]
mincurvemv <- as.numeric(mv(mincurve, n=5)) # Smoothed Median Path for Plotting.
mincurvemv[1] <- mincurve[1]
mincurvemv[2] <- mincurve[2]
mincurvemv[(length(mincurvemv)-1):length(mincurvemv)] <- mincurve[(length(mincurve)-1):length(mincurve)]
maxcurvemv <- as.numeric(mv(maxcurve, n=5)) # Smoothed Median Path for Plotting.
maxcurvemv[1] <- maxcurve[1]
maxcurvemv[2] <- maxcurve[2]
maxcurvemv[(length(maxcurvemv)-1):length(maxcurvemv)] <- maxcurve[(length(maxcurve)-1):length(maxcurve)]
plot(medianpathmv~dayseries,type = "o", pch = ".", xlab = "Day Number", ylab = "RU", main = "Median Path -- Moving Average", ylim = c(-1.2,1.2), xlim = c(-50, 1200))
polygon(c(dayseries, rev(dayseries)), c(maxcurvemv, rev(mincurvemv)), col = "lightyellow", border = NA) # Candidate Colors: lightyellow, khaki1, lightskyblue1
lines(mincurvemv~dayseries, col = "turquoise")
lines(maxcurvemv~dayseries, col = "gold")
lines(medianpathmv~dayseries,col = "gray12")


############################################################
# C.2.2.6. ACF HISTOGRAM OF ORIGINAL HITS WITH SIGNIFICANT ACF WITH BOOTSTRAP DISTRIBUTION MEAN AND SIGMA.
############################################################
meanautocorrelationrand + 1.645 * sigmarand
[1] 0.1819333
# Bootstrap Distribution.
# No Y Limit.
hist(autocorrelationrand, nclass = 200, col = rgb(0,0,1,1/2))
# Real Dataset Distribution.
hist(autocorrelation, nclass = 200, col = rgb(1,0,0,1/4), add = TRUE)
# Hits Distribution.
acfhits <- c()
for (i in 1:length(autocorrelation)) {
	if (autocorrelation[i] > meanautocorrelationrand + 1.645 * sigmarand) {
		acfhits <- c(acfhits, autocorrelation[i])
	}
}
hist(acfhits, nclass = 100, col = rgb(0,1,0,1/4), add = TRUE)
abline(v = meanautocorrelationrand + 1.645 * sigmarand, col = "red", lty = 2, lwd = 2)

# With Y Limit.
hist(autocorrelationrand, nclass = 200, col = rgb(0,0,1,1/2), ylim = c(0, 150))
# Real Dataset Distribution.
hist(autocorrelation, nclass = 200, col = rgb(1,0,0,1/4), add = TRUE)
# Hits Distribution.
acfhits <- c()
for (i in 1:length(autocorrelation)) {
	if (autocorrelation[i] > meanautocorrelationrand + 1.645 * sigmarand) {
		acfhits <- c(acfhits, autocorrelation[i])
	}
}
hist(acfhits, nclass = 100, col = rgb(0,1,0,1/4), add = TRUE)
abline(v = meanautocorrelationrand + 1.645 * sigmarand, col = "red", lty = 2, lwd = 2)




