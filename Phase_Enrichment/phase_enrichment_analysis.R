
prep_quantiles <- function(
    x, 
    probs = c(seq(0,0.5,by=0.1),
              seq(0.55,0.90,by=0.05), 
              seq(0.925,0.975,by=0.025),
              seq(0.98, 0.99, by=0.005), 
              0.9925, 0.9950, 0.9975, 0.9980, 0.9985,
              0.9990, 0.99925, 0.99950, 0.99975, 1),
    precision = 5) {
  
  qs <- sort( quantile(round(x * 10^precision)/(10^precision), prob=probs, na.rm=T ) )
  qs <- qs[ !duplicated(qs, fromLast=T) ]
  qs
  
}


coerce_binary <- function(X, na.as.false=T) {
  if (na.as.false) {
    X[is.na(X)] <- "FALSE"
    factor(X, levels=c("FALSE","TRUE"))
  } else {
    factor(X, levels=c("FALSE","TRUE"))
  }
}


phase_enrichment_analysis <- function(
    bTarget,                      # boolean, which entries are targets?
    phase_dat         = all_dat,  # data.frame with phase data
    gd.phase_colname  = "Phase",  # name of column with phase
    bin_seq           = seq(0, 360, by=10),  # bin breakpoints
    bin_offset        = diff(bin_seq)[1]/2,  # adjustment so 1st bin centre aligns with 0
    MTC_method        = "none",
    bBG = !is.na(phase_dat[,gd.phase_colname]) # default BG: all osc entries
  ) {
 
  # define bins & relabel based on centre value
  phase_dat$bins <- cut(
    phase_dat[,gd.phase_colname] - bin_offset, 
    breaks=bin_seq - bin_offset, include.lowest=T)
  levels(phase_dat$bins) <- paste0(apply( stringr::str_split_fixed( 
    stringr::str_replace_all( 
      levels(phase_dat$bins), "[()\\[\\]]", "" ), ",", 2), 1, function(x) {
        mean(as.numeric(x))} ), "°")
  
  # iterate over bins
  res <- list()
  for (bin in levels(phase_dat$bins)) {
    tab <- table( 
      "bTarget" = coerce_binary(bTarget),
      "bInBin"  = coerce_binary(# NB: genes w/o bins will NOT be discarded here
        !is.na(phase_dat$bins) & phase_dat$bins == bin)
    )
    # print(tab)
    ft  <- fisher.test(tab, alt="two.sided")
    res[[bin]] <- c(
      "nTargetsInBin" = tab["TRUE","TRUE"],
      "OR"       = unname(ft$estimate),
      "FT_pval"  = unname(ft$p.value)
    )
  }
  res <- as.data.frame(do.call(rbind, res))
  res <- tibble::rownames_to_column(res, "bins")
  res$bins <- factor(res$bins, levels=levels(phase_dat$bins))
  
  # derived values: adj p-value and two-sided signif level
  res$FT_pval <- p.adjust(res$FT_pval, method=MTC_method)
  res$signif <- symnum(
    res$FT_pval * ifelse(res$OR>1,1,-1), 
    cutpoints = c(-1, -0.1, -0.05, -0.01, -0.001, 0, 0.001, 0.01, 0.05, 0.1, 1), 
    symbols = LETTERS[1:10])
  res$signif <- c("A"=" ", "B"=".", "C"="-", "D"="--", "E"="---",
                  "F"="+++", "G"="++", "H"="+", "I"=".", "J"=" ")[
    res$signif ]
  res$signif <- factor(
    res$signif,
    levels = c(" ", ".", "-", "--", "---", "+", "++", "+++"))
  
  # return
  return(res)
  
   
}


plot_PEA_barplot <- function(
    res, drop_unused_signif=T, drop_borderline_signif=T) {

  signif_pal <- c("white", "grey80", "lightblue", "blue", "darkblue", "pink", "red", "darkred")
  names(signif_pal) <- levels(res$signif)
  if (drop_borderline_signif) {
    res$signif[res$signif=="."] <- " "
  }
  if (drop_unused_signif==T) {
    res$signif <- droplevels(res$signif)
    signif_pal <- signif_pal[levels(res$signif)]
  }
  
  ggplot(res, aes_string(x="bins", y="OR")) +
    geom_bar(aes_string(fill="signif"), stat="identity") + 
    scale_fill_manual(values=signif_pal) +
    geom_hline(yintercept = 1, col="black", linetype="dashed") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}


plot_PEA_radar <- function(
    res, maxy = max(ceiling(res$OR), na.rm=T), 
    drop_unused_signif=T, drop_borderline_signif=T) {
  
  bin_width <- 2 * pi / nrow(res)
  
  signif_pal <- c("white", "grey80", "lightblue", "blue", "darkblue", "pink", "red", "darkred")
  names(signif_pal) <- levels(res$signif)
  if (drop_borderline_signif) {
    res$signif[res$signif=="."] <- " "
  }
  if (drop_unused_signif==T) {
    res$signif <- droplevels(res$signif)
    signif_pal <- signif_pal[levels(res$signif)]
  }
  
  ggplot(res, aes_string(x="bins", y="OR")) +
    
    # convert to polar
    coord_polar(theta="x", start=-5*bin_width, direction = -1) +
    scale_y_continuous(
      limits = c(-0.5, maxy),
      breaks = seq(0, maxy, by=1),
      minor_breaks = seq(0.5, 2.5, by=1)
    ) +
    
    geom_col(aes_string(fill="signif"), color="darkgrey") + 
    scale_fill_manual(values=signif_pal) +
    
    geom_hline(yintercept = 1, col="black", linetype="dashed") +
    
    # Annotate custom scale inside plot
    annotate(
      x = nrow(res) / 4 + 1, 
      y = 1:maxy-0.11, 
      label = 1:maxy, 
      geom = "text", color = "gray12"
    ) +
    # add "axis" labels
    annotate(
      x = nrow(res) / 4 + 1,  y = maxy-0.5, 
      label = c("odds ratio"), 
      geom = "text", color = "gray12", size=5
    ) +
    annotate(
      x = 1,  y = maxy-0.5, 
      label = "Phase", 
      geom = "text", color = "gray12", size=5
    ) +
    theme(
      # Remove axis ticks and text
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      # Use gray text for the region names
      axis.text.x = element_text(color = "gray12", size = 10),
      # Move the legend to the bottom
      legend.position = "bottom",
      # background colour and gridlines
      panel.background = element_rect(
        fill = "white", color = "white"),
      panel.grid.major = element_line(
        size = 0.5, linetype = 'solid', colour = "grey"), 
      panel.grid.minor = element_line(
        size = 0.25, linetype = 'solid', colour = "lightgrey")
    )
  
}



phase_densities <- function(
    
  target_gIDs,
  
  oscDat = AllOsc_info,
  geneDat = as.data.frame(geneBody_gr),
  
  flag_gIDs = NULL, 
  adjust = 0.25,
  
  do_plot = T,
  
  ...
  
) {
  
  # establish phases of genes w/ peaks and w/o
  target_phases <- na.omit(oscDat[
    rownames(oscDat) %in% target_gIDs,"Phase"])
  nontarget_phases <- na.omit(oscDat[
    ! rownames(oscDat) %in% target_gIDs,"Phase"])
  
  # helper: circular density
  tri_density <- function(x, width=360, adjust=adjust, n=512) {
    xw <- c(x-width, x, x+width)
    d <- stats::density(xw, adjust=adjust, n=n, 
                        from=0-width, to=2*width)
    ok <- (d$x>=0 & d$x<=width)
    d$x <- d$x[ok]
    d$y <- d$y[ok]
    d$n <- d$n / 3
    d
  }
  
  
  # helper: plot density polygon
  dens_poly <- function(
    d,
    poly.border=NA, poly.fill="red", 
    line.w=2, line.col="darkred") {
    
    d2 <- d
    d2$x <- c(min(d2$x), d2$x, max(d2$x))
    d2$y <- c(0, d2$y, 0)
    polygon(d2, border=poly.border, col=poly.fill)
    lines(d, lwd=line.w, col=line.col)
    
  }
  
  
  # establish densities
  d_T  <- tri_density(target_phases, adjust=adjust)
  d_nT <- tri_density(nontarget_phases, adjust=adjust)
  
  
  # plot
  if (do_plot) {
    ylim <- range(0, d_T$y, d_nT$y)
    
    plot( c(0, 360), ylim, type="n", bty="n", 
          xlab="Phase (degrees)", ylab="density", main="" )
    abline(h=0, col="grey10")
    dens_poly(d_nT, poly.fill = alpha("grey", 0.5), line.col = "darkgrey")
    dens_poly(d_T,  poly.fill = alpha("pink", 0.5), line.col = "red")
    
    
    # flag some genes
    if (!is.null(flag_gIDs)) {
      if (any(! flag_gIDs %in% rownames(oscDat))) {
        warning(paste0(setdiff(flag_gIDs, rownames(oscDat)), collapse=", "), 
                " were not oscillating(s).\n")
      }
      flag_gIDs <- intersect(flag_gIDs, rownames(oscDat))
      if (any(! flag_gIDs %in% target_gIDs)) {
        warning(paste0(setdiff(flag_gIDs, target_gIDs), collapse=", "), 
                " were not target(s).\n")
      }
      flag_gIDs <- intersect(flag_gIDs, target_gIDs)
      if (length(flag_gIDs)>0) {
        x <- oscDat[ flag_gIDs,"Phase"]
        flags <- geneDat[flag_gIDs,"gene_name"]
        abline(v=x, col="purple", lty="dashed")
        text(x, rep(ylim[2], length(x)), flags, col="purple", 
             srt=90, adj=c(1.2, -0.2))
      }
    }
    
    
    # legend
    legend(
      "bottomright", inset=0.05,
      leg=c(
        paste0("target (", d_T$n,")"), 
        paste0("not target (", d_nT$n,")")),
      bg="white", pch=c(22, 22), 
      text.col = c("red","darkgrey"),
      pt.bg=c(alpha("pink", 0.5), alpha("grey", 0.5)),
      col=c("red","darkgrey"))
  }
  
  
  # return
  invisible(list(d_T=d_T, d_nT=d_nT))
  
}



phase_densities_polar <- function(
    obj,
    spoke_ext = 1.05,
    spoke_col = "darkgrey",
    circle_col = "darkgrey",
    plot_ext  = 1.075 ) {
  
  
  
  polar2cart <- function(df){
    data.frame(
      x = df$y * cos(df$x * 2 * pi / 360) ,
      y = df$y * sin(df$x * 2 * pi / 360)
    )
  }
  
  dat4gg <- rbind(
    cbind(polar2cart(obj$d_T),  "gene" = "target osc"),
    cbind(polar2cart(obj$d_nT), "gene" = "non-target osc")
  )
  
  # ggplot(dat4gg, aes_string(x="x", y="y")) +
  #   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  #   geom_path(aes_string(colour="gene"), size=2)
  
  add_circle <- function(
    r=1, x=0, y=0, nv=100, 
    border="black", col=NA,
    lty="solid", lwd=1,
    do_label=T) {
    # inspired by plotrix::draw.circle, but removing the link to aspect ratio
    
    angle.inc <- 2 * pi/nv
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    for (i in 1:length(r)) {
      xv <- cos(angles) * r[i] + x
      yv <- sin(angles) * r[i] + y
      polygon(xv, yv, border = border, col = col, lty = lty, lwd = lwd)
      if (do_label) {
        xv <- cos(pi/2) * r[i] + x
        yv <- sin(pi/2) * r[i] + y
        text(xv, yv, label=r[i], pos=1)
      }
    }
  }
  
  add_spokes <- function(
    r=1, x=0, y=0, angle.inc=pi/4,
    col="black", lty="solid", lwd=1,
    do_label = T
    ) {
    
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    for (i in 1:length(angles)) {
      xv <- cos(angles) * r + x
      yv <- sin(angles) * r + y
      segments(
        x0 = x, y0 = y, x1 = xv, y1 = yv, 
        col = col, lty = lty, lwd = lwd)
      if (do_label) {
        xv <- cos(angles + angle.inc / 10) * r + x
        yv <- sin(angles + angle.inc / 10) * r + y
        text( xv, yv, label = paste0(round(360 * angles / (2 * pi)),"°") )
      }
    }
  }
  
  grid_circles <- setdiff( pretty(
    abs(c(obj$d_T$y, obj$d_nT$y, dat4gg$x, dat4gg$y))), 0)
  m <- max(grid_circles) 
  
  plot(
    c(-m, m) * plot_ext, c(-m, m) * plot_ext, 
    type="n", main="", sub="", xlab="", ylab="",
    xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n" ) 
  add_spokes(m * spoke_ext, col=spoke_col)
  add_circle(grid_circles, border=circle_col)
  
  pal <- c(
    "non-target osc" = "darkgrey",
    "target osc"     = "darkblue")
  ns <- c(
    "non-target osc" = obj$d_nT$n,
    "target osc"     = obj$d_T$n)
  
  for (lev in names(pal)) {
    polygon(
      dat4gg[ dat4gg$gene==lev,c("x","y")],
      col=alpha(pal[lev],0.5), border=pal[lev])
  }
  
  legend("bottomright", 
         leg=paste0(names(ns), " (",ns,")"), 
         pch=22, pt.bg=alpha(pal),
         text.col=pal, cex = 0.85)
  
  invisible(dat4gg)
}


# EOF
