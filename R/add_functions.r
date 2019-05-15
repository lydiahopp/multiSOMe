r <- unclass(lsf.str(envir = asNamespace("oposSOM"), all = T))

for(name in r) eval(parse(text=paste0(name, '<-oposSOM:::', name)))


summarySheetsModules <- function()
{
  environment(modules.CSV.sheets) <- environment()
  environment(modules.report.sheets) <- environment()
  environment(modules.profiles) <- environment()
  environment(modules.chromosomes) <- environment()


  dirname <- paste(files.name, "- Results/Summary Sheets - Modules")
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  #### Overexpression Spots ####

  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Overexpression Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Chromosomes.pdf") )


  #### Underexpression Spots ####

  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Underexpression Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.underexpression, main="Underexpression Spots", path=file.path(dirname,"Report.pdf") )


  #### Correlation Cluster ####

  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Correlation Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.correlation, main="Correlation Cluster", path=file.path(dirname,"Report.pdf") )


  #### K-Means Cluster ####

  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","K-Means Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Chromosomes.pdf") )


  #### D-Clusters ####

  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","D-Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Chromosomes.pdf") )


  #### Group Overexpression Spots ####

  if (length(unique(group.labels)) > 1)
  {
    dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Group Overexpression Spots" )
    util.info("Writing:", file.path(dirname, "*.pdf"))
    dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

    modules.report.sheets(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Report.pdf") )
    modules.profiles(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Profiles.pdf") )
    modules.chromosomes(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Chromosomes.pdf") )
  }


  #### module gene lists CSV sheets ####

  dirname <- file.path(output.paths["CSV"], "Spot Lists")
  util.info("Writing:", file.path(dirname, "*.csv"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.CSV.sheets(spot.list=spot.list.overexpression, main="Overexpression Spots", path=dirname )
  modules.CSV.sheets(spot.list=spot.list.kmeans, main="K-Means Cluster", path=dirname )
  modules.CSV.sheets(spot.list=spot.list.dmap, main="D-Cluster", path=dirname )

  if (length(unique(group.labels)) > 1)
  {
    modules.CSV.sheets(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=dirname )
  }

}

summarySheetsModules2 <- function()
{
  environment(modules.CSV.sheets) <- environment()
  environment(modules.report.sheets) <- environment()
  environment(modules.profiles) <- environment()
  environment(modules.chromosomes) <- environment()


  dirname <- paste(files.name, "- Results/Summary Sheets - Modules",files.name)
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  #### Overexpression Spots ####

  dirname <- file.path( paste(files.name, "- Results/Summary Sheets - Modules",files.name),"Overexpression Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Chromosomes.pdf") )


  #### Underexpression Spots ####

  dirname <- file.path( paste(files.name, "- Results/Summary Sheets - Modules",files.name),"Underexpression Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.underexpression, main="Underexpression Spots", path=file.path(dirname,"Report.pdf") )


  #### Correlation Cluster ####

  dirname <- file.path( paste(files.name, "- Results/Summary Sheets - Modules",files.name),"Correlation Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.correlation, main="Correlation Cluster", path=file.path(dirname,"Report.pdf") )


  #### K-Means Cluster ####

  dirname <- file.path( paste(files.name, "- Results/Summary Sheets - Modules",files.name),"K-Means Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Chromosomes.pdf") )


  #### D-Clusters ####

  dirname <- file.path( paste(files.name, "- Results/Summary Sheets - Modules",files.name),"D-Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Chromosomes.pdf") )


  #### Group Overexpression Spots ####

  if (length(unique(group.labels)) > 1)
  {
    dirname <- file.path( paste(files.name, "- Results/Summary Sheets - Modules",files.name),"Group Overexpression Spots" )
    util.info("Writing:", file.path(dirname, "*.pdf"))
    dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

    modules.report.sheets(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Report.pdf") )
    modules.profiles(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Profiles.pdf") )
    modules.chromosomes(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Chromosomes.pdf") )
  }


  #### module gene lists CSV sheets ####

  dirname <- file.path(output.paths["CSV"], "Spot Lists")
  util.info("Writing:", file.path(dirname, "*.csv"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.CSV.sheets(spot.list=spot.list.overexpression, main="Overexpression Spots", path=dirname )
  modules.CSV.sheets(spot.list=spot.list.kmeans, main="K-Means Cluster", path=dirname )
  modules.CSV.sheets(spot.list=spot.list.dmap, main="D-Cluster", path=dirname )

  if (length(unique(group.labels)) > 1)
  {
    modules.CSV.sheets(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=dirname )
  }

}

topologyProfiles <- function()
{
  metadata.scaled <- apply(metadata, 2, function(x) { (x-min(x)) / (max(x)-min(x)) })

  filename <- file.path(paste(files.name, "- Results"), "Supporting Maps&Profiles", "Topology Profiles.pdf")
  util.info("Writing:", filename)

  if (length(unique(group.labels)) > 1)
  {
    pdf(filename, 42/2.54, 21/2.54, useDingbats=FALSE)
    par(mar=c(10, 6, 4, 5), mfrow=c(1, 2))
  } else
  {
    pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)
    par(mar=c(10, 6, 4, 5))
  }

  ### Number of overexpression spots ###
  spotdata <- get(paste("spot.list.",preferences$standard.spot.modules,sep=""))$spotdata
  n.spots <- colSums( spotdata>sd(spotdata) )

  barplot(n.spots, col=group.colors, main="Number of activated spot modules",
          names.arg="", las=2, cex.main=2,
          border=if (ncol(indata) < 80) "black" else NA)
  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes <- by(n.spots, group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, las=2, xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  ### Spot number distribution ###
  n.spots.groups <- sapply(tapply(n.spots, group.labels, c), function(x)
  {
    ret <- table(x)[as.character(0:max(n.spots))]
    ret[which(is.na(ret))] <- 0
    names(ret) <- as.character(0:max(n.spots))
    return(ret)
  })

  if (is.vector(n.spots.groups))
  {
    n.spots.groups <- matrix(n.spots.groups,nrow=1)
    colnames(n.spots.groups) <- unique(group.labels)
  }

  n.spots.groups <- sapply(unique(group.labels), function(x)
  {
    return(n.spots.groups[,x] / table(group.labels)[x])
  })

  par(mfrow=c(1, 1))

  barplot(as.vector(n.spots.groups), col=rep(groupwise.group.colors,each=max(n.spots)+1),
          names.arg=rep(c(0:max(n.spots)), length(unique(group.labels))), las=1,
          main=if(files.name=="only exp") "Fraction of samples showing respective number of overexpression spots" else  "Fraction of samples showing respective number of hypermethylated spots" ,
          cex.main=2, border=if (ncol(indata) < 80) "black" else NA, ylim=c(0, 1))

  par(mfrow=c(1, 2))

  ### Fraction of red metagenes ###
  K.red <- apply(metadata.scaled, 2, function(x) { length(which(x > 0.9)) }) / preferences$dim.1stLvlSom^2

  barplot(K.red, col=group.colors, main=if(files.name=="only exp") "Fraction of overexpressed metagenes" else "Fraction of hypermethylated metagenes", cex.main=2,
          names.arg="", las=2, border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes <- by(K.red, group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, las=2, main="", xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  ### Length of borderline ###
  K.border <- apply(metadata.scaled, 2, function(x)
  {
    m <- matrix(x, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
    m[which(m < 0.9)] <- NA
    count.border.metagenes <- 0

    for (i in 1:preferences$dim.1stLvlSom)
    {
      for (j in 1:preferences$dim.1stLvlSom)
      {
        if (!is.na(m[i,j]))
        {
          neighbours <- sapply(get.neighbors(i, j, preferences$dim.1stLvlSom),
                               function(x) { m[x[1],x[2]] })

          if (length(which(is.na(neighbours))))
          {
            count.border.metagenes <- count.border.metagenes + 1
          } else if (i %in% c(1, preferences$dim.1stLvlSom) || j %in% c(1, preferences$dim.1stLvlSom))
          {
            count.border.metagenes <- count.border.metagenes + 1
          }
        }
      }
    }

    return(count.border.metagenes)
  })

  barplot(K.border, col=group.colors, main= if(files.name=="only exp") "Length of borderline along overexpressed metagenes"else "Length of borderline along hypermethylated metagenes",
          names.arg="", las=2, cex.main=1.8,
          border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(K.border, group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, las=2, main="", xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  ### Compactness of spots ###
  C <- K.red / K.border

  barplot(C, col=group.colors, main="Compactness of spots",
          names.arg="", las=2, cex.main=2,
          border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes <- by(C, group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, las=2, main="", xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  ### Shape of spots ###
  C <- (K.red * preferences$dim.1stLvlSom^2) / K.border^2

  barplot(C, col=group.colors, main="Shape of spots",
          names.arg="", las=2, cex.main=2,
          border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes <- by(C, group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, las=2, main="",xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  dev.off()
}

entropyProfiles <- function()
{

  dirname <- paste(files.name, "- Results/Supporting Maps&Profiles")
  dir.create(dirname, showWarnings=FALSE)

  filename <- file.path(paste(files.name, "- Results"), "Supporting Maps&Profiles", "Entropy Profiles.pdf")
  util.info("Writing:", filename)

  ### Metagene Mean Expresion + Variance ###
  pdf(filename, 42/2.54, 21/2.54, useDingbats=FALSE)
  par(mar=c(10, 6, 4, 5))

  barplot(apply(metadata, 2, mean), col=group.colors,
          main=bquote("Mean metagene expression: <e"[m]^meta~">"), names.arg="", las=2,
          cex.main=2.5,border=if (ncol(indata) < 80) "black" else NA)
  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes <- by(apply(metadata, 2, mean), group.labels, c)[unique(group.labels)]

    boxplot(mean.boxes, col=groupwise.group.colors, las=2,
            main=bquote("Mean metagene expression: <e"[m]^meta~">"), cex.main=2.5, xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  barplot(apply(metadata, 2, var), col=group.colors, main=bquote("Metagene variance: var(e"[m]^meta~")"),
          names.arg="", las=2, cex.main=2.5, border=if (ncol(indata) < 80) "black" else NA)
  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.sd.boxes <- by(apply(metadata, 2, var), group.labels, c)[unique(group.labels)]

    boxplot(mean.sd.boxes, col=groupwise.group.colors, las=2,
            main=bquote("Metagene variance: var(e"[m]^meta~")"), cex.main=2.5, xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  ### Entropy ###
  q25 <- quantile(metadata,0.25)
  q75 <- quantile(metadata,0.75)

  p.metadata <- apply(metadata, 2, function(x)
  {
    hist(x, breaks=c(min(x), q25, q75, max(x)), plot=FALSE)$counts / preferences$dim.1stLvlSom^2
  })

  q25 <- quantile(metadata * som.result$node.summary[,"n.features"],0.25)
  q75 <- quantile(metadata * som.result$node.summary[,"n.features"],0.75)

  p.metadata.weighted <- apply(metadata * som.result$node.summary[,"n.features"], 2, function(x)
  {
    hist(x, breaks=c(min(x), q25, q75, max(x)), plot=FALSE)$counts / preferences$dim.1stLvlSom^2
  })

  ### Standard sample-related metagene entropy
  H <- apply(p.metadata, 2, function(p)
  {
    p <- p[p != 0]
    -sum(p * log2(p))
  })

  ylim <- c(min(H)-0.1*(max(H)-min(H)),max(H)+0.1*(max(H)-min(H)))

  barplot(H, col=group.colors, main=bquote("Standard metagene entropy: h"[m]),
          names.arg="", las=2, cex.main=2.5,
          ylim=ylim, xpd=FALSE, border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.sd.boxes <- by(H, group.labels, c)[unique(group.labels)]

    boxplot(mean.sd.boxes, col=groupwise.group.colors, las=2,
            main=bquote("Standard metagene entropy: h"[m]~""), cex.main=2.5, xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  ### Weighted sample-related metagene entropy
  H <- apply(p.metadata.weighted, 2, function(p)
  {
    p <- p[p != 0]
    -sum(p * log2(p))
  })

  ylim <- c(min(H)-0.1*(max(H)-min(H)),max(H)+0.1*(max(H)-min(H)))

  barplot(H, col=group.colors, main=bquote("Metagene entropy weighted for metagene population: h"[m]^weighted),
          names.arg="", las=2, cex.main=2.5, ylim=ylim,
          xpd=FALSE, border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.sd.boxes <- by(H, group.labels, c)[unique(group.labels)]

    boxplot(mean.sd.boxes, col=groupwise.group.colors, las=2,
            main=bquote("Metagene entropy weighted for metagene population: h"[m]^weighted~""), cex.main=2.5, xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  ### Tsallis Entropy
  q <- 0.25

  H <- apply(p.metadata, 2, function(p)
  {
    (1/(q-1)) * (1-sum(p^q))
  })

  ylim <- c(min(H)-0.1*(max(H)-min(H)),max(H)+0.1*(max(H)-min(H)))

  barplot(H, col=group.colors, main=bquote("Tsallis metagene entropy: h"[m]^Tsallis~""),
          names.arg="", las=2, cex.main=2.5, ylim=ylim,
          xpd=FALSE, border=if (ncol(indata) < 80) "black" else NA)

  box()

  if (length(unique(group.labels)) > 1)
  {
    mean.sd.boxes <- by(H, group.labels, c)[unique(group.labels)]

    boxplot(mean.sd.boxes, col=groupwise.group.colors, las=2,
            main=bquote("Tsallis metagene entropy: h"[m]^Tsallis~""), cex.main=2.5, xaxt="n")

    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2)
  }

  dev.off()
}

