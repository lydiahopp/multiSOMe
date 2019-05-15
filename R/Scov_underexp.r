ScoV.underexp=function()
{
  environment(modules.CSV.sheets) <- environment()
  environment(modules.report.sheets) <- environment()
  environment(modules.profiles) <- environment()
  environment(modules.chromosomes) <- environment()
  environment(modules.relations) <- environment()

  group.metadata <- do.call(cbind, by(t(metadata), group.labels, colMeans))[,unique(group.labels)]
  group.metadata.scaled <- apply(group.metadata, 2, function(x)   (x-min(x)) / (max(x)-min(x))   )

  ## extract sample modules ##
  sample.spot.list <- list()
  sample.spot.core.list <- list()
  n.sample.modules <- 0

  for (m in 1:ncol(group.metadata))
  {
    # define bigger core regions
    core <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
    core[which(group.metadata.scaled[,m] < 1- preferences$spot.threshold.groupmap)] <- -1

    spot.i <- 0

    while (nrow(which(core == -1, arr.ind=TRUE)) > 0)
    {
      start.pix <- which(core == -1, arr.ind=TRUE)[1,]
      spot.i <- spot.i + 1
      core <- col.pix(core, start.pix[1], start.pix[2], spot.i, preferences$dim.1stLvlSom)
    }

    # shrink each separate region to core size
    for (s.i in 1:max(core,na.rm=TRUE))
    {
      if (sum(core == s.i, na.rm=TRUE) > preferences$spot.coresize.groupmap)
      {
        core.metagenes <- which(core == s.i)
        o <- order(group.metadata[core.metagenes,m], decreasing=FALSE)[1:preferences$spot.coresize.groupmap]
        core[setdiff(core.metagenes, core.metagenes[o])] <- NA
      }
    }

    core[which(!is.na(core))] <- -1
    spot.i <- 0

    while (nrow(which(core == -1, arr.ind=TRUE)) > 0)
    {
      start.pix <- which(core == -1, arr.ind=TRUE)[1,]
      spot.i <- spot.i + 1
      core <- col.pix(core, start.pix[1], start.pix[2], spot.i, preferences$dim.1stLvlSom)
    }

    # define spot area around cores
    for (s.i in 1:max(core,na.rm=TRUE))
    {
      n.sample.modules <- n.sample.modules + 1

      sample.spot.core.list[[n.sample.modules]] <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
      sample.spot.core.list[[n.sample.modules]][which(core == s.i)] <- metadata[which(core == s.i),m]

      spot <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
      spot[which(group.metadata.scaled[,m] < 1 - preferences$spot.threshold.groupmap)] <- -1

      start.pix <- which(!is.na(sample.spot.core.list[[n.sample.modules]]), arr.ind=TRUE)
      start.pix <- start.pix[which.max(sample.spot.core.list[[n.sample.modules]][start.pix]),]

      spot <- col.pix(spot, start.pix[1], start.pix[2], 1, preferences$dim.1stLvlSom)

      sample.spot.list[[n.sample.modules]] <- matrix(NA, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom)
      sample.spot.list[[n.sample.modules]][which(spot == 1)] <- group.metadata.scaled[which(spot == 1),m]
    }
  }

  # filter
  remove <- c()

  for (i in 1:n.sample.modules)
  {
    if (sum(!is.na(sample.spot.list[[i]])) <= 1)
    {
      # empty, i.e. does not exceed threshold -> remove
      remove <- c(remove, i)
    } else if (sum(!is.na(sample.spot.list[[i]])) < sum(!is.na(sample.spot.core.list[[i]])))
    {
      # core larger than spot -> shrink core
      sample.spot.core.list[[i]] <- sample.spot.list[[i]]
    }
  }

  if (length(remove) > 0)
  {
    sample.spot.list <- sample.spot.list[-remove]
    sample.spot.core.list <- sample.spot.core.list[-remove]
    n.sample.modules <- length(sample.spot.list)
  }

  o <- order(sapply(sample.spot.core.list,function(x) mean(x,na.rm=TRUE)), decreasing=FALSE)
  sample.spot.list <- sample.spot.list[o]
  sample.spot.core.list <- sample.spot.core.list[o]

  ## merge overlapping sample cores ##
  merged <- TRUE

  while (merged)
  {
    merged <- FALSE
    i <- 1

    while (i < length(sample.spot.list))
    {
      j <- i + 1

      while (j <= length(sample.spot.list))
      {
        core1 <- which(!is.na(sample.spot.core.list[[i]]))
        core2 <- which(!is.na(sample.spot.core.list[[j]]))

        if (any(core1 %in% core2))
        {
          # merge cores
          if (length(setdiff(core1,core2)) > 0)
          {
            sample.spot.core.list[[i]][setdiff(core2,core1)] <-
              sample.spot.core.list[[j]][setdiff(core2,core1)]
          }

          # merge spots
          spot1 <- which(!is.na(sample.spot.list[[i]]))
          spot2 <- which(!is.na(sample.spot.list[[j]]))

          if (length(setdiff(spot2,spot1)) > 0)
          {
            sample.spot.list[[i]][setdiff(spot2,spot1)] <-
              sample.spot.list[[j]][setdiff(spot2,spot1)]
          }

          # remove j
          sample.spot.list <- sample.spot.list[-j]
          sample.spot.core.list <- sample.spot.core.list[-j]

          merged <- TRUE
        } else
        {
          j <- j + 1
        }
      }

      i <- i + 1
    }
  }

  o <- order(sapply(sample.spot.core.list, function(x) { mean(x, na.rm=TRUE) }), decreasing=FALSE)
  sample.spot.list <- sample.spot.list[o]
  sample.spot.core.list <- sample.spot.core.list[o]

  ## shrinking overlapping spots ##
  if (length(sample.spot.list) > 1)
  {
    for (i in seq_len(length(sample.spot.list)-1))
    {
      for (j in seq(i+1, length(sample.spot.list)))
      {
        spot1 <- which(!is.na(sample.spot.list[[i]]))
        spot2 <- which(!is.na(sample.spot.list[[j]]))

        if (any(spot1 %in% spot2))
        {
          spot12intersect <- which(spot1 %in% spot2)

          spot1 <- which(!is.na(sample.spot.list[[i]]), arr.ind=TRUE)
          spot2 <- which(!is.na(sample.spot.list[[j]]), arr.ind=TRUE)

          spot12intersect <- spot1[spot12intersect, ,drop=FALSE]

          core1.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[i]]), arr.ind=TRUE), 2, range))
          core2.center <- colMeans(apply(which(!is.na(sample.spot.core.list[[j]]), arr.ind=TRUE), 2, range))

          spot12assoc <- apply(spot12intersect, 1, function(x)
          {
            which.min(c(sum((core1.center - x) ^ 2), sum((core2.center - x) ^ 2)))
          })

          sample.spot.list[[j]][spot12intersect[which(spot12assoc==1),,drop=FALSE]] <- NA
          sample.spot.list[[i]][spot12intersect[which(spot12assoc==2),,drop=FALSE]] <- NA
        }
      }
    }
  }

  ## define overexpression spots ##


  spot.metagenes = lapply( sample.spot.list,function(x){which(!is.na(x))  }     )

  positions= sapply(spot.metagenes,function(y){colMeans(apply(som.result$node.summary[y, 1:2]+1, 2, range)) })
  o=order(positions[1,]*positions[2,],decreasing=T)
  sample.spot.list=sample.spot.list[o]
  spot.metagenes = lapply( sample.spot.list,function(x){which(!is.na(x))  }     )


  spot.list.group.underexpression = list()

  spot.list.group.underexpression$overview.map =
    apply(apply(group.metadata, 2, function(x){ (x - min(x)) / (max(x) - min(x)) }), 1, min)

  spot.list.group.underexpression$overview.mask = rep(NA, preferences$dim.1stLvlSom ^ 2)
  spot.list.group.underexpression$spots = list()


  for (i in seq_along(sample.spot.list))
  {
    spot.metagenes <- which(!is.na(sample.spot.list[[i]]))
    spot.genes <- rownames(indata)[which(som.result$feature.BMU %in% spot.metagenes)]

    if (length(spot.genes) > 0)
    {
      spot.list.group.underexpression$overview.mask[which(!is.na(sample.spot.list[[i]]))] = i
      spot.list.group.underexpression$spots[[LETTERS[i]]] = list()
      spot.list.group.underexpression$spots[[LETTERS[i]]]$metagenes = spot.metagenes
      spot.list.group.underexpression$spots[[LETTERS[i]]]$genes = spot.genes
      spot.list.group.underexpression$spots[[LETTERS[i]]]$mask = rep(NA, preferences$dim.1stLvlSom * preferences$dim.1stLvlSom)
      spot.list.group.underexpression$spots[[LETTERS[i]]]$mask[spot.metagenes] = 1

      spot.list.group.underexpression$spots[[LETTERS[i]]]$position =
        colMeans(apply(som.result$node.summary[spot.metagenes, 1:2]+1, 2, range))

      spot.list.group.underexpression$spots[[LETTERS[i]]]$beta.statistic =
        get.beta.statistic(set.data=metadata[spot.list.group.underexpression$spots[[LETTERS[i]]]$metagenes,,drop=FALSE],
                           weights=som.result$node.summary[spot.list.group.underexpression$spots[[LETTERS[i]]]$metagenes,]$nobs)
    }
  }

  spot.list.group.underexpression$spotdata =
    t(sapply(spot.list.group.underexpression$spots, function(x)
    {
      if (length(x$genes > 0))
      {
        colMeans(indata[x$genes,,drop=FALSE])
      } else
      {
        rep(0, ncol(indata))
      }
    }))

  colnames(spot.list.group.underexpression$spotdata) = colnames(indata)



  spot.fisher.p <- function(spot)
  {
    spot$Fisher.p <- GeneSet.Fisher(unique(na.omit(gene.info$ids[spot$genes])),
                                    unique.protein.ids, gs.def.list, sort=TRUE)

    return(spot)
  }

  spot.list.group.underexpression$spots =lapply( spot.list.group.underexpression$spots, spot.fisher.p)

  #spot.list.group.underexpression2=spot.list.group.underexpression

  #rm(list=ls()[-which(ls()=="spot.list.group.underexpression2")])
  #detach(env)
  # load("ScoV.RData")
  #env$spot.list.group.underexpression=NULL
  #env$spot.list.group.underexpression=spot.list.group.underexpression2
  # save(env, file="ScoV.RData")

  # rm(list=ls())
  # load("ScoV.RData")
  #attach(env)
  #library(oposSOM)

  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Group neg Correlation Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.group.underexpression, main="Group neg Correlation Spots", path=file.path(dirname,"Report.pdf") )
  modules.profiles(spot.list=spot.list.group.underexpression,main="Group neg Correlation Spots", path=file.path(dirname,"Profiles.pdf") )
  modules.chromosomes(spot.list=spot.list.group.underexpression, main="Group neg Correlation Spots", path=file.path(dirname,"Chromosomes.pdf") )

}
