# Creates a new opossom environment
opossom.new <- function(preferences=NULL)
{
  # Init the environment
  env <- new.env()
  env$color.palette.portraits <- NULL
  env$color.palette.heatmaps <- NULL
  env$t.ensID.m <- NULL
  env$Fdr.g.m <- NULL
  env$fdr.g.m <- NULL
  env$files.name <- NULL
  env$gene.info <- NULL
  env$chromosome.list <- NULL
  env$group.silhouette.coef <- NULL
  env$group.colors <- NULL
  env$group.labels <- NULL
  env$gs.def.list <- NULL
  env$samples.GSZ.scores <- NULL
  env$spot.list.correlation <- NULL
  env$spot.list.dmap <- NULL
  env$spot.list.group.overexpression <- NULL
  env$spot.list.kmeans <- NULL
  env$spot.list.overexpression <- NULL
  env$spot.list.samples <- NULL
  env$spot.list.underexpression <- NULL
  env$indata <- NULL
  env$indata.gene.mean <- NULL
  env$indata.sample.mean <- NULL
  env$metadata <- NULL
  env$n.0.m <- NULL
  env$output.paths <- NULL
  env$pat.labels <- NULL
  env$p.g.m <- NULL
  env$p.m <- NULL
  env$perc.DE.m <- NULL
  env$som.result <- NULL
  env$t.g.m <- NULL
  env$t.m <- NULL
  env$groupwise.group.colors <- NULL
  env$unique.protein.ids <- NULL
  env$WAD.g.m <- NULL

  # Generate some additional letters
  env$LETTERS <- c(LETTERS, as.vector(sapply(1:10, function(x) {
    return(paste(LETTERS, x, sep=""))
  })))

  env$letters <- c(letters, as.vector(sapply(1:10, function(x) {
    return(paste(letters, x, sep=""))
  })))

  # Set default preferences
  env$preferences <- list(dataset.name = "Unnamed",
                          dim.1stLvlSom = "auto",
                          dim.2ndLvlSom = 20,
                          training.extension = 1,
                          rotate.SOM.portraits = 0,
                          flip.SOM.portraits = FALSE,
                          activated.modules = list( "reporting" = TRUE,
                                                    "primary.analysis" = TRUE,
                                                    "sample.similarity.analysis" = TRUE,
                                                    "geneset.analysis" = TRUE,
                                                    "geneset.analysis.exact" = FALSE,
                                                    "group.analysis" = TRUE,
                                                    "difference.analysis" = TRUE ),
                          database.biomart = "ENSEMBL_MART_ENSEMBL",
                          database.host = "www.ensembl.org",
                          database.dataset = "auto",
                          database.id.type = "",
                          standard.spot.modules = "dmap",
                          spot.coresize.modules = 3,
                          spot.threshold.modules = 0.95,
                          spot.coresize.groupmap = 5,
                          spot.threshold.groupmap = 0.75,
                          adjust.autogroup.number = 0,
                          feature.centralization = TRUE,
                          sample.quantile.normalization = TRUE,
                          pairwise.comparison.list = NULL)

  # Merge user supplied information
  if (!is.null(preferences))
  {
    env$preferences <-
      modifyList(env$preferences, preferences[names(env$preferences)])
  }
  if(!is.null(preferences$indata))
  {
    env$indata <- preferences$indata
  }
  if(!is.null(preferences$group.labels))
  {
    env$group.labels <- preferences$group.labels
  }
  if(!is.null(preferences$group.colors))
  {
    env$group.colors <- preferences$group.colors
  }

  return(env)
}

# Executes the oposSOM pipeline.
multiSOMe.run <- function(env)
{
  util.info("Started:", env$preferences$started)
  util.info("Name:", env$preferences$dataset.name)

  #### Preparation & Calculation part ####

  if (!util.call(oposSOM:::pipeline.checkInputParameters, env)) {
    return()
  }

  if(env$preferences$activated.modules$primary.analysis)
  {
    env$preferences$system.info <- Sys.info()
    env$preferences$session.info <- sessionInfo()
    env$preferences$started <- format(Sys.time(), "%a %d %b %Y %X")
  }

  if(env$preferences$activated.modules$reporting)
  {
    # create output dirs
    dir.create(paste(env$files.name, "- Results"), showWarnings=FALSE)
    #    dir.create(paste(env$files.name, "- Results/CSV Sheets"), showWarnings=FALSE)

    # if(env$preferences$activated.modules$primary.analysis)
    #  {
    #    util.call(oposSOM:::pipeline.qualityCheck, env)
    #  }
  }

  if(env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
  {
    util.info("Loading gene annotation data. This may take several minutes until next notification.")
    biomart.available <- oposSOM:::biomart.available
    util.call(oposSOM:::pipeline.prepareAnnotation, env)
  }

  if(env$preferences$activated.modules$primary.analysis)
  {
    util.info("Processing SOM. This may take several time until next notification.")
    util.call(oposSOM:::pipeline.prepareIndata, env)
    util.call(oposSOM:::pipeline.generateSOM, env)

    filename <- paste(env$files.name, "pre.RData")
    util.info("Saving environment image:", filename)
    save(env, file=filename)

    util.info("Processing Differential Expression Statistics")
    #  Get.Running.Average <- oposSOM:::Get.Running.Average
    util.call(oposSOM:::pipeline.calcStatistics, env)

    util.info("Detecting Spots")
    util.call(oposSOM:::pipeline.detectSpotsSamples, env)
    util.call(oposSOM:::pipeline.detectSpotsIntegral, env)
    util.call(oposSOM:::pipeline.patAssignment, env)
    util.call(oposSOM:::pipeline.groupAssignment, env)
  }

  if (env$preferences$activated.modules$geneset.analysis)
  {
    util.info("Calculating Geneset Enrichment")
    util.call(oposSOM:::pipeline.genesetStatisticSamples, env)
    util.call(oposSOM:::pipeline.genesetStatisticIntegral, env)
  }


  if(env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
  {
    filename <- paste(env$files.name, ".RData", sep="")
    util.info("Saving environment image:", filename)
    save(list=c("env","modsom"), file=filename)

    if (file.exists(paste(env$files.name, "pre.RData")) && file.exists(filename))
    {
      file.remove(paste(env$files.name, "pre.RData"))
    }
  }

  #### Reporting part ####


  util.info("Plotting Supporting Information")
  util.call(oposSOM:::pipeline.supportingMaps, env)
  #util.call(oposSOM:::pipeline.entropyProfiles, env)
  #util.call(oposSOM:::pipeline.topologyProfiles, env)



  util.info("Finished:", format(Sys.time(), "%a %b %d %X"))
}

eachOME.run <- function(env)
{
  ids=grep(modsom$ome,colnames(env$indata))

  dir.create(paste(env$files.name, "- Results/only exp - Results/"), showWarnings=FALSE)
  dir.create(paste(env$files.name, "- Results/only meth - Results/"), showWarnings=FALSE)
  dir.create(paste(env$files.name, "- Results/ScoV - Results/"), showWarnings=FALSE)


  setwd(paste(env$files.name, "- Results/") )

  #library(doParallel)

  cl2 <- makeCluster(3)

  cat( "\n\n\nStarted:", format(Sys.time(), "%a %b %d %X\n" ) )

  registerDoParallel(cl2)

#  clusterCall(cl2,function()
#  {
#    library(igraph)
#    library(ape)
#    library(tsne)
#    library(fastICA)
#    library(scatterplot3d)
#    library(pixmap)
#    library(fdrtool)
#    library(Biobase)
#    library(biomaRt)
#    library(oposSOM)})


  w=modsom$w
  w_ome=  modsom$w_ome
  w_exp=  modsom$w_exp

  indata.exp=env$indata[,-ids]/(w)*w_exp
  metadata.exp=env$metadata[,-ids]/(w)*w_exp
  metadata.ome=env$metadata[,ids]/(1-w)*w_ome
  indata.ome=env$indata[,ids]/(1-w)*w_ome


  group.labels.exp=env$group.labels[-ids]
  group.labels.meth=env$group.labels[ids]

  group.metadata.exp = do.call(cbind, by(t(metadata.exp), group.labels.exp, colMeans))[,unique(group.labels.exp)]
  group.metadata.meth = do.call(cbind, by(t(metadata.ome), group.labels.meth, colMeans))[,unique(group.labels.meth)]


  if(group.labels.exp==group.labels.meth)
    {util.call(meanportraits, env)}


  foreach(iterat=1:3) %dopar%
  {

    if(iterat==1)
    {
      env$files.name="only exp"
      env$group.labels=env$group.labels[-ids]
      env$indata=indata.exp
      env$metadata=metadata.exp
      env$group.colors=env$group.colors[-ids]

      env$samples.GSZ.scores=env$samples.GSZ.scores[,-ids]

      env$spot.list.overexpression$spotdata=env$spot.list.overexpression$spotdata[,-ids]/(w)*w_exp
      env$spot.list.underexpression$spotdata=env$spot.list.underexpression$spotdata[,-ids]/(w)*w_exp
      env$spot.list.correlation$spotdata=env$spot.list.correlation$spotdata[,-ids]/(w)*w_exp
      env$spot.list.kmeans$spotdata=env$spot.list.kmeans$spotdata[,-ids]/(w)*w_exp
      env$spot.list.dmap$spotdata=env$spot.list.dmap$spotdata[,-ids]/(w)*w_exp
      env$spot.list.group.overexpression $spotdata=env$spot.list.group.overexpression $spotdata[,-ids]/(w)*w_exp
      env$spot.list.samples=env$spot.list.samples [-ids]
      env$pat.labels= env$pat.labels [-ids]
    }
    if(iterat==2)
    {
      env$files.name="only meth"
      env$group.labels=env$group.labels[ids]
      env$metadata=metadata.ome
      env$indata=indata.ome
      env$group.colors=env$group.colors[ids]

      env$samples.GSZ.scores=env$samples.GSZ.scores[,ids]

      env$spot.list.overexpression$spotdata=env$spot.list.overexpression$spotdata[,ids]/(1-w)*w_ome
      env$spot.list.underexpression$spotdata=env$spot.list.underexpression$spotdata[,ids]/(1-w)*w_ome
      env$spot.list.correlation$spotdata=env$spot.list.correlation$spotdata[,ids]/(1-w)*w_ome
      env$spot.list.kmeans$spotdata=env$spot.list.kmeans$spotdata[,ids]/(1-w)*w_ome
      env$spot.list.dmap$spotdata=env$spot.list.dmap$spotdata[,ids]/(1-w)*w_ome
      env$spot.list.group.overexpression $spotdata=env$spot.list.group.overexpression $spotdata[,ids]/(1-w)*w_ome
      env$spot.list.samples=env$spot.list.samples [ids]
      env$pat.labels= env$pat.labels[ids]
    }
    if(iterat==3 && group.labels.exp==group.labels.meth)
    {
      env$files.name="ScoV"
      env$metadata=sign(metadata.ome*metadata.exp)*sqrt(abs(metadata.ome*metadata.exp))
      env$indata=sign(indata.ome*indata.exp)*sqrt(abs(indata.ome*indata.exp))

      env$group.labels=env$group.labels[ids]
      env$group.colors=env$group.colors[ids]
      env$groupwise.group.colors = env$group.colors[match(unique(env$group.labels), env$group.labels)]
      #		env$samples.GSZ.scores=env$samples.GSZ.scores[,-ids]

      #		env$spot.list.overexpression$spotdata=env$spot.list.overexpression$spotdata[,-ids]
      #		env$spot.list.underexpression$spotdata=env$spot.list.underexpression$spotdata[,-ids]
      #		env$spot.list.correlation$spotdata=env$spot.list.correlation$spotdata[,-ids]
      #		env$spot.list.kmeans$spotdata=env$spot.list.kmeans$spotdata[,-ids]
      #		env$spot.list.dmap$spotdata=env$spot.list.dmap$spotdata[,-ids]
      #		env$spot.list.group.overexpression $spotdata=env$spot.list.group.overexpression $spotdata[,-ids]
      #		env$spot.list.samples=env$spot.list.samples [-ids]
      #		env$pat.labels= env$pat.labels [-ids]
      env$preferences$spot.coresize.modules = 7
      env$preferences$spot.threshold.modules = 0.93
      env$preferences$spot.coresize.groupmap = 3
      env$preferences$spot.threshold.groupmap = 0.75
      env$indata.gene.mean =rowMeans(env$indata)
      env$indata.sample.mean =colMeans(env$indata)

    }
    #library(oposSOM)

    dir.create(paste(env$files.name, "- Results/Sample Similarity Analysis"), showWarnings=FALSE)
    dir.create(paste(env$files.name, "- Results/Summary Sheets - Groups"), showWarnings=FALSE)
    dir.create(paste(env$files.name, "- Results/Summary Sheets - Modules"), showWarnings=FALSE)
    dir.create(paste(env$files.name, "- Results/CSV Sheets"), showWarnings=FALSE)

    #for(i in 1: length(dir("/homes/biertruck/hlydia/skripte/opossom/src/R/")))
    #	source(paste("/homes/biertruck/hlydia/skripte/opossom/src/R/",dir("/homes/biertruck/hlydia/skripte/opossom/src/R/")[i],sep=""),local=F)

    env$output.paths =
      c("LPE"=paste(env$files.name, "- Results/LPE"),
        "CSV"=paste(env$files.name, "- Results/CSV Sheets"),
        "Summary Sheets Samples"=paste(env$files.name, "- Results/Summary Sheets - Samples"),
        "Summary Sheets Modules"=paste(env$files.name, "- Results/Summary Sheets - Modules"))

 #   pipeline.moduleCorrelationMap <- oposSOM:::pipeline.moduleCorrelationMap
#    modules.report.sheets <- oposSOM:::modules.report.sheets

    if(iterat==3 && group.labels.exp==group.labels.meth)
    {
      util.call(oposSOM:::pipeline.calcStatistics, env)

      util.info("Detecting Spots")
      util.call(oposSOM:::pipeline.detectSpotsSamples, env)
      util.call(oposSOM:::pipeline.detectSpotsIntegral, env)
      util.call(oposSOM:::pipeline.patAssignment, env)
      util.call(oposSOM:::pipeline.groupAssignment, env)

      if (env$preferences$activated.modules$geneset.analysis)
      {
        util.call(oposSOM:::pipeline.genesetStatisticSamples, env)
        util.call(oposSOM:::pipeline.genesetStatisticIntegral, env)
      }

    }
   # util.call(oposSOM:::pipeline.supportingMaps, env)
    if(iterat!=3)
    {
      util.call(oposSOM:::pipeline.entropyProfiles, env)
      util.call(oposSOM:::pipeline.topologyProfiles, env)
    }

    util.call(oposSOM:::pipeline.sampleExpressionPortraits, env)
    util.call(oposSOM:::pipeline.patAssignment, env)
    util.call(oposSOM:::pipeline.groupAssignment, env)
    #util.call(oposSOM:::pipeline.groupAnalysis, env)


   # modules.CSV.sheets<-oposSOM:::modules.CSV.sheets
  #  groupSpecificGenesets<-  oposSOM:::pipeline.groupSpecificGenesets


    util.call(oposSOM:::pipeline.sampleSimilarityAnalysisED , env)
    util.call(oposSOM:::pipeline.sampleSimilarityAnalysisSOM, env)
    util.call(oposSOM:::pipeline.sampleSimilarityAnalysisCor, env)
    util.call(oposSOM:::pipeline.sampleSimilarityAnalysisICA, env)

    if (env$preferences$activated.modules$geneset.analysis)
    {
      dir.create(paste(env$files.name, "- Results/Geneset Analysis"), showWarnings=FALSE)

      util.info("Plotting Geneset Enrichment Heatmaps")
      util.call(oposSOM:::pipeline.genesetOverviews, env)

      util.info("Plotting Geneset Profiles and Maps")
      util.call(oposSOM:::pipeline.genesetProfilesAndMaps, env)
    }

    util.info("Calculating Cancer Hallmark Enrichment")
    util.call(oposSOM:::pipeline.cancerHallmarks, env)
    util.info("Writing Gene Lists")
    util.call(oposSOM:::pipeline.geneLists, env)

    util.info("Plotting Summary Sheets (Samples)")
    util.call(oposSOM:::pipeline.summarySheetsSamples, env)

    util.info("Plotting Summary Sheets (Modules & PATs)")

    if(iterat==3 && group.labels.exp==group.labels.meth)
    {
      util.call(ScoV.underexp,env)
    }
   # util.call(oposSOM:::pipeline.summarySheetsModules, env)
    util.call(oposSOM:::pipeline.summarySheetsModules, env)
    util.call(oposSOM:::pipeline.summarySheetsPATs, env)


    #   if(env$preferences$activated.modules$group.analysis && length(unique(env$group.labels)) >= 2)
    #{
            util.info("Processing Group-centered Analyses")
            util.call(oposSOM:::pipeline.groupAnalysis, env)
    #}

            util.info("Processing Difference Analyses")
            util.call(oposSOM:::pipeline.differenceAnalyses, env)


    filename <- paste(env$files.name, ".RData", sep="")
    util.info("Saving environment image:", filename)
    save(env, file=filename)

    #Spotsegmentation for each ome seperately

    if(iterat!=3)
    {
    util.call(oposSOM:::pipeline.calcStatistics, env)
    util.info("Detecting Spots")
    util.call(oposSOM:::pipeline.detectSpotsSamples, env)
    util.call(oposSOM:::pipeline.detectSpotsIntegral, env)

    if (env$preferences$activated.modules$geneset.analysis)
    {
      util.info("Calculating Geneset Enrichment")
      util.call(oposSOM:::pipeline.genesetStatisticIntegral, env)
    }


    env$output.paths =
      c("LPE"=paste(env$files.name, "- Results/LPE"),
        "CSV"=paste(env$files.name, "- Results/CSV Sheets",env$files.name),
        "Summary Sheets Samples"=paste(env$files.name, "- Results/Summary Sheets - Samples"),
        "Summary Sheets Modules"=paste(env$files.name, "- Results/Summary Sheets - Modules",env$files.name))

    util.call(oposSOM:::pipeline.supportingMaps, env)

    util.info("Plotting Summary Sheets (Modules & PATs)")
    util.call(summarySheetsModules2, env)
    }

  }

  util.info("Calculation for each OME finished")
  stopCluster(cl2)


}
