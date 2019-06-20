r <- unclass(lsf.str(envir = asNamespace("oposSOM"), all = T))

for(name in r) eval(parse(text=paste0(name, '<-oposSOM:::', name)))



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


