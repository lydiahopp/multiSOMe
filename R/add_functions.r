r <- unclass(lsf.str(envir = asNamespace("oposSOM"), all = T))

for(name in r) eval(parse(text=paste0(name, '<-oposSOM:::', name)))



meanportraits <- function()
{

  pdf("Mean_portraits_ScoVs_all_omes.pdf",15,(1.5*length(unique(group.labels))+1))
  layout(t(matrix(c(1:(10*(1+length(unique(group.labels))))),ncol=(1+length(unique(group.labels))))),width=c(2,2,2,2,2),height=c(1,rep(2,length(unique(group.labels)))))
  #par(mfrow=c(9,7))
  options("scipen"=-100,digits=5)
  par(mar=c(1,1,1,1))
  plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
  for(j in 1:4)
  {
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
    text(0.5,0.5,c("Mean Exp","Mean Meth", "ScoV Meth Exp", "CoV Meth Exp")[j],cex=1.5)
  }
  for(i in 1:length(unique(group.labels)))
  {

    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
    text(0.5,0.5,unique(group.labels)[i],cex=2)
    image( matrix( group.metadata.exp[,i], preferences$dim.1stLvlSom,preferences$dim.1stLvlSom ), axes=F, col = color.palette.portraits(1000),zlim=max(abs(group.metadata.exp[,i]))*c(-1,1)  )
    box()
    image( matrix( group.metadata.meth[,i], preferences$dim.1stLvlSom,preferences$dim.1stLvlSom ), axes=F, col = color.palette.portraits(1000),zlim=max(abs(group.metadata.meth[,i]))*c(-1,1)  )
    box()

    s=sign(group.metadata.meth*group.metadata.exp)*sqrt(abs(group.metadata.meth*group.metadata.exp))
    #s=do.call(cbind, by(t(s),group.labels[grep("cnv",names(group.labels))], colMeans)[unique(group.labels)])

    image( matrix( s[,i], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom ), axes=F, col = color.palette.portraits(1000),zlim=max(abs(s[,i]))*c(-1,1)  )
    #	mtext(text=paste0("(",as.character(round(min(s[,i]),3)),",",as.character(round(max(s[,i]),3)),")"),side=1,line=0.9,cex=1.4)
    box()

    s=(sapply(group.metadata.meth[,i],function(x){x-mean(group.metadata.meth[,i])}))*(sapply(group.metadata.exp[,i],function(x){x-mean(group.metadata.exp[,i])}))
    #s=do.call(cbind, by(t(s),group.labels[grep("cnv",names(group.labels))], colMeans)[unique(group.labels)])

    image( matrix( s, preferences$dim.1stLvlSom, preferences$dim.1stLvlSom ), axes=F, col = color.palette.portraits(1000),zlim=max(abs(s))*c(-1,1)  )
    #	mtext(text=paste0("(",as.character(round(min(s[,i]),3)),",",as.character(round(max(s[,i]),3)),")"),side=1,line=0.9,cex=1.4)
    box()

  }
  dev.off()

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


