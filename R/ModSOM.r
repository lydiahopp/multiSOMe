#' Preparation of data for multiSOMe analysis
#'
#' @param indata.exp  gene expression matrix, Should have the same colnames as indata.meth, suffix "exp" will be added automatically
#' @param indata.meth  DNA methylation matrix, Should have the same colnames as indata.exp. Suffix of ome will be added to the colnames.
#' @param group.labels character vector of length ncol(indata.exp). Should be in the same order as samples in expression matrix.
#' @param modpar  see manual of modpar function
#'
#' @return modified matrix for input, harmonization weights for expression and methylation (to make them comparable)
#' @export
#'

ModSOM <- function(indata.exp,indata.meth,group.labels,modpar)
  {
 # indata.meth=indata.meth[,colnames(indata.exp)]

  meth.type=modpar$meth.type
  w=modpar$weight
  mean.norm.exp=modpar$mean.norm.exp
  mean.norm.meth=modpar$mean.norm.meth
  quant.exp=modpar$quant.exp
  quant.meth=modpar$quant.meth
  M=modpar$M
  row.ids.exp=modpar$row.ids.exp
  row.ids.meth=modpar$row.ids.meth
  return.ids=modpar$return.ids
  ome="meth"
  database.dataset=modpar$database.dataset
  database.biomart=modpar$database.biomart


  if(meth.type=="27K")
  {
    data(gene_annotation_27K)
    annotation=annotation_27K
  }
  if(meth.type=="450K")
    {
    data(gene_annotation_450K)
    annotation=annotation_450K
    }
  if(meth.type!="")
  {
    indata.meth=indata.meth[na.omit(match(annotation[,2],rownames(indata.meth))),]

    if(any(is.na(match(annotation[,2],rownames(indata.meth)))) )
      annotation=annotation[-which(is.na(match(annotation[,2],rownames(indata.meth)))),]

    annotation_vec=annotation[,1]
    names(annotation_vec)=annotation[,2]
    indata.meth=do.call(rbind, by(indata.meth,annotation_vec, colMeans))[unique(annotation_vec),]
  }

  if(any(indata.meth==0))
    indata.meth[which(indata.meth==0)]=0.00001


  if(row.ids.meth!=row.ids.exp)
  {
    library("biomaRt" )

    mart<-useMart(biomart = database.biomart,host = "aug2017.archive.ensembl.org" ) #  "www.ensembl.org"
    mart<-useDataset(database.dataset ,mart=mart)

    biomart.table = getBM( c( row.ids.meth, row.ids.exp, return.ids ) ,row.ids.meth ,rownames(indata.meth), mart, checkFilters=F )

    if(length(which(apply(biomart.table,1,function(x){any(x=="")})==TRUE))>0)
      biomart.table= biomart.table[-which(apply(biomart.table,1,function(x){any(x=="")})==TRUE),]

    if(any(is.na(match(biomart.table[,2],rownames(indata.exp)))))
      biomart.table=biomart.table[-which(is.na(match(biomart.table[,2],rownames(indata.exp)))),]

    indata.meth=  indata.meth[biomart.table[,1],]
    indata.exp=  indata.exp[biomart.table[,2],]

    rn=biomart.table[,3]
    names(rn)=biomart.table[,2]
    indata.exp=do.call(rbind,by(indata.exp,rn,colMeans))[unique(rn),]

    names(rn)=biomart.table[,1]
    indata.meth=do.call(rbind,by(indata.meth,rn,colMeans))[unique(rn),]
  }else
  {
    rn=intersect(rownames(indata.exp),rownames(indata.meth))

    indata.meth= indata.meth[rn,]
    indata.exp=indata.exp[rn,]
  }


  if(ome=="meth")
  {
    library("biomaRt" )

    mart<-useMart(biomart = database.biomart, host ="aug2017.archive.ensembl.org" ) #   "jul2015.archive.ensembl.org"
    mart<-useDataset(database.dataset ,mart=mart)

    biomart.table = getBM( c( row.ids.meth,  "chromosome_name","strand" ) ,row.ids.meth ,rn, mart, checkFilters=F )

    if((length(which(biomart.table$chromosome_name=="X"))+length(which(biomart.table$chromosome_name=="Y")))>0)
    {
      indata.meth=indata.meth[- match(biomart.table[c(which(biomart.table$chromosome_name=="X"),which(biomart.table$chromosome_name=="Y")),1],rn),]

      rn=intersect(rownames(indata.exp),rownames(indata.meth))

      indata.meth= indata.meth[rn,]
      indata.exp=indata.exp[rn,]
    }
  }


  if(ome=="meth" && M==T)
  {
    indata.meth[which(indata.meth==0)]=0.00001
    M=log10(indata.meth/(1-indata.meth))
    indata.meth=M

    cat( "\nTransform Beta into M-values\n" ); flush.console()
  }

  if(quant.exp)
  {
    indata.exp = oposSOM:::Quantile.Normalization( indata.exp )

    cat( "\nQuantile normalization indata.exp\n" ); flush.console()
  }

  if(quant.meth)
  {
    indata.meth = oposSOM:::Quantile.Normalization( indata.meth )
    cat( "\nQuantile normalization indata.meth\n" ); flush.console()

  }


  if(mean.norm.exp)
  {
    indata.exp = indata.exp - rowMeans( indata.exp )
    cat( "\nZentralization indata.exp\n" ); flush.console()

  }


  if(mean.norm.meth)
  {
    indata.meth = indata.meth - rowMeans( indata.meth )
    cat( "\nZentralization indata.meth\n" ); flush.console()
  }



  w_ome=mean(unlist(abs(indata.meth)))
  indata.meth= indata.meth/w_ome
  w_exp=mean(unlist(abs(indata.exp)))
  indata.exp=indata.exp/w_exp

  cat( "\nHarmonization\n" ); flush.console()

  if(w==0)
    w=0.001
  if(w==1)
    w=0.999

  indata.meth=indata.meth*(1-w)

  indata.exp=indata.exp*(w)
  colnames(indata.meth)=paste(colnames(indata.meth),ome,sep="_")
  colnames(indata.exp)=paste(colnames(indata.exp),"exp",sep="_")
  indata = cbind(indata.exp,indata.meth)


  group.labels = c(group.labels,group.labels)
  names(group.labels)=colnames(indata)

  indata=as.matrix(indata[,order(group.labels)])
  group.labels=sort(group.labels)


  mod_par=list(indata,w_ome,w_exp,w,ome,indata.meth,indata.exp,group.labels)
  names(mod_par)=c("indata","w_ome","w_exp","w","ome","indata.meth","indata.exp","group.labels")

  return(mod_par)
}

