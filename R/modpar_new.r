#' Settings for multiSOMe analysis
#'
#' @param dataset.name all output will be stored in a folder with the name given
#' @param meth.type  type of Methylation Assay used. Either '27K' for Infinium HumanMethylation27K or 'K450' for Infinium HumanMethylation450K. If methylation data is already preprocessed and the data matrix contains methylation values per gene per sample set meth.type=""
#' @param w  weight used to either focus on gene expression (w=0), DNA methylation (w=1) or on both (w=0.5)
#' @param mean.norm.exp  enables centralization of the expression features (boolean)
#' @param mean.norm.meth  enables centralization of the methylation features (boolean)
#' @param quant.exp  enables quantile normalization of the expression samples (boolean)
#' @param quant.meth  enables quantile normalization of the expression samples (boolean)
#' @param M  enables calculation of M values out of methylation beta values (boolean)
#' @param row.ids.exp  type of rowname identifier in biomaRt database of gene expression matrix, for example 'ensembl_gene_id', "refseq_mrna" or 'external_gene_name'. For details see https://www.bioconductor.org/packages//2.7/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
#' @param row.ids.meth  type of rowname identifier in biomaRt database of DNA methylation matrix, for example 'ensembl_gene_id', "refseq_mrna" or 'external_gene_name'. For details see https://www.bioconductor.org/packages//2.7/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
#' @param return.ids  type of rowname identifier in biomaRt database to return, for example 'ensembl_gene_id' or 'external_gene_name'. For details see https://www.bioconductor.org/packages//2.7/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
#'
#' @return modpar
#' @export
modpar.new <- function(modpar=NULL)
{
  return(modpar)
}
