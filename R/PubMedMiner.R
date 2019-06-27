#' Generates a list of Pubmed results based on a UniProt ID and a user-defined keyword.
#'
#' @param df R dataframe having 5 attributes, UniProtID, IDType, TaxID, Keyword, KeywordInTitleOnly.
#' @param output.file Path to the output  file with the results of the Pubmed query.
#' @return Generates an Excel file \code{output.file} with Pubmed query results using the UniProt identifers and keyword search in \code{query.file}.
#' @examples
#' pubmedMiner_entry(potentialmarker, output.file = "potential_marker_pubmed_results.xlsx")
#' dat.input <- openxlsx::readWorkbook("~/Input.xlsx")
#' pubmedMiner_entry(dat.input, output.file = "~/Pubmed_results.xlsx")
#' @export
pubmedMiner_entry <- function(df, output.file = "pubmed_results.xlsx") {
  list.datquery = list()
  list.datpubmed = list()

  # create progress bar
  pb <- txtProgressBar(min = 0, max = nrow(df), style = 3)

  for(query.idx in 1:nrow(dat.input)) {
    Sys.sleep(0.1)
    UniProtID = df[query.idx, 'UniProtID']
    IDType = df[query.idx, "IDType"]
    taxid = df[query.idx, 'TaxID']
    keyword =  df[query.idx, 'Keyword']
    ti.only = df[query.idx,5]

    pd.res <- try(pubmed_miner(UniProtID, IDType, taxid, keyword, ti.only, query.idx=query.idx))

    if(!inherits(pd.res, 'try-error')) {

      list.datquery[[query.idx]] = pd.res$dat.query
      list.datpubmed[[query.idx]] = pd.res$dat.pubmed

    } else {
      list.datquery[[query.idx]] = data.frame('UniProtID'=UniProtID, 'GeneID' = NA, 'TaxID'=taxid, 'Synonyms' = paste(synonyms, collapse=','),
        'Keywords' = keyword, 'KeywordInTitleOnly'=ti.only, 'TotalResults'=0,
        'Category' = 3, 'False' = 0, 'PubmedQuery'=NA)

      list.datpubmed[[query.idx]] = 'No result'

    }
    # update progress bar
    setTxtProgressBar(pb, query.idx)
  }

  all.datquery = do.call(rbind, list.datquery)
  all.datquery$NQuery = 1:nrow(all.datquery)
  all.datquery = all.datquery[,c(ncol(all.datquery),1:(ncol(all.datquery)-1))]

  all.datquery[is.na(all.datquery)] = ''
  # output

  wb = openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, sheetName="query")

  openxlsx::writeData(wb,"query", all.datquery)

  for(ii in 1:length(list.datpubmed) ) {

    list.datquery[[ii]][is.na(list.datquery[[ii]])] = ''

    openxlsx::addWorksheet(wb, sheetName=paste("pubmed result", ii) )
    openxlsx::writeData(wb, paste("pubmed result", ii), list.datquery[[ii]])
    openxlsx::writeData(wb, paste("pubmed result", ii), list.datpubmed[[ii]], startRow=4)

    if(file.exists(paste('barplotNwordcloud', ii, '.png', sep='') ) )
      openxlsx::insertImage(wb, paste("pubmed result", ii), paste('barplotNwordcloud', ii, '.png', sep=''),
        width=5, height=8, startRow = 3, startCol=12)

    if(file.exists(paste('plot_dist_mesh', ii, '.png', sep='') ) )
      openxlsx::insertImage(wb, paste("pubmed result", ii), paste('plot_dist_mesh', ii, '.png', sep=''),
        width=5, height=5, startRow = 3, startCol=20)
  }

  openxlsx::saveWorkbook(wb, output.file, overwrite=T)

  list(all.datquery=all.datquery, list.datpubmed=list.datpubmed, list.datquery=list.datquery)
}






