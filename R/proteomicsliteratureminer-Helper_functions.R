

pubmed_miner <- function(UniProtID, IDType, taxid, keyword, ti.only, query.idx=1) {
  fields = "TIAB"

  if(ti.only=="Yes") fields = "TI"

  fterm = paste('[', fields, ']', sep='')


  dat.query = data.frame('UniProtID'=UniProtID, 'IDType' = IDType, 'TaxID'=taxid, 'Synonyms' = NA,
    'Keywords' = keyword, 'KeywordInTitleOnly'=ti.only, 'TotalResults'=0,
    'Category' = 3, 'False' = 0, 'PubmedQuery'=NA)

  dat.pubmed = "No pubmed results returned!"

  if(!is.null(UniProtID)) {


    synonyms <- try(getUniprotSynonyms(UniProtID, IDType, taxid))


    if(!inherits(synonyms, "try-error")) {
      if(!is.null(synonyms)) {


        vec.keyword=NA
        if(!is.na(keyword))
          vec.keyword = strsplit(keyword, split=',')[[1]]

        pubmed.result <- try(pubmed_summary(synonyms, vec.keyword, fields) )

        # search pubmed

        pubres = pubmed.result$pubres

        res.count = pubmed.result$res.count

        pub.query = pubmed.result$pub.query

        #dat.query$GeneID = geneid
        dat.query$Synonyms = paste(synonyms, collapse=',')
        dat.query$PubmedQuery = pub.query

        if(res.count > 0) {

          dat.pubmed <- try(pubmed_record(pubres, vec.keyword, synonyms, fields))

          # mining and analysis

          plot_stats(dat.pubmed, file=paste('barplotNwordcloud', query.idx, '.png', sep='') )

          # TO assess

          dat.pubmed$Cluster.byMeSH = mesh_clustering(as.character(dat.pubmed$MeSH),
            file=paste('plot_dist_mesh', query.idx, '.png', sep='') )

          dat.query$TotalResults = nrow(dat.pubmed)


          # category 0: no synonyms; 1: hits that have reviews has review; 2: hits but do not have reviews; 3: no hits
          dat.query$Category = 0
          if(length(synonyms) > 0) {
            dat.query$Category = 2
            if (sum(dat.pubmed$IsReview) > 0) dat.query$Category = 1
          }

          dat.query$False = sum(!(dat.pubmed$HasSynonyms))
        }
      }
    }
  }


  list(dat.query=dat.query, dat.pubmed=dat.pubmed)

}


# get gene synonyms using Uniprot

# modified instead use UniProtID and only look for reviewed

getUniprotSynonyms <- function(UniProtID, IDType="Accession", taxid=9606) {

  if(tolower(IDType) == "accession") {
    query = paste("https://www.uniprot.org/uniprot/?query=accession:", UniProtID, '+organism:', taxid,
      '+reviewed:yes&columns=genes&format=tab', sep='') #added +reviewed:yes for only reviewed hits and accession instead of gene
  } else {
    query = paste("https://www.uniprot.org/uniprot/?query=gene:", UniProtID, '+organism:', taxid,
      '+reviewed:yes&columns=genes&format=tab', sep='') #added +reviewed:yes for only reviewed hits and accession instead of gene
  }


  uniprot.res = try(httr::GET(query))

  list.syn = NULL
  if(!inherits(uniprot.res, "try-error"))
    if(!is.null(httr::content(uniprot.res)))
      list.syn = unique( strsplit(httr::content(uniprot.res), split="\n")[[1]][-1])
  #
  # get the first synonym in each

  res = NULL
  if(!is.null(list.syn))
    res = as.vector(na.omit(unlist(sapply(list.syn, function(x) strsplit(x, split=' ')[[1]]))))

  # remove ; in some uniprot synonyms

  res = gsub("\\;", "", res)


  as.vector(unique(res))
}



pubmed_summary <- function(synonyms, vec.keyword=NA, fields="TIAB") {


  fterm = paste('[', fields, ']', sep='')
  qkeyword = sapply(vec.keyword, function(x) paste0("\"",x,"\"") )

  qsynonym = sapply(synonyms, function(x) paste0("\"",x,"\"") )

  if(!is.na(vec.keyword)) {
    keyword4query = paste(qkeyword, collapse=paste(' AND', fterm, sep='') )

    pub_query = paste0(keyword4query, fterm, ' AND ( ', paste0(qsynonym, collapse=paste0(fterm, ' OR ') ), fterm,')' )

  } else {
    pub_query = paste0(qsynonym, collapse=paste0(fterm, ' OR ') )
  }

  pubres = RISmed::EUtilsSummary(pub_query, datetype="pdat")

  res.count = RISmed::QueryCount(pubres)


  list(pubres=pubres, res.count=res.count, pub.query=pub_query)
}


pubmed_record <- function(pubres, vec.keyword=NA, synonyms=NULL, fields="TI") {

  pubrecords = RISmed::EUtilsGet(pubres)

  titles =  RISmed::ArticleTitle(pubrecords)
  abstracts = RISmed::AbstractText(pubrecords)
  years = RISmed::YearPubmed(pubrecords)
  list.authors = RISmed::Author(pubrecords)	 # list
  list.aff = RISmed::Affiliation(pubrecords) # list
  country = RISmed::Country(pubrecords)
  pmid = RISmed::PMID(pubrecords)
  ptyp = RISmed::PublicationType(pubrecords) # list
  is.review = sapply(1:length(ptyp), function(i) "Review" %in% ptyp[[i]])

  # Check false discovery
  in.title = FALSE
  has.synm = FALSE

  if(!is.na(vec.keyword))
    in.title = sapply(titles, function(x)
      any(sapply(vec.keyword, function(y) grepl(tolower(y), tolower(x)))), USE.NAMES=FALSE )

  if(!is.null(synonyms))
    has.synm = sapply(titles, function(x)
      any(sapply(synonyms, function(y) grepl(tolower(y), tolower(x)))), USE.NAMES=FALSE )

  if(fields == "TIAB")
    has.synm = (has.synm + sapply(abstracts, function(x)
      any(sapply(synonyms, function(y) grepl(tolower(y), tolower(x)))), USE.NAMES=FALSE ) ) > 0

  list.mesh = RISmed::Mesh(pubrecords) # list

  authors = sapply(list.authors, function(a) 	paste(paste(a$Initials, a$LastName), collapse=',') )

  affs = sapply(list.aff, function(a) paste(a, collapse=';;') )


  meshs = sapply(list.mesh, function(x) {d=data.frame(x); paste(d$Heading, collapse=',')} )


  # check

  stopifnot(all(c(length(titles)==length(abstracts),
    length(years) == length(list.authors), length(pmid) == length(ptyp) ) ))


  dat.pubmed <- data.frame('Title'=titles, 'Abstract'=abstracts, 'Year'=years,
    'Author'=authors, 'Country' = country, 'Link'=paste0('https://www.ncbi.nlm.nih.gov/pubmed/',pmid),
    'ConditionInTitle'=in.title, 'IsReview'=is.review, 'MeSH'=meshs, 'HasSynonyms' = has.synm)
  # clean up
  dat.pubmed = dat.pubmed[!grepl('^\\[', dat.pubmed$Title),]

  dat.pubmed
}


plot_stats <- function(dat.pubmed, file='barplotNwordcloud.png') {

  png(file, 2000, 3000, res=300)

  layout(matrix(1:2, ncol=1))

  par(mar=c(7,4,4,4)+.1)

  yy = dat.pubmed$Year
  if(length(yy) > 0) {
    count <- table(yy)
    cumcount <- cumsum(count)

    axScale <- ceiling(max(cumcount)/max(count))

    bp <- barplot(count, main="Number of publications by year" )

    lines(bp, cumcount/axScale, lwd=2, col='red')

    yaxp = par('yaxp')

    at4 <- seq(from=yaxp[1], to=yaxp[2], length.out=yaxp[3]+1)
    labels4 <- at4*axScale

    axis(4, at=at4, labels=labels4)
  }


  all.abstract = paste(dat.pubmed$Abstract, collapse=' ')

  # word cloud of abstracts
  suppressWarnings({ wcl = try(wordcloud::wordcloud(all.abstract, max.words=200)) })

  # barplots of top 20 MeSH


  dev.off()



}

abstract_clustering <- function(abstracts, method=c('hierarchical', 'kmeans'), k=4) {

  res <- rep(NA, length(abstracts))

  idx.na = nchar(abstracts)==0

  abstract.corpus <- tm::Corpus(VectorSource(abstracts[!idx.na]))

  # stop-word removing, stemming
  suppressWarnings({
    abstract.corpus <- tm::tm_map(abstract.corpus, stripWhitespace)
    abstract.corpus <- tm::tm_map(abstract.corpus, content_transformer(tolower))
    abstract.corpus <- tm::tm_map(abstract.corpus, removeWords, stopwords("english") )
    abstract.corpus <- tm::tm_map(abstract.corpus, stemDocument)
  })


  dtm.tfidf <- tm::DocumentTermMatrix(abstract.corpus, control=list(removePunctuation=TRUE,
    removeNumbers=TRUE,
    weighting=function(x) weightTfIdf(x, normalize=FALSE)) )

  inspect(dtm.tfidf)

  mat.tfidf <- as.matrix(	dtm.tfidf)

  if(method == 'hierarchical') {
    # pairwise distance matrix
    d.cosine <- cosineDist(mat.tfidf)


    # hierarchical clustering

    hc.cosine <- hclust(d.cosine , method="ward.D")


    ct.cosine <- cutree(hc.cosine, k=k)

    clusterID.cosine <- as.vector(ct.cosine)

    res[!idx.na] <- clusterID.cosine
  } else { #kmeans

  }
  # res

}


mesh_matrix <- function(v, u) { # v and u are two mesh terms

  1*(u %in% v)

}

# meshs is a vector of MeSH terms delimited by ','
mesh_clustering <- function(meshs,k=4, file='plot_dist_mesh.png') {


  res <- rep(NA, length(meshs))

  idx.na = nchar(meshs)==0

  nna.meshs = meshs[!idx.na]

  list.meshs = lapply(nna.meshs, function(x) strsplit(x, split=',')[[1]])

  all.meshs = unlist(list.meshs)

  all.uni.meshs = unique(all.meshs)

  NumTopMeash = min(c(20, length(all.uni.meshs)))

  if(length(list.meshs) > 0) {

    mat.mesh = matrix(NA, nrow=length(list.meshs), ncol=length(all.uni.meshs))

    for(ii in 1:length(	list.meshs) )
      mat.mesh[ii,] = mesh_matrix(list.meshs[[ii]], all.uni.meshs)

    stopifnot(all(rowSums(mat.mesh) > 0))

    png(file, 2000, 2000, res=300)

    par(mar=c(12,5,4,2)+.1)
    bp = barplot(sort(table(all.meshs), decreasing=TRUE)[1:NumTopMeash],

      main=paste('Top', NumTopMeash, 'MeSH categories'),
      las=2, cex.lab=.8,
      ylab='Frequency')


    dev.off()

  }

  if(length(nna.meshs) > 10) {

    # pairwise distance matrix
    d.cosine <- cosineDist(mat.mesh)


    hc.cosine <- hclust(d.cosine , method="complete")


    ct.cosine <- cutree(hc.cosine, k=k)

    clusterID.cosine <- as.vector(ct.cosine)

    res[!idx.na] <- clusterID.cosine
  }
  res

}


# cosine similarity, cos(a,b) = a.b/norm(a)*norm(b)
cosineDist <- function(x) {
  as.dist(1-x%*%t(x)/(sqrt(rowSums(x^2) %*% t(sqrt(rowSums(x^2)))) ))

}
