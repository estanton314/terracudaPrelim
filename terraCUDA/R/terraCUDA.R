#' @title wordcloudfunc
#' @description creates a word cloud from description of workspace
#' @examples
#' wordcloudfunc("broad-firecloud-tcga","TCGA_HNSC_hg38_OpenAccess_GDCDR-12-0_DATA")
#' @param WorkspaceNamespace the workspace namespace
#' @param Workspace the workspace
#' @export
wordcloudfunc <- function(WorkspaceNamespace,Workspace){

  myspace = terra$getWorkspace("WorkspaceNamespace","Workspace")

  mycontent = httr::content(myspace) ##help

  myatts=mycontent$workspace$attributes

  mydesc = myatts$description

  descdoc <- Corpus(VectorSource(mydesc))

  toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  descdoc <- tm_map(descdoc, toSpace, "/")
  descdoc <- tm_map(descdoc, toSpace, "@")
  descdoc <- tm_map(descdoc, toSpace, "\\|")

  # Convert the text to lower case
  descdoc <- tm_map(descdoc, content_transformer(tolower))
  # Remove numbers
  descdoc <- tm_map(descdoc, removeNumbers)
  # Remove english common stopwords
  descdoc <- tm_map(descdoc, removeWords, stopwords("english"))
  # Remove your own stop word
  # specify your stopwords as a character vector
  descdoc <- tm_map(descdoc, removeWords, c("blabla1", "blabla2"))
  # Remove punctuations
  descdoc <- tm_map(descdoc, removePunctuation)
  # Eliminate extra white spaces
  descdoc <- tm_map(descdoc, stripWhitespace)
  # Text stemming
  # descdoc <- tm_map(descdoc, stemDocument)

  dtm <- TermDocumentMatrix(descdoc)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)

  wordcloud(words = d$word, freq = d$freq, min.freq = 1,
            max.words=200, random.order=FALSE, rot.per=0.35,
            colors=brewer.pal(8, "Dark2"))

}

#' @title datamaidreport
#' @description makes a data maid report based on a given table
#' @examples
#' datamaidreport("broad-firecloud-tcga","TCGA_HNSC_hg38_OpenAccess_GDCDR-12-0_DATA","sample")
#' @param WorkspaceNamespace the workspace namespace
#' @param Workspace the workspace
#' @param TableName the table name
#' @export
datamaidreport <- function(WorkspaceNamespace,Workspace,TableName){
  avtables(WorkspaceNamespace,Workspace)

  samps = avtable(TableName,WorkspaceNamespace,Workspace)

  library(dataMaid)
  makeDataReport(samps,replace=T)
}

#' @title getColumnNames
#' @description returns all columns in a given table
#' @examples
#' getColumnNames("broad-firecloud-tcga","TCGA_HNSC_hg38_OpenAccess_GDCDR-12-0_DATA","sample")
#' @param WorkspaceNamespace the workspace namespace
#' @param Workspace the workspace
#' @param TableName the table name
#' @export
getColumnNames <- function(WorkspaceNamespace,Workspace,TableName){
  yourdata=avtable(TableName,WorkspaceNamespace,Workspace)
  colNames = names(yourdata)
  return(colNames)
}

#' @title getUniqueColNames
#' @description returns all unique columns in all tables of a given workspace
#' @examples
#' getUniqueColNames("broad-firecloud-tcga","TCGA_HNSC_hg38_OpenAccess_GDCDR-12-0_DATA")
#' @param WorkspaceNamespace the workspace namespace
#' @param Workspace the workspace
#' @export
getUniqueColNames <- function(WorkspaceNamespace,Workspace){
  tabTib = avtables(WorkspaceNamespace,Workspace)
  tableNames=tabTib$table
  tempCol = c()
  for (tabName in tableNames) {
    tempCol = c(tempCol,getColumnNames(WorkspaceNamespace,Workspace,tabName))
    class(tempCol)
    #print(tempCol)
  }
  uniqCol = unique(tempCol)
  return(uniqCol)
}

#' @title workspaceContains
#' @description returns boolean for whether or not workspace contains a given column name
#' @examples
#' workspaceContains("broad-firecloud-tcga","TCGA_HNSC_hg38_OpenAccess_GDCDR-12-0_DATA","participants.items")
#' @param WorkspaceNamespace the workspace namespace
#' @param Workspace the workspace
#' @param colName the column of interest
#' @export
workspaceContains <- function(WorkspaceNamespace, Workspace, colName){
  colNames = getUniqueColNames(WorkspaceNamespace,Workspace)
  here=FALSE
  for(name in colNames){
    if(name == colName){
      here = TRUE
      return(TRUE)
    }
  }
  if(here != TRUE){
    return(FALSE)
  }
}

#' @title anyWorkspaceContains
#' @description returns list of the workspace/namespaces that contain a given column name
#' @examples
#' anyWorkspaceContains("gender")
#' @param colName the column of interest
#' @export
anyWorkspaceContains<-function(colName){
  require(AnVIL)
  terra=Terra()
  ws <- httr::content(terra$listWorkspaces())
  realnames = c()
  realnamespc = c()
  i=1
  langth = length(ws) + 1
  while(i<langth){
    tempws <- ws[[i]]
    tempname = tempws$workspace$name
    tempnamespc = tempws$workspace$namespace
    if(workspaceContains(tempnamespc,tempname,colName)){
      realnames = c(realnames, tempname)
      realnamespc = c(realnamespc, tempnamespc)
    }
    i = i+1
  }
  return(list(realnames,realnamespc))
}

