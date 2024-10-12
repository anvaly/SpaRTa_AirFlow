# Author: Anna Lyubetskaya. Date: 19-10-04
# Various tibble operations


df_transpose_my <- function(df){
  ## Transpose a tibble
  ## Tibble should have both colnames and rownames
  
  return(tibble::as_tibble(cbind(X1 = names(df), t(df))))
}


df_edit_colnames_my <- function(df, regex=" (.+)", replacement=""){
  ## Use a regex to cleanup the column names
  
  names(df) <- gsub(x=names(df), pattern=regex, replacement=replacement)
  
  return(df)
}


df_wide2long_my <- function(wide_df, key="Hugo_Symbol", val="CERES", start_col=2){
  ## Transform wide matrix to long
  
  # Gather column names of the wide matrix
  cols <- colnames(wide_df)
  
  # Reshape data from wide format to long
  long_df <- wide_df %>% 
    tidyr::gather(key=!!rlang::sym(key), value=!!rlang::sym(val), cols[start_col]:cols[length(cols)])
  
  return(long_df)
}


df_long2wide_my <- function(long_df, rows="DepMap_ID", cols="Hugo_Symbol", value="Expression"){
  ## Transform long matrix to wide
  ## spread() and gather() are not mirror functions
  ## Be careful to prep your input if necessary
  ## https://rstudio-pubs-static.s3.amazonaws.com/282405_e280f5f0073544d7be417cde893d78d0.html
  ## https://github.com/tidyverse/tidyr/issues/426
  ## https://community.rstudio.com/t/spread-why-errors/2076/3
  
  # Select only relevant columns
  long_df <- long_df %>% 
    dplyr::select(all_of(c(cols, rows, value)))
  # Clean up duplicates
  long_df <- long_df[!duplicated(long_df),]
  
  # Format data from long to wide
  wide_df <- long_df %>% 
    tidyr::spread(key=cols, value=value)
  
  # This is a solution when spread() keeps complaining about "keys shared for rows"
  nothing_else_works <- FALSE
  if(nothing_else_works == TRUE){
    wide_df <- long_df %>% 
      dplyr::group_by(!!rlang::sym(cols), !!rlang::sym(rows)) %>%  # group by everything other than the value column
      dplyr::mutate(row_id=1:n()) %>% # add an index to make rows unique
      multiplyr::ungroup() %>%  # build group index
      tidyr::spread(key=cols, value=value) # spread
    dplyr::select(-row_id)
  }
  
  return(wide_df)
}


df_zscore_my <- function(df, col_by="Hugo_Symbol", value="Expression"){
  ## Calculate z-score of given values within given groups
  ## Z-scoring of a matrix: https://stackoverflow.com/questions/46185816/how-do-i-calculate-a-grouped-z-score-in-r-using-dplyr
  
  col_name <- paste0(value, "_zscore")
  
  zscore_df <- df %>% 
    dplyr::group_by(!!rlang::sym(col_by)) %>% 
    dplyr::mutate(ZSCORE = c(scale(!!rlang::sym(value))))
  
  colnames(zscore_df) <- gsub("^ZSCORE$", paste0(value, "_zscore"), colnames(zscore_df))
    
  return(zscore_df)
}


df_remove_empty_columns_my <- function(df){
  ## Remove all empty columns
  
  return(df %>% 
           purrr::discard(~all(is.na(.x))) %>%
           purrr::map_df(~.x))
}


list2tibble_my <- function(data_list, rownames="id", transpose=FALSE){
  ## Transform a named list of lists into a tibble
  
  data_df <- tibble::as_tibble(as.data.frame(data.table::rbindlist(data_list, use.names=TRUE, fill=TRUE)), rownames=rownames)

  if(transpose == TRUE){
    # Make row columns before transposing
    data_df <- data_df %>%
      tibble::column_to_rownames(rownames)
    
    # Transpose
    data_df <- tibble::as_tibble(cbind(nms = names(data_df), t(data_df)))
    
    # Fix column names using information from a master meta data file
    colnames(data_df) <- c(rownames, names(data_list))
  }
  
  return(data_df)
}


tibble2eset_my <- function(data_df, rownames="Symbol", normalize=TRUE){
  ## Create a normalized eSet object from a tibble
  
  # Create an eSet object of raw CCLE counts
  expr_matrix <- as.matrix(data_df %>% 
                             tibble::column_to_rownames(rownames))
  
  # Make all values numeric
  class(expr_matrix) <- "numeric"
  
  # Create an eset_object
  eset_object <- Biobase::ExpressionSet(assayData = expr_matrix)
  
  if(normalize == TRUE){
    # Create a DGEList Object for edgeR
    dge <- edgeR::DGEList(counts=Biobase::exprs(eset_object), genes = Biobase::fData(eset_object))
    
    # Calculate normalization factors
    dge <- edgeR::calcNormFactors(dge, method="TMM")
    
    # Get counts per million.  
    Biobase::exprs(eset_object) <- edgeR::cpm(dge, log = T)
  }
  
  return(eset_object)
}
