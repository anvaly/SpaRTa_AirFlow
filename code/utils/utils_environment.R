# Author: Anna Lyubetskaya. Date: 20-03-04


render_setup_my <- function(){
  ## Adjust environment for HTMP rendering

  cat("Report created on:", format(Sys.Date(), format="%Y-%B-%d"))
  # Switch off warnings if rendering
  options(warn=-1)
  # Switch of messages from readr::read_delim()
  options(readr.num_columns = 0)
  # Always print all columns in a tibble
  options(tibble.width = Inf) 
  pander::panderOptions("table.split.table", Inf)
  # Setup figures
  knitr::opts_chunk$set(dev="png", out.width="100%", fig.width=8, fig.height=8, fig.show="hold")
}


run_parallel_wrap_my <- function(){
  ## Run a function using parallel library
  
  result_list <- list()
  
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  
  doParallel::registerDoParallel(cl)
  
  result_list <- foreach::foreach(input = input_file_paths, .combine='c') %dopar% {
    # Function of interest
  };
  
  parallel::stopCluster(cl)
}
