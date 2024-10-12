# Author: Anna Lyubetskaya. Date: 19-10-23
# Fitting a generalized linear model via penalized maximum likelihood using glmnet and caret

# https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
# https://daviddalpiaz.github.io/r4sl/elastic-net.html
# http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/153-penalized-regression-essentials-ridge-lasso-elastic-net/
# https://stats.stackexchange.com/questions/209009/how-to-treat-categorical-predictors-in-lasso



set.seed(7)
source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_specialized_plots.R")


tibble_clean_for_regression <- function(wide_df, outcome_col=NULL, rownames_col=NULL, factor_fields=NULL){
  ## Perform necessary checks on a wide tibble before regression
  
  # # Clean up the tibble
  # wide_df <- wide_df %>% 
  #   unique %>%
  #   tidyr::drop_na()
  
  # Find columns that have more than one value across samples
  col_names <- wide_df %>%
    df_wide2long_my(key="Columns", val="Count") %>%
    # tidyr::drop_na() %>%
    dplyr::group_by(Columns) %>%
    dplyr::summarise(Count_unique = dplyr::n_distinct(Count)) %>%
    dplyr::filter(Count_unique > 1) %>%
    dplyr::pull(Columns)
  
  # Define the outcome and select columns that have range
  wide_df <- wide_df %>% 
    dplyr::select(all_of(c(rownames_col, col_names))) 

  # Separate out outcome column
  if(!is.null(outcome_col)){
    wide_df <- wide_df %>%
      dplyr::rename("outcome" = !!sym(outcome_col))
  }
  
  # Factorize factor fields
  factor_fields <- intersect(colnames(wide_df), factor_fields)
  if(!is.null(factor_fields)){
    wide_df[factor_fields] <- lapply(wide_df[factor_fields], as.factor) 
  }
  
  # Explicitly define row names
  if(!is.null(outcome_col)){
    wide_df <- wide_df %>%
      tibble::column_to_rownames(rownames_col)
  }
  
  return(wide_df)
}


train_model_caret_my <- function(df, formula=NULL, cv_fold=5, tune_length=10){
  ## Train an Elastic Net regularized regression model using caret wraper of glmnet
  
  if(is.null(formula)){
    formula <- "outcome ~ ."
  }
  
  # Alpha: balances ridge (alpha = 0) and lasso (alpha = 1) components
  # Lambda: amount of ridge and lasso penalization against overfitting (lambda = 0 - no shrinkage)
  model <- caret::train(as.formula(formula), 
                        data = df,
                        method = "glmnet",
                        trControl = caret::trainControl(method = "cv", number = cv_fold),
                        tuneLength = tune_length)
  
  # print("Best Elastic Net model has the following parameters")
  # print(model$bestTune)
  # print("Best Elastic Net model proposes the following coefficients")
  # print(coef(model$finalModel, model$bestTune$lambda))
  
  return(model)
}


train_model_glmnet_my <- function(x, y, alpha=0, output_file=""){
  ## Train a ridge or a lasso regularized regression model using glmnet
  
  # Fit the model and explore coefficients
  fit <- glmnet::glmnet(x, y, alpha = alpha)
  
  # Visualize coefficients
  jpeg(paste0(output_file, "glmnet_coef_a", alpha,".png"))
  plot(fit)
  dev.off()
  
  # Find the best lambda using cross-validation
  cv <- glmnet::cv.glmnet(x, y, alpha = alpha)
  
  cat("Regularized regression model with alpha =", alpha, "lambda = ", cv$lambda.min, "\n")
  
  # Visualize lambda search
  jpeg(paste0(output_file, "glmnet_lambda_a", alpha, ".png"))
  plot(cv)
  dev.off()
  
  # Fit the final model on the training data
  model <- glmnet::glmnet(x, y, alpha = 0, lambda = cv$lambda.min)
  
  return(model)
}


assess_regression_coef_my <- function(fit, lambda="lambda.1se"){
  ## Report coefficients, their sum, and sum of their squares for a given regression model
  ## Lambda = "lambda.min" or "lambda.1se"
  
  coef_fit <- stats::coef(fit, s = lambda)
  
  return(list("coef" = coef_fit, 
              "coef_sum" = sum(abs(coef_fit[-1])), 
              "coef_sum_sq" = sum(coef_fit[-1] ^ 2)))
}


df2model_my <- function(df){
  ## Transform tibble into predictor variables matrix and outcome matrix
  
  return(list("data" = model.matrix(outcome ~ ., df)[,-1], 
              "outcome" = df$outcome))
}


analyze_model_my <- function(model, analysis_type="elnet", lambda_type="lambda.1se"){
  ## Analyze performance of various linear regressions
  ## Analysis types: norm, ridge, lasso, elnet
  ## Strategy to select lambda in regularized linear regression: "lambda.1se" or "lambda.min"
  
  # Declare named lists to store model parameters: alpha, lambda, coefficients
  alpha <- NULL
  lambda <- NULL
  coef <- list()
  
  # Find coefficients from various linear regression fits
  if(analysis_type == "elnet"){
    alpha <- round(model$bestTune$alpha, 3)
    lambda <- round(model$bestTune$lambda, 3)
    coef <- assess_regression_coef_my(model$finalModel, model$bestTune$lambda)
  } else{
    coef <- assess_regression_coef_my(model, lambda_type)
    
    if(analysis_type == "ridge"){
      alpha <- 0
      lambda <- round(model$lambda, 3)
    } else if(analysis_type == "lasso"){
      alpha <- 1
      lambda <- round(model$lambda, 3)
    }
  }
  
  # Gather all coefficients in a structure for plotting
  coef_matrix <- as.matrix(coef$coef)
  colnames(coef_matrix) <- "coefficients"
  
  cat(analysis_type, ": alpha =", alpha, "; lambda =", lambda, "\n")
  
  coef_df <- tibble::as_tibble(coef_matrix, rownames = "id") %>%
    dplyr::rename("variable" = "id") %>%
    dplyr::mutate_if(is.numeric, round, 3) %>%
    dplyr::arrange(!!rlang::sym(colnames(coef_matrix)[1]))
  
  return(coef_df)
}


visualize_model_my <- function(coef_df, filename, feature_list, n=30){
  ## Write model results to file and visualize
  
  # Clean up coefficient vector
  coef_df <- coef_df %>%
    dplyr::filter(coefficients != 0) %>%
    dplyr::mutate(coef_abs = abs(coefficients),
                  IsLineage = grepl("Lineage|lineage", variable)) %>%
    dplyr::arrange(dplyr::desc(coef_abs))

  # Write coefficients and their sum to file
  readr::write_delim(coef_df, file=paste0(filename, "_coefficients.txt"), delim = "\t")
  
  # Select top 20 coefficients to plot
  coef_plot_df <- coef_df %>%
    dplyr::filter(variable %in% feature_list) %>%
    dplyr::top_n(n, coef_abs)

  cat("The result of linear regression is here:", paste0(filename, "_coefficients.txt"), "\n")
  pander::pandoc.table(coef_plot_df %>%
                         dplyr::select(-c(coef_abs)), style="simple")
  
  # Compare coefficients proposed by various linear regression fits
  p <- create_bar_plot_my(coef_plot_df, 
                          x_label="variable", y_label="coefficients", 
                          fill_label="IsLineage", position="dodge", 
                          filename=paste0(filename, "_bar_coefficients"), 
                          labels=c("Model", "Coefficients", "Linear regression coefficients"), 
                          reorder_x=TRUE)
  
  # Gather coefficient sums
  cat("Model totals:\n",
      "Coefficient total =", round(sum(abs(coef_df$coefficients)), 3), "\n",
      "Coefficient total, squared =", round(sqrt(sum(abs(coef_df$coefficients))), 3), "\n")
  
  return(coef_df %>%
           dplyr::filter(IsLineage == FALSE & variable != "(Intercept)"))
}


model_split_data_my <- function(df, row_name_col="DepMap_ID", fraction=0.7){
  ## Split a tibble randomly into a training and testing set
  
  # Create training set
  train_df <- df %>% 
    dplyr::sample_frac(fraction)
  
  # Create test set
  test_df <- dplyr::anti_join(df, train_df, by = row_name_col) 
  
  return(list("train" = train_df %>% tibble::column_to_rownames(row_name_col), 
              "test" = test_df %>% tibble::column_to_rownames(row_name_col)))
}


run_all_models_my <- function(data_model, analysis_types, output_loc=NULL, formula=NULL){
  ## Fit linear regression with various parameters to input data
  ## Analysis types: norm, ridge, lasso, elnet
  ## Input is a tibble where every column is a variable and every row is a feature.
  ## One the input columns should be called "outcome" which will be the default outcome variable
  
  # Variable to store all models
  model <- list()
  
  # Prepare data for glmnet
  data_glmnet <- df2model_my(data_model)
  
  for(a in analysis_types){
    if(a == "norm"){
      # Fit regular multivariate regression and report coefficients
      # If regular regression returns NA coefficients, variable in question is linearly related to the other variables
      model[[a]] <- stats::lm(outcome ~ ., data_model)
    } else if(a == "elnet"){
      # Regularized regression using caret
      model[[a]] <- train_model_caret_my(data_model, formula, cv_fold=5, tune_length=10)
    } else if(a == "ridge"){
      # Regularized ridge regression using glmnet
      model[[a]] <- train_model_glmnet_my(data_glmnet$data, data_glmnet$outcome, alpha=0, output_loc)
    } else if(a == "lasso"){  
      # Regularized lasso regression using glmnet
      model[[a]] <- train_model_glmnet_my(data_glmnet$data, data_glmnet$outcome, alpha=1, output_loc)
    }
  }
  
  return(model)
}
