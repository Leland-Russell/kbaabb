#' Validate initial parameters
#' @noRd
validate_parameters <- function(survey_data,
                                population_data,
                                formula,
                                k,
                                strata,
                                center_scale,
                                keep_survey_id,
                                id_name) {
  if (!is.numeric(k)) {
    stop(paste0("Must provide a numeric value for k."))
  }
  if (!(is.character(strata) || is.null(strata))) {
    stop(paste0("Must provide character or null for strata"))
  } else if (!is.null(strata)) {
    if (length(strata) > 1) {
      stop(
        paste0("Can only supply one stratifying variable.",
               " If you want to use multiple stratifying variables, consider",
               " using pivot_longer().")
      )
    } else if (!(strata %in% colnames(survey_data))) {
      stop(
        paste0("Supplied strata is not in your donor dataset. Please ensure that",
               " your data include the stratifying variable.")
      )
    } else if (!(strata %in% colnames(population_data))) {
      stop(
        paste0("Supplied strata not in your recieving dataset. Please ensure that",
               " your data include the stratifying variable.")
      )
    }
  }
  if (!center_scale) {
    warning(
      paste0("You have chosen not to center or scale your variables. This may",
             " lead to biased KNN draws.")
    )
  }
  if (!inherits(formula, "formula")) {
    tryCatch({formula = stats::as.formula(formula)},
             error = function(e) {
               stop("formula ", formula, " is not a formula and ",
                    "cannot be coerced into a formula.\n", e)
             })
  }
  
  if (!is.character(id_name)) {
    stop("id_name must be a character")
  }
}

#' Validate covariates/response
#' @noRd
validate_variables = function(
    x_vars,
    y_vars,
    survey_data,
    population_data,
    keep_survey_id,
    id_name
) {
  if (length(x_vars) == 0) {
    stop(paste("Supplied no auxiliaries variables. Please supply auxiliary",
               "variables."))
  }
  if (length(y_vars) == 0) {
    stop(paste("Supplied no response variables. Cannot impute data without",
               "a response variable."))
  } 

  temp_setdiff_varnames = setdiff(y_vars, colnames(survey_data))
  if (length(setdiff(y_vars, colnames(survey_data))) != 0) {
    stop("Response variable(s) not in survey data.")
  }
  
  for (vars in list(
    x_vars,
    y_vars
  )) {
    temp_setdiff_varnames = setdiff(vars, colnames(survey_data))
    if (length(temp_setdiff_varnames) > 0) {
      stop(paste("Variable", temp_setdiff_varnames, "is not in survey data."))
    }
  }
  
  temp_intersect_varnames = intersect(y_vars, colnames(population_data))
  if (length(temp_intersect_varnames) != 0) {
    warning("Variable", ifelse(length(y_vars) > 1, "s", ""), " ",
            temp_intersect_varnames,
            ifelse(length(y_vars) > 1, " are", " is"),
            " already present in the population data. Imputing ",
            ifelse(length(y_vars) > 1, "them", "it"), " anyway.")
  }
  
  if (keep_survey_id && (id_name %in% union(x_vars, y_vars))) {
    stop("id_name must be unique, but you already have a column named ",
         id_name, ". Please rename id_name, or rename this column.")
  }
}

#' Validate the parameters for just longitudinal KBAABB creation
#' @noRd
validate_longitudinal_params <- function (
    population_data,
    longitudinal_idx,
    covariates,
    formula
) {
  if (!longitudinal_idx %in% colnames(population_data)) {
    stop("Variable `longitudinal_idx` must be a column in population data.")
  }
  
  NA_mask = apply(
    population_data[,covariates], 1, 
    function (row) {
      any(is.na(row))
    }
  )
  
  widening_attempt = population_data
  widening_attempt = widening_attempt[!NA_mask,]
  widening_attempt = dplyr::summarize(
    widening_attempt, `count` = dplyr::n(), .by = colnames(widening_attempt)
  )
  if (max(widening_attempt$count) > 1) {
    stop(
      "Observations are not be uniquely identified by auxiliary variables. ",
      "Use the following code to identify duplicate rows: \n",
      "{population_data} |>
  dplyr::summarise(n = dplyr::n(), .by = c(", paste0(covariates, collapse = ", "), 
      ")) |>
  dplyr::filter(n > 1L) "
    )
  }
  
  if (length(formula[[2]]) > 1) {
    stop("Cannot perform multiple data imputation longitudinally at this time.")
  }
}

#' Compute temporal autocorrelations of KBAABB
#' @noRd
compute_temporal_acf = function(wide_data, forecast = 20) {
  forecast = min(forecast, ncol(wide_data) + 1)
  acfs = apply(wide_data, 1, function (row) {
    row = as.numeric(row)
    if (length(unique(row)) == 1) {
      return (rep(1, forecast))
    }
    return (stats::acf(row, plot = F, lag.max = forecast-1)$acf[,,1])
  })
  
  acf_array = t(array(unlist(acfs), dim = c(forecast, length(acfs))))
  return (colMeans(acf_array))
}

#' Temporal smoother
#' @noRd
knn_smoother = function(wide_data, 
                        penalty = function (x1, x2) {1/(1 + abs(x1 - x2))},
                        window_radius = 10) {
  # nearest neighbor smoother :
  # y[i] gets replaced with its (window_radius) nearest neighbors,
  # and the relative weight nearby points are assigned in the value of y[i]
  # is dictated by penalty
  
  # This might be only slightly faster as the case below, but
  # it might be faster on longer observations (we tested on 26, could be better
  # on, say, 100 though)
  if (window_radius > ceiling((ncol(wide_data)-1)/2)) {
    smooth_kbaabb = apply(wide_data, 1, function (row) {
      y = as.numeric(row)
      n = length(y)
      y = t(array(rep(y, 16), dim = c(n,n)))
      penalties = sapply(1:n, function(i) {penalty(i, 1:n)})
      weights = y * penalties
      return (rowSums(weights)/rowSums(penalties))
    })
  } else {
    smooth_kbaabb = apply(wide_data, 1, function (row) {
      y = as.numeric(row)
      x = as.numeric(names(row))
      x = x - min(x) + 1 # so we have proper R indices here.
      
      new_y = y
      
      for (i in x) {
        window_idx = max(1,i-window_radius):min(length(y), i+window_radius)
        window = y[window_idx]
        penalties = penalty(i, window_idx)
        window = window * penalties
        new_y[i] = sum(window)/sum(penalties)
      }
      return (new_y)
    })
  }
  return (t(smooth_kbaabb))
}

#' Replace smoothed values with observed values in equal relative proportion
#' @noRd
replace_values = function(data, real_val, preds_name = "preds") {
  preds = data[[preds_name]]
  ecdf_smooth_data = stats::ecdf(preds)
  data[[preds_name]] = length(real_val) * ecdf_smooth_data(preds)
  
  impute_biomass = function(x) {
    x[x < 1] = 1
    x[x > length(real_val)] = length(real_val)
    return(real_val[round(x)])
  }
  
  data[[preds_name]] = impute_biomass(data[[preds_name]])
  return(data)
}

#' Check the inputs to the longitudinal smoothing function
#' @noRd
check_smooth_inputs = function(
  longitudinal_kbaabb,
  penalty_fn, # penalty function for spatial smoothing
  smooth_radius, # radius for spatial smoothing
  dry_run,
  dry_run_subset,
  compute_acf,
  forecast
) {
  if (!inherits(longitudinal_kbaabb, "longitudinal_kbaabb")) {
    stop("longitudinal_kbaabb must be an object of class longitudinal_kbaabb")
  }
  if (!inherits(penalty_fn, "function")) {
    stop("penalty_fn must be a function")
  }
  t = try(penalty_fn(1,1:10))
  if (inherits(t, "try-error")) {
    stop("Penalty function cannot compute weights with inputs: ", t)
  }
  
  t = try(as.logical(dry_run))
  if (inherits(t, "try-error") || is.na(t)) {
    stop(
      "`dry_run` must be a logical, or be easily coerced into a logical. Not", dry_run
    )
  }
  
  t = try(as.logical(compute_acf)) 
  if (inherits(t, "try-error") || is.na(t)) {
    stop(
      "`compute_acf` must be a logical, or easily coerced into a logical. Not ", compute_acf
    )
  }
  
  if (dry_run) {
    I = 2
  } else {
    I = 1
  }
  
  for (i in 1:I) {
    var = c(smooth_radius, dry_run_subset)[i]
    name = c("smooth_radius", "dry_run_subset")[i]
    if (!is.numeric(var)) {
      stop("`", name, "` must be a numeric. Not ", var)
    } else if (!is.infinite(var) && !(var %% 1)) {
      stop("`", name, "` must be an integer or Inf. Not ", var)
    } else if (var <= 0) {
      stop("`", name, "` must be a positive integer. Not ", var)
    }
  }
}
