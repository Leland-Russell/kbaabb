#' @title Longitudinal KBAABB Populations
#'
#' @description \code{longitudinal_kbaabb()} generates an imputed population 
#' dataset based on methodology introduced in White et al. (2025) for each value
#' of some indexing variable, intending to mimic repeated observations of plots
#'
#' @param survey_data A dataframe containing the survey data to be used for the 
#' imputation (i.e. "the sample"). This dataframe must contain all variables 
#' listed in \code{formula} and \code{strata}.
#' @param population_data A dataframe containing the population data for
#' imputation to occur on (i.e. "the population"). This dataframe must contain 
#' all auxiliary variables listed in \code{formula} and \code{strata}. 
#' @param formula The formula specified for imputation, taking the form 
#' \code{y ~ x_1 + x_2 + ... + x_n} where \{x_1, x_2, ..., x_n\} is the set of 
#' auxiliary variables used for imputation and y is the response. 
#' @param longitudinal_idx name of variable indexing multiple observations of plots. A 
#' separate KBAABB population will be made for each unique value of this variable.
#' @param k Integer. The number of neighbors used in the \code{k} nearest neighbors
#' imputation
#' @param strata Character. The name of a variable to be used for 
#' stratification. If \code{NULL} (default), imputation is performed without
#' stratification. Otherwise, stratification occurs based on the variable
#' specified in \code{strata}. It is advised to provide strata.
#' @param center_scale Logical. If \code{TRUE} (default), auxiliary variables are
#' centered and scaled (mean = 0, variance = 1) based on the population data. 
#' Otherwise, the original sample and population dataframes supplied by the user
#' are used in an unmodified form. 
#' @param seed A seed to be set for reproducibility.
#' @param keep_survey_id Logical, should the rownames of the survey data which
#' correspond to the imputed values be kept? Allows for a lookup of which rows
#' were used for which imputed observations.
#' @param id_name The name for the column containing the rownames in 
#' \code{keep_survey_id}.
#' @param progress_bar logical, should a progress bar be displayed
#' @param ... Currently ignored. For extendability.
#'
#' @return A \code{longitudinal_kbaabb} object including:
#' \itemize{
#'  \item \code{imputed_population_data}: The imputed population
#'  \item \code{imputed_population_data_wide}: The imputed population in wide form
#'  \item \code{survey_data}: Original sample data
#'  \item \code{longitudinal_idx}: The variable used to index longitudinal populations
#'  \item \code{k}: Value for \code{k}
#'  \item \code{stratified}: Indicator of if stratifying occured, as well as
#'  \item \code{strata}: The stratifying variable
#'  \item \code{center_scale}: If centering and scaling occured
#'  \item \code{formula}: The formula used for population imputation
#' }
#' 
#' @seealso \code{\link[kbaabb]{kbaabb}}, \code{\link[kbaabb]{smooth_longitudinal_kbaabb}}
#' 
#' @examples
#' 
#' pop = SJC_population[1:1000,c("island", "tcc", "elev")]
#' pop = dplyr::summarize(
#'   pop, `count` = dplyr::n(), .by = colnames(pop)
#' )
#' pop = pop[,c("island", "tcc", "elev")]
#' for (i in 1:4) {
#'   pop = rbind(pop, pop)
#' }
#' pop$id = rep(1:(2^4), each = nrow(pop)/2^4)
#' 
#' longitudinal_kbaabb(
#'   SJC_sample,
#'   pop,
#'   biomass ~ tcc + elev,
#'   longitudinal_idx = "id",
#'   strata = "island"
#' )
#' 
#' @export
longitudinal_kbaabb <- function(
    survey_data, # dataframe (to be coerced into a matrix)
    population_data, # dataframe (to be coerced into a matrix)
    formula, # formula
    longitudinal_idx, # variable name indexing the creation of multiple KBAABB populations
    k = 10, # positive integer; nearest neighbors
    strata = NULL, # NULL or character 
    center_scale = TRUE, # logical
    seed = NULL,  # numeric
    keep_survey_id = TRUE,
    id_name = "survey_row_id",
    progress_bar = TRUE, # logical; do you want a progress bar?
    ...
) {
  validate_parameters(
    survey_data,
    population_data,
    formula,
    k,
    strata,
    center_scale,
    keep_survey_id,
    id_name
  )
  
  if (!is.null(seed)) {
    # we set the seed here instead of in the KBAABB populations
    set.seed(seed)
  }
  
  covariates = all.vars(formula[[3]])
  if (!is.null(strata)) {
    covariates = union(covariates, strata)
  }
  covariates = union(covariates, longitudinal_idx)
  
  validate_longitudinal_params(
    population_data,
    longitudinal_idx,
    covariates,
    formula
  )
  
  pop_vals = unique(population_data[[longitudinal_idx]])
  message("Making Populations")
  if (progress_bar) {
    pb = utils::txtProgressBar(
      min = 0, max = length(pop_vals), title = "Making KBAABB populations"
    )
  }
  
  NA_mask = apply(
    population_data[,covariates], 1, 
    function (row) {
      any(is.na(row))
    }
  )
  
  pops = list()
  for (i in 1:length(pop_vals)) {
    theval = pop_vals[i]
    
    tmp_population_data = population_data 
    tmp_mask = tmp_population_data[[longitudinal_idx]] == theval
    tmp_mask = tmp_mask & !NA_mask
    tmp_population_data = tmp_population_data[tmp_mask,]
    
    # not ideal warning handling here, but kbaabb already gives a lot of warnings
    # as is, if we generate multiple populations, it will be even more
    # think more about this...
    if (i == 1) {
      wrapper = function (x) {x}
    } else {
      wrapper = function (x) {suppressWarnings(x)}
    }
    
    imputed_data = wrapper(kbaabb(
      survey_data = survey_data,
      population_data = tmp_population_data,
      formula = formula,
      strata = strata,
      k = k,
      center_scale = center_scale,
      seed = NULL,
      keep_survey_id = keep_survey_id,
      id_name = id_name
    ))
    pops[[i]] = imputed_data$imputed_population_data
    rm(imputed_data) # these can get kind of large for smaller computers
    if (progress_bar) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  # we will do a little divide and conquer for merging these, just because
  # for many large populations, it can take quite a while to do it naively
  # although I'm just going to break them into bins of four to do it. This 
  # could become an argument, but is not necessary
  
  tree_depth = 0
  n_pops = length(pop_vals)
  while (n_pops > 1) {
    n_pops = ceiling(n_pops/4)
    tree_depth = tree_depth + 1
  }
  
  num_rows = 0
  for (i in 1:length(pops)) {
    num_rows = num_rows + nrow(pops[[i]])
  }

  if (progress_bar) {
    close(pb)
  }
  message("Merging Populations")
  # since we are doing D&C, it won't be entirely accurate to display what percentage
  # of the way done we are with merges, but rather how much work we do. And we do 
  # proportionally more work merging larger things, so the progress bar will show the
  # percentage of work we've done
  if (progress_bar) {
    pb = utils::txtProgressBar(
      min = 0, max = tree_depth * num_rows, title = "Merging KBAABB populations"
    )
    pb_counter = 0
  }
  n_pops = length(pop_vals)
  
  merge_window = 4
  while (length(pops) > 1) {
    old_pops = pops
    pops = list()
    n_pops = ceiling(n_pops/4)
    for (i in 1:n_pops) {
      obj_name_list = c()
      idx = merge_window * (i-1) + 1 # starting index of the list
      for (j in idx:min(idx + merge_window - 1, length(old_pops))) {
        obj_name_list = union(obj_name_list, paste0("old_pops[[", j, "]]"))
      }
      
      args = paste0(obj_name_list, collapse = ", ")
      pops[[i]] = eval(parse(text = paste0(
        "rbind(", args, ")"
      )))
      if (progress_bar) {
        pb_counter = pb_counter + 
          nrow(pops[[i]])
        
        utils::setTxtProgressBar(pb, pb_counter)
      }
    }
  }
  if (progress_bar) {
    close(pb)
  }
  
  merged_data = pops[[1]]
  rm(pops, old_pops)
  
  warning(
    "These are unsmoothed populations. ",
    "Use smooth_longitudinal_kbaabb() to smooth them"
  )
  
  if (keep_survey_id) {
    # the survey ids will mess with the join, so remove the id name for wide form
    # it will still exist in long form
    merged_data = merged_data[,-which(colnames(merged_data) == id_name)]
  } 
  wide_pop_data = tidyr::pivot_wider(
    merged_data, names_from = dplyr::all_of(longitudinal_idx), values_from = formula[[2]]
  )
  
  ret = list(
    imputed_population_data = merged_data,
    imputed_population_data_wide = wide_pop_data,
    survey_data = survey_data,
    longitudinal_idx = longitudinal_idx,
    k = k,
    strata = strata,
    center_scale = center_scale,
    formula = formula,
    keep_survey_id = keep_survey_id,
    id_name = id_name
  )
  
  class(ret) <- "longitudinal_kbaabb"
  return (ret)
}

#' @title Smooth longitudinal KBAABB populations
#' 
#' @description
#' this can only be used if the variable indexing multiple populations is 
#' numeric, or can be cast to a numeric in a sensible way
#' 
#' @param longitudinal_kbaabb an object of class \code{longitudinal_kbaabb}
#' @param penalty_fn a function of two numeric arguments, the second of which can
#' accept vector arguments, defining the weights used for a weighted average smooth
#' @param smooth_radius numeric, radius about a given point to compute weighted averages
#' @param dry_run logical, should the number of unique units considered be subset? 
#' This may help with computation times while tuning \code{penalty_fn} or \code{smooth_radius}
#' @param dry_run_subset how many unique units should be considered? All longitudinal observations
#' from this unit will be preserved
#' @param compute_acf logical, should a plot showing the autocorrelation of the 
#' smoothed and unsmoothed observations be created?
#' @param forecast numeric, forecast distance for computing autocorrelation
#' 
#' @returns If \code{compute_acf} is \code{TRUE}, then a list including
#' \itemize{
#'  \item \code{kbaabb}: A \code{longitudinal_kbaabb} object with smoothed observations
#'  replacing unsmoothed observations in \code{kbaabb$imputed_population_data} and
#'  \code{kbaabb$imputed_population_data_wide}
#'  \item \code{plot}: A \code{ggplot2} plot showing the estimated autocorrelation
#'  function values lags dictated by \code{forecast}.
#' }
#' otherwise an object of class \code{longitudinal_kbaabb}, as described above
#' 
#' @examples
#' pop = SJC_population[1:1000,c("island", "tcc", "elev")]
#' pop = dplyr::summarize(
#'   pop, `count` = dplyr::n(), .by = colnames(pop)
#' )
#' pop = pop[,c("island", "tcc", "elev")]
#' for (i in 1:4) {
#'   pop = rbind(pop, pop)
#' }
#' pop$id = rep(1:(2^4), each = nrow(pop)/2^4)
#' 
#' long_kbaabb = longitudinal_kbaabb(
#'   SJC_sample,
#'   pop,
#'   biomass ~ tcc + elev,
#'   longitudinal_idx = "id",
#'   strata = "island"
#' )
#' 
#' smooth_longitudinal_kbaabb(
#'   long_kbaabb
#' )
#' 
#' @export
smooth_longitudinal_kbaabb <- function (
  longitudinal_kbaabb,
  penalty_fn = function (x1, x2) {1/(1 + abs(x1 - x2))}, # penalty function for spatial smoothing
  smooth_radius = Inf, # radius for spatial smoothing
  dry_run = FALSE, # subset the data for tuning?
  dry_run_subset = 1e4, # how much to subset data
  compute_acf = TRUE, # compute autocorrelation?
  forecast = Inf # maximum order of autocorrelation
) {
  check_smooth_inputs(
    longitudinal_kbaabb,
    penalty_fn, # penalty function for spatial smoothing
    smooth_radius, # radius for spatial smoothing
    dry_run,
    dry_run_subset,
    compute_acf,
    forecast
  )
  
  wide_data = longitudinal_kbaabb$imputed_population_data_wide
  og_long_data = longitudinal_kbaabb$imputed_population_data
  idx = longitudinal_kbaabb$longitudinal_idx
  
  common_aux = intersect(
    colnames(og_long_data),
    colnames(wide_data)
  )
  
  remaining_vars = setdiff(
    colnames(wide_data), common_aux
  )
  
  res_name = longitudinal_kbaabb$formula[[2]]
  smooth_var_name = paste0(res_name, "_smooth")
  
  smooth_obs = knn_smoother(
    wide_data[,remaining_vars],
    penalty = penalty_fn,
    window_radius = smooth_radius
  )
  
  wide_data[,remaining_vars] = smooth_obs
  long_data = tidyr::pivot_longer(
    wide_data, tidyr::all_of(remaining_vars),
    names_to = idx,
    values_to = smooth_var_name
  )
  
  long_data[[idx]] = as.numeric(
    long_data[[idx]]
  )
  
  if (!is.numeric(
    og_long_data[[idx]]
  )) {
    og_long_data[[idx]] = as.numeric(
      og_long_data[[idx]]
    )
  }
  
  long_data = dplyr::left_join(
    og_long_data,
    long_data,
    by = c(common_aux, idx)
  )
  
  response_vec = longitudinal_kbaabb$survey_data[[res_name]]
  response_vec = sort(response_vec)
  print(response_vec)
  
  long_data = 
    replace_values(long_data, response_vec, preds_name = res_name)

  wide_data = tidyr::pivot_wider(
    long_data[,-which(colnames(long_data) == res_name)], 
    names_from = tidyr::all_of(idx), 
    values_from = tidyr::all_of(smooth_var_name)
  )
  
  if (compute_acf) {
    acf_smooth_kbaabb = compute_temporal_acf(
      wide_data[,remaining_vars], forecast = forecast
    )
    og_acf_smooth_kbaabb = compute_temporal_acf(
      longitudinal_kbaabb$imputed_population_data_wide[,remaining_vars], forecast = forecast
    )
    
    acf_smooth_kbaabb = data.frame(
      Autocorrelation = acf_smooth_kbaabb, 
      lag = 1:length(acf_smooth_kbaabb)-1
    )

    og_acf_smooth_kbaabb = data.frame(
      Autocorrelation = og_acf_smooth_kbaabb, 
      lag = 1:length(og_acf_smooth_kbaabb)-1
    )
    
    plot = ggplot2::ggplot() + 
      ggplot2::geom_hline(yintercept = 0, color = "grey") +
      ggplot2::geom_point(
        ggplot2::aes(x = acf_smooth_kbaabb$lag, y = acf_smooth_kbaabb$Autocorrelation, 
                     color = "Smoothed Population")
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = og_acf_smooth_kbaabb$lag, y = og_acf_smooth_kbaabb$Autocorrelation, 
                     color = "Unsmoothed Population")
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = paste0("Lags (", idx, "s)"), 
        y = "Autocorrelation",
        title = "Smoothed and Unsmoothed autocorrelation of repeated KBAABB populations"
      ) +
      ggplot2::scale_color_manual(
        values = c("red", "black"),
        name = "Data Type"
      ) +
      ggplot2::ylim(0,1)
    
    ret = list(kbaabb = longitudinal_kbaabb)
    ret$kbaabb$imputed_population_data = long_data
    ret$kbaabb$imputed_population_data_wide = wide_data
    ret$plot = plot
    return (ret)
  } else {
    ret = longitudinal_kbaabb
    ret$imputed_population_data = long_data
    ret$imputed_population_data_wide = wide_data
    return (ret)
  }
}
