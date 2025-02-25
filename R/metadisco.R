#' @export
setClass("Metadisco", slots = list(
  details = "list",
  group1_spomics = "list",
  group2_spomics = "list",
  results = "list"))

#' @export
create_metadisco <- function(group1_name, group1_spomics, group2_name, group2_spomics, process_spomic_lists = FALSE) {
  # Later, come back and add options for different types of inputs for group1_spomics and group2_spomics
  # For example, instead of a list of spomics, it would be nice to provide a vector of file paths

  if(process_spomic_lists) {

  }

  metadisco <- new("Metadisco",
                   details = list(),
                   group1_spomics = group1_spomics,
                   group2_spomics = group2_spomics,
                   results = list()
  )

  metadisco@details$group1_name <- group1_name
  metadisco@details$group2_name <- group2_name
  metadisco@details$hyperparameters <- list()

  return(metadisco)
}

#' @export
set_metadisco_hyperparameters <- function(metadisco, r, colocalization_type = "crossK.inhom", tau_estimator = "SJ") {
  metadisco@details$hyperparameters$r <- r
  metadisco@details$hyperparameters$colocalization_type <- colocalization_type
  metadisco@details$hyperparameters$tau_estimator <- tau_estimator

  for(i in seq_along(metadisco@group1_spomics)) {
    metadisco@group1_spomics[[i]] <- spomic::set_spomic_hyperparameters(spomic = metadisco@group1_spomics[[i]],
                                                                r = r,
                                                                colocalization_type = colocalization_type)
  }
  for(i in seq_along(metadisco@group2_spomics)) {
    metadisco@group2_spomics[[i]] <- spomic::set_spomic_hyperparameters(spomic = metadisco@group2_spomics[[i]],
                                                                r = r,
                                                                colocalization_type = colocalization_type)
  }

  return(metadisco)
}

#' @export
determine_pairs_to_test <- function(metadisco, nonrandom_consistency) {
  # nonrandom_consistency is a threshold representing the percentage of samples in
  # a group that must have nonrandom crossK values for a given cell pair

  group1_nonrandom_pairs <- c()
  for(i in seq_along(metadisco@group1_spomics)) {
    group1_nonrandom_pairs <- c(group1_nonrandom_pairs, metadisco@group1_spomics[[i]]@results$nonrandom_pairs)
  }
  group1_consistency_table <- table(group1_nonrandom_pairs)/length(metadisco@group1_spomics)
  group1_consistency_table <- group1_consistency_table[group1_consistency_table >= nonrandom_consistency]
  group1_pairs <- names(group1_consistency_table)

  group2_nonrandom_pairs <- c()
  for(i in seq_along(metadisco@group2_spomics)) {
    group2_nonrandom_pairs <- c(group2_nonrandom_pairs, metadisco@group2_spomics[[i]]@results$nonrandom_pairs)
  }
  group2_consistency_table <- table(group2_nonrandom_pairs)/length(metadisco@group2_spomics)
  group2_consistency_table <- group2_consistency_table[group2_consistency_table >= nonrandom_consistency]
  group2_pairs <- names(group2_consistency_table)

  pairs_to_test <- intersect(group1_pairs, group2_pairs)

  metadisco@details$hyperparameters$nonrandom_consistency <- nonrandom_consistency
  metadisco@results$pairs_to_test <- pairs_to_test

  return(metadisco)
}

#' @export
run_random_effects_meta_analysis <- function(kcross_distributions, method) {
  cell_pairs <- unique(kcross_distributions$i_j)
  pooled_estimates <- list()
  for(pair in cell_pairs) {
    yi <- kcross_distributions |> dplyr::filter(i_j == pair) |> dplyr::pull(kcross)
    vi <- kcross_distributions |> dplyr::filter(i_j == pair) |> dplyr::pull(kcross_var)
    slab <- kcross_distributions |> dplyr::filter(i_j == pair) |> dplyr::pull(sample)
    model <- metafor::rma.uni(yi = yi, vi = vi, slab = slab, method = method)
    # metafor::forest(model)

    coefficients <- coef(summary(model))
    pooled_estimates[[pair]] <- data.frame(i_j = pair,
                                           colocalization_estimate = coefficients$estimate,
                                           # colocalization_variance = length(yi)*coefficients$se^2,
                                           colocalization_variance = coefficients$se^2,
                                           colocalization_se = coefficients$se)
  }
  return(dplyr::bind_rows(pooled_estimates))
}

#' @export
run_metadisco <- function(metadisco, nonrandom_consistency=0.5) {
  metadisco <- determine_pairs_to_test(metadisco, nonrandom_consistency = nonrandom_consistency)
  # Start by running the loh bootstrapping on the relevant cell_pairs
  for(pair in metadisco@results$pairs_to_test) {
    # print(pair)
    split_pair <- strsplit(pair, "_")[[1]]
    i <- split_pair[1]
    j <- split_pair[2]

    for(a in seq_along(metadisco@group1_spomics)) {
      metadisco@group1_spomics[[a]] <- spomic::get_kcross(spomic = metadisco@group1_spomics[[a]], i = i, j = j)
    }
    for(b in seq_along(metadisco@group2_spomics)) {
      metadisco@group2_spomics[[b]] <- spomic::get_kcross(spomic = metadisco@group2_spomics[[b]], i = i, j = j)
    }
  }

  group1_kcross_distributions <- list()
  for(a in seq_along(metadisco@group1_spomics)) {
    group1_kcross_distributions[[a]] <- spomic::get_kcross_summary(metadisco@group1_spomics[[a]])
  }
  group1_kcross_distributions <- dplyr::bind_rows(group1_kcross_distributions)

  group2_kcross_distributions <- list()
  for(b in seq_along(metadisco@group2_spomics)) {
    group2_kcross_distributions[[b]] <- spomic::get_kcross_summary(metadisco@group2_spomics[[b]])
  }
  group2_kcross_distributions <- dplyr::bind_rows(group2_kcross_distributions)

  group1_pooled_estimates <- run_random_effects_meta_analysis(kcross_distributions = group1_kcross_distributions,
                                                              method = metadisco@details$hyperparameters$tau_estimator)
  group2_pooled_estimates <- run_random_effects_meta_analysis(kcross_distributions = group2_kcross_distributions,
                                                              method = metadisco@details$hyperparameters$tau_estimator)

  pooled_estimates <- dplyr::inner_join(group1_pooled_estimates, group2_pooled_estimates, by = "i_j", suffix = c("_group1", "_group2"))

  differential_testing <- pooled_estimates |> dplyr::mutate(log2fc = log2(colocalization_estimate_group2 + 1) - log2(colocalization_estimate_group1 + 1),
                                                            z_score = (colocalization_estimate_group2 - colocalization_estimate_group1)/sqrt(colocalization_se_group1^2 + colocalization_se_group2^2),
                                                            pval = (1-pnorm(abs(z_score), mean = 0, sd = 1, lower.tail = TRUE)) * 2,
                                                            FDR = p.adjust(pval, method = "fdr"),
                                                            holm = p.adjust(pval, method = "holm"))
  metadisco@results$differential_testing <- differential_testing

  return(metadisco)
}


