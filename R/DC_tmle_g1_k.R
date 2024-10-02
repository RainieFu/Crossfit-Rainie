#' Estimate Average Treatment Effect (ATE) using the TMLE estimator with cross-fit algorithm (generalization 1)
#'
#' This function estimates the Average Treatment Effect (ATE) by applying the Targeted Maximum Likelihood Estimator (TMLE) using a cross-fit approach. The package provides tools for data simulation and model fitting.
#'
#' @param data A data frame or tibble containing the dataset.
#' @param exposure Name of the exposure variable (treatment indicator).
#' @param outcome Name of the outcome variable.
#' @param covarsT A vector of names of covariates for the treatment model.
#' @param covarsO A vector of names of covariates for the outcome model.
#' @param family.y The family for the outcome model. It can be `binomial()` (default) or `"gaussian"`.
#' @param learners Similar to \code{SL.library()} in the `SuperLearner` package, specifying the learners to be used.
#' @param control Similar to \code{cvControl()} in the `SuperLearner` package, specifying the control parameters.
#' @param num_cf Number of repetitions. The default is 5.
#' @param n_split Number of splits used. The default is `n_split = 3`.
#' @param rand_split Logical value. If `FALSE` (default), discordant splits for exposure and outcome models are chosen systematically; otherwise, they are chosen randomly.
#' @param gbounds A numeric vector of length 2 specifying the bounds for truncation of predicted probabilities of exposure. The defaults are 5/sqrt(n)log(n) and 1-5/sqrt(n)log(n). See \code{tmle::tmle()} for more information.
#' @param Qbounds A numeric vector of length 2 specifying the bounds for predicted probabilities of outcomes. The defaults are 5e-04 and 1-5e-04.
#' @param seed Numeric value for reproducibility of splits.
#' @param conf.level Confidence level for confidence intervals. Default is 0.95.
#' @param stat Name of the summary statistics, "median" (default) or "mean", to calculate the overall ATE from the repetition-specific estimates.
#'
#' @return A tibble containing risk (or mean) difference (`ATE`), standard error (`se`), lower and upper confidence intervals (`lower.ci` and `upper.ci`, respectively).
#'
#' @section Data Generation:
#' The package includes a function, \code{gen.data}, to generate synthetic datasets for simulation studies. The function creates covariates and treatment indicators, and can generate counterfactual outcomes based on user-defined parameters.
#'
#' @references
#' Kang, J. D. Y., & Schafer, J. L. (2007). Demystifying Double Robustness: A Comparison of Alternative Strategies for Estimating a Population Mean from Incomplete Data. \emph{Statistical Science}, 22(4), 523-539.
#'
#' @import dplyr tibble tidyr purrr furrr
#'
#' @export
#'
#' @examples
#'
#' # Generate a synthetic dataset with 1200 observations
#' data <- gen.data(n = 1200, do.misp = TRUE, my.transform = TRUE, seed = 234)
#' head(data)
#'
#' # Fit the TMLE model to the generated data
#' fit_tmle_g1 <- DC_tmle_g1_k(data = data, exposure = "X", outcome = "Y",
#'                             covarsT = c("C1", "C2", "C3", "C4"),
#'                             covarsO = c("C1", "C2", "C3", "C4"),
#'                             family.y = "gaussian",
#'                             learners = c("SL.glm", "SL.mean"),
#'                             control = list(V = 3, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
#'                             num_cf = 5, n_split = 3, rand_split = FALSE,
#'                             Qbounds = 0.0025, gbounds = 0.0025, seed = 236,
#'                             conf.level = 0.95, stat = "median")
#'
#' # Display the results
#' print(fit_tmle_g1)
#'
DC_tmle_g1_k <- function(data,
                         exposure,
                         outcome,
                         covarsT,
                         covarsO,
                         family.y="binomial",
                         learners,
                         control,
                         n_split,
                         num_cf,
                         rand_split=FALSE,
                         gbounds = NULL,
                         Qbounds = 5e-04,
                         seed=146,
                         conf.level=0.95,
                         stat = "median"){

  tryCatch({
    print("Function started!")

    runs <- list()
    print("Setting seed and initializing variables")
    set.seed(seed)
    cf_seed = sample(num_cf)
    n = nrow(data)
    if(is.null(gbounds)) {
      print(paste("n =", n))
      print(paste("class of n:", class(n)))
      gbounds = 5/sqrt(n)/log(n)
      print(paste("gbounds =", gbounds))
    }

    print("Starting main loop")
    debug_info <- list()
    for(cf in 1:num_cf){
      print(paste("Processing split", cf, "of", num_cf))
      seed1 = cf_seed[cf]
      fit_result = tryCatch({
        print(nrow(data))
        tmle_single_g1_p(data,
                         exposure,
                         outcome,
                         covarsT,
                         covarsO,
                         family.y,
                         learners,
                         control,
                         n_split,
                         rand_split,
                         gbounds,
                         Qbounds,
                         seed=seed1)
      }, error = function(e) {
        print(paste("Error in tmle_single_g1_p:", e$message))
        print("Stack trace:")
        print(sys.calls())
        return(NULL)
      })

      if (is.null(fit_result)) {
        fit_sngle <-  data.frame(rd=NA, var = NA)
      } else {
        fit_sngle <- fit_result
      }
      runs[[cf]] <- fit_sngle
      debug_info[[cf]] <- list(seed = seed1, fit_result = fit_result)

      print(paste("Completed split", cf, "Result:", paste(names(fit_sngle), fit_sngle, sep="=", collapse=", ")))
    }

    print("Binding rows")
    res = dplyr::bind_rows(runs)

    print("Calculating statistics")
    if(stat == "mean"){
      medians <- apply(res, 2, mean, na.rm = TRUE)
      res <- res %>% mutate(var0 = var + (rd - medians[1])^2)
      results <- apply(res, 2, mean, na.rm = TRUE)
    }
    if(stat == "median"){
      medians <- apply(res, 2, median, na.rm = TRUE)
      res <- res %>% mutate(var0 = var + (rd - medians[1])^2)
      results <- apply(res, 2, median, na.rm = TRUE)
    }

    print("Calculating t-value")
    print(paste("conf.level:", conf.level))
    print(paste("nrow(data):", nrow(data)))
    print(paste("Class of nrow(data):", class(nrow(data))))

    if (!is.numeric(conf.level) || !is.numeric(nrow(data))) {
      stop("conf.level or nrow(data) is not numeric")
    }

    t.value = qt((1-conf.level)/2, nrow(data), lower.tail = F)
    print(paste("t.value:", t.value))

    print("Calculating confidence intervals")
    print(paste("results[1]:", results[1]))
    print(paste("results[3]:", results[3]))

    if (!is.numeric(results[1]) || !is.numeric(results[3])) {
      stop("results[1] or results[3] is not numeric")
    }

    l_ci = results[1] - t.value*sqrt(results[3])
    u_ci = results[1] + t.value*sqrt(results[3])

    res1 = tibble(ATE=results[1], se = sqrt(results[3]), lower.ci = l_ci, upper.ci = u_ci)

    print("Returning results")
    return(list(result = res1, debug_info = debug_info))
  }, error = function(e) {
    message <- paste("Error in DC_tmle_g1_k:", e$message)
    print(message)
    print(paste("Error occurred at:", e$call))
    # print("Stack trace:")
    # print(sys.calls())
    stop(message)
  })
}
