library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(furrr)
library(gam)
library(xgboost)
library(SuperLearner)

fit_tmle_g1 <- DC_tmle_g1_k(data=ObsData,
                           exposure="X",
                           outcome="Y",
                           covarsT=c("C1", "C2", "C3", "C4"),
                           covarsO=c("C1", "C2", "C3", "C4"),
                           family.y="gaussian",
                           learners = c("SL.glm", "SL.mean"),
                           # learners=c("SL.glm.dcTMLE", "SL.mean.dcTMLE", "SL.gam4.dcTMLE", "SL.xgboost.dcTMLE"),
                           control=list(V = 3, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                           num_cf=5,
                           n_split=3,
                           rand_split=FALSE,
                           Qbounds = 0.0025,
                           gbounds = 0.0025,
                           seed=236,
                           conf.level=0.95,
                           stat = "median")

fit_tmle_g1

fit_tmle_g2 <- DC_tmle_g2_k(data=ObsData,
                            exposure="X",
                            outcome="Y",
                            covarsT=c("C1", "C2", "C3", "C4"),
                            covarsO=c("C1", "C2", "C3", "C4"),
                            family.y="gaussian",
                            learners = c("SL.glm", "SL.mean"),
                            # learners=c("SL.glm.dcTMLE", "SL.mean.dcTMLE", "SL.gam4.dcTMLE", "SL.xgboost.dcTMLE"),
                            control=list(V = 3, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                            num_cf=5,
                            n_split=3,
                            rand_split=FALSE,
                            Qbounds = 0.0025,
                            gbounds = 0.0025,
                            seed=236,
                            conf.level=0.95,
                            stat = "median")

fit_tmle_g2

