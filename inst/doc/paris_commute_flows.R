## ---- include = FALSE---------------------------------------------------------
old_opt <- options(max.print = 200)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=TRUE, results=FALSE, message=FALSE---------------------------
library("spflow")
library("sf")
library("spdep")

## -----------------------------------------------------------------------------
data("paris10km_municipalities")
data("paris10km_commuteflows")

## -----------------------------------------------------------------------------
drop_area <- names(paris10km_municipalities) != "AREA"
plot(paris10km_municipalities[drop_area])

## -----------------------------------------------------------------------------
old_par <- par(mfrow = c(1, 3), mar = c(0,0,1,0))

mid_points <- suppressWarnings({
    st_point_on_surface(st_geometry(paris10km_municipalities))})

paris10km_nb <- list(
  "by_contiguity" = spdep::poly2nb(paris10km_municipalities),
  "by_distance" = spdep::dnearneigh(mid_points,d1 = 0, d2 = 5),
  "by_knn" = spdep::knn2nb(knearneigh(mid_points,3))
)

plot(st_geometry(paris10km_municipalities))
plot(paris10km_nb$by_contiguity, mid_points, add = T, col = rgb(0,0,0,alpha=0.5))
title("Contiguity") 

plot(st_geometry(paris10km_municipalities))
plot(paris10km_nb$by_distance,mid_points, add = T, col = rgb(0,0,0,alpha=0.5)) 
title("Distance") 

plot(st_geometry(paris10km_municipalities))
plot(paris10km_nb$by_knn, mid_points, add = T, col = rgb(0,0,0,alpha=0.5))
title("3 Nearest Neighbors") 

par(old_par)

## -----------------------------------------------------------------------------
paris10km_commuteflows 

## -----------------------------------------------------------------------------
paris10km_net <- 
  sp_network_nodes(
    network_id = "paris",
    node_neighborhood = nb2mat(paris10km_nb$by_contiguity),
    node_data = st_drop_geometry(paris10km_municipalities),
    node_key_column = "ID_MUN")

paris10km_net

## -----------------------------------------------------------------------------
paris10km_net_pairs <- sp_network_pair(
  orig_net_id = "paris",
  dest_net_id = "paris",
  pair_data = paris10km_commuteflows,
  orig_key_column = "ID_ORIG",
  dest_key_column = "ID_DEST")

paris10km_net_pairs

## -----------------------------------------------------------------------------
paris10km_multi_net <- sp_multi_network(paris10km_net,paris10km_net_pairs)
paris10km_multi_net

## -----------------------------------------------------------------------------
results_default <- spflow(
  flow_formula = log(1 + COMMUTE_FLOW) ~ . + G_(log( 1 + DISTANCE)),
  sp_multi_network = paris10km_multi_net)

results_default

## -----------------------------------------------------------------------------
clog <- function(x) {
  log_x <- log(x)
  log_x - mean(log_x)
}

flow_formula  <- 
  log(COMMUTE_FLOW + 1) ~
  D_(log(NB_COMPANY) + clog(MED_INCOME)) +
  O_(log(POPULATION) + clog(MED_INCOME)) +
  I_(log(POPULATION)) +
  G_(log(DISTANCE + 1))

results_mle  <- spflow(
  flow_formula,
  paris10km_multi_net)
results_mle

## -----------------------------------------------------------------------------
sdm_formula <- ~
  O_(log(POPULATION) + clog(MED_INCOME)) +
  D_(log(NB_COMPANY) + clog(MED_INCOME))

cntrl <- spflow_control(
  estimation_method = "mcmc",
  sdm_variables = sdm_formula,
  model = "model_7") # restricts \rho_w = 0

results_mcmc  <- spflow(
  flow_formula,
  paris10km_multi_net,
  flow_control = cntrl)

results_mcmc

## ----"cleanup", include=FALSE-------------------------------------------------
options(old_opt)

