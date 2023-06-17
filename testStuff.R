n=120; relevance=0.15
set.seed(20)
data = StroblData(n, relevance, yFactor = FALSE)
suppressWarnings(rf_model <- ranger(
  formula = y ~ ., 
  data = data,
  num.trees = 200,
  importance = "impurity",
  #mtry = 1,
  keep.inbag = T,
  classification = T
))

unified = ranger.unify(rf_model,data)
#inbag/oob treeshap calculation for smoothing
treeshap_obj = treewise_shap(unified,data,rf_inbag = rf_model$inbag.counts)
shap_vals_in_sim = treeshap_obj[[1]]$shaps[,-6]
shap_vals_oob_sim = treeshap_obj[[2]]$shaps[,-6]
shap_avs_oob_sim <- shap_vals_oob_sim %>% abs() %>% colMeans() 
shap_avs_in_sim <- shap_vals_in_sim %>% abs() %>% colMeans() 
###Raw treeshap calculation with original treeshap::treeshap function
shap_obj = treeshap(unified,data)
shap_vals_sim = shap_obj$shaps[,-6]
shap_avs_sim = shap_vals_sim %>% abs() %>% colMeans()

par(mfrow=c(3,2))
for (j in 1:5) plot(shap_vals_in_sim[,j], shap_vals_oob_sim[,j], pch=20, col = rgb(0,0,1,0.5))
for (j in 1:5) hist(shap_vals_sim[,j])
