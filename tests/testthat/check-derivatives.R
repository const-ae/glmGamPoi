# ## This file checks all the loglikelihood, score, derivative of score functions and makes
# ## sure that they are correct
#
# samples <- distraltparam::raltnbinom(n = 30, mean = 4, dispersion = 0.7)
# mu <- rnorm(n = 30, mean = 4)
# X <- matrix(rnorm(n = 30 * 4), nrow = 30, ncol = 4)
#
# # Loglikelihood vs conventional score
# nloptr::check.derivatives(.x = 2.1,
# func = function(log_theta){
#   - conventional_loglikelihood_fast(samples, mu = mu, log_theta = log_theta,
#                                     model_matrix = X, do_cr_adj = TRUE)
# }, func_grad = function(log_theta){
#   - conventional_score_function_fast(samples, mu = mu, log_theta = log_theta,
#                                      model_matrix = X, do_cr_adj = TRUE)
# })
#
# # Conventional score vs conventional score derivative
# nloptr::check.derivatives(.x = 2.1,
# func = function(log_theta){
#   - conventional_score_function_fast(samples, mu = mu, log_theta = log_theta,
#                                     model_matrix = X, do_cr_adj = TRUE)
# }, func_grad = function(log_theta){
#   - conventional_deriv_score_function_fast(samples, mu = mu, log_theta = log_theta,
#                                      model_matrix = X, do_cr_adj = TRUE)
# })
#
#
