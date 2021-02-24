# When log(0) need to use inverse hyperbolic sine transformation (REF NEEDED)
#https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions#Inverse_hyperbolic_sine
ihs = function(x) {
  return(log(x + sqrt(x^2+1)))
}

# run this function - critical for workflow for plotting results:
tidy = function(model, ...) {
  if (attributes(model)$class == "felm") {
    data = cbind.data.frame(term = names(coef(model)),
                            coef(summary(model, robust= T)),
                            confint(model))
    names(data) = c("term","estimate","std.error", "statistic","p.value","conf.low","conf.high")
    return(data)
  } else if (attributes(model)$class == "fixest") {
    data = cbind.data.frame(term = names(coef(model)),
                            coef(summary(model)),
                            se(model),
                            confint(model))
    names(data) = c("term","estimate","std.error","conf.low","conf.high")
    return(data)
  }
}

# run this function - critical for workflow for plotting results:
# tidy = function(model, ...) {
#   data = cbind.data.frame(term = names(coef(model)),
#                           coef(summary(model, robust= T)),
#                           confint(model))
#   names(data) = c("term","estimate","std.error", "statistic","p.value","conf.low","conf.high")
#   return(data)
# }