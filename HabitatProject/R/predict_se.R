predict_se <- function(model, data) {
  v <- predict(model, data, se.fit=TRUE, type = 'response',na.action = na.pass)
  cbind(p=as.vector(v$fit), se=as.vector(v$se.fit))
}