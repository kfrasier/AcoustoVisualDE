weighted_logloss <- function (y_pred, y_true, y_weights) 
{
  eps <- 1e-15
  y_pred <- pmax(pmin(y_pred, 1 - eps), eps)
  LogLoss <- -mean((y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred))
                   * (y_weights),na.rm = TRUE)
  return(LogLoss)
}