w_KS_Stat <- function (y_pred, y_true) 
{
  n_pos <- sum(y_true == 0)
  n_neg <- sum(y_true == 1)
  rpp_vec <- numeric(length = 0)
  rnp_vec <- numeric(length = 0)
  for (i in sort(y_pred)) {
    pos_neg <- y_true[y_pred <= i]
    rpp_vec <- append(rpp_vec, sum(pos_neg == 0)/n_pos)
    rnp_vec <- append(rnp_vec, sum(pos_neg == 1)/n_neg)
  }
  KS_Stat <- max(rpp_vec - rnp_vec) * 100
  return(KS_Stat)
}