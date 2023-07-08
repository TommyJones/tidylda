# below is an experiment. Looking to fix "recover_counts_from_probs"


prob_matrix

prior_matrix

total_vector

Cd_test <- prob_matrix * total_vector + prob_matrix * rowSums(prior_matrix) - prior_matrix

Cd_test <- 
  Cd_test |> apply(1, function(x) {
    # no negative probabilities
    p <- x + abs(x) + .Machine$double.eps 
    
    result <- rmultinom(
      n = 1,
      size = sum(x),
      prob = p
    )
    
    result[, 1]
  }) |> 
  t()
