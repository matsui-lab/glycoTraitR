cal_entropy <- function(vec) {
  if(length(vec) == 1){
    entropy <- 0
  } else {
    entropy <- -sum(vec * log(vec)) / log(length(vec))
  }
  entropy
}

add_trait_entropy <- function(sub_mat, trait) {
  trait_mean <- sum(sub_mat$Area * sub_mat[[trait]])
  trait_entropy <- sub_mat %>%
    group_by(!!sym(trait)) %>%
    reframe(new_area = sum(Area, na.rm = TRUE)) %>%
    pull(new_area) %>%
    cal_entropy()
  add_vec <- c(trait_mean, trait_entropy)
  names(add_vec) <- paste0(trait, "_" , c("mean", "entropy"))
  add_vec
}
