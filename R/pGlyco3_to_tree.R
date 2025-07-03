pGlyco3_to_tree <- function(expr) {
  expr <- strsplit(expr, "")[[1]]
  seq_raw <- character()
  lin_vec <- character()
  stack <- c()
  parent <- NULL
  idx <- 1
  for (ch in expr) {
    if (ch == "(") {
      stack <- c(parent, stack)
    } else if (ch == ")") {
      parent <- stack[1]
      stack <- stack[-1]
    } else if (grepl("[A-Za-z]", ch)) {
      label <- letters[idx]
      idx <- idx + 1
      seq_raw <- c(seq_raw, ch)
      if (!is.null(parent)) {
        lin_vec <- c(lin_vec, paste0(parent, "-", label))
      }
      parent <- label
    }
  }
  list(seq_raw = seq_raw, lin_vec = lin_vec)
}


tree_to_pGlyco3 <- function(list) {
  seq_raw <- list$seq_raw
  lin_vec <- list$lin_vec

  n_nodes  <- length(seq_raw)
  labels   <- letters[seq_len(n_nodes)]
  mono_map <- setNames(seq_raw, labels)

  adj <- lapply(labels, function(x) character())
  names(adj) <- labels
  for (edge in lin_vec) {
    parts <- strsplit(edge, "-", fixed = TRUE)[[1]]
    a <- parts[1];  b <- parts[2]
    adj[[a]] <- c(adj[[a]], b)
    adj[[b]] <- c(adj[[b]], a)
  }

  parent <- setNames(rep(NA_character_, n_nodes), labels)
  parent["a"] <- ""                     # sentinel
  queue <- "a"
  while (length(queue)) {
    node  <- queue[1]; queue <- queue[-1]
    for (nb in adj[[node]]) {
      if (is.na(parent[nb])) {
        parent[nb] <- node
        queue      <- c(queue, nb)
      }
    }
  }

  children <- lapply(labels, function(x) character())
  names(children) <- labels
  for (ch in labels[parent != "" & !is.na(parent)]) {
    par <- parent[ch]
    children[[par]] <- c(children[[par]], ch)
  }
  children <- lapply(children, sort)

  recurse <- function(node) {
    expr <- mono_map[[node]]
    for (ch in children[[node]]) {
      expr <- paste0(expr, "(", recurse(ch), ")")
    }
    expr
  }
  paste0("(", recurse("a"), ")")
}

# expr <- "(N(F)(N(H(N)(H(N)(N(H)(F)))(H(N)(N(F)(H(A)))))))"
# expr %>% pGlyco3_to_tree() %>% tree_to_pGlyco3() == expr

