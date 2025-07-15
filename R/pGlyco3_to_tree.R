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
      label <- c(letters, LETTERS)[idx]
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
  labels   <- c(letters, LETTERS)[seq_len(n_nodes)]
  mono_map <- setNames(seq_raw, labels)

  ## -------- build undirected adjacency ---------------------------------
  adj <- lapply(labels, function(x) character())
  names(adj) <- labels
  for (edge in lin_vec) {
    parts <- strsplit(edge, "-", fixed = TRUE)[[1]]
    a <- parts[1];  b <- parts[2]
    adj[[a]] <- c(adj[[a]], b)
    adj[[b]] <- c(adj[[b]], a)
  }

  ## -------- BFS to get parent map --------------------------------------
  parent <- setNames(rep(NA_character_, n_nodes), labels)
  parent["a"] <- ""                      # sentinel root
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

  ## -------- children list (unordered) ----------------------------------
  children <- lapply(labels, function(x) character())
  names(children) <- labels
  for (ch in labels[parent != "" & !is.na(parent)]) {
    par <- parent[ch]
    children[[par]] <- c(children[[par]], ch)
  }

  ## -------- recursive canonical builder --------------------------------
  recurse <- function(node) {
    child_exprs <- vapply(children[[node]], recurse, character(1))

    if (length(child_exprs)) {
      ord <- order(nchar(child_exprs), child_exprs)
      child_exprs <- child_exprs[ord]
      paste0(
        mono_map[[node]],
        paste0("(", child_exprs, ")", collapse = "")
      )
    } else {
      mono_map[[node]]
    }
  }

  paste0("(", recurse("a"), ")")
}


# expr <- gd_pglyco[sample(length(gd_pglyco), 1)]
# expr %>% pGlyco3_to_tree() %>% tree_to_pGlyco3() == expr

