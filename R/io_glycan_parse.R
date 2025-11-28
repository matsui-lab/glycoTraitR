res_cat <- c(
  "[a2122h-1x_1-5_2*NCC/3=O]"      = "N",
  "[a1221m-1x_1-5]"                = "F",
  "[a1122h-1x_1-5]"                = "H",
  "[a2122h-1b_1-5_2*NCC/3=O]"      = "N",
  "[a1122h-1b_1-5]"                = "H",
  "[a1122h-1a_1-5]"                = "H",
  "[a2122h-1a_1-5_2*NCC/3=O]"      = "N",
  "[u2122h_2*NCC/3=O]"             = "N",
  "[a2112h-1b_1-5]"                = "H",
  "[a1221m-1a_1-5]"                = "F",
  "[Aad21122h-2a_2-6_5*NCC/3=O]"   = "A",
  "[a2112h-1a_1-5]"                = "H",
  "[a2112h-1b_1-5_2*NCC/3=O]"      = "N",
  "[Aad21122h-2a_2-6_5*NCCO/3=O]"  = "G",
  "[a1221m-1b_1-5]"                = "F",
  "[a2112m-1a_1-5]"                = "F",
  "[o2122h_2*NCC/3=O]"             = "N",
  "[a2112h-1b_1-4]"                = "H",
  "[a4334h-1b_1-4]"                = "H",
  "[a4344h-1x_1-5_2*NCC/3=O]"      = "N",
  "[axxxxh-1x_1-5_2*NCC/3=O]"      = "N",
  "[axxxxh-1x_1-?_2*NCC/3=O]"      = "N",
  "[a2122h-1b_1-5]"                = "H",
  "[a2112h-1x_1-5]"                = "H",
  "[axxxxh-1b_1-5_2*NCC/3=O]"      = "N",
  "[a2112h-1x_1-5_2*NCC/3=O]"      = "N",
  "[axxxxh-1b_1-?_2*NCC/3=O]"      = "N",
  "[a2112h-1a_1-4]"                = "H",
  "[a2122h-1a_1-5]"                = "H",
  "[a4334h-1a_1-4]"                = "H",
  "[a4344h-1a_1-?]"                = "H",
  "[a2112m-1b_1-5]"                = "F",
  "[Aad21122h-2x_2-6_5*NCC/3=O]"   = "A",
  "[a2112h-1x_1-4]"                = "H",
  "[a4334h-1x_1-4]"                = "H",
  "[axxxxh-1x_1-5]"                = "H",
  "[a1121h-1a_1-5]"                = "H",
  "[a2122h-1x_1-5]"                = "H",
  "[a4344h-1x_1-?]"                = "H",
  "[Aad21122h-2x_2-6_5*NCCO/3=O]"  = "G",
  "[a2112h-1a_1-5_2*NCC/3=O]"      = "N",
  "[axxxxh-1x_1-?]"                = "H",
  "[a4334h-1b_1-5]"                = "H",
  "[a4344h-1b_1-5_2*NCC/3=O]"      = "N",
  "[a2122h-1b_1-?]"                = "H",
  "[a2122h-1b_1-?_2*NCC/3=O]"      = "N"
)

# WURCS string into tree structure
wurcs_to_tree <- function(w){
  core <- sub("^WURCS=[^/]+/", "", w)

  parts <- character()
  cur <- ""
  depth <- 0L
  for(ch in strsplit(core, "")[[1]]){
    if (ch == "[") depth <- depth + 1L
    if (ch == "]") depth <- depth - 1L
    if (ch == "/" && depth == 0L){
      parts <- c(parts, cur)
      cur <- ""
    } else {
      cur <- paste0(cur, ch)
    }
  }
  parts <- c(parts, cur)

  ## UniqueRES --------
  res_raw  <- regmatches(parts[2],
                         gregexpr("\\[[^]]+\\]", parts[2], perl = TRUE))[[1]]
  res_raw <- res_cat[res_raw]

  ## RES-Sequence -------
  res_idx  <- as.integer(strsplit(parts[3], "-", fixed = TRUE)[[1]])
  node  <- res_raw[res_idx]

  ## LIN --------
  edge <- unlist(strsplit(parts[4], "_", fixed = TRUE))
  edge <- unlist(strsplit(edge, "_", fixed = TRUE))
  edge <- unlist(strsplit(edge, "\\|"))
  edge <- grep("^[A-Za-z][0-9?]+-[A-Za-z][0-9?]+$", edge, value = TRUE)
  edge <- sub("^([A-Za-z])[0-9?]+-([A-Za-z])[0-9?]+$", "\\1-\\2", edge, perl = TRUE)

  list(node = node, edge = edge)
}

# pGlyco3 format into tree structure
pGlyco3_to_tree <- function(expr) {
  expr <- strsplit(expr, "")[[1]]
  node <- character()
  edge <- character()
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
      node <- c(node, ch)
      if (!is.null(parent)) {
        edge <- c(edge, paste0(parent, "-", label))
      }
      parent <- label
    }
  }
  list(node = node, edge = edge)
}


