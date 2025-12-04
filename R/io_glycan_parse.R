# Mapping of UniqueRES tokens (WURCS) to simplified residue types
# Internal constant used by `wurcs_to_tree()`
WURCS_RES_MAP <- c(
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

#' Convert a WURCS string into a glycan tree structure
#'
#' Parse a WURCS (WURCS 2.0) glycan annotation and extract:
#' \itemize{
#'   \item residue sequence (as a character vector)
#'   \item edge list (parent–child relationships)
#' }
#'
#' This function is used internally to convert WURCS
#' strings into a tree representation that can support structural
#' trait computation and graph construction.
#'
#' @details
#' The function performs several parsing steps:
#' \enumerate{
#'   \item Remove the WURCS prefix and split the string into components.
#'   \item Extract \emph{UniqueRES} entries and map them to residue symbols
#'         using the internal \code{WURCS_RES_MAP} table.
#'   \item Follow the \emph{RES sequence} index to reconstruct the residue
#'         vector.
#'   \item Parse \emph{LIN} entries and normalize them into simple "X-Y" edges.
#' }
#'
#' @param w A character string containing a WURCS 2.0 glycan annotation.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{node}: character vector of residue types in order
#'   \item \code{edge}: character vector of edges in "A-B" format
#' }
#'
#' @examples
#' # Example WURCS glycan string
#' w <- "WURCS=2.0/4,9,8/[u2122h_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-2-3-4-4-4-4-4-4/a4-b1_b4-c1_c3-d1_c6-f1_d2-e1_f3-g1_f6-h1_h2-i1"
#' tree <- wurcs_to_tree(w)
#' tree
#'
#' @keywords internal
wurcs_to_tree <- function(w) {
  core <- sub("^WURCS=[^/]+/", "", w)

  parts <- character()
  cur <- ""
  depth <- 0L
  for (ch in strsplit(core, "")[[1]]) {
    if (ch == "[") depth <- depth + 1L
    if (ch == "]") depth <- depth - 1L
    if (ch == "/" && depth == 0L) {
      parts <- c(parts, cur)
      cur <- ""
    } else {
      cur <- paste0(cur, ch)
    }
  }
  parts <- c(parts, cur)

  ## UniqueRES --------
  res_raw <- regmatches(
    parts[2],
    gregexpr("\\[[^]]+\\]", parts[2], perl = TRUE)
  )[[1]]
  res_raw <- WURCS_RES_MAP[res_raw]

  ## RES-Sequence -------
  res_idx <- as.integer(strsplit(parts[3], "-", fixed = TRUE)[[1]])
  node <- res_raw[res_idx]

  ## LIN --------
  edge <- unlist(strsplit(parts[4], "_", fixed = TRUE))
  edge <- unlist(strsplit(edge, "_", fixed = TRUE))
  edge <- unlist(strsplit(edge, "\\|"))
  edge <- grep("^[A-Za-z][0-9?]+-[A-Za-z][0-9?]+$", edge, value = TRUE)
  edge <- sub("^([A-Za-z])[0-9?]+-([A-Za-z])[0-9?]+$", "\\1-\\2", edge, perl = TRUE)

  list(node = node, edge = edge)
}

#' Convert a pGlyco3 glycan string into a glycan tree structure
#'
#' Parse a pGlyco3-style glycan expression (e.g. \code{"N(H(H))"})
#' and reconstruct the residue sequence and edge relationships
#' as a tree suitable for downstream structural analysis. This parser assumes
#' simple pGlyco3 monosaccharide symbols (e.g. "N", "H", "A", "F").
#'
#' @details
#' This function interprets parentheses as branch delimiters and assigns:
#' \enumerate{
#'   \item one residue per character (e.g. \code{N}, \code{H}, \code{A})
#'   \item parent–child edges based on bracket nesting
#' }
#'
#' Each residue is assigned a synthetic node label
#' (\code{a}, \code{b}, \code{c}, …), ensuring compatibility with
#' graph-based trait extraction.
#'
#' @param expr A character string representing the pGlyco3 glycan structure.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{node}: character vector of residue types
#'   \item \code{edge}: character vector of edges in "a-b" format
#' }
#'
#' @examples
#' # Example: parse a pGlyco3-style glycan expression into a tree
#' pGlyco_expr <- "(N(N(H(H(H))(H(H)(H)(H(H))))))"
#' # Convert to glycan tree structure
#' tree <- pGlyco3_to_tree(pGlyco_expr)
#' tree
#'
#' @keywords internal
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
