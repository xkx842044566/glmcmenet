#' Build sibling/cousin group index lists from CME column names
#'
#' @description
#' Creates index groups for **siblings** and **cousins** based on a design
#' matrix whose columns are main effects (e.g., `"A"`) and conditional main
#' effects (CMEs) named like `"A|B+"` or `"A|B-"`.
#' For each effect name `E`, this returns:
#' \itemize{
#'   \item \code{S<E>} — indices of columns where \code{E} is the **main** effect
#'         (the ME column \code{"E"} itself, and all CMEs \code{"E|*+/-"}).
#'   \item \code{C<E>} — indices of columns where \code{E} appears as the
#'         **conditioning** effect (all CMEs \code{"*|E+/-"} and the ME \code{"E"}).
#' }
#'
#' @param X A numeric (or logical) matrix whose column names encode main effects
#'   and CMEs using the convention \code{"<Main>|<Cond>+"} and
#'   \code{"<Main>|<Cond>-"} for conditional effects, and \code{"<Main>"}
#'   for standalone main effects.
#'
#' @details
#' Group names are constructed as \code{"S<effect>"} and \code{"C<effect>"}.
#' Effects are taken from the column-name substrings before/after \code{"|"} with
#' any trailing sign \code{"+"} or \code{"-"} removed for the conditioning effect.
#' If an effect appears only as a conditioner (never as a standalone ME column),
#' its \code{C<effect>} group is still created and populated.
#'
#' @return
#' A named list of integer vectors. Each element contains 1-based column indices
#' into \code{X} for the corresponding group.
#'
#'
#' @seealso \code{\link{full.model.mtx}} for constructing ME/CME design matrices.
#' @export



grouplist <- function(X) {

  names <- colnames(X)
  mainEffects <- unique(gsub("\\|.*$", "", names))

  groups <- list()

  # Initialize groups for each main effect
  for (effect in mainEffects) {
    groups[[paste0("S", effect)]] <- integer()
    groups[[paste0("C", effect)]] <- integer()
  }

  for (col in names) {
    if (grepl("\\|", col)) {
      # This is a conditional effect
      parts <- strsplit(col, "\\|")[[1]]
      mainEffect <- parts[1]
      conditionEffect <- gsub("\\+|-", "", parts[2])

      # Add to S(main) and C(condition)
      groups[[paste0("S", mainEffect)]] <- c(groups[[paste0("S", mainEffect)]], which(colnames(X) == col))
      groups[[paste0("C", conditionEffect)]] <- c(groups[[paste0("C", conditionEffect)]], which(colnames(X) == col))
    } else {
      # This is a main effect
      groups[[paste0("S", col)]] <- c(groups[[paste0("S", col)]], which(colnames(X) == col))
      groups[[paste0("C", col)]] <- c(groups[[paste0("C", col)]], which(colnames(X) == col))
    }
  }
  return(groups)
}



