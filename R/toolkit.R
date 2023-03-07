#' Merge lists
#'
#'
#' @export
#'
tl_merge_lists <- function(lst_to, lst_from) {

    for (i in seq(length(lst_from))) {
        cur_name <- names(lst_from)[i]

        if (cur_name %in% names(lst_to))
            next

        lst_to[[cur_name]] <- lst_from[[cur_name]]
    }

    lst_to
}
