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


#' Draw sample size for each arm
#'
#' @param sample_size total sample size
#' @param ratio_by_arm randomization ratio by arm
#' @param fixed Fixed ratio or multinomial ratio
#' @param seed Random seed
#'
#' @export
#'
tl_draw_arm_size <- function(sample_size,
                             ratio_by_arm,
                             fixed = TRUE,
                             seed = NULL) {

    if (fixed) {
        n_by_arm <- sample_size
        n_by_arm <- n_by_arm * ratio_by_arm / sum(ratio_by_arm)
        n_by_arm <- floor(n_by_arm)
    } else {
        if (!is.null(seed))
            old_seed <- set.seed(seed)

        n_by_arm <- rmultinom(1, sample_size, ratio_by_arm)

        if (!is.null(seed))
            set.seed(old_seed)

    }

    n_by_arm
}
