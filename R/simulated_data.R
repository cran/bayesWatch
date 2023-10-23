#' Simulated Data with Imposed Change-points.
#'
#' Data simulated using the BDgraph package.  A change-point is imposed between days 5
#' and 6.  The change only occurs in variables 3 and 4.
#'
#' `full_data` is a data frame with 1,000 rows and 5 columns.
#' `day_of_observations`; is a timestamp of each of `full_data`'s 1,000 rows.
#' `day_dts`; is a vector of unique elements from `day_of_observations`..
#'
#' @format `full_data` is a matrix, the latter two are vectors.
#' @examples
#' full_data
#' day_of_observations
#' day_dts
"full_data"

#' @rdname full_data
#' @format NULL
"day_of_observations"

#' @rdname full_data
#' @format NULL
"day_dts"
