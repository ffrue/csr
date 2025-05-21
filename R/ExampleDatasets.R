#' atp_2017
#'
#' This datasets contains the name of the top male tennis players in 2017 as well as their tour wins and ranking points during that year.
#'
#' @format A data frame with 47 rows and 3 variables:
#' \describe{
#'   \item{player_slug}{string. player names}
#'   \item{ranking_points}{numeric. gained ATP points that year}
#'   \item{wins}{numeric. ATP matches won that year}
#' }
#' @source https://github.com/datasets/atp-world-tour-tennis-data
#'
#' @examples
#' head(atp_2017)
"atp_2017"


#' movies
#'
#' This datasets contains information about movies releases for big franchises since 2000.
#'
#' @format A data frame with 60 rows and 5 variables:
#' \describe{
#'   \item{Title}{string. movie title}
#'   \item{Year}{numeric. year of release}
#'   \item{Lifetime.Gross}{numeric. cumulated earnings over time in current US-dollars}
#'   \item{Budget}{numeric. budget in current US-dollars}
#'   \item{Runtime}{numeric. runtime in minutes}
#' }
#' @source https://www.kaggle.com/datasets/thedevastator/global-movie-franchise-revenue-and-budget-data
#'
#' @examples
#' head(movies)
"movies"
