\name{tartan}
\alias{tartan}
\title{Create graphics for association of mutation pairs and features from two (comparable) analyses.}
\description{Creates a plot for two results of co-mutation analysis.}
\usage{tartan(path_to_file=NULL,
 path_to_file2=NULL, save_name, space, colors, name_positions, names,
 ticks, first_position_1, second_position_1, value_1, first_position_2,
 second_position_2, value_2, with_distance_matrix, path_to_distance_matrix)}
\arguments{
  \item{path_to_file}{a csv file with results from co mutation without HLA types. For reference please look in example file.}
  \item{path_to_file2}{a csv file with results from co mutation without HLA types. For reference please look in example file.}
  \item{save_name}{the file name of the result file in pdf format}
  \item{space}{the space between blocks}
  \item{colors}{the colors of the tartan plot}
  \item{name_positions}{the positions on x and y axis where the labels should be inserted}
  \item{names}{the labels to be inserted on x and y axis. Beware! name_position and names must have the same length!}
  \item{ticks}{the ticks which should mark intresting spots and are also the breakpoints on which the blocks are seperated.}
  \item{first_position_1}{column in which the first position is located in the first csv file}
  \item{second_position_1}{column in which the second position is located in the first csv file}
  \item{value_1}{column in which the value is located in the first csv file}
  \item{first_position_2}{column in which the first position is located in the second csv file}
  \item{second_position_2}{column in which the second position is located in the second csv file}
  \item{value_2}{column in which the value is located in the csv second file}
  \item{with_distance_matrix}{if there is a distance matrix available and should be used}
  \item{path_to_distance_matrix}{the path to the distance matrix. Can be zero.}
}
\details{A page in an pdf is created with a sequence x sequences graphic. In this graphic every dot (one position with the other position) above the thr.sig.fi is based on it's value marked with a colour. The colour code is on the right side. Please see that the values in you're result csv are allready in a comparebale range! To achieve this you may use logarithm or similar mechanics.}

\author{Bettina Budeus}
\note{Only use files generated without allel usage!}


\seealso{\code{\link{assocpair}}}
\examples{
ex <- system.file("extdata", "shared_mutations_result_for_tartan.csv", package="SeqFeatR")
	ex2 <- system.file("extdata", "co_mutation_results_wo_allels.csv", package="SeqFeatR")
	tartan(
 ex,
 ex2,
 "tartan_plot.pdf",
 5,
 c("wheat", "darkblue", "black", "green"),
 c(1,30),
 c("S","F"),
 c(10,30),
 2,
 3,
 4,
 2,
 3,
 11,
 FALSE,
 NULL)
}
