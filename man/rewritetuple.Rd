\name{rewritetuple}
\alias{rewritetuple}
\title{Rewrite result from tuple/feature association analysis}
\description{This is a function to put the raw data from the shared mutations program into an usable form for tartan plot}
\usage{rewritetuple(input_file = NULL, save_name, 
first_position, second_position, value_position, sep, cutoff)}
\arguments{
  \item{input_file}{the result file from get_shared_mutations}
  \item{save_name}{the name of the result file from this program}
  \item{first_position}{the column position of the first sequence position}
  \item{second_position}{the column position of the second sequence position}
  \item{value_position}{the column position of the p-value}
  \item{sep}{the seperator of the input file}
  \item{cutoff}{the cutoff below which p-value the data should be included in the result}
}
\details{It just extracts the two sequence positions and the p-value if it is below the cut off.}

\author{Bettina Budeus}


\seealso{\code{\link{assoctuple}}}
\examples{
ex <- system.file("extdata", "shared_mutations_result.csv", package="SeqFeatR")
	rewritetuple(
	input_file = ex, 
save_name = "shared_mutations_result_for_tartan.csv", 
first_position = 3, 
second_position = 4, 
value_position = 10, 
sep = "\t", 
cutoff = 0.1
)
}
