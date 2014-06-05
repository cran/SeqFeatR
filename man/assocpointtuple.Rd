\name{assocpointtuple}
\alias{assocpointtuple}
\title{Compare results of co-mutation and point mutation association with feature(s)}
\description{Checks which positions in the results of epitope finder are also in the results of given co-mutation. Therefore just compares certain columns of the results, so they have to be in the output format from point mutation with feature association and co-mutation.}
\usage{assocpointtuple(path_to_file_s = NULL, 
    path_to_file_e = NULL, path_to_file_m = NULL,
    p.value = 0.05, save_name_result, save_name_pos)}
\arguments{
  \item{path_to_file_s}{file with sequence data in FASTA format.  For reference please look in example file.}
  \item{path_to_file_e}{file with result from point mutation association with feature(s). For reference please look in example file.}
  \item{path_to_file_m}{file with results from co-mutation analysis. For reference please look in example file.}
  \item{p.value}{p-value to be defined as significant.}
  \item{save_name_result}{the file name of the result file.}
  \item{save_name_pos}{the file name of the possible compensatory file. Only written if there are such possible compensatory mutations.}
}
\details{This function takes the results from point mutation association with feature(s) and co-mutation analysis and tries to combine them into one result. It checks if a possible epitope, as it is an output of the first, is also one of the positions which are found in the co-mutation. If this is the case, they may be a point of compensatory mutation.}
\value{the corresponding table with an added column with q-values}

\author{Bettina Budeus}



\seealso{\code{\link{assocpoint}}

\code{\link{assocpairfeat}}

\code{\link{assocpair}}}
\examples{
ex <- system.file("extdata", "Example_aa.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "epitope_results.csv", package="SeqFeatR")
	co <- system.file("extdata", "co_mutation_results.csv", package="SeqFeatR")
	assocpointtuple(ex, ep, co, 0.05,
	 "co_mut_in_epitopes_result.csv",
	 "possible_compensatory_mutation.csv")
}
