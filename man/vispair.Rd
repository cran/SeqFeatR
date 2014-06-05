\name{vispair}
\alias{vispair}
\title{Create co-mutation plot without allels.}
\description{creates a plot for the results of co-mutation analysis}
\usage{vispair(path_to_file = NULL, 
    save_name, thr.sig.fi = 0.01)}
\arguments{
  \item{path_to_file}{a csv file with results from co mutation without HLA types. For reference please look in example file.}
  \item{save_name}{the file name of the result file in pdf format}
  \item{thr.sig.fi}{the significance value below which the results of the co-mutation analysis are considered to be relevant enough to be plotted}
}
\details{A page in an pdf is created with a sequence x sequences graphic. In this graphic every dot (one position with the other position) above the thr.sig.fi is based on it's value marked with a colour. The colour code is on the right side.}


\author{Bettina Budeus}
\note{Only use files generated without allel usage!}


\seealso{\code{\link{assocpair}}}
\examples{
ex <- system.file("extdata", "co_mutation_results_wo_allels.csv", package="SeqFeatR")
	vispair(
 ex,
 "co_mut_g_results_wo_allels.pdf",
 0.05)
}
