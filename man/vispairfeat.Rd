\name{vispairfeat}
\alias{vispairfeat}
\title{Create co-mutation plot with allels.}
\description{creates a plot for the results of co-mutation analysis}
\usage{vispairfeat(path_to_file = NULL, save_name, 
    thr.sig.fi = 0.01)}
\arguments{
  \item{path_to_file}{a csv file with results from co mutation with HLA types. For reference please look in example file.}
  \item{save_name}{the file name of the result file in pdf format}
  \item{thr.sig.fi}{the significance value below which the results of the co-mutation analysis are considered to be relevant enough to be plotted}
}
\details{For every HLA type (if only one p-value is below thr.sig.fi) a page in an pdf is created with a sequence x sequences graphic. In this graphic every dots (one position with the other position) colour is based on it's value. The colour code is on the right side.}


\author{Bettina Budeus}
\note{Only use files generated without allel usage!}


\seealso{\code{\link{assocpairfeat}}}
\examples{
ex <- system.file("extdata", "co_mutation_results.csv", package="SeqFeatR")
	vispairfeat(ex,
 "co_mut_g_results.pdf", 0.05)
}
