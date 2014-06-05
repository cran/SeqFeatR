\name{founderelim}
\alias{founderelim}
\title{Tries to identify possible problems in the sequence data due to a founder effect}
\description{takes the result of co-mutation analysis and corrects them with an added column with a note if the row is in the same branch.}
\usage{founderelim(path_to_file_s = NULL,
    path_to_file_p = NULL, path_to_file_m = NULL,
    save_name, value = 7)}
\arguments{
  \item{path_to_file_s}{a FASTA file with sequence data. For reference please look in example file.}
  \item{path_to_file_p}{a csv file with the results of the co-mutation analysis. For reference please look in example file.}
  \item{path_to_file_m}{a nexus file with tree data. For reference please look in example file.}
  \item{save_name}{the file name of the result file.}
  \item{value}{number of tips so that the whole bulk is considered as in one branch.}
}
\details{Takes the positions of the co-mutation results and a tree of the sequences with which the co-mutation was created and counts how often the pair of AAs is in one branch. If there are more than 'value' times the number of found co-mutations in one generated subtree, they are considered to be because of the founder effect.}
\value{the same file as the co muation input file with an added column in which the user can read if the row is found possibly due to the founder effect.}
\references{Mayr, Ernst (1954). "Change of genetic environment and evolution". In Julian Huxley. 
Evolution as a Process. London: George Allen & Unwin. OCLC 974739}
\author{Bettina Budeus}
\note{Only use files generated without allel usage!}


\seealso{\code{\link{assocpair}}}
\examples{
ex <- system.file("extdata", "Example_aa.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "co_mutation_results_wo_allels.csv", package="SeqFeatR")
	co <- system.file("extdata", "Example_tree.nh", package="SeqFeatR")
	founderelim(ex,
 ep,
 co,
 "d.csv",
 7)
}
