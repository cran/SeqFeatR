\name{assocpair}
\alias{assocpair}
\title{Test for co-mutation without HLA types}
\description{Calculates if there are two positions in the given sequence which may have a co-mutation.}
\usage{assocpair(path_to_file_s = NULL,
    path_to_file_c = NULL, save_name, dna = FALSE,
    thr.sig.f = 0.05, statistical_correction = "bonferroni")}
\arguments{
  \item{path_to_file_s}{a FASTA file with sequence data. For reference please look in example file.}
  \item{path_to_file_c}{a FASTA file with consensus data. For reference please look in example file.}
  \item{save_name}{the file name of the result file in csv format}
  \item{dna}{if the data is in DNA or amino Acid Code}
  \item{thr.sig.f}{p-value threshold below which the results from fishers exact test should be added to output.}
  \item{statistical_correction}{the statistical correction applied to the p-values. Input can be: "holm",\cr
 "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".}
}
\details{For every position in the sequence given within the FASTA file a fishers exact test is made with every other position in the sequence and the Consensus sequence. The result p-values are collected in one big table if they are below a given value (thr.sig.f). If there are any HLA types in the FASTA description they are simply ignored. Also uses p.adjust from stats package to calculate some p-value correction additionally as an extra column in the csv output.}
\value{a table with every possible co-mutation below the given p-value.}

\author{Bettina Budeus}
\note{If you want an graphical output, please use:\cr
\code{\link{vispair}}.}


\seealso{\code{\link{vispairfeat}}}
\examples{
\dontrun{
ex <- system.file("extdata", "Example_aa.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "Example_Consensus_aa.fasta", package="SeqFeatR")
	assocpair(ex,
 ep,
 "co_mutation_results_wo_allels.csv",
 FALSE,
 0.05,
 "bonferroni")
}
}
