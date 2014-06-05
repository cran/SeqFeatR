\name{assocpairfeat}
\alias{assocpairfeat}
\title{Tries to identify sequence positions pairs in accordance to sequence feature.}
\description{Calculates if there are two positions in the given sequence which may have a co-mutation.}
\usage{assocpairfeat(path_to_file = NULL, save_name, dna = FALSE, 
    patnum.threshol = 1, thr.sig.f = 0.05, 
    A11a, A12a, A21a, A22a, B11a, B12a, B21a, B22a,
    statistical_correction = "bonferroni")}
\arguments{
  \item{path_to_file}{a FASTA file with sequence data. For reference please look in example file.}
  \item{save_name}{the file name of the result file in csv format}
  \item{dna}{if the data is in DNA or amino Acid Code}
  \item{patnum.threshol}{the minimum number of patients of one HLA type to consider in the calculation.}
  \item{thr.sig.f}{p-value threshold below which the results from fishers exact test should be added to output.}
  \item{A11a}{the position of the start of the first HLA A Allel in the description block of the FASTA file.}
  \item{A12a}{the position of the end of the first HLA A Allel in the description block of the FASTA file.}
  \item{A21a}{the position of the start of the second HLA A Allel in the description block of the FASTA file.}
  \item{A22a}{the position of the end of the second HLA A Allel in the description block of the FASTA file.}
  \item{B11a}{the position of the start of the first HLA B Allel in the description block of the FASTA file.}
  \item{B12a}{the position of the end of the first HLA B Allel in the description block of the FASTA file.}
  \item{B21a}{the position of the start of the second HLA B Allel in the description block of the FASTA file.}
  \item{B22a}{the position of the end of the second HLA B Allel in the description block of the FASTA file.}
  \item{statistical_correction}{the statistical correction applied to the p-values. Input can be: "holm",\cr
 "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".}
}
\details{For every position in the sequence given within the FASTA file a fishers exact test is made with every other position in the sequence and every HLA-type. The result p-values are collected in one big table if they are below a given value (thr.sig.f). Also uses p.adjust from stats package to calculate some p-value correction additionally as an extra column in the csv output.}
\value{a table with every possible co-mutation below the given p-value.}

\author{Bettina Budeus}
\note{If you want an graphical output, please use \code{\link{vispairfeat}}.
If a patient has a homozygot HLA Allel, then please change the second one to "00" (without ") instead!}


\seealso{\code{\link{vispairfeat}}}
\examples{
\dontrun{
ex <- system.file("extdata", "Example_aa.fasta", package="SeqFeatR")
	assocpairfeat(ex,
 "co_mutation_results.csv",
 FALSE,
 1,
 0.05,
 10,
 11,
 13,
 14,
 17,
 18,
 20,
 21,
 "bonferroni")
}
}
