\name{hlaTree}
\alias{hlaTree}
\title{Find possible escape mutations}
\description{Takes a tree and a set of sequences and searches for possible escape mutations}
\usage{hlaTree(path_to_file_s = NULL, path_to_file_t = NULL,
    save_name_csv, patnum.threshol = 1, thr.sig.f = 0.05,
    A11, A12, A21, A22, B11, B12, B21, B22)}
\arguments{
  \item{path_to_file_s}{a FASTA file with sequence data. Homzygot patiens have to have a 00 instead of the HLA typeFor reference please look in
example file.}
  \item{path_to_file_t}{a tree file of the sequences. For reference please look in example file.}
  \item{save_name_csv}{the file name of the result file in csv format}
  \item{patnum.threshol}{the minimum number of patients of one HLA type to consider in the calculation.}
  \item{thr.sig.f}{p-value threshold below which the results from fishers exact test should be added to output.}
  \item{A11}{the position of the start of the first HLA A Allel in the description block of the FASTA file.}
  \item{A12}{the position of the end of the first HLA A Allel in the description block of the FASTA file.}
  \item{A21}{the position of the start of the second HLA A Allel in the description block of the FASTA file.}
  \item{A22}{the position of the end of the second HLA A Allel in the description block of the FASTA file.}
  \item{B11}{the position of the start of the first HLA B Allel in the description block of the FASTA file.}
  \item{B12}{the position of the end of the first HLA B Allel in the description block of the FASTA file.}
  \item{B21}{the position of the start of the second HLA B Allel in the description block of the FASTA file.}
  \item{B22}{the position of the end of the second HLA B Allel in the description block of the FASTA file.}
}



\author{Bettina Budeus}

\examples{
\dontrun{
ex <- system.file("extdata", "Example_aa.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "Example_tree.nh", package="SeqFeatR")
	hlaTree(ex,
 ep,
 "ancestral_analysis_results.csv",
 1,
 0.05,
 10,
 11,
 13,
 14,
 17,
 18,
 20,
 21)
}	
}
