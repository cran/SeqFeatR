\name{assocpoint}
\alias{assocpoint}
\title{Calculate point mutations vs. feature(s)}
\description{Searches in a given sequence for possible associations between amino acid/nucleotide and sequence position.}
\usage{assocpoint(path_to_file_s = NULL,
    path_to_file_p = NULL, path_to_file_m = NULL, save_name,
    save_name_csv, dna = FALSE, 
    patnum.threshold = 1, optical.significance_level = 0.05, 
    star.significance_level = 0.001, pot.epitope_level = 0.05, 
    window_threshold = 0.05, A11, A12, A21, A22, B11, B12, B21, 
    B22, statistical_correction = "bonferroni",
    with_phylogenetic_comparison = FALSE,
    matrix_for_phylo = NULL, reference_sequence=NULL, one_ident=FALSE, window_size=9)}
\arguments{
  \item{path_to_file_s}{a FASTA file with sequence data. Homzygot patiens have to have a 00 instead of the HLA typeFor reference please look in
example file.}
  \item{path_to_file_p}{a csv file with known epitopes. For reference please look in example file.}
  \item{path_to_file_m}{a csv file with HLA binding motifs if available. For reference please look in example file.}
  \item{save_name}{the file name of the result file in pdf format}
  \item{save_name_csv}{the file name of the result file in csv format}
  \item{dna}{if the data is in DNA or amino Acid Code}
  \item{patnum.threshold}{the minimum number of patients of one HLA type to consider in the calculation.}
  \item{optical.significance_level}{the height of the optical horizontal line in the graphical output. It should be yout chosen p-value limit, like 0.05.}
  \item{star.significance_level}{the height of the invisible horizontal line above which alle points are marked as stars. Should be a high p-value limit,
like 0,001.}
  \item{pot.epitope_level}{the level below which a p-value should be consideres in the calculation of possible epitopes. 
Should be the same as optical.significance_leve.}
  \item{window_threshold}{the level below which a p-value should be consideres in the calculation of the sliding window. 
Should be the same as optical.significance_leve. }
  \item{A11}{the position of the start of the first HLA A Allel in the description block of the FASTA file.}
  \item{A12}{the position of the end of the first HLA A Allel in the description block of the FASTA file.}
  \item{A21}{the position of the start of the second HLA A Allel in the description block of the FASTA file.}
  \item{A22}{the position of the end of the second HLA A Allel in the description block of the FASTA file.}
  \item{B11}{the position of the start of the first HLA B Allel in the description block of the FASTA file.}
  \item{B12}{the position of the end of the first HLA B Allel in the description block of the FASTA file.}
  \item{B21}{the position of the start of the second HLA B Allel in the description block of the FASTA file.}
  \item{B22}{the position of the end of the second HLA B Allel in the description block of the FASTA file.}
  \item{statistical_correction}{the statistical correction applied to the p-values. Input can be: "holm",\cr
 "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".}
  \item{with_phylogenetic_comparison}{if the sequences with the certain trait and Amino Acid with the lowest p.value should be compared in sequence similarity with all sequences. See details.}
  \item{matrix_for_phylo}{the matrix for phylogenetic analysis. Has to be a matrix which is with this name available in the Biostrings package.}
  \item{reference_sequence}{Reference sequence for better comparison of results. Will be added as extra column in the results.}
  \item{one_ident}{if there is only one feature. See details.}
  \item{window_size}{sizwe for the sliding window - only used in graphical output!}
}
\details{For every position in the sequences a fisher's exact test of (Certain feature identifier, other feature identifier, aa and not this aa) is calculated and the p-value is pasted in a table.
If you have only one kind of feature, the please set one_ident to TRUE. This can be if you are to analyse tropisms or treatment with only Yes and No. If you have only one feature it has to be seperated from the rest of the comment line in the fasta file at the end with a semicolon. If you have two kinds of feature (like HLA types) you can give there exact position inside the comment line.
Every HLA type has a column, every position a row. There is also a graphical output generated in which every HLA type has it's own page. In this output, a line is inserted, which shows the user where the p-value limit is considered above to be significant. 
Integrated in this graphic are the known epitopes if there are any and the calculated possible epitopes, which are calculated with the known HLA binding motifs. As an extras information a second graphical output is added in which the user can see the level of significant muations in a 9er sliding window.
Each point in this graphic is the number of mutations in the next nine sequence positions.
Uses EasyGregExpr from http://www.r-bloggers.com/developing-a-user-friendly-regular-expression-function-easygregexpr/
Also uses p.adjust from stats package to calculate some p-value correction additionally for the graphical output and as an extra column in the csv output.
Please consider that you have to use the dna/aa switch depending on your data!
If pylogenetic comparison was choosen, the programm compares for every lowest p.value below a certain user given value the sequence similarity of the sequences with this trai with the similarity of all sequences with a wilcoxon test and adds the p.values in another column.
Please use a "B" for unknown amino acid.}


\author{Bettina Budeus}
\note{If a patient has a homozygot HLA Allel, then please change the second one to "00" (without ") instead!}


\seealso{\code{\link{assocpointtuple}}}
\examples{
ex <- system.file("extdata", "Example_aa.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "Example_epitopes_aa.csv", package="SeqFeatR")
	co <- system.file("extdata", "Example_HLA_binding_motifs_aa.csv", package="SeqFeatR")
	assocpoint(ex,
 ep,
 co,
 "epitope_results.pdf",
 "epitope_results.csv",
 FALSE,
 1,
 0.01, 
 0.5,
 0.01,
 0.01,
 10,
 11,
 13,
 14,
 17,
 18,
 20,
 21,
 "bonferroni",
 FALSE,
 BLOSUM62,
 NULL,
 FALSE,
 9)
}
