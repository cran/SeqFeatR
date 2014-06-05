\name{assoctuple}
\alias{assoctuple}
\title{Get co-mutations of tuple positions above threshold from point mutation analysis}
\description{Calculates a Fisher's exact test for a tuple of sequence positions and sequence assosiations like HLA type.}
\usage{assoctuple(original_fasta_seq,	epi_results, threshold,
	min_number_of_ele_in_tupel, max_number_of_ele_in_tuple,
	result_filename, column, column_of_values, column_of_aas,
	A11, A12, A21, A22, B11, B12, B21, B22, one_ident, idet)}
\arguments{
  \item{original_fasta_seq}{a fasta file with the sequences which created the point mutation result file. Must have an unique identifier inside every comment line to identify the sequence by type. See details. Homozygot patiens have to have a 00 instead of the HLA type. For reference please look in example file.}
  \item{epi_results}{the result file from epitope analysis through epitope finder or a similar csv file with a column with sequence positions and a column with corresponding p-values}
  \item{threshold}{the threshold for the selection which position should enter the tuple process}
  \item{min_number_of_ele_in_tupel}{the minimum number of members in the tuple}
  \item{max_number_of_ele_in_tuple}{the maximum number of members in the tuple}
  \item{result_filename}{the filename of the output file}
  \item{column}{the column in which the type is located for which the analysis should be done}
  \item{column_of_values}{the column from which the (p) values should be taken. See details.}
  \item{column_of_aas}{the column from which the amino acids should be taken. See details.}
  \item{A11}{the position of the start of the first HLA A Allel in the description block of the FASTA file.}
  \item{A12}{the position of the end of the first HLA A Allel in the description block of the FASTA file.}
  \item{A21}{the position of the start of the second HLA A Allel in the description block of the FASTA file.}
  \item{A22}{the position of the end of the second HLA A Allel in the description block of the FASTA file.}
  \item{B11}{the position of the start of the first HLA B Allel in the description block of the FASTA file.}
  \item{B12}{the position of the end of the first HLA B Allel in the description block of the FASTA file.}
  \item{B21}{the position of the start of the second HLA B Allel in the description block of the FASTA file.}
  \item{B22}{the position of the end of the second HLA B Allel in the description block of the FASTA file.}
  \item{one_ident}{if there is only one feature}
  \item{idet}{the feature identifier which should be analzed. See details.}
}
\details{For this to work properly you must have an csv file in which at least one column describes the sequence position and another the p-value at this position, as well as an fasta file from which those file originated. 
If you have only one kind of feature, the please set one_ident to TRUE. This can be if you are to analyse tropisms or treatment with only Yes and No. If you have only one feature it has to be seperated from the rest of the comment line in the fasta file at the end with a semicolon. If you have two kinds of feature (like HLA types) you can give there exact position inside the comment line.
It takes only those sequence positions which have a p-value lower than the threshold. Be aware that it takes a realy long time to calculate all positions with all other positions and that it could result in low quality data!
Please be also aware that it just uses the position in your csv file. If this is NOT the correct position, because you removed some positions, then please correct them first in the csv file or keep in mind, that the result file will also be just the relative position from your fasta file. Please use the same fasta file you had for epitope analysis! The values for the analysis to estimate which positions should be included are normaly the corrected p-values but can be anything else. Even an own added column. You have to add the column with the amino acids with the lowest p-value, that the algorithm know which amino acid to search for.
The tuples can be anything from 2 to number of rows in the csv result file, but be aware that larger tuples may take some time to calculate.
Please use a "B" for unknown amino acid. 
}

\value{A csv list of tupel positions and the p-value of there association.}

\author{Bettina Budeus}
\note{If a patient has a homozygot HLA Allel, then please change the second one to "00" (without ") instead!}


\seealso{\code{\link{assocpointtuple}}}
\examples{
\dontrun{
	ex <- system.file("extdata", "Example_aa.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "epitope_results.csv", package="SeqFeatR")
	assoctuple(ex,
	ep,
	threshold = 0.2,
	min_number_of_ele_in_tupel = 2,
	max_number_of_ele_in_tuple = 2,
	result_filename = "shared_mutations_result.csv",
	column = 9,
	column_of_values = 9,
	column_of_aas = 12,
	A11 = 10,
	A12 = 11, 
	A21 = 13, 
	A22 = 14, 
	B11 = 17, 
	B12 = 18, 
	B21 = 20, 
	B22 = 21,
	FALSE,
	"")
}
}
