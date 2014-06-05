\name{orPlot}
\alias{orPlot}
\title{Create a second graphical output for point mutation vs. feature result}
\description{This tool creates a graphic ion which you can see not only the best p.values for every position of your sequence, but also the direction in which the p.value points. See details for more infomation.}
\usage{orPlot(input_file = NULL,
	save_name, seperator = ";",
	number_of_cases = 2, odds.position = c(6),
	p.value.position = c(2), name = c(2), freq = 7,
	aminoacid.position = NULL, high_log_p = 50, intervall = 10, has_color = 0)}
\arguments{
  \item{input_file}{the input file with data.  For reference please look in example file.}
  \item{save_name}{the name of the save file. Will be processed further with numbers. See details.}
  \item{seperator}{value with which the csv is seperated.}
  \item{number_of_cases}{the number of different search types like different HLA genes or tropism.}
  \item{odds.position}{columnnumber of the odds ratio to determine if up or down. Can be a repeated value as every 4th column.}
  \item{p.value.position}{columnnumber of the p-value to determine the height. Can be a repeated value as every 4th column.}
  \item{name}{the position of the name for the right lable.}
  \item{freq}{frequency of repeat in columns with more than two cases.}
  \item{aminoacid.position}{columnnumber of the amino acid sequence to show above the seq. Can be a repeated value as every 4th column.}
  \item{high_log_p}{the estimated guess of the best p-value. Usuall the number behind the e as positiv number.}
  \item{intervall}{the intervall on the y axis}
  \item{has_color}{if there should be presented another kind of information as color. See details.}
}
\details{This tool took the results from the search for epitopes and combines them in an graphical output where each position in the sequence has its own bar. The heigt of the bar is the log10 of the p-value, the direction is from the odds ratio (below one is down, above one is up). If the user adds a column number for an amino acid or DNA base, this sequence is added above each bar. Further sequences (such as a kind of reference or consensus sequence) besides the one created with find_possible_epiopes have to be added manually to the input *.csv file.
If your bar plot should be additionally colored by another type of information - like eg entropy of the sequence - you can add the position of the column with this information in "has_color". If you do not want any color, than just insert 0. The column should contain only numbers between zero and one.}

\note{The reference sequence has to be aligned with the sequences used for possible epitope analysis to be shown correctly.}

\author{Bettina Budeus}

\seealso{\code{\link{assocpoint}}}

\examples{ep <- system.file("extdata", "epitope_results.csv", package="SeqFeatR")
	en <- system.file("extdata", "final_graphic", package="SeqFeatR")
	orPlot(ep, 
	en,
	";",
	8,
	6,
	2,
	2,
	7,
	NULL,
	4,
	0.2,
	0)
}
