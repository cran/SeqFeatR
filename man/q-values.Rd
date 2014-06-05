\name{q-values}
\alias{qvalues}
\title{Calculate q-values from p-values}
\description{calculates and adds q_values to a list of given p_values in csv format with p_value as column header for this column.}

\usage{qvalues(path_to_file = NULL, save_name)}
\arguments{
	\item{path_to_file}{file with a column with p-values. Has to be in csv and contain a column called p_value. For reference please look in example file.}
	\item{save_name}{output file name}
}
\details{Takes a csv file with a column which is called p_values and uses the qvalues package to calculate from this column the corresponding q-values. Uses qvalue package and the calculation within this package to estimate the q-values.}

\value{returns the input file as matrix with added column with q-values}


\author{Bettina Budeus}

\examples{mut <- system.file("extdata", "co_mutation_results.csv", package="SeqFeatR")
	qvalues(mut, "csv_with_q_values.csv")
}



