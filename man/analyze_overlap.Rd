% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_database.linux.R
\name{analyze_overlap}
\alias{analyze_overlap}
\title{A helper method that allows for multi-threading of calculating of overlap coefficients
This is an internal method used by process_overlap and requires
the conn list of connection information to be passed into it by way of
cluster import.  This method will loop through all organisms_id values in in_set
and populate a temporary table that stores the percentage of peptide overlap between
two organisms in the.  Before exiting this method all calculated values are then
copied to the overlap_coefficient table.  If a file_output string is passed using the
cluster import the process information is appended to the end of the file otherwise
the process information is the return value for the method.}
\usage{
analyze_overlap(in_set, conn, file_output = "")
}
\arguments{
\item{in_set}{a vector of organism_id values from the database}

\item{conn}{database connection values to create a MySQL connection}

\item{file_output}{the name of the file to output all messages to, if not provided
the return value will be the debug output}
}
\value{
T if a file_output was privided or process information messages if not provided
}
\description{
A helper method that allows for multi-threading of calculating of overlap coefficients
This is an internal method used by process_overlap and requires
the conn list of connection information to be passed into it by way of
cluster import.  This method will loop through all organisms_id values in in_set
and populate a temporary table that stores the percentage of peptide overlap between
two organisms in the.  Before exiting this method all calculated values are then
copied to the overlap_coefficient table.  If a file_output string is passed using the
cluster import the process information is appended to the end of the file otherwise
the process information is the return value for the method.
}
\author{
Dustin Crockett
}
