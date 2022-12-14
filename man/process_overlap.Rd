% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_database.linux.R
\name{process_overlap}
\alias{process_overlap}
\title{Populate the database with overlap coefficient of peptides between all organisms
This method is multi threaded.  Each call analyze_overlap creates its own temperary table to load values
into, once all overlap coefficients are calculated for that call the value are copied to the
overlap_coefficients table.  The number times analyze_overlap is called defined by the core_count variable,
the default value is 4.  If a file_output value is provided the the file
is appended to with process information messages, otherwise the messages will be sent as the return value
once all values are calculated and saved.}
\usage{
process_overlap(
  conn,
  core_count = 2,
  split_size = 10,
  start_index = 1,
  file_output = ""
)
}
\arguments{
\item{conn}{database connection values to create a MySQL connection}

\item{core_count}{This is to total number cores to be used during parallel processing
it should be no greater than 1/2 of total cores}

\item{split_size}{This is size of list that should be passed to each execution of
analyze_strong_count.}

\item{start_index}{the minimum organism_id to start processing at, this can be used
for times when this method fails while processing to allow it to continue
where it left off}

\item{file_output}{the name of the file to output all messages to, if not provided
the return value will be the debug output}
}
\value{
T or debug output if no file_output was provided
}
\description{
\code{process_overlap} calculates the overlap coefficient of peptides between all organisms.
}
\author{
Dustin Crockett
}
