% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/digest_ent_protein.R
\name{digest_ent_protein}
\alias{digest_ent_protein}
\title{Digest one protein from KEGG ENT file}
\usage{
digest_ent_protein(protein_info, hydrogen_mass = 1.00727646627)
}
\arguments{
\item{protein_info}{a character vector of length 2 containing the protein_id
and protein sequence}

\item{hydrogen_mass}{mass of hydrogen to use when calculating neutral mass.
Default is 1.00727646627 which comes from
\href{https://physics.nist.gov/cgi-bin/cuu/Category?view=html&Atomic+and+nuclear.x=123&Atomic+and+nuclear.y=21}{NIST proton molar mass}}
}
\value{
Returns a dataframe with the protein ID, peptide string, and peptide
mass.
}
\description{
This function is meant to be called as part of parse_ent() which parses the
KEGG database ENT files to prepare them for upload to the MySQL database.
}
\author{
Sarah C. Jenson
}
