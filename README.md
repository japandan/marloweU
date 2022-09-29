# marloweU
Create a new Candidate MySql database with source data from Uniprot. 
The package contains a function to extract data fields from UniProt FASTA files and build rData files which 
are imported into a complex MariaDB/MySQL database for use with MARLOWE.

This is my attempt at an extension of the MARLOWE code developed by PNNL.
Their program uses a Candidate MySql database which was constructed by digesting .ent
files in-silico from the kegg.jp database.  Basically, taking the amino acid sequences and splitting
them using an R OrgMassSpecR::digest function that digests with Trypsin.

The MARLOWE algorithm takes output data from MassSpectrometry performed on an unknown sample.
It trys to find likely organisms which are contained in the sample.

In addition, I am including documentation and R functions to build the R environment and MySQL database on Linux.
