# marloweU
What the Package Does (Title Case)

This is my attempt at an extension of the MARLOWE code developed by PNNL.
Their program uses a Candidate MySql database which was constructed by digesting .ent
files in-silico from the kegg.jp database.  Basically, taking the amino acid sequences and splitting
them using an R OrgMassSpecR::digest function that digests with Trypsin.

The MARLOWE algorithm takes output data from MassSpectrometry performed on an unknown sample.
It trys to find likely organisms which are contained in the sample.

My goal is to create a new Candidate MySql database with source data from Uniprot.
I also may digest the amino acid sequences using multiple proteases creating options for searches
when the sample has been digested something besides Trypsin.
