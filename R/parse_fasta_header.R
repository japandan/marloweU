#' parse_fasta_header.R
#'
#' Parse a UniProt uniref or uniprotkb FASTA AAStringSet
#'
#'# input variables
#' @param AA_list is an AAStringSet with headers from uniprotkb or uniref.
#'
#' @return Default is to return a 4 column matrix to the caller with id,protein name, organism, TaxonID
#' @export
#'
#' @author Daniel T. Vogel
#'
#' @importFrom plyr adply ldply
#' @importFrom stringr str_match str_trim str_replace_all str_extract_all str_match_all
#' @importFrom assertthat assert_that
#' @importFrom Biostrings readAAStringSet
#
parse_fasta_header <- function( AA_list ){

  #checking inputs
  # Add assert to verify that AA_list is an AAStringSet

  # UniProtKB Fasta header format
  #
  # These files, composed of canonical and additional sequences, are non-redundant
  # FASTA sets for the sequences of each reference proteome.
  # The additional set contains isoform/variant sequences for a given gene, and its
  # FASTA header indicates the corresponding canonical sequence ("Isoform of ...").
  # The FASTA format is the standard UniProtKB format.
  #
  # For further references about the standard UniProtKB format, please see:
  #  http://www.uniprot.org/help/fasta-headers
  #  http://www.uniprot.org/faq/38
  #
  # >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
  #
  # >tr|B9R7K7|B9R7K7_RICCO Cytochrome P450 OS=Ricinus communis OX=3988 GN=RCOM_1592680 PE=3 SV=1
  # >tr|B9R8T7|B9R8T7_RICCO Cinnamoyl-CoA reductase, putative OS=Ricinus communis OX=3988 GN=RCOM_1602080 PE=4 SV=1
  # >tr|B9R9L1|B9R9L1_RICCO 3-ketoacyl-CoA synthase OS=Ricinus communis OX=3988 GN=RCOM_1498550 PE=3 SV=1
  # >tr|B9RAM4|B9RAM4_RICCO 26S proteasome non-atpase regulatory subunit, putative OS=Ricinus communis OX=3988 GN=RCOM_1507340 PE=4 SV=1
  # >tr|B9RBT9|B9RBT9_RICCO C-4 methyl sterol oxidase, putative OS=Ricinus communis OX=3988 GN=RCOM_1681040 PE=3 SV=1
  #
  # Where:
  #
  # db is 'sp' for UniProtKB/Swiss-Prot and 'tr' for UniProtKB/TrEMBL.
  # UniqueIdentifier is the primary accession number of the UniProtKB entry.
  # EntryName is the entry name of the UniProtKB entry.
  # ProteinName is the recommended name of the UniProtKB entry as annotated in the RecName field. For UniProtKB/TrEMBL entries without a RecName field, the SubName field is used. In case of multiple SubNames, the first one is used. The 'precursor' attribute is excluded, 'Fragment' is included with the name if applicable.
  # OrganismName is the scientific name of the organism of the UniProtKB entry.
  # OrganismIdentifier is the unique identifier of the source organism, assigned by the NCBI.
  # GeneName is the first gene name of the UniProtKB entry. If there is no gene name, OrderedLocusName or ORFname, the GN field is not listed.
  # ProteinExistence is the numerical value describing the evidence for the existence of the protein.
  # SequenceVersion is the version number of the sequence.
  #
  # UniRef Fasta header format
  #
  # UniqueIdentifier is the primary accession number of the UniRef cluster.
  # ClusterName is the name of the UniRef cluster.
  # Members is the number of UniRef cluster members.
  # TaxonName is the scientific name of the lowest common taxon shared by all UniRef cluster members.
  # TaxonIdentifier is the NCBI taxonomy identifier of the lowest common taxon shared by all UniRef cluster members.
  # RepresentativeMember is the entry name of the representative member of the UniRef cluster.
  # e.g.
  # >UniqueIdentifier ClusterName n=Members Tax=TaxonName TaxID=TaxonIdentifier RepID=RepresentativeMember
  # >UniRef50_A0A5A9P0L4 Peptidylprolyl isomerase n=1 Tax=Triplophysa tibetana TaxID=1572043 RepID=A0A5A9P0L4_9TELE

  # Check if format is uniref or uniprotkb
  fastaformat <- "uniprotkb"
  if ( tolower( substr( names(AA_list)[1], 1, 6 ) ) == "uniref" ) {
    fastaformat <- "uniref"
  }

  field_matrix <- NULL

    if ( fastaformat == "uniref") {
    #Use regular expressions to pull out the needed fields.  Need to be returned as.matrix

        #protein_id. <- UniqueIdentifier (first field, can contain underlines)
        field_matrix <- stringr::str_match( names(AA_list), "(^UniRef\\S+)\\s(.+)\\sn\\=\\d+.+Tax=(.+)\\sTaxID=(.+)\\s.+")
    }


    # >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier
    if ( fastaformat == "uniprotkb") {
        #Use regular expressions to pull out the needed fields.  Need to be returned as.matrix
        #field_matrix <- stringr::str_match( names(AA_list), "\\|(.+)\\|(.+)\\s(\\w+)\\sOS=(\\w+)\\sOX=(\\w)\\s" )
        field_matrix <- stringr::str_match( names(AA_list), "\\|(\\w+)\\|\\w+\\s(.+)\\sOS=(.+)\\sOX=(\\w+)\\s.+" )
    }


  # remove the first column containing the regex match and keep this in matrix format
  field_matrix <-  field_matrix[ , -1, drop=FALSE ]

  return( field_matrix )


}





