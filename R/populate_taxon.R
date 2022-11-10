
#'This method is executed after the Kegg database is fully loaded.  This method populates
#'taxon information into the database in adds taxon ids for each organism
#'A lineage table is also loaded that provides a full parent lineage for a taxon id
#'
#'`populate_taxon` Populates taxon information into a fully populated database.
#'
#'@param conn database connection values to create a MySQL connection
#'@param all_names a dataframe of all taxon names
#'@param all_nodes a dataframe of all taxons
#'@param kegg_org_taxon a dataframe with the taxon id for each organism in kegg
#'@param lineage_table a dataframe of the entire lineage for each taxon id
#'
#'@return T
#'
#'@importFrom assertthat assert_that
#'@importFrom DBI dbWriteTable dbExecute dbConnect dbDisconnect
#'@importFrom RMariaDB MariaDB
#'
#'@author Dustin Crockett
#'
#'@export
populate_taxon <- function(conn, all_names, all_nodes, kegg_org_taxon, lineage_table) {
  cat("============Populating taxon information============\n")
  begin_time <- Sys.time()
  start_time <- begin_time

  #validate the conn list
  # verify that the conn is a list
  assertthat::assert_that(is.list(conn),
                          msg = "conn is not a list. Make sure conn is a list.")
  # verify that conn has the correct fields
  assertthat::assert_that("dbname" %in% names(conn),
                          msg = "dbname is missing from the list conn")
  assertthat::assert_that("host" %in% names(conn),
                          msg = "host is missing from the list conn")
  assertthat::assert_that("port" %in% names(conn),
                          msg = "port is missing from the list conn")
  assertthat::assert_that("user" %in% names(conn),
                          msg = "user is missing from the list conn")
  assertthat::assert_that("password" %in% names(conn),
                          msg = "password is missing from the list conn")



  # verify that input values are dataframes has the correct fields
  assertthat::assert_that(is.data.frame(all_names),
                          msg = "all_names is not a dataframe")
  assertthat::assert_that(is.data.frame(all_nodes),
                          msg = "all_nodes is not a dataframe")
  assertthat::assert_that(is.data.frame(kegg_org_taxon),
                          msg = "kegg_org_taxon is not a dataframe")
  assertthat::assert_that(is.data.frame(lineage_table),
                          msg = "lineage_table is not a dataframe")


  #verify that the data frame all_names has the correct columns
  assertthat::assert_that("taxon_id" %in% colnames(all_names),
                          msg = "taxon_id column missing from all_names")
  assertthat::assert_that("name_txt" %in% colnames(all_names),
                          msg = "name_txt column missing from all_names")
  assertthat::assert_that("unique_name" %in% colnames(all_names),
                          msg = "unique_name column missing from all_names")
  assertthat::assert_that("name_class" %in% colnames(all_names),
                          msg = "name_class column missing from all_names")

  #verify that the data frame all_nodes has the correct columns
  assertthat::assert_that("taxon_id" %in% colnames(all_nodes),
                          msg = "taxon_id column missing from all_nodes")
  assertthat::assert_that("parent_taxon_id" %in% colnames(all_nodes),
                          msg = "parent_taxon_id column missing from all_nodes")
  assertthat::assert_that("rank" %in% colnames(all_nodes),
                          msg = "rank column missing from all_nodes")
  assertthat::assert_that("embl_code" %in% colnames(all_nodes),
                          msg = "embl_code column missing from all_nodes")
  assertthat::assert_that("division_id" %in% colnames(all_nodes),
                          msg = "division_id column missing from all_nodes")
  assertthat::assert_that("inherited_div_flag" %in% colnames(all_nodes),
                          msg = "inherited_div_flag column missing from all_nodes")
  assertthat::assert_that("genetic_code_id" %in% colnames(all_nodes),
                          msg = "genetic_code_id column missing from all_nodes")
  assertthat::assert_that("inherited_GC_flag" %in% colnames(all_nodes),
                          msg = "inherited_GC_flag column missing from all_nodes")
  assertthat::assert_that("mitochondrial_genetic_code_id" %in% colnames(all_nodes),
                          msg = "mitochondrial_genetic_code_id column missing from all_nodes")
  assertthat::assert_that("inherited_MGC_flag" %in% colnames(all_nodes),
                          msg = "inherited_MGC_flag column missing from all_nodes")
  assertthat::assert_that("GenBank_hidden_flag" %in% colnames(all_nodes),
                          msg = "GenBank_hidden_flag column missing from all_nodes")
  assertthat::assert_that("hidden_subtree_root_flag" %in% colnames(all_nodes),
                          msg = "hidden_subtree_root_flag column missing from all_nodes")
  assertthat::assert_that("comments" %in% colnames(all_nodes),
                          msg = "comments column missing from all_nodes")
  assertthat::assert_that("plastid_genetic_code_id" %in% colnames(all_nodes),
                          msg = "plastid_genetic_code_id column missing from all_nodes")
  assertthat::assert_that("inherited_PGC_flag" %in% colnames(all_nodes),
                          msg = "inherited_PGC_flag column missing from all_nodes")
  assertthat::assert_that("specified_species" %in% colnames(all_nodes),
                          msg = "specified_species column missing from all_nodes")
  assertthat::assert_that("hydrogenosome_genetic_code_id" %in% colnames(all_nodes),
                          msg = "hydrogenosome_gentic_code_id column missing from all_nodes")
  assertthat::assert_that("inherited_HCG_flag" %in% colnames(all_nodes),
                          msg = "inherited_HCG_flag column missing from all_nodes")

  #verify that the data frame kegg_org_taxon has the correct columns
  assertthat::assert_that("Org_code" %in% colnames(kegg_org_taxon),
                          msg = "Org_code column missing from kegg_org_taxon")
  assertthat::assert_that("taxon_id" %in% colnames(kegg_org_taxon),
                          msg = "taxon_id column missing from create")

  #verify that the data frame lineage_table has the correct columns
  assertthat::assert_that("node_taxon_id" %in% colnames(lineage_table),
                          msg = "node_taxon_id column missing from lineage_table")
  assertthat::assert_that("node_tax_name" %in% colnames(lineage_table),
                          msg = "node_tax_name column missing from lineage_table")
  assertthat::assert_that("node_rank" %in% colnames(lineage_table),
                          msg = "node_rank column missing from lineage_table")
  assertthat::assert_that("superkingdom" %in% colnames(lineage_table),
                          msg = "superkingdom column missing from lineage_table")
  assertthat::assert_that("superkingdom_name" %in% colnames(lineage_table),
                          msg = "superkingdom_name column missing from lineage_table")
  assertthat::assert_that("kingdom" %in% colnames(lineage_table),
                          msg = "kingdom column missing from lineage_table")
  assertthat::assert_that("kingdom_name" %in% colnames(lineage_table),
                          msg = "kingdom_name column missing from lineage_table")
  assertthat::assert_that("superphylum" %in% colnames(lineage_table),
                          msg = "superphylum column missing from lineage_table")
  assertthat::assert_that("superphylum_name" %in% colnames(lineage_table),
                          msg = "superphylum_name column missing from lineage_table")
  assertthat::assert_that("phylum" %in% colnames(lineage_table),
                          msg = "phylum column missing from lineage_table")
  assertthat::assert_that("phylum_name" %in% colnames(lineage_table),
                          msg = "phylum_name column missing from lineage_table")
  assertthat::assert_that("subphylum" %in% colnames(lineage_table),
                          msg = "subphylum column missing from lineage_table")
  assertthat::assert_that("subphylum_name" %in% colnames(lineage_table),
                          msg = "subphylum_name column missing from lineage_table")
  assertthat::assert_that("superclass" %in% colnames(lineage_table),
                          msg = "superclass column missing from lineage_table")
  assertthat::assert_that("superclass_name" %in% colnames(lineage_table),
                          msg = "superclass_name column missing from lineage_table")
  assertthat::assert_that("class" %in% colnames(lineage_table),
                          msg = "class column missing from lineage_table")
  assertthat::assert_that("class_name" %in% colnames(lineage_table),
                          msg = "class_name column missing from lineage_table")
  assertthat::assert_that("subclass" %in% colnames(lineage_table),
                          msg = "subclass column missing from lineage_table")
  assertthat::assert_that("subclass_name" %in% colnames(lineage_table),
                          msg = "subclass_name column missing from lineage_table")
  assertthat::assert_that("infraclass" %in% colnames(lineage_table),
                          msg = "infraclass column missing from lineage_table")
  assertthat::assert_that("infraclass_name" %in% colnames(lineage_table),
                          msg = "infraclass_name column missing from lineage_table")
  assertthat::assert_that("cohort" %in% colnames(lineage_table),
                          msg = "cohort column missing from lineage_table")
  assertthat::assert_that("cohort_name" %in% colnames(lineage_table),
                          msg = "cohort_name column missing from lineage_table")
  assertthat::assert_that("subcohort" %in% colnames(lineage_table),
                          msg = "subcohort column missing from lineage_table")
  assertthat::assert_that("subcohort_name" %in% colnames(lineage_table),
                          msg = "subcohort_name column missing from lineage_table")
  assertthat::assert_that("superorder" %in% colnames(lineage_table),
                          msg = "superorder column missing from lineage_table")
  assertthat::assert_that("superorder_name" %in% colnames(lineage_table),
                          msg = "superorder_name column missing from lineage_table")
  assertthat::assert_that("order" %in% colnames(lineage_table),
                          msg = "order column missing from lineage_table")
  assertthat::assert_that("order_name" %in% colnames(lineage_table),
                          msg = "order_name column missing from lineage_table")
  assertthat::assert_that("suborder" %in% colnames(lineage_table),
                          msg = "suborder column missing from lineage_table")
  assertthat::assert_that("suborder_name" %in% colnames(lineage_table),
                          msg = "suborder_name column missing from lineage_table")
  assertthat::assert_that("infraorder" %in% colnames(lineage_table),
                          msg = "infraorder column missing from lineage_table")
  assertthat::assert_that("infraorder_name" %in% colnames(lineage_table),
                          msg = "infraorder_name column missing from lineage_table")
  assertthat::assert_that("superfamily" %in% colnames(lineage_table),
                          msg = "superfamily column missing from lineage_table")
  assertthat::assert_that("superfamily_name" %in% colnames(lineage_table),
                          msg = "superfamily_name column missing from lineage_table")
  assertthat::assert_that("family" %in% colnames(lineage_table),
                          msg = "family column missing from lineage_table")
  assertthat::assert_that("family_name" %in% colnames(lineage_table),
                          msg = "family_name column missing from lineage_table")
  assertthat::assert_that("subfamily" %in% colnames(lineage_table),
                          msg = "subfamily column missing from lineage_table")
  assertthat::assert_that("subfamily_name" %in% colnames(lineage_table),
                          msg = "subfamily_name column missing from lineage_table")
  assertthat::assert_that("tribe" %in% colnames(lineage_table),
                          msg = "tribe column missing from lineage_table")
  assertthat::assert_that("tribe_name" %in% colnames(lineage_table),
                          msg = "tribe_name column missing from lineage_table")
  assertthat::assert_that("subtribe" %in% colnames(lineage_table),
                          msg = "subtribe column missing from lineage_table")
  assertthat::assert_that("subtribe_name" %in% colnames(lineage_table),
                          msg = "subtribe_name column missing from lineage_table")
  assertthat::assert_that("genus" %in% colnames(lineage_table),
                          msg = "genus column missing from lineage_table")
  assertthat::assert_that("genus_name" %in% colnames(lineage_table),
                          msg = "genus_name column missing from lineage_table")
  assertthat::assert_that("subgenus" %in% colnames(lineage_table),
                          msg = "subgenus column missing from lineage_table")
  assertthat::assert_that("subgenus_name" %in% colnames(lineage_table),
                          msg = "subgenus_name column missing from lineage_table")
  assertthat::assert_that("section" %in% colnames(lineage_table),
                          msg = "section column missing from lineage_table")
  assertthat::assert_that("section_name" %in% colnames(lineage_table),
                          msg = "section_name column missing from lineage_table")
  assertthat::assert_that("subsection" %in% colnames(lineage_table),
                          msg = "subsection column missing from lineage_table")
  assertthat::assert_that("subsection_name" %in% colnames(lineage_table),
                          msg = "subsection_name column missing from lineage_table")
  assertthat::assert_that("series" %in% colnames(lineage_table),
                          msg = "series column missing from lineage_table")
  assertthat::assert_that("series_name" %in% colnames(lineage_table),
                          msg = "series_name column missing from lineage_table")
  assertthat::assert_that("species_group" %in% colnames(lineage_table),
                          msg = "species_group column missing from lineage_table")
  assertthat::assert_that("species_group_name" %in% colnames(lineage_table),
                          msg = "species_group_name column missing from lineage_table")
  assertthat::assert_that("species_subgroup" %in% colnames(lineage_table),
                          msg = "species_subgroup column missing from lineage_table")
  assertthat::assert_that("species_subgroup_name" %in% colnames(lineage_table),
                          msg = "species_subgroup_name column missing from lineage_table")
  assertthat::assert_that("species" %in% colnames(lineage_table),
                          msg = "species column missing from lineage_table")
  assertthat::assert_that("species_name" %in% colnames(lineage_table),
                          msg = "species_name column missing from lineage_table")
  assertthat::assert_that("subspecies" %in% colnames(lineage_table),
                          msg = "subspecies column missing from lineage_table")
  assertthat::assert_that("subspecies_name" %in% colnames(lineage_table),
                          msg = "subspecies_name column missing from lineage_table")
  assertthat::assert_that("varietas" %in% colnames(lineage_table),
                          msg = "varietas column missing from lineage_table")
  assertthat::assert_that("varietas_name" %in% colnames(lineage_table),
                          msg = "varietas_name column missing from lineage_table")
  assertthat::assert_that("forma" %in% colnames(lineage_table),
                          msg = "forma column missing from lineage_table")
  assertthat::assert_that("forma_name" %in% colnames(lineage_table),
                          msg = "forma_name column missing from lineage_table")

  con <- DBI::dbConnect(RMariaDB::MariaDB(),
                        dbname = conn$dbname,
                        host = conn$host,
                        port = conn$port,
                        user = conn$user,
                        password = conn$password)
  cat("Loading all_names data\n")

  #With the data validated, first we upload all of the data into temporary tables
  DBI::dbWriteTable(
    conn = con,
    DBI::Id(schema="candidate",table="upload_all_names"),
    value = all_names,
    overwrite = T,
    temporary = T,
  )
  end_time <- Sys.time()
  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(start_time, units = "secs"))
  cat(paste0("Finished uploading all_names ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds, starting all_nodes\n"))
  start_time <- end_time

  DBI::dbWriteTable(
    conn = con,
    DBI::Id(schema="candidate",table="upload_all_nodes"),
    value = all_nodes,
    overwrite = T,
    temporary = T,
  )

  end_time <- Sys.time()
  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(start_time, units = "secs"))
  cat(paste0("Finished uploading all_nodes ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds, starting kegg_org_taxon\n"))
  start_time <- end_time

  DBI::dbWriteTable(
    conn = con,
    DBI::Id(schema="candidate",table="upload_kegg_org_taxon"),
    value = kegg_org_taxon,
    overwrite = T,
    temporary = T,
  )

  end_time <- Sys.time()
  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(start_time, units = "secs"))
  cat(paste0("Finished uploading kegg_org_taxon ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds, starting lineage_table\n"))
  start_time <- end_time

  DBI::dbWriteTable(
    conn = con,
    DBI::Id(schema="candidate",table="upload_lineage_table"),
    value = lineage_table,
    overwrite = T,
    temporary = T,
  )

  end_time <- Sys.time()
  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(start_time, units = "secs"))
  cat(paste0("Finished uploading lineage_table ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds, starting taxon population\n"))
  start_time <- end_time

  #next we migrate the data into properly named and indexed tables
  DBI::dbExecute(con, 'insert into candidate.taxons (taxon_id, parent_taxon_id, `rank`, embl_code, division_id , inherited_div_flag,genetic_code_id,inherited_Gc_flag,
			mitochondrial_genetic_code_id,inherited_MGC_flag,GenBank_hidden_flag,hidden_subtree_root_flag,comments,plastid_genetic_code_id,
			inherited_PGC_flag,specified_species,hydrogenosome_genetic_code_id,inherited_HCG_flag)
		select taxon_id, parent_taxon_id, `rank`, embl_code, division_id , inherited_div_flag,genetic_code_id,inherited_Gc_flag,
			mitochondrial_genetic_code_id,inherited_MGC_flag,GenBank_hidden_flag,hidden_subtree_root_flag,comments,plastid_genetic_code_id,
			inherited_PGC_flag,specified_species,hydrogenosome_genetic_code_id,inherited_HCG_flag from candidate.upload_all_nodes')

  end_time <- Sys.time()
  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(start_time, units = "secs"))
  cat(paste0("Finished populating taxon ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds, starting taxon_names population\n"))
  start_time <- end_time

  DBI::dbExecute(con, 'insert into candidate.taxon_names(taxon_id, name_txt, unique_name, name_class)
	select taxon_id, name_txt, unique_name, name_class from candidate.upload_all_names ')

  end_time <- Sys.time()
  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(start_time, units = "secs"))
  cat(paste0("Finished populating taxon_names ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds, starting taxon_lineages population\n"))
  start_time <- end_time

  DBI::dbExecute(con, 'insert into candidate.taxon_lineages (taxon_id, taxon_name, taxon_rank,superkingdom,kingdom,subkingdom,superphylum,phylum,subphylum ,
            superclass ,class , subclass ,infraclass ,cohort,superorder,`order`,suborder,infraorder,parvorder,superfamily,
            family,subfamily, tribe,subtribe,genus,subgenus,`section`,subsection,series,species_group,species_subgroup, species,subspecies,
            varietas,forma,superkingdom_name,kingdom_name,subkingdom_name,superphylum_name,phylum_name,subphylum_name ,
            superclass_name ,class_name , subclass_name ,infraclass_name ,cohort_name,superorder_name,order_name,suborder_name,infraorder_name,parvorder_name,
            superfamily_name,family_name,subfamily_name, tribe_name,subtribe_name,genus_name,subgenus_name,section_name,subsection_name,series_name,
            species_group_name,species_subgroup_name, species_name,subspecies_name,varietas_name,forma_name)
select node_taxon_id, node_tax_name, node_rank,superkingdom,kingdom,subkingdom,superphylum,phylum,subphylum ,
            superclass ,class , subclass ,infraclass ,cohort,superorder,`order`,suborder,infraorder,parvorder,superfamily,
            family,subfamily, tribe,subtribe,genus,subgenus,`section`,subsection,series,species_group,species_subgroup, species,subspecies,
            varietas,forma,superkingdom_name,kingdom_name,subkingdom_name,superphylum_name,phylum_name,subphylum_name ,
            superclass_name ,class_name , subclass_name ,infraclass_name ,cohort_name,superorder_name,order_name,suborder_name,infraorder_name,parvorder_name,
            superfamily_name,family_name,subfamily_name, tribe_name,subtribe_name,genus_name,subgenus_name,section_name,subsection_name,series_name,
            species_group_name,species_subgroup_name, species_name,subspecies_name,varietas_name,forma_name from candidate.upload_lineage_table ')

  end_time <- Sys.time()
  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(start_time, units = "secs"))
  cat(paste0("Finished populating taxon_lineages ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds, starting organism population\n"))
  start_time <- end_time

  #now we update the organism table with a foreign key to the taxon_id table
  DBI::dbExecute(con, 'update candidate.organisms o
    	join candidate.upload_kegg_org_taxon ukt on ukt.org_code = o.kegg_org_code
    	set o.taxon_id = ukt.taxon_id ')

  end_time <- Sys.time()
  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(start_time, units = "secs"))
  cat(paste0("Finished populating organisms ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds\n"))

  #finally lets clean up the database connection
  DBI::dbDisconnect(con)

  diff_time <- format(as.numeric(end_time, units = "secs") - as.numeric(begin_time, units = "secs"))
  cat(paste0("====Taxon population, started at ", start_time, ", and finished at ", end_time, ", took a total duration of ", diff_time, " seconds====\n"))

  return(T)
}

populate_taxon_id<-function(){
  # MySQL connection credentials.  These need to match your installation
  conn <- list(
    "dbname"= "candidate",
    "host" = "localhost",
    "port" = 3306,
    "user" = "msdba",
    "password" = "MassSpec2021!"
  )

  con <- DBI::dbConnect(RMariaDB::MariaDB(),
                        dbname = conn$dbname,
                        host = conn$host,
                        port = conn$port,
                        user = conn$user,
                        password = conn$password)
  library(tidyverse)
  library(magrittr)

  # get the taxons from a file
  load("data/taxon/taxons.RData")
  # update the database table taxons
  # Error: Cannot drop table 'taxons' referenced by a foreign key constraint 'organisms_ibfk_2' on table 'organisms'. [3730]
  DBI::dbWriteTable(
    conn = con,
    DBI::Id(schema="candidate",table="upload_all_nodes"),
    value = taxons,
    overwrite = T,
    temporary = F,
  )

  #next we migrate the data into properly named and indexed tables
  DBI::dbExecute(con, 'insert into candidate.taxons (taxon_id, parent_taxon_id, `rank`, embl_code, division_id , inherited_div_flag,genetic_code_id,inherited_Gc_flag,
			mitochondrial_genetic_code_id,inherited_MGC_flag,GenBank_hidden_flag,hidden_subtree_root_flag,comments,plastid_genetic_code_id,
			inherited_PGC_flag,specified_species,hydrogenosome_genetic_code_id,inherited_HCG_flag)
		select taxon_id, parent_taxon_id, `rank`, embl_code, division_id , inherited_div_flag,genetic_code_id,inherited_Gc_flag,
			mitochondrial_genetic_code_id,inherited_MGC_flag,GenBank_hidden_flag,hidden_subtree_root_flag,comments,plastid_genetic_code_id,
			inherited_PGC_flag,specified_species,hydrogenosome_genetic_code_id,inherited_HCG_flag from candidate.upload_all_nodes')


  #get the kegg_id from the organisms in the database
  organisms <- DBI::dbGetQuery(con, "select * from candidate.organisms")
  organisms

  #fill in the kegg_org_code with the kegg_id removed "U" prefix character
  sqlupdate<-"update organisms set kegg_org_code = substring(kegg_id,2)"
  DBI::dbExecute(con, sqlupdate)

  # The first time we do this, the taxon_id is NA.
  # Set the taxon_id to be the same as the kegg_id without the "U"
  kegg_org_taxon <- organisms[c("kegg_id","taxon_id")]
  colnames( kegg_org_taxon )[1]<-"org_code"
  kegg_org_taxon$taxon_id <- substring( kegg_org_taxon$org_code, 2)
  kegg_org_taxon$org_code <- kegg_org_taxon$taxon_id
  kegg_org_taxon

  # we need a data structure called kegg_org_taxon which has 2 columns
  #verify that the data frame kegg_org_taxon has the correct columns
  assertthat::assert_that("org_code" %in% colnames(kegg_org_taxon),
                          msg = "Org_code column missing from kegg_org_taxon")
  assertthat::assert_that("taxon_id" %in% colnames(kegg_org_taxon),
                          msg = "taxon_id column missing from create")

  # build a temporary table with kegg_id, taxonid
  DBI::dbWriteTable(
    conn = con,
    DBI::Id(schema="candidate",table="upload_kegg_org_taxon"),
    value = kegg_org_taxon,
    overwrite = T,
    temporary = F,
  )


  #now we update the organism table with a foreign key to the taxon_id table
  DBI::dbExecute(con, 'update candidate.organisms o
    	join candidate.upload_kegg_org_taxon ukt on ukt.org_code = o.kegg_org_code
    	set o.taxon_id = ukt.taxon_id ')

  organisms_after <- DBI::dbGetQuery(con, "select * from candidate.organisms")
  organisms_after

  DBI::dbDisconnect(con)
}

#populate_taxon_id()
