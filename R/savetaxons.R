
conn<- list(
  "dbname"= "candidate",
  "host" = "localhost",
  "port" = 3306,
  "user" = "root",
  "password" = "MassSpec2021!"
)

con <- DBI::dbConnect(RMariaDB::MariaDB(),
                      dbname = conn$dbname,
                      host = conn$host,
                      port = conn$port,
                      user = conn$user,
                      password = conn$password)


# taxon_lineages <- DBI::dbGetQuery(con, "select * from taxon_lineages")
# taxon_names <- DBI::dbGetQuery(con, "select * from taxon_names")
# taxons <- DBI::dbGetQuery(con, "select * from taxons" )
# taxon_lineages_ids <- DBI::dbGetQuery(con, "select * from taxon_lineages_ids")
#
# save( taxon_lineages, file="D:\\taxon_lineages.RData" )
# save( taxon_lineages_ids, file="D:\\taxon_lineages_ids.RData" )
# save( taxons, file="D:\\taxons.RData" )
# save( taxon_names, file="D:\\taxon_names.RData" )
#
# table_taxon_lineages<-dbReadTable( con, "taxon_lineages")
# table_taxon_lineages_ids<-dbReadTable( con, "taxon_lineages_ids")
# table_taxons<-dbReadTable( con, "taxons")
# table_taxon_names<-dbReadTable( con, "taxon_names")
#
# save( table_taxon_lineages, file="D:\\table_taxon_lineages.RData" )
# save( table_taxon_lineages_ids, file="D:\\table_taxon_lineages_ids.RData" )
# save( table_taxons, file="D:\\table_taxons.RData" )
# save( table_taxon_names, file="D:\\table_taxon_names.RData" )
#
# dbWriteTable(con, "taxon_lineages", taxon_lineages, overwrite=TRUE )
# dbWriteTable(con, "taxon_lineages_ids", taxon_lineages_ids, overwrite=TRUE )
# dbWriteTable(con, "taxon_names", taxon_names, overwrite=TRUE )
# dbWriteTable(con, "taxons", taxons, overwrite=TRUE)

# original code to update taxon_id
#now we update the organism table with a foreign key to the taxon_id table
# DBI::dbExecute(con, 'update candidate.organisms o
#     	join candidate.upload_kegg_org_taxon ukt on ukt.org_code = o.kegg_org_code
#     	set o.taxon_id = ukt.taxon_id ')


DBI::dbExecute(con, 'update candidate.organisms o
    	set o.taxon_id = substring( o.kegg_id, 2 )')
# Error: Cannot add or update a child row: a foreign key constraint fails (`candidate`.`organisms`, CONSTRAINT `organisms_ibfk_2` FOREIGN KEY (`taxon_id`) REFERENCES `taxons` (`taxon_id`)) [1452]

# > dbWriteTable(con, "taxon_lineages_ids", taxon_lineages_ids, overwrite=TRUE )
# Error: Table 'taxon_lineages_ids' already exists [1050]
# > dbWriteTable(con, "taxon_names", taxon_names, overwrite=TRUE )
# > dbWriteTable(con, "taxons", taxons, overwrite=TRUE)
# Error: Cannot drop table 'taxons' referenced by a foreign key constraint 'organisms_ibfk_2' on table 'organisms'. [3730]
