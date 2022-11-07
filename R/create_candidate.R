# this will create a database
library(CandidateSearchDatabase)
conn<- list("dbname" = "candidate", "host" = "localhost", "port" = 3306, "user" = "msdba", "password" = "MassSpec2021!")

## creates the database ##
create_datamodel( conn )

delete_database( conn, T )
create_datamodel( conn= conn, large_storage_found = T,  large_storage = "K:\\mysql_data" )
## "C:\Program Files\MySQL\MySQL Server 8.0\bin\mysql" -e "show tables;" candidate -u root -p

## test with show tables;
## "C:\Program Files\MySQL\MySQL Server 8.0\bin\mysql" -e "show tables;" candidate -u root -p > dumpoutput 2>&1
## Execute the database restore .sql script
## "C:\Program Files\MySQL\MySQL Server 8.0\bin\mysql" -e "source dump.sql" candidate -u root -p > dumpoutput 2>&1


library( DBI )
#download taxon datasets
conn
con <- DBI::dbConnect(RMariaDB::MariaDB(),
                      dbname = conn$dbname,
                      host = conn$host,
                      port = conn$port,
                      user = conn$user,
                      password = conn$password)


#for the query to work with thousands of organisms we first need to update group_concat_max_len
DBI::dbExecute(con, 'SET SESSION group_concat_max_len = 1000000')

#for the query to work with thousands of organisms we first need to update group_concat_max_len
DBI::dbExecute(con, 'SET SESSION collation_connection = latin1_swedish_ci')


taxon_names<-DBI::dbGetQuery(con, "select * from taxon_names")
## 2256409 obs of 18 variables
save( taxon_names, file="D:\\marlowe\\taxon_names.RData")


all_nodes<-DBI::dbGetQuery(con, "select * from taxons")
## 3077899 obs of 4 variables
save( all_nodes, file="D:\\marlowe\\all_nodes.RData")

kegg_org_taxons<-DBI::dbGetQuery(con, "select kegg_id,kegg_org_code,taxon_id from organisms")
## 5851 obs of 3 variables
save( kegg_org_taxons, file="D:\\marlowe\\kegg_org_taxons.RData")

#now we update the organism table with a foreign key to the taxon_id table
#DBI::dbExecute(con, 'update candidate.organisms o
#    	join candidate.upload_kegg_org_taxon ukt on ukt.org_code = o.kegg_org_code
#    	set o.taxon_id = ukt.taxon_id ')

#verify that the data frame kegg_org_taxon has the correct columns
#assertthat::assert_that("Org_code" %in% colnames(kegg_org_taxon),
#                        msg = "Org_code column missing from kegg_org_taxon")
#assertthat::assert_that("taxon_id" %in% colnames(kegg_org_taxon),
#                        msg = "taxon_id column missing from kegg_org_taxon")


#clean up connection
DBI::dbDisconnect(con)
