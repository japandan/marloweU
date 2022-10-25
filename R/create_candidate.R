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


