
#'Populate the database with a given sample
#'
#'`upload_sample` populates the database with the provided sample.
#'
#'
#'@param conn database connection values to create a MySQL connection
#'@param sample_id The ID of the sample to be processed
#'@param sample A list of data.frames
#'
#'@return true when successful
#'@importFrom assertthat assert_that
#'@importFrom DBI sqlInterpolate dbGetQuery dbWriteTable dbExecute dbDisconnect dbConnect
#'@importFrom RMariaDB MariaDB
#'
#' @author Dustin Crockett
#'
#'@export
upload_sample <- function(conn, sample_id, sample) {
  start_time = Sys.time()
  #validate that the sample_id is valid
  assertthat::assert_that(is.character(sample_id) | is.numeric(sample_id),
                          msg = "The sample_id must be either a number or Character ")

  print(paste0("The sample_id is: ", sample_id))
  #validate the format of the sample object

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

  con <- DBI::dbConnect(RMariaDB::MariaDB(),
                   dbname = conn$dbname,
                   host = conn$host,
                   port = conn$port,
                   user = conn$user,
                   password = conn$password)

  #check to see if the sample id already exists
  sqlcmd <- DBI::sqlInterpolate(con,"select count(*) from candidate.sample_sets where name = ?name", name = sample_id)
  count = DBI::dbGetQuery(con, sqlcmd)
  if(count > 0) {
    #delete the sample id
    delete_sample(conn,sample_id,T)
    #warning(paste0("The sample ", sample_id," was already in the database.  Deleting the old one and replacing with a new one"))
    print(paste0("The sample ", sample_id," was already in the database.  Deleting the old one and replacing with a new one"))
  }
  stop_time <- Sys.time()
  diff_time <- formatC(as.numeric(stop_time, units = "secs") - as.numeric(start_time, units = "secs"), digits=1, format= "f")
  print(paste0("Starting upload of sample ", sample_id ," Started at ", start_time, ", and finished at ", stop_time, " took ", diff_time, " seconds with ", nrow(sample), " tag samples."))
  start_time <- stop_time

  #upload the sample object into the sample_upload table
  DBI::dbWriteTable(
    conn = con,
    DBI::Id(schema="candidate",table="upload_samples"),
    value = sample,
    field.types = c(original_peptide="varchar(4000)",pep_mass="numeric(13,4)",min_mass="numeric(13,4)",max_mass="numeric(13,4)",
                    min_mass_ppm="int2", max_mass_ppm="int2",
                    Weight="numeric(9,5)", tag="varchar(250)", scan_num="numeric(13,4)", LC_threshold="numeric(13,4)", sample_tag_id="varchar(500)"),
    overwrite = T
  )

  print(paste0("creating new sample set ", sample_id ," with ", nrow(sample), " tag samples."))

  #execute the SQL function to process the sample into the final tables
  sqlcmd <- DBI::sqlInterpolate(con,"insert into candidate.sample_sets(name)
	values(?name);", name = sample_id)
  DBI::dbExecute(con, sqlcmd)

  print(paste0("creating samples for ", sample_id ," with ", nrow(sample), " tag samples."))

  sqlcmd <- DBI::sqlInterpolate(con,"insert into candidate.samples(sample_set_id, peptide, peptide_tag, characters_left, mass, weight, min_mass, max_mass, min_mass_ppm, max_mass_ppm, original_peptide , sample_tag_id, scan_num, lc_threshold)
    select ss.id, us.tag, substring(us.tag, 1,5) peptide_tag, length(us.tag) - 5 characters_left, us.Pep_Mass, us.weight, us.min_mass, us.max_mass, us.min_mass_ppm, us.max_mass_ppm, us.original_peptide, us.sample_tag_id, us.scan_num, us.lc_threshold
      from candidate.upload_samples us
        join candidate.sample_sets ss on ss.name = ?name", name = sample_id)
  DBI::dbExecute(con, sqlcmd)

  stop_time <- Sys.time()
  diff_time <- formatC(as.numeric(stop_time, units = "secs") - as.numeric(start_time, units = "secs"), digits=1, format= "f")
  print(paste0("created samples for ", sample_id ," Started at ", start_time, ", and finished at ", stop_time, " took ", diff_time, " seconds"))
  start_time <- stop_time

  #now populate the sample to peptide mappings
  sqlcmd <- DBI::sqlInterpolate(con,"insert into candidate.samples_to_peptides ( sample_id, peptide_id, first_location, count_in_peptide, organism_count)
                          select s.id sample_id, pm.peptide_id, pm.first_location, pm.count_in_peptide, p.organism_count
                        	from candidate.sample_sets ss
                           join candidate.samples s on ss.id = s.sample_set_id and ss.name = ?name
                           join candidate.peptide_tags pt on pt.peptide_tag = s.peptide_tag
                           join candidate.peptide_map pm use index(peptide_map_peptide_tag_id_idx)on pt.id = pm.peptide_tag_id and pm.mass > s.min_mass and pm.mass < s.max_mass
                           join candidate.peptides p on pm.peptide_id = p.id and p.peptide_substitution like concat('%', s.peptide, '%')", name = sample_id)
  DBI::dbExecute(con, sqlcmd)

  stop_time <- Sys.time()
  diff_time <- formatC(as.numeric(stop_time, units = "secs") - as.numeric(start_time, units = "secs"), digits=1, format= "f")
  print(paste0("populated sample_to_peptides for ", sample_id ," Started at ", start_time, ", and finished at ", stop_time, " took ", diff_time, " seconds"))
  start_time <- stop_time

  #now update the samples_peptide_to_organisms table to include the averages

  sqlcmd <- DBI::sqlInterpolate(con,"insert into candidate.samples_peptide_to_organisms(sample_id, organism_id, peptide_id, is_strong_peptide, total_count)
                            select s.id, o.id organism_id, stp.peptide_id,case when sp.strong_percent is null then 0 else 1 end is_strong, otp.total_count
                              from candidate.samples s
                                join candidate.sample_sets ss on ss.id = s.sample_set_id and ss.name = ?name
                              	join candidate.samples_to_peptides stp use index (Primary) on stp.sample_id = s.id
                              	join candidate.organisms_to_peptides otp on stp.peptide_id = otp.peptide_id
                								join candidate.organisms o on o.id = otp.organism_id
                								join candidate.species spe on spe.id = o.species_id
                                left join candidate.strong_peptides sp on sp.peptide_id = stp.peptide_id and sp.strong_genus_id = spe.genus_id", name = sample_id)
  DBI::dbExecute(con, sqlcmd)

  stop_time <- Sys.time()
  diff_time <- formatC(as.numeric(stop_time, units = "secs") - as.numeric(start_time, units = "secs"), digits=1, format= "f")
  print(paste0("creating sample_peptide_to_organisms for ", sample_id ," Started at ", start_time, ", and finished at ", stop_time, " took ", diff_time, " seconds"))
  start_time <- stop_time

  #clean up connection
  DBI::dbDisconnect(con)

  print(paste0("Uploaded sample ", sample_id ," with ", nrow(sample), " tag samples."))

  return(T)
}

#'Removes all information in the database for the provided sample ID
#'
#'`delete_sample` Removes all information in the database for the provided sample ID.
#'
#'@param conn database connection values to create a MySQL connection
#'@param sample_id The ID for the sample
#'@param force When true don't prompt the user to verify if they want to delete the sample
#'
#'@return the total count of rows deleted
#'@importFrom assertthat assert_that
#'@importFrom DBI sqlInterpolate dbExecute dbDisconnect dbConnect
#'@importFrom RMariaDB MariaDB
#'
#' @author Dustin Crockett
#'
#'@export
delete_sample <- function(conn, sample_id, force=F) {
  print(c("inside delete_sample, sample_id= ", sample_id))

  #validate that the sample_id is valid
  assertthat::assert_that(is.character(sample_id) | is.numeric(sample_id),
                          msg = "The sample_id must be either a number or Character ")

  #validate the the user really wants to delete the sample
  if(!is.logical(force) | !isTRUE(force)) {
    response <- readline(prompt=paste0("Are you sure you want to delete sample ",sample_id, "?"))
    if(!toupper(response) == "Y" & !toupper(response) == "YES" ) {
      return(0)
    }
  }

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

  con <- DBI::dbConnect(RMariaDB::MariaDB(),
                   dbname = conn$dbname,
                   host = conn$host,
                   port = conn$port,
                   user = conn$user,
                   password = conn$password)

  #delete the sample_to_organisms  - this was commented out by Dustin
  # sql <- sqlInterpolate(con, "delete from candidate.samples_to_organisms so
  #                       where so.sample_id in (select s.id from candidate.sample_sets ss
  #                       join candidate.samples s on s.sample_set_id = ss.id
  #                       where ss.name = ?sample)",
  #                       sample = sample_id)
  #rowcount <- dbExecute(con, sql)



  #delete the sample_peptide_to_organisms
  print( "get ready for the error we seem to get with delete samples now...")
  # sample_id <- "Rcom_1_M0_AM_R1_7Mar16_Samwise_15-08-55"

  # SELECT FROM candidate.samples_peptide_to_organisms AS sp
  # WHERE sp.sample_id in
  # (SELECT s.id
  # FROM candidate.sample_sets AS ss
  # JOIN candidate.samples AS s on s.sample_set_id = ss.id
  # WHERE ss.name =  "Rcom_1_M0_AM_R1_7Mar16_Samwise_15-08-55");


  # sql <- DBI::sqlInterpolate(con,  "DELETE FROM candidate.samples_peptide_to_organisms sp
  #                                   WHERE sp.sample_id in ( SELECT s.id
  #                                                           FROM candidate.sample_sets ss
  #                                                           JOIN candidate.samples s on s.sample_set_id = ss.id
  #                                                           WHERE ss.name = ?sample)",
  #                                                           sample = sample_id)



  sql <- DBI::sqlInterpolate(con, "delete from candidate.samples_peptide_to_organisms sp
                        where sp.sample_id in (select s.id from candidate.sample_sets ss
                        join candidate.samples s on s.sample_set_id = ss.id
                        where ss.name = ?sample)",
                        sample = sample_id)

  print( "delete from candidate.samples_peptide_to_organisms")
  print(sql)
  rowcount <- DBI::dbExecute(con, sql)

  print( "delete from candidate.samples_peptides")
  #delete the sample_to_peptides
  sql <- DBI::sqlInterpolate(con, "delete from candidate.samples_to_peptides sp
                        where sp.sample_id in (select s.id from candidate.sample_sets ss
                        join candidate.samples s on s.sample_set_id = ss.id
                        where ss.name = ?sample)",
                        sample = sample_id)
  rowcount <- rowcount + DBI::dbExecute(con, sql)
  #delete all samples for the sample set
  sql <- DBI::sqlInterpolate(con, "delete from candidate.samples
      where sample_set_id in (select id from candidate.sample_sets ss where ss.name = ?sample)",
                        sample = sample_id)
  rowcount <- rowcount + DBI::dbExecute(con, sql)
  #delete the sample set
  sql <- DBI::sqlInterpolate(con, "delete from candidate.sample_sets ss
                        where ss.name = ?sample",
                        sample = sample_id)
  rowcount <- rowcount + DBI::dbExecute(con, sql)

  #clean up connection
  DBI::dbDisconnect(con)

  return(rowcount)
}

#'Returns Candidate ratings for the sample ID provided
#'
#'`candidate_ratings` Returns Candidate ratings for the sample ID provided.
#'
#'@param conn database connection values to create a MySQL connection
#'@param sample_id The ID for the sample
#'
#'@return a data.frame of ratings for all organisms in the KEGG
#'@importFrom assertthat assert_that
#'@importFrom DBI sqlInterpolate dbGetQuery dbDisconnect dbConnect
#'@importFrom RMariaDB MariaDB
#'
#' @author Dustin Crockett
#'
#'@export
candidate_ratings <- function(conn, sample_id) {
  #validate that the sample_id is valid
  assertthat::assert_that(is.character(sample_id) | is.numeric(sample_id),
                          msg = "The sample_id must be either a number or Character ")

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

  con <- DBI::dbConnect(RMariaDB::MariaDB(),
                   dbname = conn$dbname,
                   host = conn$host,
                   port = conn$port,
                   user = conn$user,
                   password = conn$password)

  #retrieve the candidate ratings for the given sample_id
  sqlcmd <- DBI::sqlInterpolate(con,"select  o.kegg_id, g.name genus, s.name species, o.name organism, o.peptide_count total_peptides_unique, pc.pep_unique_count peptide_unique_hits, pc.pep_total_count peptide_total_hits, sample_unique_count,
	mean_org_count, median.median median_org_count,weighted_hit_unique_count,weighted_hit_total_count, strong_count
from
candidate.organisms o
join candidate.species s on o.species_id = s.id
join candidate.genus g on s.genus_id = g.id
left outer join
(select coalesce(sum(weight),0) weighted_hit_unique_count,coalesce(sum(total_weight),0) weighted_hit_total_count, coalesce( sum(is_strong_peptide ),0) strong_count, kegg_id from
(select s.id sample_id, max(s.weight) weight, sum(s.weight) total_weight, o.kegg_id, ss.id set_id, max(spto.is_strong_peptide) is_strong_peptide
	from candidate.sample_sets ss
	join candidate.samples s on ss.id = s.sample_set_id and ss.name = ?code
	join candidate.samples_to_peptides sp on s.id = sp.sample_id
	join candidate.organisms_to_peptides otp on sp.peptide_id = otp.peptide_id
	join candidate.organisms o on otp.organism_id = o.id
     join candidate.samples_peptide_to_organisms spto on spto.organism_id = o.id and spto.peptide_id = sp.peptide_id and spto.sample_id = s.id
	group by s.id, o.kegg_id) as innerQ
group by kegg_id ) as k on o.kegg_id = k.kegg_id
left outer join
(select  o.kegg_id, count(distinct sp.peptide_id ) pep_unique_count, count(sp.peptide_id ) pep_total_count, count(distinct sp.sample_id ) sample_unique_count, avg(sp.organism_count) mean_org_count
	from candidate.sample_sets ss
	join candidate.samples s on ss.id = s.sample_set_id and ss.name = ?code
	join candidate.samples_to_peptides sp on s.id = sp.sample_id
	join candidate.organisms_to_peptides otp on sp.peptide_id = otp.peptide_id
	join candidate.organisms o on otp.organism_id = o.id
    join candidate.samples_peptide_to_organisms spto on spto.organism_id = o.id and spto.peptide_id = sp.peptide_id and spto.sample_id = s.id
	group by o.kegg_id) as pc on o.kegg_id = pc.kegg_id
left outer join
(select kegg_id, avg(1.0 * organism_count ) as median
from(
select d.kegg_id, d.organism_count, count(*) over (partition by kegg_id) as cnt, row_number() over(partition by kegg_id order by organism_count asc) as rn
from
(select  o.kegg_id, sp.organism_count
	from candidate.sample_sets ss
	join candidate.samples s on ss.id = s.sample_set_id and ss.name = ?code
	join candidate.samples_to_peptides sp on s.id = sp.sample_id
	join candidate.organisms_to_peptides otp on sp.peptide_id = otp.peptide_id
	join candidate.organisms o on otp.organism_id = o.id
    join candidate.samples_peptide_to_organisms spto on spto.organism_id = o.id and spto.peptide_id = sp.peptide_id and spto.sample_id = s.id) d
) as x where rn in ((cnt + 1)/2, (cnt+2)/2) group by kegg_id) as median  on o.kegg_id = median.kegg_id
order by o.kegg_id asc", code = sample_id)
  results <- DBI::dbGetQuery(con, sqlcmd)

  #clean up connection
  DBI::dbDisconnect(con)

  return(results)
}

#'Returns all details about the sample peptide matches to candidates
#'
#'`peptide_matches` Returns all details about the sample peptide matches to candidates.
#'
#'@param conn database connection values to create a MySQL connection
#'@param sample_id The ID for the sample
#'
#'@return a data.frame of all details
#'@importFrom assertthat assert_that
#'@importFrom DBI sqlInterpolate dbGetQuery dbDisconnect dbConnect
#'@importFrom RMariaDB MariaDB
#'
#' @author Dustin Crockett
#'
#'@export
peptide_matches <- function(conn, sample_id) {
  #validate that the sample_id is valid
  assertthat::assert_that(is.character(sample_id) | is.numeric(sample_id),
                          msg = "The sample_id must be either a number or Character ")

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

  con <- DBI::dbConnect(RMariaDB::MariaDB(),
                   dbname = conn$dbname,
                   host = conn$host,
                   port = conn$port,
                   user = conn$user,
                   password = conn$password)

  #retrieve all details about the sample peptide matches to candidates
  sqlcmd <- DBI::sqlInterpolate(con,"select s.peptide sample_peptide , sp.first_location, sp.count_in_peptide, pr.id protein_id, pr.protein_kegg_id, unipro.key unipro, NCBI_prot.key as `NCBI-prot`, NCBI_Gene.key as `NCBI-Gene`,
                            o.kegg_id, o.kegg_org_code, o.name organism, spe.name species, g.name genus, p.mass recorded_mass, s.mass measured_mass, p.peptide_original, sp.organism_count
                           , so.total_count organism_total_count, stp.total_count protein_total_count, sample_hits.total_count sample_total_count, so.is_strong_peptide
                           , s.original_peptide novor_peptide, s.sample_tag_id, s.scan_num, tn.taxon_id, tn.name_txt taxon_name
                           from candidate.sample_sets ss
                           join candidate.samples s on ss.id = s.sample_set_id and ss.name = ?code
                           join candidate.samples_to_peptides sp on s.id = sp.sample_id
                           join candidate.peptides p on sp.peptide_id = p.id
                           join candidate.sequence_to_peptides stp on p.id = stp.peptide_id
                           join candidate.proteins pr on stp.sequence_id = pr.sequence_id
                           join candidate.organisms o on pr.organism_id = o.id
                           join candidate.species spe on o.species_id = spe.id
                           join candidate.genus g on spe.genus_id = g.id
                           join candidate.samples_peptide_to_organisms so on so.sample_id = s.id and so.organism_id = o.id and so.peptide_id = p.id
                           join (select s.peptide, count(s.id) total_count from candidate.samples s
            									join candidate.sample_sets ss on ss.id = s.sample_set_id and ss.name = ?code
            								group by s.peptide)	sample_hits on sample_hits.peptide = s.peptide
                           left outer join candidate.protein_to_dblinks unipro on unipro.protein_id = pr.id and unipro.dblink_id in (select id from candidate.dblinks d where `database` = 'UniProt' )
                           left outer join candidate.protein_to_dblinks NCBI_prot on NCBI_prot.protein_id = pr.id and NCBI_prot.dblink_id in (select id from candidate.dblinks d where `database` = 'NCBI-ProteinID' )
                           left outer join candidate.protein_to_dblinks NCBI_Gene on NCBI_Gene.protein_id = pr.id and NCBI_Gene.dblink_id in (select id from candidate.dblinks d where `database` = 'NCBI-GeneID')
                           left outer join candidate.taxon_names tn on tn.taxon_id = o.taxon_id", code = sample_id)
  results <- DBI::dbGetQuery(con, sqlcmd)

  #clean up connection
  DBI::dbDisconnect(con)

  return(results)
}


#'Returns all details about the sample peptide matches to candidates
#'
#'`sample_overlap` Returns all samples and organism pairs with hit counts and is_strong.
#'
#'@param conn database connection values to create a MySQL connection
#'@param sample_id The ID for the sample
#'
#'@return a data.frame of sample and organism pairs with hit counts and is_strong, this can be pivoted to make a matrix
#'@importFrom assertthat assert_that
#'@importFrom DBI sqlInterpolate dbGetQuery dbDisconnect dbConnect
#'@importFrom RMariaDB MariaDB
#'
#' @author Dustin Crockett
#'
#'@export
sample_overlap <- function(conn, sample_id) {
  #validate that the sample_id is valid
  assertthat::assert_that(is.character(sample_id) | is.numeric(sample_id),
                          msg = "The sample_id must be either a number or Character ")

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

  con <- DBI::dbConnect(RMariaDB::MariaDB(),
                   dbname = conn$dbname,
                   host = conn$host,
                   port = conn$port,
                   user = conn$user,
                   password = conn$password)

  #retrieve all details about the sample peptide matches to candidates
  sqlcmd <- DBI::sqlInterpolate(con,"select main.kegg_id, main.peptide, coalesce(samp_data.sample_found, 0) sample_found, coalesce(samp_data.is_strong_peptide, 0) is_strong_peptide, coalesce(samp_data.hit_count, 0) hit_count
		from (select o.id organism_id, o.kegg_id, samples.peptide, samples.id sample_id from candidate.organisms o,
                           (select s.id, s.peptide from candidate.sample_sets ss
                           join candidate.samples s on s.sample_set_id = ss.id
                           where ss.name = ?code) samples) main
                           left outer join (select s.id, spo.organism_id, max(spo.is_strong_peptide) is_strong_peptide, count(spo.peptide_id) hit_count, 1 sample_found
                           from candidate.sample_sets ss
                           join candidate.samples s on s.sample_set_id = ss.id
                           join candidate.samples_peptide_to_organisms spo on spo.sample_id = s.id
                           where ss.name = ?code
                           group by s.id, spo.organism_id) samp_data on  main.organism_id = samp_data.organism_id and main.sample_id = samp_data.id", code = sample_id)
  results <- DBI::dbGetQuery(con, sqlcmd)

  #clean up connection
  DBI::dbDisconnect(con)

  return(results)
}



#'Returns the count of peptide hits for each sample and organism pair even when there are no hits
#'
#'`sample_count_by_organism` Returns the count of peptide hits for each sample and organism pair even when there are no hits.
#'
#'@param conn database connection values to create a MySQL connection
#'@param sample_id The ID for the sample
#'
#'@return a data.frame of sample and organism pairs with hit counts in matrix form
#'@importFrom assertthat assert_that
#'@importFrom DBI sqlInterpolate dbGetQuery dbDisconnect dbConnect
#'@importFrom RMariaDB MariaDB
#'@importFrom utils read.csv
#' @author Dustin Crockett
#'
#'@export
sample_count_by_organism <- function(conn, sample_id) {
  #validate that the sample_id is valid
  assertthat::assert_that(is.character(sample_id) | is.numeric(sample_id),
                          msg = "The sample_id must be either a number or Character ")

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

  #retrieve all hit counts
  sqlcmd <- DBI::sqlInterpolate(con,"select 'peptide_tag' peptide_tag, GROUP_CONCAT( o.kegg_id
                                      ORDER BY o.id ASC
                                      SEPARATOR ',') matrix
                                      from organisms o
                                union
                                select concat(peptide,'_',mass),
                                	   GROUP_CONCAT( hit_count
                                        ORDER BY organism_id ASC
                                        SEPARATOR ',')
                                from
                                (select main.organism_id, main.peptide, main.mass, coalesce(samp_data.hit_count, 0) hit_count, sample_id
                                		from (select o.id organism_id, o.kegg_id, samples.peptide, samples.mass, samples.id sample_id from candidate.organisms o,
                                                           (select s.id, s.peptide, s.mass from candidate.sample_sets ss
                                                           join candidate.samples s on s.sample_set_id = ss.id
                                                           where ss.name = ?code) samples) main
                                                           left outer join (select s.id, spo.organism_id, count(spo.peptide_id) hit_count
                                                           from candidate.sample_sets ss
                                                           join candidate.samples s on s.sample_set_id = ss.id
                                                           join candidate.samples_peptide_to_organisms spo on spo.sample_id = s.id
                                                           where ss.name = ?code
                                                           group by s.id, spo.organism_id) samp_data on  main.organism_id = samp_data.organism_id and main.sample_id = samp_data.id) main
                                group by sample_id", code = sample_id)
  results <- DBI::dbGetQuery(con, sqlcmd)

  #clean up connection
  DBI::dbDisconnect(con)

  #let R read database results like a csv file
  #results$peptide_tag is a vector of strings with the first string being the column header, and becomes column 1
  #results$matrix will supply the rest of the columns
  output_df <- cbind(
    read.csv(text = results$peptide_tag, stringsAsFactors = F),
    read.csv(text = results$matrix, stringsAsFactors = F)
  )

  return(output_df)

}

#'Returns the count of strong peptide hits for each sample and organism pair even when there are no hits
#'
#'`sample_strong_count_by_organism` Returns the count of strong peptide hits for each sample and organism pair even when there are no hits.
#'
#'@param conn database connection values to create a MySQL connection
#'@param sample_id The ID for the sample
#'
#'@return a data.frame of sample and organism pairs with hit counts in matrix form
#'@importFrom assertthat assert_that
#'@importFrom DBI sqlInterpolate dbGetQuery dbDisconnect dbConnect
#'@importFrom RMariaDB MariaDB
#'@importFrom utils read.csv
#' @author Dustin Crockett
#'
#'@export
sample_strong_count_by_organism <- function(conn, sample_id) {
  #validate that the sample_id is valid
  assertthat::assert_that(is.character(sample_id) | is.numeric(sample_id),
                          msg = "The sample_id must be either a number or Character ")

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

  #retrieve all hit counts
  sqlcmd <- DBI::sqlInterpolate(con,"select 'peptide_tag' peptide_tag, GROUP_CONCAT( o.kegg_id
                                      ORDER BY o.id ASC
                                      SEPARATOR ',') matrix
                                      from organisms o
                                union
                                select concat(peptide,'_',mass),
                                	   GROUP_CONCAT( hit_count
                                        ORDER BY organism_id ASC
                                        SEPARATOR ',')
                                from
                                (select main.organism_id, main.peptide, main.mass, coalesce(samp_data.hit_count, 0) hit_count, sample_id
                                		from (select o.id organism_id, o.kegg_id, samples.peptide, samples.mass, samples.id sample_id from candidate.organisms o,
                                                           (select s.id, s.peptide, s.mass from candidate.sample_sets ss
                                                           join candidate.samples s on s.sample_set_id = ss.id
                                                           where ss.name = ?code) samples) main
                                                           left outer join (select s.id, spo.organism_id, count(spo.peptide_id) hit_count
                                                           from candidate.sample_sets ss
                                                           join candidate.samples s on s.sample_set_id = ss.id
                                                           join candidate.samples_peptide_to_organisms spo on spo.sample_id = s.id and spo.is_strong_peptide = 1
                                                           where ss.name = ?code
                                                           group by s.id, spo.organism_id) samp_data on  main.organism_id = samp_data.organism_id and main.sample_id = samp_data.id) main
                                group by sample_id", code = sample_id)
  results <- DBI::dbGetQuery(con, sqlcmd)

  #clean up connection
  DBI::dbDisconnect(con)

  #let R read database results like a csv file
  #results$peptide_tag is a vector of strings with the first string being the column header, and becomes column 1
  #results$matrix will supply the rest of the columns
  output_df <- cbind(
    read.csv(text = results$peptide_tag, stringsAsFactors = F),
    read.csv(text = results$matrix, stringsAsFactors = F)
  )

  return(output_df)

}

#'Returns the sample_Id and peptide_tag for each sample in a given sample set
#'
#'`sample_info` Returns the sample_Id and peptide_tag for each sample in a given sample set.
#'
#'@param conn database connection values to create a MySQL connection
#'@param sample_id The ID for the sample
#'
#'@return a data.frame of sample_id and peptide_tag, this is to be used to pivot other dataframes into a matrix
#'@importFrom assertthat assert_that
#'@importFrom DBI sqlInterpolate dbGetQuery dbDisconnect dbConnect
#'@importFrom RMariaDB MariaDB
#'
#' @author Dustin Crockett
#'
#'@export
sample_info <- function(conn, sample_id) {
  #validate that the sample_id is valid
  assertthat::assert_that(is.character(sample_id) | is.numeric(sample_id),
                          msg = "The sample_id must be either a number or Character ")

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

  con <- DBI::dbConnect(RMariaDB::MariaDB(),
                   dbname = conn$dbname,
                   host = conn$host,
                   port = conn$port,
                   user = conn$user,
                   password = conn$password)

  #retrieve all details about the sample peptide matches to candidates
  sqlcmd <- DBI::sqlInterpolate(con,"select s.id sample_id, s.peptide peptide_tag
                                	from candidate.sample_sets ss
                                			join samples s on ss.id = s.sample_set_id
                                			where ss.name = ?code", code = sample_id)
  results <- DBI::dbGetQuery(con, sqlcmd)

  #clean up connection
  DBI::dbDisconnect(con)

  return(results)
}



#' Get sample ids in database
#'
#' @param conn database connection values to create a MySQL connection
#'
#' @return Returns a single column dataframe with the sample_ids of the samples
#'   currently loaded in the database.
#' @export
#'
#'@importFrom assertthat assert_that
#'@importFrom DBI sqlInterpolate dbGetQuery dbDisconnect dbConnect
#'@importFrom RMariaDB MariaDB
get_sample_ids <- function(conn) {

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

  #making the db connection object
  con <- DBI::dbConnect(RMariaDB::MariaDB(),
                        dbname = conn$dbname,
                        host = conn$host,
                        port = conn$port,
                        user = conn$user,
                        password = conn$password)

  #retrieve all sample ids currently in database
  results <- DBI::dbGetQuery(con, "select ss.name from candidate.sample_sets ss")

  #clean up connection
  DBI::dbDisconnect(con)

  return(results)

}
