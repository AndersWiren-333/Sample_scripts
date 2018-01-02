# EDIT ROWS 101 AND 106! FILL IN THE CORRECT SAMPLE SIZE.

# Prepare the database by adding indices to some columns in some tables (speeds up the script very much)
#alter table catalog_index add index cat_id_index (cat_id asc);
#alter table samples change column group_id sex varchar(32);
#alter table catalog_snps add index col_index (col asc);
#alter table catalog_snps add index rank_1_index (rank_1 asc);
#alter table catalog_snps add index rank_2_index (rank_2 asc);
#alter table catalog_snps add index rank_3_index (rank_3 asc);
#alter table catalog_snps add index rank_4_index (rank_4 asc);

# List all loci (catalog ids, called "cat_id") that have 1 SNP and 2 alleles (required for assay design to work)
create table yx01_loci (
select cat_id from catalog_index
where snps = 1
and alleles = 2 and batch_id = 1);
alter table yx01_loci add index cat_id_index (cat_id asc);

# List all loci (cat_ids) that have their SNP between positions 35 and 105 in the read (required for assay design to work)
create table yx01_column (
select distinct tag_id as cat_id from sumstats
where col between 35 and 105 and batch_id = 1);
alter table yx01_column add index cat_id_index (cat_id asc);

create table yx01_proper_loci (select t1.cat_id from yx01_loci as t1
inner join yx01_column as t2 on t1.cat_id=t2.cat_id);
alter table yx01_proper_loci add index cat_id_index (cat_id asc);

drop table yx01_loci, yx01_column;

# Gather onformation about which samples have which alleles at which loci, and which sex those samples have
create table yx02_samples (select t1.catalog_id as cat_id, t1.sample_id, t1.allele, t2.sex from matches as t1
inner join samples as t2 on t1.sample_id = t2.sample_id where t1.batch_id = 1);
alter table yx02_samples add index cat_id_index (cat_id asc);
alter table yx02_samples add index sample_id_index (sample_id asc);
alter table yx02_samples add index allele_index (allele asc);
alter table yx02_samples add index sex_index (sex asc);

# List all loci that that occur in males, and all that occur in females (these may overlap)
create table yx03_loci_M (select distinct cat_id from yx02_samples where sex='m'); alter table yx03_loci_M add index cat_id_index (cat_id asc);
create table yx03_loci_F (select distinct cat_id from yx02_samples where sex='f'); alter table yx03_loci_F add index cat_id_index (cat_id asc);

# Make a list of the loci in the male list that have the required no. of SNPs, alleles and a suitable position in the read
create table yx04_proper_loci_M (select t1.cat_id from yx03_loci_M as t1 inner join yx01_proper_loci as t2 on t1.cat_id=t2.cat_id);
alter table yx04_proper_loci_M add index cat_id_index (cat_id asc);

# Make a list of the loci in the female list that have the required no. of SNPs, alleles and a suitable position in the read
create table yx04_proper_loci_F (select t1.cat_id from yx03_loci_F as t1 inner join yx01_proper_loci as t2 on t1.cat_id=t2.cat_id);
alter table yx04_proper_loci_F add index cat_id_index (cat_id asc);

drop table yx03_loci_F, yx03_loci_M;

# Make a list of suitable loci for which both males and females have genotypes (= potential X-loci)
create table yx05_loci_inclusive_F (select f.cat_id as cat_id_F, m.cat_id as cat_id_M from yx04_proper_loci_F as f
left join yx04_proper_loci_M as m on f.cat_id=m.cat_id);
alter table yx05_loci_inclusive_F add index cat_id_F_index (cat_id_F asc);

create table yx06_loci_MF (select cat_id_F as cat_id from yx05_loci_inclusive_F where cat_id_M is not null order by cat_id_F);
alter table yx06_loci_MF add index cat_id_index (cat_id asc);

drop table yx05_loci_inclusive_F;

# Sort out loci from yx04_proper_loci_M where females don't have a genotype (= presumed Y-loci)
create table yx06_loci_M_1 (select t1.cat_id as cat_id_M, t2.cat_id as cat_id_F from yx04_proper_loci_M as t1
left join yx04_proper_loci_F as t2 on t1.cat_id = t2.cat_id);
alter table yx06_loci_M_1 add index cat_id_M_index (cat_id_M asc);

create table yx06_loci_M (select cat_id_M as cat_id from yx06_loci_M_1
where cat_id_F is null);
alter table yx06_loci_M add index cat_id_index (cat_id asc);

drop table yx06_loci_M_1, yx01_proper_loci, yx04_proper_loci_F, yx04_proper_loci_M;

# Integrate genotype information for each individual into the tables with male and female loci
create table yx07_samples_MF (select t1.* from yx02_samples as t1 inner join yx06_loci_MF as t2 on t1.cat_id = t2.cat_id);
alter table yx07_samples_MF add index cat_id_index (cat_id asc);
alter table yx07_samples_MF add index sample_id_index (sample_id asc);
alter table yx07_samples_MF add index sex_index (sex asc);
alter table yx07_samples_MF add index allele_index (allele asc);

create table yx07_samples_M (select t1.* from yx02_samples as t1 inner join yx06_loci_M as t2 on t1.cat_id = t2.cat_id);
alter table yx07_samples_M add index cat_id_index (cat_id asc);
alter table yx07_samples_M add index sample_id_index (sample_id asc);
alter table yx07_samples_M add index sex_index (sex asc);
alter table yx07_samples_M add index allele_index (allele asc);

# Add zygosity column to the samplesM/F tables
alter table yx07_samples_M add column zygosity varchar(40) after sample_id;
alter table yx07_samples_MF add column zygosity varchar(40) after sample_id;

# Add information about zygosity to each row of the samples_M/F tables, using the setZygosity stored procedures
# (they work by counting the number of rows for the "allele" column of these tables, on a per locus/sample basis. If a certain sample
# has two allele-rows in the table for a given locus, it means that the individual is heterozygous for that locus. 
# If it has only one row it is homozygous. 'HET' or 'HOM' is filled into the 'zygosity' column depending on which.
call setZygosity_M();
call setZygosity_MF();

# Remove those loci from samples_MF and samples_M where males are heterozygous (males can only be hemizygous for x-linked loci)
delete from yx07_samples_MF where cat_id in (select distinct cat_id from (select cat_id from yx07_samples_MF as t2 where sex='m' and zygosity='HET') as t1); 
delete from yx07_samples_M where cat_id in (select distinct cat_id from (select cat_id from yx07_samples_M as t2 where sex='m' and zygosity='HET') as t1);
delete from yx07_samples_MF where zygosity is null;

# Choose those X-loci that have data for all individuals
create table yx08_presumed_X (select cat_id from yx07_samples_MF group by cat_id
having count(distinct file) = 16);
alter table yx08_presumed_X add index cat_id_index (cat_id asc);

# Choose those Y-loci that have data for all (male) individuals
create table yx08_presumed_Y_hiRep (select cat_id from yx07_samples_M group by cat_id
having count(distinct file) < 6); 
alter table yx08_presumed_Y_hiRep add index cat_id_index (cat_id asc);

drop table yx07_samples_M, yx07_samples_MF, yx06_loci_M, yx06_loci_MF, yx02_samples;

alter table yx08_presumed_X rename to z_Potential_X_SNPs;
alter table yx08_presumed_Y_hiRep rename to z_Potential_Y_SNPs;

# -------------------------------------------------------

create table z01_X_seq (select t1.cat_id, t2.seq from z_Potential_X_SNPs as t1 inner join catalog_tags as t2 on t1.cat_id = t2.tag_id);
alter table z01_X_seq add index cat_id_index (cat_id asc);
create table z01_Y_seq (select t1.cat_id, t2.seq from z_Potential_Y_SNPs as t1 inner join catalog_tags as t2 on t1.cat_id = t2.tag_id);
alter table z01_Y_seq add index cat_id_index (cat_id asc);

create table z02_X_info (select t1.*, t2.col, t2.rank_1 as p_nuc, t2.rank_2 as q_nuc from z01_X_seq as t1 inner join catalog_snps as t2 on t1.cat_id = t2.tag_id);
alter table z02_X_info add index cat_id_index (cat_id asc);

create table z02_Y_info (select t1.*, t2.col, t2.rank_1 as p_nuc, t2.rank_2 as q_nuc from z01_Y_seq as t1 inner join catalog_snps as t2 on t1.cat_id = t2.tag_id);
alter table z02_Y_info add index cat_id_index (cat_id asc);

drop table z_Potential_X_SNPs, z_Potential_Y_SNPs, z01_X_seq, z01_Y_seq;

alter table z02_X_info rename to z_Potential_X_SNPs;
alter table z02_Y_info rename to z_Potential_Y_SNPs;

select distinct cat_id from z_Potential_X_SNPs into outfile '/Users/ided0001/Desktop/whitelistX' 
fields terminated by ';' lines terminated by '\n';

select distinct cat_id from z_Potential_Y_SNPs into outfile '/Users/ided0001/Desktop/whitelistY' 
fields terminated by ';' lines terminated by '\n';



# ---------------------------------------------------------

create table z04_X (select cat_id, concat('[',p_nuc,'/',q_nuc,']') as snp, col, seq  from z_Potential_X_SNPs);
alter table z04_X add index cat_id_index (cat_id asc);

create table z04_Y (select cat_id, concat('[',p_nuc,'/',q_nuc,']') as snp, col, seq  from z_Potential_Y_SNPs);
alter table z04_Y add index cat_id_index (cat_id asc);

drop table z_Potential_X_SNPs, z_Potential_Y_SNPs;

alter table z04_X rename to z_potential_X_SNPs;
alter table z04_Y rename to z_potential_Y_SNPs;

select cat_id from z_Potential_X_SNPs into outfile '/Users/ided0001/Desktop/IdasX_SNP' 
fields terminated by ';' lines terminated by '\n';

select cat_id from z_Potential_Y_SNPs into outfile '/Users/ided0001/Desktop/IdasY_SNP' 
fields terminated by ';' lines terminated by '\n';