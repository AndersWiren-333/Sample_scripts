# This SQL script is designed to work with a MySQL database which is the result of

create table catalog_genotypes_Orig (select * from catalog_genotypes);

# Add population names to the populations table
alter table populations 
ADD COLUMN pop_name VARCHAR(128) NULL DEFAULT NULL AFTER pop_id;

update populations 
set pop_name='Capreolus'
where pop_id=1;

update populations 
set pop_name='Cervus'
where pop_id=2;

update populations 
set pop_name='Dama'
where pop_id=3;

update populations 
set pop_name='Alces'
where pop_id=4;

update populations 
set pop_name='Rangifer'
where pop_id=5;



# Add columns to the samples table
alter table samples add column pop_name char(100);
alter table samples add column sex char(1);
alter table samples add column sample_name char(100);
update samples set sample_name = file;

# Manually add sex information for individuals in the samples table


# Add indices to tables
alter table catalog_index add index cat_id_index (cat_id asc);
alter table sumstats add index col_index (col asc);

alter table matches add index allele_index (allele asc);
alter table matches add index depth_index (depth asc);
