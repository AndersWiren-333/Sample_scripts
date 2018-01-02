# Add population names to the populations table
ALTER TABLE populations 
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
alter table samples add column sex varchar(32) after pop_name;

# Manually add sex information for individuals in the samples table


# Add indices to tables
alter table catalog_index add index cat_id_index (cat_id asc);
alter table catalog_index add index snps_index (snps asc);
alter table catalog_index add index alleles_index (alleles asc);

alter table sumstats add index tag_id_index (tag_id asc);
alter table sumstats add index col_index (col asc);

