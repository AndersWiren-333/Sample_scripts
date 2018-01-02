alter table matches add index sample_id_index (sample_id asc);
alter table matches add index catalog_id_index (catalog_id asc);

create table uw01_noSNP_loci (select distinct cat_id from catalog_index where snps=0);
alter table uw01_noSNP_loci add index cat_id_index (cat_id asc);

create table uw02_samples (select t1.catalog_id as cat_id, t1.sample_id, t1.allele, t2.group_id as sex, t2.file as sample_name, t2.pop_id 
from matches as t1 inner join samples as t2 on t1.sample_id = t2.sample_id);
alter table uw02_samples add index cat_id_index (cat_id asc);

create table uw03_noSNP_samples (
select t1.* from uw02_samples as t1
inner join uw01_noSNP_loci as t2 on t1.cat_id = t2.cat_id);
alter table uw03_noSNP_samples add index cat_id_index (cat_id asc);

create table temp1 (select cat_id from uw03_noSNP_samples where sex = 'f');
create table temp2 (select distinct cat_id from temp1);
drop table temp1;
alter table temp2 add index cat_id_index (cat_id asc);

delete from uw03_noSNP_samples where cat_id in (select cat_id from temp2);
drop table temp2;
alter table uw03_noSNP_samples rename to uw03_noSNP_samples_M;

create table uw_Potential_Y_conserved (select cat_id from uw03_noSNP_samples_M group by cat_id
having count(distinct sample_name) = 4); #between 4 and 8);
alter table uw_Potential_Y_conserved add index cat_id_index (cat_id asc);