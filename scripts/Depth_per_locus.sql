create table zz401_Depth_master (select catalog_id, sample_id, depth from matches);
alter table zz401_Depth_master add index catalog_id_index (catalog_id asc);
alter table zz401_Depth_master add index sample_id_index (sample_id asc);
alter table zz401_Depth_master add index depth_index (depth asc);

create table zz402_Depth_per_locus (select catalog_id, sum(depth) as depth from zz401_Depth_master 
group by catalog_id);

select avg(depth) from zz402_Depth_per_locus;
select stddev_samp(depth) from zz402_Depth_per_locus;

call depth_per_sample();

select avg(mean), stddev_samp(mean) from zz403_Depth_per_sample_per_locus;

create table zz405_Depth_per_sample_per_locus_120 (select t1.* from zz403_Depth_per_sample_per_locus as t1
inner join x6_Fluidigm_order as t2 on t1.cat_id = t2.namn);

create table zz404_Depth_master_120 (select t1.* from zz401_Depth_master as t1
inner join x6_Fluidigm_order as t2 on t1.catalog_id = t2.namn);

create table zz406_Depth_per_locus_120 (select catalog_id, sum(depth) as depth from zz404_Depth_master_120 group by catalog_id);

select * from zz406_Depth_per_locus_120;
select avg(depth), stddev_samp(depth) from zz406_Depth_per_locus_120;

select avg(mean), stddev_samp(mean) from zz405_Depth_per_sample_per_locus_120;