-- --------------------------------------------------------------------------------
-- Routine DDL
-- Note: comments before and after the routine body will not be stored by the server
-- --------------------------------------------------------------------------------
DELIMITER $$

CREATE DEFINER=`aswn0002`@`localhost` PROCEDURE `depth_per_sample`()
BEGIN
block1: begin
	declare locus int;
	declare medel float;
	declare savv float;
	declare finished_loci int default 0;
	declare locus_cursor cursor for (select distinct catalog_id from zz401_Depth_master);
	declare continue handler for not found set finished_loci=1;

	create table zz403_Depth_per_sample_per_locus (cat_id int, mean float, sdev float);
	alter table zz403_Depth_per_sample_per_locus add index cat_id_index (cat_id asc);
	alter table zz403_Depth_per_sample_per_locus add index mean_index (mean asc);
	alter table zz403_Depth_per_sample_per_locus add index sdev_index (sdev asc);

	open locus_cursor;
	
	loop_loci: loop
		fetch locus_cursor into locus;
		if finished_loci=1 then leave loop_loci; end if;

		create temporary table vala (select * from zz401_Depth_master where catalog_id = locus);
		alter table vala add index catalog_id_index (catalog_id asc);
		alter table vala add index sample_id_index (sample_id asc);
		alter table vala add index depth_index (depth asc);
		create temporary table arda (select sample_id, sum(depth) as depth from vala group by sample_id);
		alter table arda add index sample_id_index (sample_id asc);
		alter table arda add index depth_index (depth asc);

		select avg(depth) into medel from arda;
		select stddev_samp(depth) into savv from arda;

		insert into zz403_Depth_per_sample_per_locus (cat_id, mean, sdev)
		values (locus, medel, savv);

		drop temporary table vala;
		drop temporary table arda;

	end loop;
	close locus_cursor;
end block1;
END