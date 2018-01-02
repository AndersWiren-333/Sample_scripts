create table zz01_seq (select t1.cat_id, t2.seq from XXX_All_Loci_final as t1
inner join catalog_tags as t2 on t1.cat_id = t2.tag_id
where pop_name != 'Alces');

create table zz02_seq_col_nucs (select t1.*, t2.col, t2.rank_1 as p_nuc, t2.rank_2 as q_nuc 
from zz01_seq as t1 left join catalog_snps as t2 on t1.cat_id = t2.tag_id);

drop table zz01_seq;

create table zz03_concat (select cat_id, concat('[',p_nuc,'/',q_nuc,']') as snp, col, seq from zz02_seq_col_nucs as t1);
alter table zz03_concat add column namn varchar(45) after cat_id;
ALTER TABLE zz03_concat CHANGE COLUMN cat_id cat_id VARCHAR(45);
update zz03_concat set namn = concat('Anders_',cat_id);

drop table zz02_seq_col_nucs;

# Add info for Ida's moose SNPs
# Begin by saving her three worksheets separately and then saving the excel files as csv on the desktop

create table zz04_Idas (cat_id int(11), snp varchar(5), col int(10), seq tinytext);

load data local infile '/Users/aswn0002/Desktop/alces_Y.csv'
into table zz04_Idas
fields terminated by ';'
lines terminated by '\r';

load data local infile '/Users/aswn0002/Desktop/alces_X.csv'
into table zz04_Idas
fields terminated by ';'
lines terminated by '\r';

load data local infile '/Users/aswn0002/Desktop/alces_Y_cons.csv'
into table zz04_Idas
fields terminated by ';'
lines terminated by '\r';

update zz04_Idas set snp = NULL where col=0;
update zz04_Idas set col = NULL where col=0;

alter table zz04_Idas add column namn varchar(45) after cat_id;
ALTER TABLE zz04_Idas CHANGE COLUMN cat_id cat_id VARCHAR(45);
update zz04_Idas set namn = concat('Ida_',cat_id);

alter table zz03_concat rename to zz04_Anders;

insert into zz04_Anders (select * from zz04_Idas);
drop table zz04_Idas;
ALTER TABLE zz04_Anders CHANGE COLUMN cat_id cat_id INT(11) UNSIGNED NOT NULL;
alter table zz04_Anders rename to zz01_Fluidigm_details;

alter table zz01_Fluidigm_details add column category varchar(45) after namn;
update zz01_Fluidigm_details set category = 'SNP' where col is not null;
update zz01_Fluidigm_details set category = 'dummy' where col is null;

create table zz01_Fluidigm_details2 (select t1.*, (round((length(seq)/2), 0)) as col2 from zz01_Fluidigm_details as t1);
drop table zz01_Fluidigm_details;

create table zz01_Fluidigm_details3 (select t1.*, (col2-1) as col3 from zz01_Fluidigm_details2 as t1);
drop table zz01_Fluidigm_details2;
alter table zz01_Fluidigm_details3 rename to zz01_Fluidigm_details2;
alter table zz01_Fluidigm_details2 drop column col2;
alter table zz01_Fluidigm_details2 change col3 col2 int(11);

create table zz01_Fluidigm_details3 (select t1.*, mid(seq, col2+1, 1) as p_nuc from zz01_Fluidigm_details2 as t1);
drop table zz01_Fluidigm_details2;
alter table zz01_Fluidigm_details3 rename to zz01_Fluidigm_details2;

alter table zz01_Fluidigm_details2 add column q_nuc varchar(1) after p_nuc;
update zz01_Fluidigm_details2 set q_nuc = 'G' where p_nuc = 'A';
update zz01_Fluidigm_details2 set q_nuc = 'C' where p_nuc = 'T';
update zz01_Fluidigm_details2 set q_nuc = 'A' where p_nuc = 'G';
update zz01_Fluidigm_details2 set q_nuc = 'T' where p_nuc = 'C';

create table zz01_Fluidigm_details3 (select t1.*, concat('[',p_nuc,'/',q_nuc,']') as snp2 from zz01_Fluidigm_details2 as t1);
drop table zz01_Fluidigm_details2;

update zz01_Fluidigm_details3 set snp = snp2 where category = 'dummy';
update zz01_Fluidigm_details3 set col = col2 where category = 'dummy';
alter table zz01_Fluidigm_details3 drop column col2;
alter table zz01_Fluidigm_details3 drop column snp2;
alter table zz01_Fluidigm_details3 drop column p_nuc;
alter table zz01_Fluidigm_details3 drop column q_nuc;
alter table zz01_Fluidigm_details3 rename to zz01_Fluidigm_details;

create table zz02_Fluidigm_order (
select namn, (select insert(seq, col+1, 1, snp)) as seq from zz01_Fluidigm_details);
#where cat_id = ANY (select cat_id from x5));

select * from zz02_Fluidigm_order
into outfile '/Users/aswn0002/Desktop/Ungulate_Fluidigm_order.txt'
fields terminated by '\t'
lines terminated by '\n';

# Check that allele depths are not unreasonably high
#create table zz05_depth (select catalog_id as cat_id, avg(depth) as Mean_al_depth, min(depth) as Min_al_depth, max(depth) as Max_al_depth from matches
#group by catalog_id
#having catalog_id in (select distinct cat_id from XXX_All_Loci_final));