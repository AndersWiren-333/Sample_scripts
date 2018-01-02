create table x2_Final_list (
cat_id int);
load data local infile '/Users/aswn0002/Desktop/whitelist2'
into table x2_Final_list
fields terminated by ';'
lines terminated by '\n';
create table x3 (
select x1_First_SNP_list.cat_id, col, p_nuc, q_nuc, seq from x1_First_SNP_list
inner join x2_Final_list on x1_First_SNP_list.cat_id = x2_Final_list.cat_id);
create table x4 (select cat_id, concat('[',p_nuc,'/',q_nuc,']') snp from x3);
create table x5 (
select x3.cat_id, col, seq, snp from x3
inner join x4 on x3.cat_id = x4.cat_id);
drop table x3;
drop table x4;
create table x6_Fluidigm_order (
select cat_id as namn, (select insert(seq, col+1, 1, snp)) as seq from x5
where cat_id = ANY (select cat_id from x5));
drop table x5;
select * from x6_Fluidigm_order
into outfile '/Users/aswn0002/Desktop/Fluidigm_order.txt'
fields terminated by '\t'
lines terminated by '\n';