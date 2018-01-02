create table z05_Fluidigm_order_X (select cat_id as namn, (select insert(seq, col+1, 1, snp)) as seq from z04_X);
alter table z05_Fluidigm_order_X add index namn_index (namn asc);
create table z05_Fluidigm_order_Y (select cat_id as namn, (select insert(seq, col+1, 1, snp)) as seq from z04_Y);
alter table z05_Fluidigm_order_Y add index namn_index (namn asc);


select * from z05_Fluidigm_order_X into outfile '/Users/aswn0002/Desktop/Fluidigm_order_X.txt'
fields terminated by '\t' lines terminated by '\n';

select * from z05_Fluidigm_order_Y into outfile '/Users/aswn0002/Desktop/Fluidigm_order_Y.txt'
fields terminated by '\t' lines terminated by '\n';