create table x1 (
select sumstats.tag_id, col, pop_id, p_nuc, q_nuc, n, p, obs_het, obs_hom, exp_het, exp_hom, pi, fis, seq from sumstats
inner join catalog_tags on sumstats.tag_id = catalog_tags.tag_id
where n between 50 and 55
and p between 0 and 0.7
and obs_het between 0 and 0.5
and col between 40 and 100
and fis > 0);
create table x2 (
select cat_id, tag_id, snps, max_pct, alleles, ratio from catalog_index
where snps between 1 and 1 
and alleles between 2 and 2
and max_pct between 0 and 100);
create table x3 (
select x2.cat_id, x1.tag_id, col, pop_id, p_nuc, q_nuc, n, p, obs_het, obs_hom, exp_het, exp_hom, pi, fis, seq, snps, max_pct, alleles, ratio from x1
inner join x2 on x1.tag_id = x2.tag_id);
drop table x1;
drop table x2;
create table x4 (
select catalog_id, geno_map from markers);
create table x1_First_SNP_list (
select cat_id, tag_id, col, pop_id, p_nuc, q_nuc, n, p, obs_het, obs_hom, exp_het, exp_hom, pi, fis, seq, snps, max_pct, alleles, ratio, geno_map from x3
inner join x4 on x3.tag_id = x4.catalog_id);
drop table x3;
drop table x4;
select cat_id from x1_First_SNP_list
into outfile '/Users/aswn0002/Desktop/whitelist1'
fields terminated by ';'
lines terminated by '\n';