--Running a script with parameters at the command line (or from e.g. a perl system call)

mysql -h hostname -u username -ppassword < script.sql
mysql -h hostname -u username -ppassword -e "set @var1:='Elsa Brandstrom'; set @var2:=3; source path/to/script.sql";
system("mysql -h hostname -u username -ppassword -e \"set \@parameter1:='value1'; set \@parameter2:=2; source C:/path/to/script.sql;\"");

--  Databases
create database tolkien;		-- Create a database
create database if not exists tolkien;		-- Create a database if it doesn't already exist
use tolkien;					-- Select a database to work with
drop database tolkien;			-- Delete a database
show databases;					-- List all databases on current server

-- Tables 
create table if not exists tablename (
	colname1 coltype1 colAttribute1,
	colname2 coltype2, colAttribute2,
	primary key(colname),
	foreign key (colname) references tablename (colname)
)
drop table tableName;			--  Delete a table
show tables;					-- List all tables in current (or default) database
alter table add colDefinition;
alter table drop colname;
alter table add foreign key (colname) references tableName (colname);
alter table drop foreign key constraintName;

-- Rows
insert into tablename values (1,"Greta",3);		-- Insert a new row
insert into tablename values (1,"Greta",3), (2,"Hans",4);	-- Insert two rows
insert into tablename (colname1, colname3)		-- Insert values into specific columns
	values ("Anna", "Åke");
delete from tablename where col1="Anders";
update tablename set col3 = "Excellent" where col1 = "Karin";

-- Write a MySQL specific statement
/*!select "This is a MySQL specific statement" as ''*/;

-- Write a query
select * from `table1` where `column2` = "string" and `column3` = 3;

-- Write a parameterized query (a prepared statement, safer against injection attacks when the
-- parameters (questionmarks) in the query are supplied by a user)
prepare query1 from 'select * from `inventarier` where `kategori` = ? and `inventarienummer` = ?';	-- You can't use a parameter to specify which table to work with
execute query1 using 'Elektronik', 1;
deallocate prepare query1;

-- Declare and initialize a (user defined) variable
set @var1 = "envanliggronskasrikadrakt";

-- Use a variable
select * from `table1` where `column2` = @var1;

-- Print something (Mysql doesn't have a print function, but this is a workaround)
select "Print this string";
select "Print this string" as '';

-- End the script
exit