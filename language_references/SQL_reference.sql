--Running a script with parameters at the command line (or from e.g. a perl system call)

mysql -h hostname -u username -ppassword < script.sql
mysql -h hostname -u username -ppassword -e "set @var1:='Elsa Brandstrom'; set @var2:=3; source path/to/script.sql";
system("mysql -h hostname -u username -ppassword -e \"set \@parameter1:='value1'; set \@parameter2:=2; source C:/path/to/script.sql;\"");

-- Create a database
create database tolkien;

-- Select a database to work with
use tolkien;

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