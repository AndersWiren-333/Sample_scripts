--Running a script with parameters at the command line

mysql -e "set @firstname:='$1'; set @lastname='$2'; source run.sql;"