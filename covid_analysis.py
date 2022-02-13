#!/usr/bin/env python3
from data_analyzer import *
import mysql.connector

# Create connection to database
con = mysql.connector.connect(**config)
cursor = con.cursor(buffered=True)

# Delete any existing tables from database
run_sql_command(commands, 'DROP', cursor)

# Create testing_data and case_info tables
run_sql_command(TABLES, 'CREATE', cursor)

# Import DOH csv into the created tables
run_sql_command(imports, 'LOAD', cursor)

# Do data analysis
print('Analyzing data....')
for location in PLACES:
    analyzer(location, con, cursor)
    print('{} done'.format(location))

# Map data
print('Creating maps....')
data_mapper()
# Close connection to database
cursor.close()
con.close()
