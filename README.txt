You are now in SpatialHydrology directory. This is home directory to run SpatHy
with Django www-server.

SpatHy model with sample data can be found in Spathy directory.
Vihti arcGIS data is in slauniainen directory.

To fetch FMI data from 'weather' database at Luke
see ildata.py in weather/ildata directory.
Run 'python ildata.py --help' for usage.

To fetch weather forecast from Norway database
see weather/yrforcast/yrforecast.py and the first
function forecast_yr (experimantal only). To write
dataframe returned to excel see the last function
write_excel. Run 'python yrforecast.py -h' for usage.

All python files require now python3.4 or higher.  It is recommended to create
python virtual environment (name it e.g. 'spathy'), activate it and install the following
packages with pip (pip install <package.name>):

1. scipy
2. matplotlib
3. cython
4. netcdf4
5. pandas
6. pillow
7. psycopg2
8. seaborn
9. xlswriter

In 'spathy/lib64/python3.4/site-packages' you should then have (in alphabetcal order):

1. cftime
2. cycler
3. Cython
4. dateutil
5. django
6. kiwisolver
7. matplotlib
8. netCDF4
9. numpy
10. pandas
11. Pillow
12. psycopg2
13. pymysql
14. pyparsing
15. python_dateutil
16. pytz
17. pyximport
18. scipy
19. seaborn
20. spotpy
21. xlsxwriter

