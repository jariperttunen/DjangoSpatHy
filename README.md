You are now in DjangoSpaFHy directory. This is home directory to 
run SpaFHy with Django www-server.

DjangoSpaFHy assumes that the required modules (from git repository)
are cloned (git clone) side by side as follows:

    DjangoSpaFHy
    SpaFHy
    SpaFHyData
    SpaFHyWeather

SpaFHy is the SpatialHydrology model itself.
SpaFHyData shall contain catchment data for
the areas to be simulated. (Data is not
necessary in git, contact Samuli Launiainen
for details).
SpaFHyWeather contain python scripts to fetch
weather data from Luke 'weather' database 
and from Norway (forecasts, at the moment experimental)
See the README files in these modules for deatils.

Python files in DjangoSpaFHy, SpaFHyWeather and in SpaFHy
require now python3.4 or higher and python packages
they depend on. It is recommended (read: mandatory)
to create python virtual environment  common to 
these three modules (name it e.g. 'spathy') 
and activate it. To create a virtual environment 
in python see 'https://docs.python.org/3/library/venv.html'
Then  install the following packages with 
pip (pip install <package name>):

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

