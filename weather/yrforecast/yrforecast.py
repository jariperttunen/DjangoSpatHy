# -*- coding: utf-8 -*-
"""
Fetching weather forecast from yr.no

!In order to use the free weather data from yr no, you HAVE to display 
the following text clearly visible on your web page. The text should be a 
link to the specified URL.
"Weather forecast from Yr, delivered by the Norwegian Meteorological Institute and the NRK"
!Please read more about our conditions and guidelines at http://om.yr.no/verdata/
English explanation at http://om.yr.no/verdata/free-weather-data/

@author: Kersti Haahti, Luke 18.6.2018
"""
import os
import xml.etree.ElementTree as ET
import urllib
from matplotlib import pyplot as plt
import pandas as pd
import xlsxwriter
import numpy as np
import argparse

http_proxy='http://10.88.2.10:8080/'
https_proxy='https://10.88.2.10:8080/'

def forecast_yr(lat, lon, hourly=False, plot=False, excel=True,proxy=True):
    """
    Fetch weather forecast from yr.no based on coordinates.
    ! Works only for Finland at the moment
    Args:
        lat (float): latitude in decimal degrees (wgs84)
        lon (float): longitude in decimal degrees (wgs84)
        hourly (boolean):
            True fetches hourly forecast for 48 hours
            False fetches 10 day forecast at 6 hour interval
        plot (boolean): plot fetched data
    Returns:
        df (DataFrame):
            precipitation (mm per time interval)
            temperature (celsius)
            pressure (hPa)
            wind speed (mps)
            symbol_id, symbol_description (see https://om.yr.no/symbol/)
    """
    # Find closest location to coordinates from geonames database (http://www.geonames.org/)
    geoname = nearest_geoname_location(lat, lon)
    # Request xml from yr.no for geoname
    root, url = request_forecast(place=geoname['name'],
                                 region=geoname['region'],
                                 country=geoname['country'],
                                 hourly=hourly,proxy=proxy)

    # Parse data and save to dataframe
    df = forecast_to_df(root)

    # Plotting
    if plot:
        print("Plotting dataframe")
        df.plot(subplots=True,
                title='Forecast for %s, %s, %s (%.2f, %.2f) \nat %.2f km distance from searched point'
                %(geoname['name'], geoname['region'], geoname['country'], geoname['latitude'], geoname['longitude'], geoname['distance']))
        plt.xlabel('Weather forecast from Yr, delivered by the Norwegian Meteorological Institute and the NRK \n%s' % url, fontsize=8)
    if excel:
        write_excel('Forecast'+geoname['name']+'.xlsx',geoname['name'],df)
    return df

def forecast_to_df(root):
    """
    Parse data and return to dataframe.
    """
    # List with parsed data from each time step
    structure_data = [parse_timestep(child) for child in root.findall("./forecast/tabular/time")]
    # To dataframe
    df = pd.DataFrame(structure_data)
    # Set datetime index
    df.index =  pd.to_datetime(df['time (from)'].values, yearfirst=True)

    return df

def parse_timestep(child):
    """
    Parse timestep data.
    """
    # Initialize dict for parsed data
    parsed = dict()
    # Keys include strat and end time for time timestep
    for key in child.keys():
        parsed['time (' + key +')'] = child.attrib.get(key)
    # Variables for timestep are in childs
    for subchild in list(child):
        # precipitation, temperature, pressure
        if 'value' in subchild.attrib.keys():
            key = '%s (%s)' % (subchild.tag, subchild.attrib.get('unit'))
            parsed[key] = float(subchild.attrib.get('value'))
        # wind speed
        if 'mps' in subchild.attrib.keys():
            key = '%s (mps)' % subchild.tag
            parsed[key] = float(subchild.attrib.get('mps'))
        # symbol describing weather conditions
        if subchild.tag == 'symbol':
            parsed['symbol_code'] = int(subchild.attrib.get('numberEx'))
            parsed['symbol_description'] = subchild.attrib.get('name')

    return parsed

def request_forecast(place, region, country, hourly=False, proxy=True):
    """
    Request data from yr.no for given place, region, country.
    """
    # Form url
    url = "http://www.yr.no/place/%s/%s/%s/" %(country, region.replace(" ", "_"), place.replace(" ", "_"))
    url1 = url + "forecast.xml"
    if hourly:
        url1 = url + "forecast_hour_by_hour.xml"
    print("URL1",url1)
    response=None
    if proxy:
        #print("Luke proxies required")
        proxy_handler=urllib.request.ProxyHandler({'http':http_proxy,'https':https_proxy})
        # Request url and parse it to element tree
        opener = urllib.request.build_opener(proxy_handler)
        response = opener.open(url1)
    else:
        response = urllib.request.urlopen(url1)
    tree = ET.parse(response)

    return tree.getroot(), url

def nearest_geoname_location(lat, lon):
    """
    Finds nearest place in geonames database (FI.txt) based on coordinates.
    """
    # Read geonames for Finland from file
    #Location of this file, absolute path
    cwd=os.path.abspath(__file__)
    cwd=os.path.dirname(cwd)+'/'
    geonames = pd.read_csv(cwd+"FI/FI.txt", sep='\t',
                           header=None, usecols=[0,1,4,5,10],
                           names=['geoname_id','name','latitude','longitude','maakunta_id'])
    # Calculate distance to given coordinate point
    geonames['distance'] = distance(lat, lon, geonames['latitude'].values, geonames['longitude'].values)
    # Select nearest location
    nearest_geoname = geonames[geonames['distance']==min(geonames['distance'])]
    # Read information of regions (läänit vs maakunnat) - lääni is what is used in yr.no url as region
    regions = pd.read_csv(cwd+"FI/regions.txt", sep='\t',
                          header=None, names=['id','maakunta', 'laani'])
    # Add region name and country
    nearest_geoname.loc[:,'region']=regions[regions['id']==int(nearest_geoname['maakunta_id'])]['laani'].values
    nearest_geoname.loc[:,'country']='Finland'
    # return in dict format
    return {col: nearest_geoname[col].values[0] for col in nearest_geoname}

def distance(lat1, lon1, lat2, lon2):
    """
    Calculate distance, Haversine formula, between two coordinate points.
    """
    # Radius of earth in km
    R = 6373.0
    # Coordinates in radians
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    # Calculate distance
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (np.sin(dlat/2))**2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon/2))**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    distance = R * c

    return distance

def write_excel(file_name,sheet_name,df):
    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')
    df.to_excel(writer,sheet_name=sheet_name)
    writer.save()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-x',type=float,dest='x',help='Latitude, decimal degrees (wgs84)')
    parser.add_argument('-y',type=float,dest='y',help='Longitude, decimal degrees (wgs84)')
    parser.add_argument('-w','--hourly',dest='h',action="store_true",help='48 hour forecast')
    args = parser.parse_args()
    df = forecast_yr(args.x,args.y,hourly=args.h)
