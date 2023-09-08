import certifi
import ssl
from datetime import datetime
from tzwhere import tzwhere
from pytz import timezone, utc

import geopy.geocoders
from geopy.geocoders import Photon
ctx = ssl.create_default_context(cafile=certifi.where())
geopy.geocoders.options.default_ssl_context = ctx

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle

from skyfield.api import Star, load, wgs84
from skyfield.data import hipparcos, mpc, stellarium
from skyfield.projections import build_stereographic_projection
from skyfield.constants import GM_SUN_Pitjeva_2005_km3_s2 as GM_SUN

import warnings
warnings.filterwarnings("ignore")

class TelescopeLocation:
    def __init__(self):
        self.name = None
        self.latitude = 0
        self.longitude = 0
        self.timezone = "UTC"

    def set_name(self, name):
        self.name = name

    def set_latitude(self, lat):
        self.latitude = lat

    def set_longitude(self, lon):
        self.longitude = lon
    
    def set_timezone(self, tz):
        self.timezone = tz

class GBTLocation(TelescopeLocation):
    def __init__(self):
        super().__init__()
        self.assign_values()

    def assign_values(self):
        self.set_name("GBT")
        self.set_latitude(38.433129)
        self.set_longitude(-79.839839)

def load_data():
    eph = load('de421.bsp')
    
    with load.open(hipparcos.URL) as f:
        stars = hipparcos.load_dataframe(f)
    
    url = ('https://raw.githubusercontent.com/Stellarium/stellarium/master'
           '/skycultures/modern_st/constellationship.fab')

    with load.open(url) as f:
        constellations = stellarium.parse_constellations(f)
        
    return eph, stars, constellations

def collect_celestial_data(location, when):
    # Get site information
    lat, long = location.latitude, location.longitude
    local = timezone(location.timezone)
    
    # Get current time as datetime object
    dt = datetime.strptime(when, '%Y-%m-%d %H:%M')

    # Get UTC time from local timezone and datetime
    local_dt = local.localize(dt, is_dst=None)
    utc_dt = local_dt.astimezone(utc)
    ts = load.timescale()
    t = ts.from_datetime(utc_dt)

    # Load ephemerides
    sun = eph['sun']
    earth = eph['earth']

    # Define an observer on the Earth's surface
    observer = wgs84.latlon(latitude_degrees=lat, longitude_degrees=long).at(t)
    position = observer.from_altaz(alt_degrees=90, az_degrees=0)
    
    # Center the observation point in the middle of the sky
    ra, dec, distance = observer.radec()
    center_object = Star(ra=ra, dec=dec)

    # Build a sky projection with 180-degree view
    center = earth.at(t).observe(center_object)
    projection = build_stereographic_projection(center)
    field_of_view_degrees = 180.0

    # Calculate star positions and project them onto a plain space
    star_positions = earth.at(t).observe(Star.from_dataframe(stars))
    stars['x'], stars['y'] = projection(star_positions)
    
    edges = [edge for name, edges in constellations for edge in edges]
    edges_star1 = [star1 for star1, star2 in edges]
    edges_star2 = [star2 for star1, star2 in edges]

    return stars, edges_star1, edges_star2

def create_star_chart(location, when, chart_size, max_star_size):
    stars, edges_star1, edges_star2 = collect_celestial_data(location, when)
    limiting_magnitude = 10
    bright_stars = (stars.magnitude <= limiting_magnitude)
    magnitude = stars['magnitude'][bright_stars]
    
    fig, ax = plt.subplots(figsize=(chart_size, chart_size),facecolor='#041A40')

    marker_size = max_star_size * 10 ** (magnitude / -2.5)
    ax.scatter(stars['x'][bright_stars], stars['y'][bright_stars],
               s=marker_size, color='white', marker='.', linewidths=0,
               zorder=2)
    
    # Draw the constellation lines.
    xy1 = stars[['x', 'y']].loc[edges_star1].values
    xy2 = stars[['x', 'y']].loc[edges_star2].values
    lines_xy = np.rollaxis(np.array([xy1, xy2]), 1)

    ax.add_collection(LineCollection(lines_xy, colors='#ffff', linewidths=0.15))


    gbt_pointing = Circle((0,0), 0.1)
    ax.add_artist(gbt_pointing)
    # set the aspect ratio of the plot to be equal
    ax.set_aspect('equal')
    
    # other settings
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    plt.axis('off')
    when_datetime = datetime.strptime(when, '%Y-%m-%d %H:%M')
    plt.title(f"Observation Location: {location.name}, Time: {when_datetime.strftime('%Y-%m-%d %H:%M')}", loc='right',color = 'white', fontsize=10)
    
    plt.show()
    plt.close()

def generate_star_chart(locations, whens):
    for location in locations:
        for when in whens:
            chart_size = 12
            max_star_size = 100

            # generate the plot
            create_star_chart(location, when, chart_size, max_star_size)

            # save the plot
            when_datetime = datetime.strptime(when, '%Y-%m-%d %H:%M')
            filename = f"{location}_{when_datetime.strftime('%Y%m%d_%H%M')}.png"
            #plt.savefig(filename, dpi=300)
            #plt.savefig(filename, format='svg', dpi=1200)
            plt.close()

            # print confirmation message
            print(f"Plot saved for location {location} and time {when}")

eph, stars, constellations = load_data()
location = GBTLocation() #'Green Bank, WV'
when = '2023-04-21 00:00'
chart_size=16
max_star_size=400
create_star_chart(location, when, chart_size, max_star_size)