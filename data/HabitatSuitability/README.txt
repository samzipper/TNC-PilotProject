These data are intrinsic potential habitat for various salmonid species, downloaded from the NOAA website:
https://www.westcoast.fisheries.noaa.gov/maps_data/maps_and_gis_data.html
Dataset name: "Intrinsic Potential GIS Data - North-central California Coast updated January 2017 56MB"
Downloaded on November 6, 2018

For each of the three species, I clipped the shapefile for the entire Northern CA Region to the Navarro River Watershed
using the QGIS Vector->Research tools->Select by Location tool, and reprojected to EPSG26910 (NAD83/UTM Zone 10N).

The files are described in detail with the Word docs and PDFs in the 'docs' subdirectory. These are the fields most likely to be relevant:

For Coho:
CO_IP_CURVE: use for mapping - shows the quality of the IP (from 0.0 to 1.0). 
COIPINT = Integrated IP: use if comparing to other watersheds. The area over which IP was calculated is taken into account (values are > 1). 
To get total IP miles for a population, sum all the COIPINT records in that population. 
CO_IP215_CURV: use for mapping - shows the quality of the IP (from 0.0 to 1.0) with a 21.5 degree celsius mean August air temperature mask. 
COIPINT_215 = Integrated IP: use if comparing to other watersheds. The area over which IP was calculated is taken into account (values are > 1) with a 21.5 degree celsius mask. 
To get total IP miles for a population with a 21.5 degree celsius mask, sum all the COIPINT_215 records in that population

For Chinook: 
CHK_IP_CURVE: use for mapping - shows the quality of the IP (from 0.0 to 1.0).
CHKIPINT = Integrated IP: use if comparing to other watersheds. The area over which IP was calculated is taken into account (values are > 1).

For Steelhead:
ST_IP_CURVE12_RE: use for mapping - shows the quality of the IP (from 0.0 to 1.0). 
STIPINT_RE = Integrated IP: use if comparing to other watersheds. The area over which IP was calculated is taken into account (values are > 1). 

Email from Jen, 11/6/2018: 
For coho, you want to use the CO_IP215_CURV, with the temperature mask. That is the standard field that NMFS and others used for mapping and analysis. For steelhead you want the field called ST_IP_CURVE12_RE, and for Chinook you want CHK_IP_CURVE – they didn’t use temperature masks for them.