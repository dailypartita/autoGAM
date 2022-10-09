README release HYDE version 3.2.1, d.d. September 12th 2017.

This HYDE version 3.2.1 replaces former 2017_beta_releases 3.2.000 (from august 2016) and 3.2.001 (from october 2016). 
Major updates concern use of new FAO data until the year 2014 (instead of 2012), improved computations of (physical) Rice areas, and an adjustment of the cropland area for India before 1000 CE.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Disclaimer: The Creative Commons License (CC BY 3.0) applies to all of the HYDE data. (see https://creativecommons.org/licenses/by/3.0/)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Reference:  Klein Goldewijk et al. (2017), Earth Syst. Sci. Data, 9, 927-953. 
or in full: Klein Goldewijk, K., A. Beusen, J.Doelman and E. Stehfest (2017), Anthropogenic land use estimates for the Holocene; HYDE 3.2, Earth System Science Data, 9, 927-953.

This database presents an update and expansion of the History Database of the Global Environment (HYDE, v 3.2.1). 
HYDE is and internally consistent combination of updated historical population estimates and land use. Categories include cropland, 
with a new distinction into irrigated and rain fed crops (other than rice) and irrigated and rain fed rice. 
Also grazing lands are provided, divided into more intensively used pasture, converted rangeland and non-converted natural (less intensively used) rangeland. 
Population is represented by maps of total, urban, rural population and population density as well as built-up area. 
The period covered is 10 000 BCE to 2015 CE. Spatial resolution is 5 arc minutes (approx. 85 km2 at the equator), the files are in Arcmap asciigrid format.

- Each header of the gridded .asc files looks like:

ncols 4320
nrows 2160
xllcorner -180
yllcorner -90
cellsize 0.0833333
NODATA_value -9999

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
1) The HYDE 3.2.1 release contains 3 scenarios (Baseline estimate, Lower estimate, Upper estimate) and a General Files directory.

### Directory /General files ###

This directory contains general files:

garea_cr.asc	map with total gridcell area in km2, spherical Earth
landlake_cr.asc	map with (land-sea mask on 5'), including lakes (land = 1, lakes = 0, rest = -9999)
maxln_cr.asc	map with maximum landarea available per gridcell in km2

iso_cr.asc	map with national iso numerical country codes 
sub_iso_cr.asc	map with sub-national codes (provinces, states)

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
2) Each scenario contains 4 directories: png, zip, txt, anthromes

### Directory /png ### 

Contains all figures for each category for all years in graphical .png format.

### Directory /zip ###

Conatains zipped files with the following formats and content:
yearBC_pop.zip (population) and yearBC_lu.zip (landuse) for the BCE period ; yearAD_pop.zip (population) and yearAD_lu.zip (landuse) for the CE period.
Time intervals are 1000 yr for the BCE period, then 100 yr till 1700, 10 yr till 2000, and from 2000 - 2015 1 year timesteps.

- Each xxxx_pop.zip contains: 
popc_<yr>.asc (population counts, in inhabitants/gridcell)
popd_<yr>.asc (population density, in inhabitants/km2 per gridcell)
rurc_<yr>.asc (rural population counts, in inh/gridcell)
urb_<yr>.asc  (urban population counts, in inh/gridcell)
uopp_<yr>.asc (total built-up area, such as towns, cities, etc, in km2 per grid cell)

- Each xxxx_lu.zip contains: 
cropland<yr>.asc   	(total cropland area, in km2 per grid cell), after 1960 identical to FAO's category 'Arable land and permanent crops'.
grazing<yr>.asc    	(total land used for grazing, in km2 per grid cell), after 1960 identical to FAO's category 'Permanent Pasture'.
pasture<yr>.asc    	(total pasture area, in km2 per grid cell), defined as Grazing land with an aridity index > 0.5, assumed to be more intensively managed.(converted in climate models)
rangeland<yr>.asc  	(total rangeland area, in km2 per grid cell), defined as Grazing land with an aridity index < 0.5, assumed to be less or not managed.(not converted in climate models)
conv_rangeland<yr>.asc (total converted rangeland area, in km2 per grid cell), defined as rangeland being placed on a potential forest area, and thus being converted in climate models.
ir_norice<yr>.asc  	(irrigated other crops area (no rice) area, in km2 per grid cell).
rf_norice<yr>.asc  	(rainfed other crops area (no rice) area, in km2 per grid cell).
ir_rice<yr>.asc    	(irrigated rice area (no rice), in km2 per grid cell).
rf_rice<yr>.asc     (rainfed rice area (no rice), in km2 per grid cell).
tot_irri<yr>.asc    (total actual irrigated area, in km2 per grid cell).
tot_rainfed<yr>.asc (total rainfed other crops (no rice),area, in km2 per grid cell).
tot_rice<yr>.asc    (total rice area, in km2 per grid cell).

### Directory /txt
- For each scenario summary files (plain ascii) are presented with the country or IMAGE 3.0 regional breakdown totals (in km2). See also www.pbl.nl/image for regional breakdown.
pop_c.txt        	(total population per country over time, in 1000)
pop_r.txt        	(total population per region over time, in 1000)
urb_c.txt	 	(total urban population per country over time, in 1000)
urb_r.txt	 	(total urban population per region over time, in 1000)
rur_c.txt	 	(total rural population per country over time, in 1000)
rur_r.txt	 	(total rural population per region over time, in 1000)
uopp_c.txt 	 	(total built-up area per country over time, in km2 per grid cell)
uopp_r.txt	 	(total built-up area per region over time, in km2 per grid cell)
croplandc.txt	 	(total cropland area for each seperate country over time, in km2)
croplandr.txt	 	(total cropland area for each region over time, in km2)
rf_rice_r.txt	   	(total rainfed rice for each region over time, in km2)
rf_rice_c.txt	   	(total rainfed rice for each country over time, in km2)
rf_norice_r.txt	   	(total rainfed other crops (no rice)for each region over time, in km2)
rf_norice_c.txt	   	(total rainfed other crops (no rice)for each country over time, in km2)
ir_rice_c.txt 	   	(total irrigated rice for each country over time, in km2)
ir_rice_r.txt 	   	(total irrigated rice for each region over time, in km2)
ir_norice_r.txt	   	(total irrigated other crops (no rice)for each region over time, in km2)
ir_norice_c.txt	   	(total irrigated other crops (no rice)for each country over time, in km2)

tot_irri_r.txt	   	(total irrigated area for each region over time, in km2)
tot_irri_c.txt	   	(total irrigated area for each country over time, in km2)
tot_rainfed_r.txt  	(total rainfed area for each region over time, in km2)
tot_rainfed_c.txt  	(total rainfed area for each country over time, in km2)

grazingc.txt	 	(total grazingland area for each seperate country over time, in km2)
grazingr.txt 	 	(total grazingland area for each region over time, in km2)
pasture_c.txt		(total pasture area for each country over time, in km2)
pasture_r.txt		(total pasture area for each region over time, in km2)
conv_rangeland_c.txt  	(total converted rangeland for each seperate country over time, in km2)
conv_rangeland_r.txt  	(total converted rangeland for each region over time, in km2)
rangeland_c.txt  	(total rangeland for each seperate country over time, in km2)
rangeland_r.txt  	(total rangeland for each region over time, in km2)

### Directory /Anthromes ###

Contains a re-classification of HYDE 3.2 into Anthropogenic Biomes (Anthromes v2), using a classification system based on Ellis et al. (2010), 
also with /zip, /png and /txt subdirectories with the ascii gridded files, the graphic maps and the summaries txt files for each anthrome and a summary.

The Anthromes classes are:

"11 Urban"
"12 Dense settlements"

"21 Village, Rice"
"22 Village, Irrigated"
"23 Village, Rainfed"
"24 Village, Pastoral"

"31 Croplands, residential irrigated"
"32 Croplands, residential rainfed"
"33 Croplands, populated"
"34 Croplands, pastoral"

"41 Rangeland, residential"
"42 Rangeland, populated"
"43 Rangeland, remote"

"51 Semi-natural woodlands, residential"
"52 Semi-natural woodlands, populated"
"53 Semi-natural woodlands, remote"
"54 Semi-natural treeless and barren lands"

"61 Wild, remote - woodlands"
"62 Wild, remote - treeless & barren"
"63 Wild, remote - ice"

"70 No definition"

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Kees Klein Goldewijk (c.g.m.kleingoldewijk@uu.nl)
September 12, 2017

Utrecht University
Faculty of Geosciences
Copernicus Institute of Sustainable Development
Heidelberglaan 2
Utrecht, The Netherlands

and

PBL Netherlands Environmental Assessment Agency
Bezuidenhoudseweg 30
The Hague, The Netherlands
