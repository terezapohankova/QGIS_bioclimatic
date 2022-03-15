# CalcBioClimVariables

The script calculated selected bioclimatic variables (NDVI, Land-surface Temperature, Albedo, Ground Heat Flux) and their corresponding Zonal Statistic.

**The provided data for testing are covering Olomouc Area in CZE**

## Requirements
* QGIS (tested for v. 3.20 and 3.22) with Python (tested for v. 3.9.5.)
* numpy v. 1.21.2
* tifffile v. 2022.2.2

## Instalation

To intregrate the script with other user scripts in QGIS, place the downloaded file in folder **scripts** in path _root\QGIS\QGIS3\profiles\default\processing\scripts_. 
After starting QGIS the scripts should be visible in Processing Toolbox under the Example scripts in Python Scripts:
<p align="center">
<img src="https://user-images.githubusercontent.com/60270092/158354957-985b1ff0-b0e2-4308-8925-8620a484ecb8.png">
</p>
  
## Data

The data required are 
* **CSV with air temperature** [°C]

  - the file contains two columns: sensing date RRRRMMDD and temperature value 

  - test dataset is provided
```
20180615,15.5
20170212,5.45
```

* **Zones for zonal statistcs** (optional)
  
  - as vector layer, ideally a geopackage
 
* **Clipping mask**
  
  - a vector layer for clipping the satellite image by ROI
 
* **Landsat 8 or Landsat 9 images**
  
  - older Landsats are also posibble, but a lot of the calculations are tested only ony L8 
  - the folder should be unzipped but with original name 
 
  ## Running the tool
  After the script is placed into the scripts folder, you should be able to run it:

<p align="center">
<video src="https://user-images.githubusercontent.com/60270092/158362055-85e2bb64-4af4-4a60-9f72-49221dd11036.mp4"></video>
</p>

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.
  
