# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""


from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessingAlgorithm,
                       QgsProcessingParameterFolderDestination,
                       QgsProcessingParameterFile,
                       QgsProcessingParameterVectorLayer,
                       QgsProcessingParameterBoolean,
                       QgsProcessingParameterDefinition
                       )
from qgis import processing
from osgeo import gdal, osr
import os, sys
import numpy as np
import tifffile as tf

    
   
class ExampleProcessingAlgorithm(QgsProcessingAlgorithm):
    INPUT_FOLDER = 'INPUT_FOLDER'
    MASK_LAYER = 'MASK_LAYER'
    OUTPUT = 'OUTPUT'
    AIR_TEMP = 'AIR_TEMP'
    STAT_ZONES = 'STAT_ZONES'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return ExampleProcessingAlgorithm()

    def name(self):
        return 'Calculation of selected bioclimatic variables'

    def displayName(self):
        return self.tr('Calculation of selected bioclimatic variables')

    def group(self):
        return self.tr('Example scripts')

    def groupId(self):
        return 'examplescripts'

    def shortHelpString(self):
        return self.tr("The script calculates Albedo, NDVI, Land Surface Temperature and Ground Heat Flux from Landsat images. \n\
        An optional calculation of Zonal Statistic is possible, if the zones are provided. \n\
        Author: Mgr. Tereza Nováková, Palacky University Olomouc, 2022")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(self.INPUT_FOLDER, self.tr('Input Folder with Landsat 8 or 9 Images'), behavior=QgsProcessingParameterFile.Folder))
        self.addParameter(QgsProcessingParameterVectorLayer(self.MASK_LAYER, self.tr('Polygon clipping mask (GPKG)'), defaultValue=None))
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT, self.tr('Output Folder')))
        self.addParameter(QgsProcessingParameterFile(self.AIR_TEMP, self.tr('Air Temperature (CSV)'), behavior=QgsProcessingParameterFile.File))
        
        #advanced
        ZONAL_STAT = QgsProcessingParameterBoolean('calc_zonal_stat', self.tr('Calculate Zonal Statistics '), defaultValue=None)
        ZONAL_STAT.setFlags(ZONAL_STAT.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(ZONAL_STAT)
        
        ZONES = QgsProcessingParameterVectorLayer(self.STAT_ZONES, self.tr('Input Zones'), defaultValue=None, optional=True)
        ZONES.setFlags(ZONES.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(ZONES)
        
    
    def rasterToFile(self, input_array, src_dataset_path, output_path): #https://gist.github.com/jkatagi/a1207eee32463efd06fb57676dcf86c8        
        cols = input_array.shape[1]
        rows = input_array.shape[0]
        
        dataset = gdal.Open(src_dataset_path, gdal.GA_ReadOnly)
        originX, pixelWidth, b, originY, d, pixelHeight = dataset.GetGeoTransform() 
        driver = gdal.GetDriverByName('GTiff')
        band_num = 1
        GDT_dtype = gdal.GDT_Float32
        outRaster = driver.Create(output_path, cols, rows, band_num, GDT_dtype)
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
        outband = outRaster.GetRasterBand(band_num)
        outband.WriteArray(input_array)
        prj=dataset.GetProjection()
        outRasterSRS = osr.SpatialReference(wkt=prj)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        
        return True
        
    def getImageName(self, input_folder):
        self.path_list_folder_all = []
        self.name_clipped_tif = []
        self.input_folder = input_folder
        
        for root, dirs, files in os.walk(self.input_folder, topdown=False):
            self.path_list_folder = []
            for name in files:
                if name.endswith('.TIF') and name.startswith('LC'):
                    self.path_list_folder.append(os.path.join(root, name))
                    #'C:\\Users\\Tereza\\Documents\\PhD\\zkousky\\2_vyvoj_pro_FOSS\\snimky\\LC08_L1TP_190025_20170528_20200903_02_T1\\LC08_L1TP_190025_20170528_20200903_02_T1_B5.TIF'
            if self.path_list_folder:
                self.path_list_folder_all.append(self.path_list_folder)          
            
        return self.path_list_folder_all

    def clipLandsatBand(self, input_band, out_path, mask):
        self.input_band = input_band
        self.out_path = out_path
        self.mask = mask
        
        alg_params = {
            'ALPHA_BAND': False,
            'CROP_TO_CUTLINE': True,
            'DATA_TYPE': 0,  # Use Input Layer Data Type
            'EXTRA': '',
            'INPUT': self.input_band,
            'KEEP_RESOLUTION': True,
            'MASK': self.mask,
            'MULTITHREADING': False,
            'NODATA': None,
            'OPTIONS': '',
            'SET_RESOLUTION': True,
            #'SOURCE_CRS': self.input_band.crs(),
            #'TARGET_CRS': self.input_band.crs(),
            'X_RESOLUTION': None,
            'Y_RESOLUTION': None,
            'OUTPUT': self.out_path
        }
        
        #self.feedback.pushInfo(str(self.out_path))
        processing.run('gdal:cliprasterbymasklayer', alg_params, is_child_algorithm = True)
        return

    def readCSV(self, input_csv):
        items = []

        csv_open = open(input_csv, 'r')
        for line in csv_open:
            items.append(line.strip('\n').split(','))
            items_tuple = [tuple(x) for x in items]
        csv_open.close()
        
        return items_tuple
    
    def NDVI(self, red_band, nir_band):
        
        np.seterr(all = "ignore")
        red_band = np.array(tf.imread(red_band)).astype(np.float32)
        nir_band = np.array(tf.imread(nir_band)).astype(np.float32)

        return np.divide((nir_band - red_band), (nir_band + red_band))

    def readMTL(self, input_folder, constant):
        MTL_file = []
        constant_value = 0
        for root, dirs, files in os.walk(input_folder):
            if 'LC' in root:
                for file in files:
                    if file.endswith("MTL.txt"):
                        MTL_file.append(file)
                        MTL_open = open(os.path.join(root, file), 'r')
                        for line in MTL_open:
                            if constant in line:
                                constant_value = float(line.strip().split('=')[1])
                        MTL_open.close()
                        return constant_value
    
    def radTOA(self, Ml, thermal_band, Al):
        thermal_band = np.array(tf.imread(thermal_band)).astype(np.float32)

        return  (Ml * thermal_band) + Al

    def reflTOA(self, Ap, Mp, band):
        band = np.array(tf.imread(band)).astype(np.float32)
        return (Mp * band) + Ap

    def brightTemp(self, TOA, K1, K2):
        np.seterr(all = "ignore")
        return K2/np.log(((K1/TOA+1))) - 273.15
                
    def vegCover(self, ndvi): #brom
        #np.seterr(all = "ignore")
        return np.divide(np.power(ndvi, 2), 0.3) 

    def emis(self, ndvi, vegcover, red_band):
        red_band = np.array(tf.imread(red_band)).astype(np.float32)

        result_e = (0.004 * vegcover) + 0.986
        result_e = np.where(ndvi < 0.2, 1 - red_band, result_e) 
        result_e = np.where(ndvi > 0.5, 0.99, result_e)
        result_e[result_e > 1] = 0.99
        result_e[result_e < 0.8] = 0.8
        return result_e

    def LST(self, BT, emis):
        return (BT / (1 + ((0.0015 * BT)/1.4488) * np.log(emis)))

    def Albedo(self, refTOA1, reflTOA3, reflTOA4, reflTOA5, reflTOA7):
        return (((0.356 * refTOA1) + (0.130 * reflTOA3) + (0.373 * reflTOA4) + (0.085 * reflTOA5) + (0.072 * reflTOA7))- 0.0018) / 1.016

    def radLong(self, temp, emis): #for LongOut -> LST, emis_surf ; for LonngIn -> emis_atm, Air temp in C
        return emis * 5.67 * 10 ** (-8.0) * (temp + 273.15) ** 4

    def radShortOut(self, albedo, radShortIn):
        return albedo * radShortIn
   
    def netRad(self, radShortIn, radShortOut, radLongIn, radLongout):
        return (radShortIn - radShortOut) + (radLongIn - radLongout)

    def groundHeatFlux(self, albedo, lst, ndvi, netRadiation):
        return (lst /albedo ) * ((0.0038 * albedo) + 0.0074 * albedo)**2 * (1 - (0.98 * ndvi ** 4)) * netRadiation
    
    def zonalStatistics(self, input_zones, input_tif, name):
        alg_params = {
            'COLUMN_PREFIX': 'ZS_',
            'INPUT': input_zones,
            'INPUT_RASTER': input_tif,
            'RASTER_BAND': 1,
            'STATISTICS': [0,1,2,3,5,6,7,8,9],
            'OUTPUT': os.path.join(self.parameters['OUTPUT'], name + '_ZonalStat.gpkg')
        }
        processing.run('native:zonalstatisticsfb', alg_params, is_child_algorithm=True)
        return
    
    # Run zonalStatistcs for images in folder
    def RunZonalStatOnTif(self, input_folder, output_suffix):
        for root, dirs, files in os.walk(input_folder):
            for file in files:
                if file.endswith('.TIF'):
                    self.zonalStatistics(self.parameters['STAT_ZONES'], os.path.join(input_folder, file), file.split('_')[3] + output_suffix)
        return


    def processAlgorithm(self, parameters, context, feedback):
        self.parameters = parameters
        self.context = context
        self.feedback = feedback
    
        result_folders = ['clipped_bands','NDVI', 'LandSurfaceTemperature', 'Albedo', 'GroundHeatFlux']
        
        self.clipped_bands_folder = os.path.join(self.parameters['OUTPUT'], result_folders[0]) #create folder for clippeed tiffs
        #↨self.TOA_folder = os.path.join(self.parameters['OUTPUT'], result_folders[1])
        self.NDVI_folder = os.path.join(self.parameters['OUTPUT'], result_folders[1]) #create folder for NDVI
        #self.BT_folder = os.path.join(self.parameters['OUTPUT'], result_folders[3])
        #self.VC_folder = os.path.join(self.parameters['OUTPUT'], result_folders[4])
        #self.E_folder = os.path.join(self.parameters['OUTPUT'], result_folders[5])
        self.LST_folder = os.path.join(self.parameters['OUTPUT'], result_folders[2])
        self.Albedo_folder = os.path.join(self.parameters['OUTPUT'], result_folders[3])
        #self.Radiation_folder = os.path.join(self.parameters['OUTPUT'], result_folders[8])
        self.G_folder = os.path.join(self.parameters['OUTPUT'], result_folders[4])

        self.get_original_name = self.getImageName(self.parameters['INPUT_FOLDER']) #get complete path to original tiff in original folder
        #'root\\input_folder\\LC08_L1TP_190025_20170731_20200903_02_T1\\LC08_L1TP_190025_20170731_20200903_02_T1_B9.TIF'
        
        self.read_temp = self.readCSV(self.parameters['AIR_TEMP'])
        
        self.RADIANCE_MULT_BAND_10 = self.readMTL(self.parameters['INPUT_FOLDER'], 'RADIANCE_MULT_BAND_10') 
        self.RADIANCE_ADD_BAND_10 = self.readMTL(self.parameters['INPUT_FOLDER'], 'RADIANCE_ADD_BAND_10')

        self.REFLECTANCE_MULT_BAND_X = self.readMTL(self.parameters['INPUT_FOLDER'], 'REFLECTANCE_MULT_BAND_1') 
        self.REFLECTANCE_ADD_BAND_X = self.readMTL(self.parameters['INPUT_FOLDER'], 'REFLECTANCE_ADD_BAND_1') 
        
        self.KELVIN_CONSTANT_1 = self.readMTL(self.parameters['INPUT_FOLDER'], 'K1_CONSTANT_BAND_10')
        self.KELVIN_CONSTANT_2 = self.readMTL(self.parameters['INPUT_FOLDER'], 'K2_CONSTANT_BAND_10')

        for folder in result_folders:
            if not os.path.exists(os.path.join(self.parameters['OUTPUT'], folder)):
                os.mkdir((os.path.join(self.parameters['OUTPUT'], folder)))
                
        for input_band in self.get_original_name: 
            #'root\\LC08_L1TP_190025_20170731_20200903_02_T1\\LC08_L1TP_190025_20170731_20200903_02_T1_B8.TIF'
            #for each tier to original tiff
           
            self.image_name = os.path.basename(input_band[0]) 
            #'LC08_L1TP_190025_20170528_20200903_02_T1_B1.TIF'
            
            self.sensing_date = self.image_name.split('_')[3] 
            #extract the sensing date from the tiff name
            
            for item in self.read_temp:
                if item[0] == self.sensing_date:
                    air_temp = float(item[1])
            
            #self.read_temp = self.readCSV(self.parameters['AIR_TEMP'])
            self.calcRadLongIn = self.radLong(0.8, air_temp) #emis_atm, Ta_C
            
                       
            self.root_ouput_path = os.path.join(self.clipped_bands_folder, self.sensing_date) 
            
            #path to clipped tiffs + sensing date
            
            if not os.path.exists(self.root_ouput_path): 
                    os.makedirs(self.root_ouput_path)
            #create path to clipped tiffs if it does not exist        
            
            for file in input_band: #for one path in list
                #root\\LC08_L1TP_190025_20170731_20200903_02_T1\\LC08_L1TP_190025_20170731_20200903_02_T1_B9.TIF
                
                self.image_name = os.path.basename(file) 
                #get base name to each day in clipped tiffs folder
                #'root\\output_folder\\clipped_bands\\20170528'
                
                self.clipped_output_path = os.path.join(self.root_ouput_path, self.image_name.replace('.TIF', '_clipped.TIF')) 
                #create output path for clipped tiffs
                #'root\\output_folder\\\clipped_bands\\20170731\\LC08_L1TP_190025_20170731_20200903_02_T1_B1_clipped.TIF'

                self.clipLandsatBand(file, self.clipped_output_path, self.parameters['MASK_LAYER'])
                # clip the tiffs woth mask and save them into output folder

        self.clipped_paths_list =  self.getImageName(self.clipped_bands_folder)
        # make a list from paths to clipped tiffs

        for clipped_path in self.clipped_paths_list:
            #'root\\output_folder\\clipped_bands\\20170731\\LC08_L1TP_190025_20170731_20200903_02_T1_B9_clipped.TIF'
            
            self.red_band = ''
            self.nir_band = ''
            self.thermal_band = ''
            for band in clipped_path:
               
                if 'B1' in band:
                    self.aerosol_band = band
                if 'B3' in band:
                    self.blue_band = band
                if 'B4' in band:
                    self.red_band = band
                if 'B5' in band:
                    self.nir_band = band
                if 'B7' in band:
                    self.swir2_band = band
                if 'B10' in band:
                    self.thermal_band = band
               

            self.calcReflTOA_BX = self.reflTOA(self.REFLECTANCE_ADD_BAND_X, self.REFLECTANCE_MULT_BAND_X, self.aerosol_band)
            #self.rasterToFile(self.calcReflTOA_BX, self.thermal_band, os.path.join(self.TOA_folder, os.path.basename(self.thermal_band).replace('B10_clipped.TIF', 'BX_REFL_TOA.TIF')))
            self.calcRadTOA = self.radTOA(self.RADIANCE_MULT_BAND_10, self.thermal_band, self.RADIANCE_ADD_BAND_10)
            #self.rasterToFile(self.calcRadTOA, self.thermal_band, os.path.join(self.TOA_folder, os.path.basename(self.thermal_band).replace('B10_clipped.TIF', 'TOA.TIF')))

            try:
                self.calcAlbedo = self.Albedo(self.calcReflTOA_BX, self.calcReflTOA_BX, self.calcReflTOA_BX, self.calcReflTOA_BX, self.calcReflTOA_BX)
                self.rasterToFile(self.calcAlbedo, self.nir_band, os.path.join(self.Albedo_folder, os.path.basename(self.nir_band).replace('B5_clipped.TIF', 'Albedo.TIF')))
                self.feedback.pushInfo(str('[SUCCESS] Albedo was generated'))
            except:
                self.feedback.pushInfo(str('[FAIL] Albedo was not generated'))
            
            try:
                self.calcNDVI = self.NDVI(self.red_band, self.nir_band)
                self.rasterToFile(self.calcNDVI, self.nir_band, os.path.join(self.NDVI_folder, os.path.basename(self.nir_band).replace('B5_clipped.TIF', 'NDVI.TIF')))
                self.feedback.pushInfo(str('[SUCCESS] NDVI was generated'))
            except:
                self.feedback.pushInfo(str('[FAIL] NDVI was generated'))


            self.calcBT = self.brightTemp(self.calcRadTOA, self.KELVIN_CONSTANT_1, self.KELVIN_CONSTANT_2)
            #self.rasterToFile(self.calcBT, self.nir_band, os.path.join(self.BT_folder, os.path.basename(self.nir_band).replace('B5_clipped.TIF', 'BT.TIF')))
            
            self.calcVegCover= self.vegCover(self.calcNDVI)
            #self.rasterToFile(self.calcVegCover, self.nir_band, os.path.join(self.VC_folder, os.path.basename(self.nir_band).replace('B5_clipped.TIF', 'VEG.TIF')))
            
            self.calcEmis = self.emis(self.calcNDVI, self.calcVegCover, self.red_band)
            #self.rasterToFile(self.calcEmis, self.nir_band, os.path.join(self.E_folder, os.path.basename(self.nir_band).replace('B5_clipped.TIF', 'EMIS.TIF')))
            
            try:
                self.calcLST = self.LST(self.calcBT, self.calcEmis)
                self.rasterToFile(self.calcLST, self.nir_band, os.path.join(self.LST_folder, os.path.basename(self.nir_band).replace('B5_clipped.TIF', 'LST.TIF')))
                self.feedback.pushInfo(str('[SUCCESS] Land Surface Temperature was generated'))
            except:
                self.feedback.pushInfo(str('[FAIL] Land Surface Temperature was generated'))

            
            self.calcRadLongOut = self.radLong(self.calcLST, self.calcEmis)
            #self.rasterToFile(self.calcRadLongOut, self.nir_band, os.path.join(self.Radiation_folder, os.path.basename(self.nir_band).replace('B5_clipped.TIF', 'RadLongOut.TIF')))
            
            #Ollila, Antero. (2013). Earth's energy balance for clear, cloudy and all-sky conditions. Development in Earth Science. 1. 
            self.caclRadShortIn = 248.9

            self.calcRadShortOut = self.radShortOut(self.calcAlbedo, self.caclRadShortIn)
            #self.rasterToFile(self.calcRadShortOut, self.nir_band, os.path.join(self.Radiation_folder, os.path.basename(self.nir_band).replace('B5_clipped.TIF', 'RadShortOut.TIF')))

            self.calcNetRadiation = self.netRad(self.caclRadShortIn, self.calcRadShortOut, self.calcRadLongIn, self.calcRadLongOut)
            #self.rasterToFile(self.calcNetRadiation, self.nir_band, os.path.join(self.Radiation_folder, os.path.basename(self.nir_band).replace('B5_clipped.TIF', 'NetRadiation.TIF')))

            try:
                self.calcGroundHeatFlux = self.groundHeatFlux(self.calcAlbedo, self.calcLST, self.calcNDVI, self.calcNetRadiation)
                self.rasterToFile(self.calcGroundHeatFlux, self.nir_band, os.path.join(self.G_folder, os.path.basename(self.nir_band).replace('B5_clipped.TIF', 'G.TIF')))
                self.feedback.pushInfo(str('[SUCCESS] Ground Heat Flux was generated'))
            except:
                self.feedback.pushInfo(str('[FAIL] Ground Heat Flux was generated'))


        if self.parameters['calc_zonal_stat']:
            self.RunZonalStatOnTif(self.Albedo_folder, '_Albedo')
            self.RunZonalStatOnTif(self.NDVI_folder, '_NDVI')
            self.RunZonalStatOnTif(self.LST_folder, '_LST')
            self.RunZonalStatOnTif(self.G_folder, '_GroundHeatFlux')
            self.feedback.pushInfo(str('[SUCCESS] ZONAL STATISTICS WAS CALCULATED'))
        
        
        self.feedback.pushInfo(str('[SUCCESS]'))
        
        return {}
