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
from qgis.core import (QgsProcessing,
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterFolderDestination,
                       QgsProcessingParameterFile,
                       QgsProcessingParameterVectorLayer)
from qgis import processing
from osgeo import gdal, osr
import os

try:
    import numpy as np
    import tifffile as tf

except:
    os.system('python -m pip install numpy --user')
    os.system('python -m pip install tifffile --user')
    import numpy as np
    import tifffile as tf
    
   
class ExampleProcessingAlgorithm(QgsProcessingAlgorithm):
    INPUT_FOLDER = 'INPUT_FOLDER'
    MASK_LAYER = 'MASK_LAYER'
    OUTPUT = 'OUTPUT'

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return ExampleProcessingAlgorithm()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return '2_vyvoj_pro_foss_cas_B'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('2_vyvoj_pro_foss_cas_B')

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr('Example scripts')

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'examplescripts'

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it..
        """
        return self.tr("Example algorithm short description")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile(self.INPUT_FOLDER, self.tr('Folder containing directories with Landsat images'), behavior=QgsProcessingParameterFile.Folder))
        self.addParameter(QgsProcessingParameterVectorLayer(self.MASK_LAYER, self.tr('Polygon clipping mask'), defaultValue=None))
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT, self.tr('Output Folder')))
    
      
    def NDVI(self, red_band, nir_band):
        self.red_band = red_band
        self.nir_band = nir_band
        
        np.seterr(all = "ignore")
        red_band = np.array(tf.imread(red_band)).astype(np.float32)
        nir_band = np.array(tf.imread(nir_band)).astype(np.float32)
        
        return (nir_band - red_band) / (nir_band + red_band)
      
      
    def rasterToFile(self, input_array, src_dataset_path, output_path): #https://gist.github.com/jkatagi/a1207eee32463efd06fb57676dcf86c8
        self.inpu_array = input_array
        self.src_dataset_path = src_dataset_path
        self.output_path = output_path
        
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
        
    def findPathToImage(self, input_folder):
        self.path_list_folder = []
        self.input_folder = input_folder
        self.landsat_name_dict = {}
        
        for root, dirs, files in os.walk(input_folder, topdown=False):
            for name in files:
                if name.endswith('.TIF') and name.startswith('LC') :
                    self.path_list_folder.append(os.path.join(root, name))
        return self.path_list_folder
        
        """for landsat_path in self.path_list_folder:
            self.metadata = {}
            self.landsat_name = os.path.basename(landsat_path)
            
            if 'LC' in self.landsat_name:
                date = self.landsat_name.split('_')[3]
                if date in self.landsat_name_dict:
                    self.landsat_name_dict[date].append(landsat_path)
                else:
                    self.landsat_name_dict[date] = [landsat_path]
        return self.landsat_name_dict"""
    
        
    def clipLandsatBand(self, input_band, out_path ,mask):
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
        
        self.feedback.pushInfo(str(self.out_path))
        processing.run('gdal:cliprasterbymasklayer', alg_params, is_child_algorithm=True)
        return
        
    def getLandsatIdentifier(self, input_folder):
        self.input_folder = input_folder
        self.path_list_metadata = []
        self.landsat_metadata = []
       
        for root, dirs, files in os.walk(input_folder, topdown=False):
            for name in files:
                if name.endswith('MTL.txt') and name.startswith('LC') :
                    self.path_list_metadata.append(os.path.join(root, name))
                    
                    
        for path in self.path_list_metadata:
            self.metadata_file = open(path, 'r')
            
            for line in self.metadata_file:
                if 'LANDSAT_PRODUCT_ID ' in line:
                    line_split = line.strip().split(' = ')[1].split('_')[3] #"LC08_L1TP_190025_20170528_20200903_02_T1"
                    #self.feedback.pushInfo(str((line_split)))
                    self.landsat_metadata.append(line_split)
                    break
                    #self.feedback.pushInfo(str((self.landsat_metadata)))
            
        return self.landsat_metadata
        
    def processAlgorithm(self, parameters, context, feedback):
        self.parameters = parameters
        self.context = context
        self.feedback = feedback
        
        self.clipped_bands_folder = os.path.join(self.parameters['OUTPUT'], 'clipped_bands')
        
        self.find_paths_to_images = self.findPathToImage(self.parameters['INPUT_FOLDER'])
        self.list_of_dates = self.getLandsatIdentifier(self.parameters['INPUT_FOLDER'])
        
        for identifier in self.list_of_dates:
            #self.feedback.pushInfo(str(identifier))
            os.makedirs(os.path.join(self.clipped_bands_folder, identifier))
            
            self.clipLandsatBand(self.find_paths_to_images, os.path.join(self.clipped_bands_folder, identifier), self.parameters['MASK_LAYER'])
        
        return {}