#Marconi Python Code for Raster Data Manipulation 
#Used for binning raster data from lidar scans into classes for histogram generation
#Useful in counting house groups at archaeological sites, or in preliminary laboratory surveying of new field areas where remotely sensed products are available 

#initial setup for .las file, so data can be manipulated as raster
import arcpy


#Setting up local variables for conversion into raster format 
try:
    inLas = arcpy.GetParameterAsText(0)
    recursion = arcpy.GetParameterAsText(1)
    surfCons = arcpy.GetParameterAsText(2)
    classCode = arcpy.GetParameterAsText(3)
    returnValue = arcpy.GetParameterAsText(4)
    spatialRef = arcpy.GetParameterAsText(5)
    lasD = arcpy.GetParameterAsText(6)
    outRaster = arcpy.GetParameterAsText(7)
    cellSize = arcpy.GetParameter(8)
    zFactor = arcpy.GetParameter(9)

#CreateLasDataset
    arcpy.management.CreateLasDataset(inLas, lasD, recursion, surfCons, sr)
#Execute MakeLasDatasetLayer for JTREE .las file 
    lasLyr = arcpy.CreateUniqueName('J_Tree')
    arcpy.management.MakeLasDatasetLayer(lasD, lasLyr, classCode, returnValue)
#Execute LasDatasetToRaster to convert file 
    arcpy.conversion.LasDatasetToRaster(lasLyr, outRaster, 'ELEVATION',
                              'TRIANGULATION LINEAR WINDOW_SIZE 10', 'FLOAT',
                              'CELLSIZE', cellSize, zFactor)
    print(arcpy.GetMessages())

except arcpy.ExecuteError:
    print(arcpy.GetMessages())

except Exception as err:
    print(err.args[0])

finally:
    arcpy.management.Delete(lasLyr)

#import modules for statistical analyses, conversion of file formats and plotting raster/array data
import numpy as np
import gdal
import matplotlib.pyplot as plt

fileName = 'C:\Users\Ben\Documents\Utah PhD\Fall 2017\Geoprocessing with Python\Project\LiDAR\points.las'
dataSet = gdal.Open(fileName)


# Display the dataset's x and y dimensions, its number of bands, and geotransform info pertaining to location of the dataset's corners and coordinate system information
columns = dataSet.RasterXSize; print('Number of Columns:',columns)
rows = dataSet.RasterYSize; print('Number of Rows:',rows)
print('Numbers of Bands:',dataSet.RasterCount)
print('Driver:',dataSet.GetDriver().LongName)

print('Projection:',dataSet.GetProjection())

print('Geotransform:',dataSet.GetGeoTransform())


mapInformation = dataSet.GetGeoTransform()
x_Minimum = mapInformation[0]
y_Maximum = mapInformation[3]


#Dividing by the pixel width to prepare to convert into singluar raster band 
x_Maximum = x_Minimum + dataSet.RasterXSize/mapInformation[1]


#Dividing by pixel height, iusing the + sign (not - because of z length direction values) 
y_Minimum = y_Maximum + dataSet.RasterYSize/mapInformation[5]
CHM_Extent = (x_Minimum,x_Maximum,y_Minimum,y_Maximum)
print('J-Tree Raster Extent:',CHM_Extent)


CHM_Raster = dataSet.GetRasterBand(1)
noDataValue = CHM_Raster.GetnoDataValueue(); print('No Data Value:',noDataValue)
scaleFactor = CHM_Raster.GetScale(); print('scale factor:',scaleFactor)
CHM_Statistics = CHM_Raster.GetStatistics(True,True)
print('J-Tree Statistics: Minimum=%.2f, Maximum=%.2f, Mean=%.3f, StDev=%.3f' %
      (CHM_Statistics[0], CHM_Statistics[1], CHM_Statistics[2], CHM_Statistics[3]))


CHM_Array = dataSet.GetRasterBand(1).ReadAsArray(0,0,columns,rows).astype(np.float)
CHM_Array[CHM_Array==int(noDataValue)]=np.nan
CHM_Array = CHM_Array/scaleFactor

#Show values in array
print('J-Tree Array:\n',CHM_Array) 


#Display stats of raster data like min, max, mean; then numpy.nanmin calculates the minimum without the NaN values.
print('J-Tree Array Statistics:')
print('Minimum:',round(np.nanmin(CHM_Array),2))
print('Maximum:',round(np.nanmax(CHM_Array),2))
print('Mean:',round(np.nanmean(CHM_Array),2))


# Calculate the % of pixels that are NaN and non-zero, so that the histrogram output is not so severely right-skewed
nonzeroPixels = np.count_nonzero(np.isnan(CHM_Array))/(rows*columns)
print('% NaN:',round(nonzeroPixels*100,2))
print('% Non-Zero:',round(100*np.count_nonzero(CHM_Array)/(rows*columns),2))


#Define the plot_band_array function to graph data (see next comment with plt.hist function below) 
def plot_band_array(band_array,refl_extent,colorlimit,ax=plt.gca(),title='',cbar ='on',cmap_title='',colormap='spectral'):
    plot = plt.imshow(band_array,extent=refl_extent,clim=colorlimit);
    if cbar == 'on':
        cbar = plt.colorbar(plot,aspect=40); plt.set_cmap(colormap);
        cbar.set_label(cmap_title,rotation=90,labelpad=20);
    plt.title(title); ax = plt.gca();
    ax.ticklabel_format(useOffset=False, style='plain'); 
    plot_band_array(CHM_Array,CHM_Extent,(0,80),title='J-Tree Analysis',cmap_title='')

import copy

CHM_NoNaN_Array = copy.copy(CHM_Array)
CHM_NoNaN_Array = CHM_NoNaN_Array[~np.isnan(CHM_Array)]
plt.hist(CHM_NoNaN_Array,weights=np.zeros_like(CHM_NoNaN_Array)+1./
         (CHM_Array.shape[0]*CHM_Array.shape[1]),bins=50);
plt.title('Distribution of J-Tree Struture Heights')
plt.xlabel(''); plt.ylabel('Relative Frequency')


CHM_NoNaN_Array = copy.copy(CHM_Array)
CHM_NoNaN_Array[CHM_Array==0]=np.nan
CHM_NonZero_Array = CHM_NoNaN_Array[~np.isnan(CHM_NoNaN_Array)]


#Use weighting to generate new realtive frequency histogram without right skewedness
plt.hist(CHM_NonZero_Array,weights=np.zeros_like(CHM_NonZero_Array)+1./
         (CHM_Array.shape[0]*CHM_Array.shape[1]),bins=50);


plt.title('Distribution of J-Tree Non-Zero Statistics')
plt.xlabel(''); plt.ylabel('Relative Frequency')


print('min:',np.amin(CHM_NonZero_Array),'m')
print('max:',round(np.amax(CHM_NonZero_Array),2),'m')
print('mean:',round(np.mean(CHM_NonZero_Array),2),'m')


#Final Plot 
plot_band_array(CHM_Array,CHM_Extent,(0,45),title='J-Tree Analysis',cmap_title='',colormap='BuGn')

