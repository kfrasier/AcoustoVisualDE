from GeoEco.DataManagement.HDFs import HDF


HDF.SDSToArcGISRaster('C:\\temp14\\hdf\\199001.s04m1pfv50-sst-16b.hdf',
'C:\\temp14\\Pathfinder_Example_Raster\\sst199001',
 'sst', -180, -90, 0.0439453125, 0,
 coordinateSystem="GEOGCS['GCS_WGS_1984',
 DATUM['D_WGS_1984',SPHEROID['WGS_1984',
 6378137.0,298.257223563]],
 PRIMEM['Greenwich',0.0],
 UNIT['Degree',0.0174532925199433]]",
 mapAlgebraExpression='float(inputRaster) * 0.075 - 3.0',
  buildPyramids=True, overwriteExisting=True)
