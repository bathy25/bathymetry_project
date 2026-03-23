// Create land mask data for Dublin

Requirements (Imports needed above)
  'POI'            = Point of interest
  'boundary'       = Feature polygon - Draw below of area of interest to discard inland water

Parameters:                                                   Defaults
  'START_DATE'            = yyyy-mm-dd                          2021-01-01
  'END_DATE'              = yyyy-mm-dd                          2021-09-20
  'CLOUD_THRES'           = threshold to mask clouds            40  (40% cloud probability)
  'NIR_thres'             = Value for near-infrared threshold   1000
                            for land / water 
Output:
  'description'           = description of output raster       'UCD_workshop_Portrane_landmask',
  'assetId'               = id of output raster                'UCD_workshop_Portrane_landmask',
*/

var START_DATE = '2021-01-01';  
var END_DATE = '2021-09-20';       
var NIR_thres = 1000 
var CLOUD_THRES = 20

// Outputs
var description = 'UCD_workshop_Portrane_landmask'
var assetId = 'UCD_workshop_Portrane_landmask'

// Create Sentinel 2 composite //
// Load Sentinel 2 raw imagery and filter it to all images in 2021.
var S2data = ee.ImageCollection('COPERNICUS/S2_SR')
    .filterDate(ee.Date(START_DATE), ee.Date(END_DATE))
    .filterBounds(POI)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_THRES))

print(S2data)

// Print the number of images in the collection
print('Number of Sentinel 2 images in median', S2data.size())

// Take median
var S2composite = S2data.median()

//  ------------------------------------------------------------------   //
// Create Land mask //
// Select near-infrared band
var b8 = S2composite.select('B8')

// Use threshold to define land and water boundary
var landmask = b8.where(b8.gt(NIR_thres),0).where(b8.lt(NIR_thres),1);

//  ------------------------------------------------------------------   //
// Add Layers to Map //
// Visualisation Parameters
var rgbVis = {min: 0, max: 3000, bands: ['B8']};

Map.addLayer(S2composite, rgbVis, 'S2 NIR Band')
Map.addLayer(landmask.clip(boundary), {bands: ['B8'], gamma: 1, max: 1, opacity: 1},'S2 Land Mask')

//  ------------------------------------------------------------------   //
// Export the image to an Earth Engine asset //
Export.image.toAsset({
  image: landmask,
  region: boundary,
  description: description,
  assetId: assetId,
  scale: 10,
  maxPixels: 551733655,
});
