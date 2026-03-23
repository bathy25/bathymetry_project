/ ==========================================================================
// Bathymetry Layer Creation for Dublin using Lyzenga 1978 using Sentinel-2
// ==========================================================================

// ---------------------- PARAMETERS -------------------------------------
var split = 0.7;
var depthName = 'Depth';
var CLOUD_THRES = 100;
var description = 'Rosses_bathy_lyzenga';
var assetId = 'Rosses_bathy_lyzenga';
var roi = landmask.geometry();

// Map settings
Map.centerObject(roi, 12);
Map.addLayer(roi, {color:"red"}, 'Rosses ROI', 0);

// ---------------------- DATE RANGES ----------------------
var rossesDates = [
  {start: '2022-07-01', end: '2022-09-30'},
];

// ---------------------- HELPER FUNCTIONS ----------------------
function maskRosses(image) {
  return image.updateMask(landmask).clip(roi);
}

// ---------------------- LOAD BATHYMETRY DATA ----------------------
var bathymetry = ee.FeatureCollection(bathy_points)
  .select(depthName)
  .filterBounds(roi);

var bathyVis = {min:-15, max:6, palette: ['darkred','red','orange','green','darkgreen']};
Map.addLayer(bathymetry, bathyVis, 'Bathymetry Heights'); 

// Split bathymetry points for training/validation
var bathyPoints = bathymetry.randomColumn('random');
var bathyTrain = bathyPoints.filter(ee.Filter.lt('random', split));
var bathyValid = bathyPoints.filter(ee.Filter.gte('random', split));

// ---------------------- REGRESSION FUNCTION ----------------------
function regressImage(image) {
  var bands = ['B2','B3','B4'];
  var logBands = bands.map(function(b){ return 'log_' + b; });
  
  image = image.addBands(image.select(bands).log().rename(logBands));
  image = image.addBands(ee.Image.constant(1).rename('constant'));
  
  var predictors = ee.List(logBands).add('constant');
  
  var training = image.select(predictors).sampleRegions({
    collection: bathyTrain,
    properties: [depthName],
    scale: 10
  });

  var linearRegression = ee.Dictionary(training.reduceColumns({
    reducer: ee.Reducer.linearRegression({numX: 4, numY: 1}),
    selectors: ['constant','log_B2','log_B3','log_B4',depthName]
  }));

  var coefList = ee.Array(linearRegression.get('coefficients')).toList();
  
  var predicted = image.select(['log_B2']).multiply(ee.Number(coefList.get(1)))
    .add(image.select(['log_B3']).multiply(ee.Number(coefList.get(2))))
    .add(image.select(['log_B4']).multiply(ee.Number(coefList.get(3))))
    .add(ee.Number(coefList.get(0))).rename('predictedDepth');

  return image.addBands(predicted).set('system:time_start', image.get('system:time_start'));
}

// ---------------------- RMSE FUNCTION ----------------------
function calculateRMSE(image) {
  var samples = image.select(['predictedDepth']).sampleRegions({
    collection: bathyValid,
    properties: [depthName],
    scale: 10
  });
  
  var rmse = ee.Number(samples.map(function(sample){
    var error = ee.Number(sample.get(depthName)).subtract(ee.Number(sample.get('predictedDepth')));
    return sample.set('squaredError', error.pow(2));
  }).reduceColumns({reducer: ee.Reducer.mean(), selectors: ['squaredError']}).get('mean')).sqrt();
  
  return image.set('rmse', rmse);
}

// ---------------------- R² FUNCTION ----------------------
function calculateRSquared(image) {
  var samples = image.select(['predictedDepth']).sampleRegions({
    collection: bathyValid,
    properties: [depthName],
    scale: 10
  });
  
  var meanDepth = samples.reduceColumns(ee.Reducer.mean(), [depthName]).get('mean');
  
  var ssTotal = samples.map(function(sample){
    var depth = ee.Number(sample.get(depthName));
    return ee.Feature(null, {'squaredTotal': depth.subtract(meanDepth).pow(2)});
  }).reduceColumns(ee.Reducer.sum(), ['squaredTotal']).get('sum');
  
  var ssResidual = samples.map(function(sample){
    var depth = ee.Number(sample.get(depthName));
    var predictedDepth = ee.Number(sample.get('predictedDepth'));
    return ee.Feature(null, {'squaredResidual': depth.subtract(predictedDepth).pow(2)});
  }).reduceColumns(ee.Reducer.sum(), ['squaredResidual']).get('sum');
  
  var rSquared = ee.Number(1).subtract(ee.Number(ssResidual).divide(ssTotal));
  return image.set('rSquared', rSquared);
}

// ---------------------- PROCESS IMAGES AND EXPORT ----------------------
rossesDates.forEach(function(dateRange){
  
  var images = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(POI)
    .filterDate(dateRange.start, dateRange.end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_THRES))
    .map(maskRosses);
  
  var predictedCollection = images.map(regressImage);
  var rmseResults = predictedCollection.map(calculateRMSE);
  var rsquaredCollection = predictedCollection.map(calculateRSquared);

  // Export RMSE and R² CSVs
  Export.table.toDrive({
    collection: rmseResults,
    description: 'rmse_export_' + dateRange.start,
    fileFormat: 'CSV'
  });
  Export.table.toDrive({
    collection: rsquaredCollection,
    description: 'r2_export_' + dateRange.start,
    fileFormat: 'CSV'
  });

  // Export Sentinel-2 images using server-side safe loop
  var imageList = images.toList(images.size());
  var nImages = imageList.size().getInfo(); // number of images

  for (var i = 0; i < nImages; i++) {
    var img = ee.Image(imageList.get(i));
    var dateStr = ee.Date(img.get('system:time_start')).format('YYYYMMdd').getInfo();

    Export.image.toDrive({
      image: img.select(['B2','B3','B4']),
      description: 'Rosses_S2_' + dateStr,
      folder: 'Rosses_Sentinel2',
      fileNamePrefix: 'Rosses_S2_' + dateStr,
      region: roi,
      scale: 10,
      maxPixels: 1e13
    });
  }
});

print('Processing complete.');

