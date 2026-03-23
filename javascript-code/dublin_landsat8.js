// ========================================================================
// Bathymetry Layer Creation for Dublin using Lyzenga 1978b using Landsat-8
// =========================================================================

// ---------------------- PARAMETERS ----------------------
var split = 0.7;
var depthName = 'Depth';
var CLOUD_THRES = 100;
var roi = landmask.geometry();

// Map settings
Map.centerObject(roi, 12);
Map.addLayer(roi, {color:"red"}, 'Rosses ROI', 0);

// ---------------------- DATE RANGE ----------------------
var startDate = '2021-04-01';
var endDate = '2021-11-30';

// ---------------------- HELPER FUNCTIONS ----------------------
function maskClouds(image) {
  var qa = image.select('QA_PIXEL');
  // Mask clouds (bit 3), cloud shadows (bit 4), snow (bit 5)
  var cloudMask = qa.bitwiseAnd(1 << 3).eq(0)
    .and(qa.bitwiseAnd(1 << 4).eq(0))
    .and(qa.bitwiseAnd(1 << 5).eq(0));
  return image.updateMask(cloudMask).updateMask(landmask).clip(roi);
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
function safeRegress(image) {
  var bands = ['SR_B2','SR_B3','SR_B4'];
  var logBands = bands.map(function(b){ return 'log_' + b; });
  var predictors = ee.List(logBands).add('constant');

  image = image.addBands(image.select(bands).log().rename(logBands));
  image = image.addBands(ee.Image.constant(1).rename('constant'));

  var training = image.select(predictors).sampleRegions({
    collection: bathyTrain,
    properties: [depthName],
    scale: 30
  });

  var count = training.size();
  return ee.Algorithms.If(count.gt(0), (function() {
    var linearRegression = ee.Dictionary(training.reduceColumns({
      reducer: ee.Reducer.linearRegression({numX: 4, numY: 1}),
      selectors: ['constant','log_SR_B2','log_SR_B3','log_SR_B4',depthName]
    }));

    var coefList = ee.Array(linearRegression.get('coefficients')).toList();

    var predicted = image.select(['log_SR_B2']).multiply(ee.Number(coefList.get(1)))
      .add(image.select(['log_SR_B3']).multiply(ee.Number(coefList.get(2))))
      .add(image.select(['log_SR_B4']).multiply(ee.Number(coefList.get(3))))
      .add(ee.Number(coefList.get(0))).rename('predictedDepth');

    return image.addBands(predicted).set('system:time_start', image.get('system:time_start'));
  })(), ee.Image(0).rename('dummy')); // return dummy image if no training
}

// ---------------------- RMSE FUNCTION ----------------------
function calculateRMSE(image) {
  var samples = image.select(['predictedDepth']).sampleRegions({
    collection: bathyValid,
    properties: [depthName],
    scale: 30
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
    scale: 30
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
var images = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .map(maskClouds);

print('\uD83D\uDDD3️ Date range:');
print('\uD83D\uDCF8 Images found:');

var imageList = images.toList(images.size());
var nImages = imageList.size().getInfo();
print('\uD83D\uDD22 Number of images:');

var metrics = [];

// client-side loop for exports + metrics collection
for (var i = 0; i < nImages; i++) {
  var img = ee.Image(imageList.get(i));
  var dateStr = ee.Date(img.get('system:time_start')).format('YYYYMMdd').getInfo();

  var regressed = ee.Image(safeRegress(img));

  // Only process if predictedDepth exists
  var hasBand = regressed.bandNames().contains('predictedDepth');
  if (hasBand.getInfo()) {
    var withRMSE = calculateRMSE(regressed);
    var withRSquared = calculateRSquared(withRMSE);

    // Export raw RGBs
    Export.image.toDrive({
      image: img.select(['SR_B2','SR_B3','SR_B4']),
      description: 'Dublin_L8_' + dateStr,
      folder: 'Dublin_Landsat8',
      fileNamePrefix: 'Dublin_L8_' + dateStr,
      region: roi,
      scale: 30,
      maxPixels: 1e13
    });

    // Store metrics
    metrics.push(ee.Feature(null, {
      'date': dateStr,
      'rmse': withRSquared.get('rmse'),
      'rSquared': withRSquared.get('rSquared')
    }));
  } else {
    print('Skipping image ' + dateStr + ' (no training points).');
  }
}


// Convert metrics to FeatureCollection
var allMetrics = ee.FeatureCollection(metrics);

// Export metrics CSV
Export.table.toDrive({
  collection: allMetrics,
  description: 'Dublin_L8_metrics_all',
  fileFormat: 'CSV'
});

print('✅ Landsat 8 bathymetry processing complete.');
