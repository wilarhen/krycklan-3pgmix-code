// ---- USER INPUT: EXTENT ----
var xmin = 19.44167;
var xmax = 19.975;
var ymin = 64.10833;
var ymax = 64.33333;

// Create rectangle geometry from extent
var aoi = ee.Geometry.Rectangle([xmin, ymin, xmax, ymax]);
Map.centerObject(aoi, 10);
Map.addLayer(aoi, {color: 'red'}, 'AOI');

// ---- TIME RANGE ----
var startYear = 2000;
var endYear   = 2020;

// ---- LOAD MODIS MCD18A1.062 ----
var dataset = ee.ImageCollection('MODIS/062/MCD18A1')
                .filterBounds(aoi)
                .map(function(img) {
                  return img.clip(aoi);
                });

// ---- 3-hour DSR bands ----
var bands = [
  'GMT_0000_DSR','GMT_0300_DSR','GMT_0600_DSR','GMT_0900_DSR',
  'GMT_1200_DSR','GMT_1500_DSR','GMT_1800_DSR','GMT_2100_DSR'
];

// ---- Compute daily total DSR + valid mask ----
var addDailyDSR = function(img) {
  // Mask: keep only good quality pixels
  var qualityMask = img.select('DSR_Quality').eq(0);
  var masked = img.updateMask(qualityMask);

  // Daily radiation from 3-hour windows (MJ m-2 d-1)
  var totalWm2 = masked.select(bands).reduce(ee.Reducer.sum());
  var dailyMJ = totalWm2.multiply(10800).divide(1e6)
                        .rename('DSR_MJ_m2_day');

  // Add valid pixel flag (1 for valid day, 0 for invalid)
  var validDay = qualityMask.rename('valid_day');

  return dailyMJ.addBands(validDay)
                .copyProperties(img, ['system:time_start']);
};

var dailyDSR = dataset.map(addDailyDSR);

// ---- Loop: Monthly statistics ----
var years = ee.List.sequence(startYear, endYear);
var months = ee.List.sequence(1, 12);

years.getInfo().forEach(function(y) {
  months.getInfo().forEach(function(m) {

    var startMonth = ee.Date.fromYMD(y, m, 1);
    var endMonth = startMonth.advance(1, 'month');

    var monthly = dailyDSR.filterDate(startMonth, endMonth);

    // Monthly mean DSR
    var monthlyMean = monthly.select('DSR_MJ_m2_day').mean();

    // Count valid days per pixel
    var validCount = monthly.select('valid_day').sum().rename('valid_days');

    // Combine into one image
    var monthlyOut = monthlyMean.addBands(validCount);

    // ---- Export ----
    Export.image.toDrive({
      image: monthlyOut,
      description: 'DSR_' + y + '_' + (m < 10 ? '0' + m : m),
      folder: 'MODIS_DSR',
      region: aoi,
      scale: 1000,
      crs: 'EPSG:4326',
      maxPixels: 1e13
    });

  });
});
