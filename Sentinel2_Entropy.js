/**
 * Sentinel-2 Vegetation Indices and Shannon Entropy Time-Series Analysis
 *
 * Author: Alessandro Mateos Sierra
 * Platform: Google Earth Engine JavaScript Code Editor
 * Dataset: COPERNICUS/S2_SR_HARMONIZED
 *
 * Purpose:
 * Calculate nine vegetation indices, normalize them using collection-wide
 * minimum and maximum values, and calculate Shannon entropy from 64-bin
 * frequency histograms for each acquisition and index.
 *
 * Required input:
 * Define or import an ee.Geometry / FeatureCollection named `geometry`
 * in the Earth Engine Code Editor before running this script. The geometry
 * should represent the study-area mask used in the analysis.
 *
 * Methodological notes:
 * - The date range is 2023-06-01 through 2023-09-30 (end date 2023-10-01 is exclusive).
 * - Scenes are filtered to CLOUDY_PIXEL_PERCENTAGE < 30.
 * - SCL classes 4, 5, and 6 are retained at pixel level.
 * - Index values are normalized across the entire filtered collection.
 * - Entropy is calculated within Earth Engine from per-image histograms.
 * - Exports are CSV tables sent to Google Drive via the Tasks panel.
 * - GRVI is defined here as NIR/GREEN (as in the source script).
 *
 * Reproducibility: supply the original study-area geometry and record the
 * acquisition IDs, processing dates, and Earth Engine dataset version.
 * The script has not been executed against the Earth Engine service here.
 */

// =====================================================
// SENTINEL-2 VEGETATION INDICES + SHANNON ENTROPY
// SINGLE STUDY-AREA POLYGON: geometry
// =====================================================


// =====================================================
// 1. PARAMETERS
// =====================================================

// Replace these dates only when intentionally defining a different analysis period.
var startDate = '2023-06-01';
var endDate   = '2023-10-01'; // Exclusive: includes 30 September 2023

var cloudThreshold = 30;

// Number of bins used for Shannon entropy
var nBins = 64;

// Vegetation indices
var indices = [
  'CLgreen',
  'GNDVI',
  'GRVI',
  'MSR',
  'MTVI',
  'NDVI',
  'NRVI',
  'RVI',
  'SAVI'
];


// =====================================================
// 2. SENTINEL-2 SCENE CLASSIFICATION MASK
// =====================================================

function maskS2clouds(image) {

  var scl = image.select('SCL');

  // Retain the following SCL classes:
  // 4 = vegetation
  // 5 = bare soil
  // 6 = water

  var mask = scl.eq(4)
    .or(scl.eq(5))
    .or(scl.eq(6));

  return image.updateMask(mask);
}


// =====================================================
// 3. SENTINEL-2 IMAGE COLLECTION
// =====================================================

var s2_raw = ee.ImageCollection(
    'COPERNICUS/S2_SR_HARMONIZED'
  )
  .filterBounds(geometry)
  .filterDate(startDate, endDate)
  .filter(
    ee.Filter.lt(
      'CLOUDY_PIXEL_PERCENTAGE',
      cloudThreshold
    )
  )
  .map(maskS2clouds);


// =====================================================
// 4. VEGETATION INDEX CALCULATION
// =====================================================

function addIndices(img) {

  // Sentinel-2 surface reflectance bands are scaled by 10,000
  var scaled = img.divide(10000);

  var nir   = scaled.select('B8');
  var red   = scaled.select('B4');
  var green = scaled.select('B3');


  // -----------------------------
  // NDVI
  // -----------------------------

  var ndvi = nir
    .subtract(red)
    .divide(
      nir.add(red)
    )
    .rename('NDVI');


  // -----------------------------
  // GNDVI
  // -----------------------------

  var gndvi = nir
    .subtract(green)
    .divide(
      nir.add(green)
    )
    .rename('GNDVI');


  // -----------------------------
  // MSR
  // -----------------------------

  var ratioNR = nir.divide(red);

  var msr = ratioNR
    .subtract(1)
    .divide(
      ratioNR
        .add(1)
        .sqrt()
    )
    .rename('MSR');


  // -----------------------------
  // SAVI
  // L = 0.5
  // -----------------------------

  var L = 0.5;

  var savi = nir
    .subtract(red)
    .divide(
      nir
        .add(red)
        .add(L)
    )
    .multiply(1 + L)
    .rename('SAVI');


  // -----------------------------
  // Chlorophyll Index Green
  // -----------------------------

  var clgreen = nir
    .divide(green)
    .subtract(1)
    .rename('CLgreen');


  // -----------------------------
  // RVI
  // NIR / RED
  // -----------------------------

  var rvi = nir
    .divide(red)
    .rename('RVI');


  // -----------------------------
  // GRVI
  // Definition used in this analysis:
  // NIR / GREEN
  // -----------------------------

  var grvi = nir
    .divide(green)
    .rename('GRVI');


  // -----------------------------
  // NRVI
  // -----------------------------

  var nrvi = rvi
    .subtract(1)
    .divide(
      rvi.add(1)
    )
    .rename('NRVI');


  // -----------------------------
  // MTVI
  // -----------------------------

  var mtvi = nir
    .subtract(green)
    .multiply(1.2)
    .subtract(
      red
        .subtract(green)
        .multiply(2.5)
    )
    .multiply(1.2)
    .rename('MTVI');


  // Return an image containing only the vegetation indices
  return ee.Image.cat([
      ndvi,
      gndvi,
      msr,
      savi,
      clgreen,
      rvi,
      grvi,
      nrvi,
      mtvi
    ])
    .copyProperties(
      img,
      ['system:time_start']
    );
}


var s2_col = s2_raw.map(addIndices);

print(
  'Sentinel-2 collection with indices:',
  s2_col
);

print(
  'Number of Sentinel-2 images:',
  s2_col.size()
);


// =====================================================
// 5. GLOBAL MIN-MAX NORMALIZATION
// =====================================================
//
// IMPORTANT:
//
// Normalization uses the entire filtered image collection.
//
// For each vegetation index:
//
// x_norm = (x - min_global) /
//          (max_global - min_global)
//
// This applies the same 0–1 scale across all
// acquisition dates.
// =====================================================


var globalMin = s2_col
  .select(indices)
  .min()
  .reduceRegion({
    reducer: ee.Reducer.min(),
    geometry: geometry,
    scale: 10,
    maxPixels: 1e9,
    tileScale: 4
  });


var globalMax = s2_col
  .select(indices)
  .max()
  .reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: geometry,
    scale: 10,
    maxPixels: 1e9,
    tileScale: 4
  });


print(
  'Global minimum by index:',
  globalMin
);

print(
  'Global maximum by index:',
  globalMax
);


// Convert dictionaries to constant multiband images
var minImage = ee.Image.constant(
  indices.map(function(name) {
    return ee.Number(
      globalMin.get(name)
    );
  })
).rename(indices);


var maxImage = ee.Image.constant(
  indices.map(function(name) {
    return ee.Number(
      globalMax.get(name)
    );
  })
).rename(indices);


// Normalization
function normalizeImage(img) {

  var normalized = img
    .select(indices)
    .subtract(minImage)
    .divide(
      maxImage.subtract(minImage)
    )
    .clamp(0, 1);

  return normalized
    .copyProperties(
      img,
      ['system:time_start']
    );
}


var s2_col_norm = s2_col.map(
  normalizeImage
);


// =====================================================
// 6. VEGETATION INDEX TIME SERIES
// =====================================================

var chartIndices = ui.Chart.image.series({
  imageCollection: s2_col_norm,
  region: geometry,
  reducer: ee.Reducer.mean(),
  scale: 10
})
.setOptions({

  title:
    'Temporal series of normalized vegetation indices',

  lineWidth: 2,

  pointSize: 4,

  hAxis: {
    title: 'Date',
    format: 'dd-MM-yyyy',
    slantedText: true,
    slantedTextAngle: 45
  },

  vAxis: {
    title: 'Normalized index value',
    minValue: 0,
    maxValue: 1
  },

  legend: {
    position: 'right'
  }
});


print(chartIndices);


// =====================================================
// 7. SHANNON ENTROPY
// =====================================================
//
// Use globally normalized vegetation indices.
//
// All normalized indices use:
// minimum = 0
// maximum = 1
//
// The same bin boundaries are therefore used
// intervalos para TODAS acquisition dates.
//
// Example with 64 bins:
//
// 0.000000 - 0.015625
// 0.015625 - 0.031250
// ...
// 0.984375 - 1.000000
//
// This makes entropy values temporally comparable.
//
// H = -sum(p_i * log2(p_i))
// =====================================================


// =====================================================
// 8. ENTROPY CALCULATION FUNCTION
// =====================================================

function calculateEntropy(image, bandName) {

  bandName = ee.String(bandName);

  var band = image
    .select([bandName])
    .clamp(0, 1);


  // ---------------------------------------------
  // Assign normalized values to integer bins 0 to (nBins - 1)
  // ---------------------------------------------

  var classified = band
    .multiply(nBins)
    .floor()
    .min(nBins - 1)
    .toInt()
    .rename('bin');


  // ---------------------------------------------
  // Histogram of fixed-width bins
  // ---------------------------------------------

  var histResult = classified.reduceRegion({

    reducer: ee.Reducer.frequencyHistogram(),

    geometry: geometry,

    scale: 10,

    maxPixels: 1e9,

    tileScale: 4
  });


  // Check whether valid pixel data are available
  var hasData = histResult.contains('bin');


  var entropy = ee.Algorithms.If(

    hasData,

    (function() {

      var histogram = ee.Dictionary(
        histResult.get('bin')
      );


      // Counts
      var counts = ee.List(
        histogram.values()
      );


      // Total number of valid pixels
      var total = ee.Number(
        counts.reduce(
          ee.Reducer.sum()
        )
      );


      // -----------------------------------------
      // Calculate:
      //
      // p_i = n_i / total
      //
      // -p_i log2(p_i)
      // -----------------------------------------

      var entropyParts = counts.map(
        function(count) {

          var p = ee.Number(count)
            .divide(total);

          return p
            .multiply(
              p.log()
                .divide(Math.log(2))
            )
            .multiply(-1);
        }
      );


      // -----------------------------------------
      // Sum contributions
      // -----------------------------------------

      return ee.Number(
        entropyParts.reduce(
          ee.Reducer.sum()
        )
      );

    })(),

    null
  );


  return entropy;
}


// =====================================================
// 9. ENTROPY BY DATE AND INDEX
// =====================================================

var entropyCollection = s2_col_norm
  .map(function(image) {

    var date = ee.Date(
      image.get('system:time_start')
    );


    var features = ee.List(indices)
      .map(function(indexName) {

        indexName = ee.String(indexName);


        var entropy =
          calculateEntropy(
            image,
            indexName
          );


        return ee.Feature(
          null,
          {

            date:
              date.format(
                'YYYY-MM-dd'
              ),

            time:
              date.millis(),

            index:
              indexName,

            entropy:
              entropy
          }
        );

      });


    return ee.FeatureCollection(
      features
    );

  })
  .flatten()
  .filter(
    ee.Filter.notNull(
      ['entropy']
    )
  );


print(
  'Shannon entropy by date and index:',
  entropyCollection
);


// =====================================================
// 10. ENTROPY TIME-SERIES CHART
// =====================================================

var entropyChart =
  ui.Chart.feature.groups({

    features:
      entropyCollection,

    xProperty:
      'time',

    yProperty:
      'entropy',

    seriesProperty:
      'index'

  })
  .setChartType(
    'LineChart'
  )
  .setOptions({

    title:
      'Temporal evolution of Shannon entropy',

    lineWidth: 2,

    pointSize: 5,

    hAxis: {

      title:
        'Date',

      format:
        'dd-MM-yyyy',

      slantedText:
        true,

      slantedTextAngle:
        45
    },

    vAxis: {

      title:
        'Shannon entropy (bits)',

      viewWindow: {
        min: 0,
        max: 5
      }
    },

    legend: {
      position: 'right'
    }

  });


print(entropyChart);


// =====================================================
// 11. EXPORT ENTROPY TABLE
// =====================================================

Export.table.toDrive({

  collection:
    entropyCollection,

  description:
    'Shannon_entropy_Sentinel2_indices',

  fileNamePrefix:
    'Shannon_entropy_Sentinel2_indices',

  fileFormat:
    'CSV',

  selectors: [
    'date',
    'index',
    'entropy'
  ]

});


// =====================================================
// 12. EXPORT NORMALIZED INDEX MEANS
// =====================================================
//
// Optional:
// Export the mean of each index for every acquisition date.
// =====================================================

var meanCollection = s2_col_norm
  .map(function(image) {

    var date = ee.Date(
      image.get('system:time_start')
    );


    var means = image.reduceRegion({

      reducer:
        ee.Reducer.mean(),

      geometry:
        geometry,

      scale:
        10,

      maxPixels:
        1e9,

      tileScale:
        4
    });


    return ee.Feature(
      null,
      means
    )
    .set(
      'date',
      date.format(
        'YYYY-MM-dd'
      )
    )
    .set(
      'time',
      date.millis()
    );

  });


print(
  'Mean normalized indices by date:',
  meanCollection
);


Export.table.toDrive({

  collection:
    meanCollection,

  description:
    'Mean_normalized_Sentinel2_indices',

  fileNamePrefix:
    'Mean_normalized_Sentinel2_indices',

  fileFormat:
    'CSV'
});


// =====================================================
// 13. MAP VISUALIZATION
// =====================================================

Map.centerObject(
  geometry,
  13
);


// Display area of interest
Map.addLayer(
  geometry,
  {
    color: 'red'
  },
  'Study area'
);


// =====================================================
// 14. SENTINEL-2 RGB COMPOSITE
// =====================================================

var rgb = s2_raw
  .median()
  .divide(10000);


Map.addLayer(
  rgb,
  {
    bands: [
      'B4',
      'B3',
      'B2'
    ],

    min: 0,

    max: 0.3
  },
  'RGB median'
);


// =====================================================
// 15. MEAN VEGETATION INDEX IMAGES
// =====================================================

var meanImg =
  s2_col.mean();


// NDVI
Map.addLayer(
  meanImg.select('NDVI'),
  {
    min: -1,
    max: 1,
    palette: [
      'blue',
      'white',
      'green'
    ]
  },
  'NDVI',
  false
);


// GNDVI
Map.addLayer(
  meanImg.select('GNDVI'),
  {
    min: -1,
    max: 1,
    palette: [
      'purple',
      'white',
      'green'
    ]
  },
  'GNDVI',
  false
);


// MSR
Map.addLayer(
  meanImg.select('MSR'),
  {
    min: 0,
    max: 3,
    palette: [
      'white',
      'green'
    ]
  },
  'MSR',
  false
);


// SAVI
Map.addLayer(
  meanImg.select('SAVI'),
  {
    min: -1,
    max: 1,
    palette: [
      'brown',
      'yellow',
      'green'
    ]
  },
  'SAVI',
  false
);


// CLgreen
Map.addLayer(
  meanImg.select('CLgreen'),
  {
    min: 0,
    max: 5,
    palette: [
      'white',
      'green'
    ]
  },
  'CLgreen',
  false
);


// RVI
Map.addLayer(
  meanImg.select('RVI'),
  {
    min: 0,
    max: 10,
    palette: [
      'white',
      'blue'
    ]
  },
  'RVI',
  false
);


// GRVI
Map.addLayer(
  meanImg.select('GRVI'),
  {
    min: 0,
    max: 10,
    palette: [
      'white',
      'green'
    ]
  },
  'GRVI',
  false
);


// NRVI
Map.addLayer(
  meanImg.select('NRVI'),
  {
    min: -1,
    max: 1,
    palette: [
      'blue',
      'white',
      'green'
    ]
  },
  'NRVI',
  false
);


// MTVI
Map.addLayer(
  meanImg.select('MTVI'),
  {
    min: -1,
    max: 1,
    palette: [
      'purple',
      'white',
      'green'
    ]
  },
  'MTVI',
  false
);

// =====================================================
// 16. EXPORT SOURCE IMAGE METADATA FOR REPRODUCIBILITY
// =====================================================
// One row per image after date, footprint, and scene cloud filtering.
// Pixel-level masking does not remove scenes from the collection.
var acquisitionMetadata = ee.FeatureCollection(s2_raw.map(function(image) {
  return ee.Feature(null, {
    acquisition_date: ee.Date(image.get('system:time_start')).format('YYYY-MM-dd'),
    system_index: image.get('system:index'),
    cloud_percentage: image.get('CLOUDY_PIXEL_PERCENTAGE')
  });
}));

Export.table.toDrive({
  collection: acquisitionMetadata,
  description: 'Sentinel2_acquisition_metadata',
  fileNamePrefix: 'Sentinel2_acquisition_metadata',
  fileFormat: 'CSV'
});
