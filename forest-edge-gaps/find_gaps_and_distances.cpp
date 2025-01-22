
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sstream>
#include <iostream>
// #include <filesystem>

#include <sys/time.h>
#include <sys/stat.h>
#include <stdlib.h>

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"

#include "vertex_grid.h"
#include "vertex.h"

// double gettime();


int main ( int argc, char* argv[] ) {

  /* 
  */

  // must be passed as argumets
  std::string height_tif = "";
  std::string gap_height_str = "";
  std::string max_distance_str = "";
  std::string border_length_str = "";
  std::string save_prefix = "";
  std::string save_dir = "";

  // optional
  std::string save_nodata_value_str = "";

  // reading argumets
  bool arguments_ok = true;
  int i = 1;
  while ( (i < argc) && (arguments_ok) ) {
    std::string argument = std::string(argv[i]);

    if (argument.compare("-height_tif") == 0) {
      if (argv[i+1] == NULL) {
	std::cout << "main: Error: command not ok." << std::endl;
	arguments_ok = false;
      } else {
	i++;
	height_tif = std::string(argv[i]);
      }

    } else if (argument.compare("-gap_height") == 0) {
      if (argv[i+1] == NULL) {
	std::cout << "main: Error: command not ok." << std::endl;
	arguments_ok = false;
      } else {
	i++;
	gap_height_str = std::string(argv[i]);
      }

    } else if (argument.compare("-max_distance") == 0) {
      if (argv[i+1] == NULL) {
	std::cout << "main: Error: command not ok." << std::endl;
	arguments_ok = false;
      } else {
	i++;
	max_distance_str = std::string(argv[i]);
      }

    } else if (argument.compare("-border_length") == 0) {
      if (argv[i+1] == NULL) {
	std::cout << "main: Error: command not ok." << std::endl;
	arguments_ok = false;
      } else {
	i++;
	border_length_str = std::string(argv[i]);
      }

    } else if (argument.compare("-save_prefix") == 0) {
      if (argv[i+1] == NULL) {
	std::cout << "main: Error: command not ok." << std::endl;
	arguments_ok = false;
      } else {
	i++;
	save_prefix = std::string(argv[i]);
      }

    } else if (argument.compare("-save_dir") == 0) {
      if (argv[i+1] == NULL) {
	std::cout << "main: Error: command not ok." << std::endl;
	arguments_ok = false;
      } else {
	i++;
	save_dir = std::string(argv[i]);
      }

    } else if (argument.compare("-save_nodata_value") == 0) {
      if (argv[i+1] == NULL) {
	std::cout << "main: Error: command not ok." << std::endl;
	arguments_ok = false;
      } else {
	i++;
	save_nodata_value_str = std::string(argv[i]);
      }

    } else {
      std::cout << "main: Error: command not ok: " << argument << std::endl;
      arguments_ok = false;
    }

    // std::cout << "main: argument: " << argument << " value: "
    //      << std::string(argv[i]) << std::endl;

    i++;
  }

  if (height_tif.compare("") == 0) {
    std::cout << "main: no height_tif." << std::endl;
    arguments_ok = false;
  }

  if (gap_height_str.compare("") == 0) {
    std::cout << "main: no gap_height." << std::endl;
    arguments_ok = false;
  }

  if (max_distance_str.compare("") == 0) {
    std::cout << "main: no max_distance." << std::endl;
    arguments_ok = false;
  }

  if (border_length_str.compare("") == 0) {
    std::cout << "main: no border_length." << std::endl;
    arguments_ok = false;
  }

  if (save_prefix.compare("") == 0) {
    std::cout << "main: no save_prefix." << std::endl;
    arguments_ok = false;
  }

  if (save_dir.compare("") == 0) {
    std::cout << "main: no save_dir." << std::endl;
    arguments_ok = false;
  }

  if (arguments_ok == false) {
    exit(1);
  }


  const double gap_height = atof(gap_height_str.c_str());
  const int max_distance = atoi(max_distance_str.c_str());

  const int default_band = 1;

  const int ew_off = 0;
  const int ns_off = 0;
  const int nLineSpace = 0;
  const int psExtraArg = 0;

  GDALAllRegister();

  CPLErr eErr = CE_None;

  printf("grabbing height...\n");

  GDALDataset* height_ds;

  height_ds = (GDALDataset*) GDALOpen(height_tif.c_str(), GA_ReadOnly);
  if( height_ds == NULL ) {
    printf("height_ds == NULL");
    exit(1);
  }

  static const int ns_size = GDALGetRasterYSize( height_ds ); 
  static const int ew_size = GDALGetRasterXSize( height_ds ); 
  // static const int n_bands = GDALGetRasterCount( height_ds );

  // possible bug: without this, the epsg seg faults
  const std::string height_projection \
    = (std::string) height_ds->GetProjectionRef();

  const OGRSpatialReference* height_srs = height_ds->GetSpatialRef();
  const char* height_epsg_char = height_srs->GetAuthorityCode(NULL);

  if (height_epsg_char == NULL) {
    printf("ERROR: height_epsg_char == NULL.\n");
    exit(1);
  }

  const int height_epsg = atoi(height_epsg_char);
  const std::string height_epsg_string = (std::string) height_epsg_char;

  double height_transform[6];
  height_ds->GetGeoTransform(height_transform);
  double ew_start = height_transform[0];
  double ew_res = height_transform[1];
  double ew_skew = height_transform[2];
  double ns_start = height_transform[3];
  double ns_skew = height_transform[4];
  double ns_res = height_transform[5];


  const double border_length = atof(border_length_str.c_str());

  // double remainder_1 = fmod(border_length, ew_res);
  // double remainder_2 = remainder(border_length, ew_res);
  // double remainder_3 = remainderf(border_length, ew_res);
  // printf("remainder_1: %f\n", remainder_1);
  // printf("remainder_2: %f\n", remainder_2);
  // printf("remainder_3: %f\n", remainder_3);

  if (remainderf(border_length, ew_res) != 0) {
    printf("ERROR: remainderf(border_length, ew_res) != 0.\n");
    printf("  border_length: %f\n", border_length);
    printf("  ew_res: %f\n", ew_res);
    exit(1);
  }

  if (remainderf(border_length, ns_res) != 0) {
    printf("ERROR: remainderf(border_length, ns_res) != 0.\n");
    printf("  border_length: %f\n", border_length);
    printf("  ns_res: %f\n", ns_res);
    exit(1);
  }

  const long int ew_border_int = abs(lround(border_length/ew_res));
  const long int ns_border_int = abs(lround(border_length/ns_res));


  printf("gap_height: %f\n", gap_height);
  printf("max_distance: %d\n", max_distance);
  printf("border_length: %f\n", border_length);

  printf("border_length: %f\n", border_length);
  printf("ew_res: %f\n", ew_res);
  printf("ns_res: %f\n", ns_res);


  GDALRasterBand* height_band;
  height_band = (GDALRasterBand*) height_ds->GetRasterBand(default_band);
  float height_nodata_value = height_band->GetNoDataValue();  

  float save_nodata_value;
  if (save_nodata_value_str.compare("") == 0) {
    save_nodata_value = float(height_nodata_value);
  } else {
    save_nodata_value = std::stof(save_nodata_value_str);
  }


  Vertex_Grid vertex_array = Vertex_Grid ( ns_size, ew_size,
                                           height_epsg, 
                                           ew_start, ew_res, ew_skew,
                                           ns_start, ns_skew, ns_res);

  float* height_array_buffer;
  height_array_buffer = (float*) CPLMalloc(sizeof(float) * ew_size * ns_size);

  eErr = height_band->RasterIO( GF_Read,
                                ew_off, ns_off, ew_size, ns_size,
                                height_array_buffer, ew_size, ns_size,
                                GDT_Float32,
                                nLineSpace, psExtraArg );

  for ( int row = 0; eErr == CE_None && row < ns_size; row++ ) {
    for ( int col = 0; col < ew_size; col++ ) {
      Vertex* current_vertex = vertex_array.GetVertex(row, col);

      current_vertex->
        SetVertexHeightPtr(&height_array_buffer[(row * ew_size) + col]);
    }
  }

  GDALClose(height_ds);


  printf("checking geotiff driver...\n");

  const char* format_char = "GTiff";
  GDALDriver* geotiff_driver;
  char** driver_metadata;

  geotiff_driver = GetGDALDriverManager()->GetDriverByName(format_char);
  if ( geotiff_driver == NULL ) {
    printf("ERROR: geotiff_driver == NULL \n");
    exit(1);
  }

  driver_metadata = geotiff_driver->GetMetadata();

  if ( ! CSLFetchBoolean( driver_metadata, GDAL_DCAP_CREATE, FALSE ) ) {
    printf( "ERROR: Driver %s doesn't support Create() method.\n",
            format_char );
    exit(1);
  }

  char** options_char = NULL;
  options_char = CSLSetNameValue( options_char, "COMPRESS", "LZW" );
  options_char = CSLSetNameValue( options_char, "TILED", "YES"  );
  options_char = CSLSetNameValue( options_char, "BIGTIFF", "YES"  );


  double result_transform[6] = {
    height_transform[0] + border_length,
    height_transform[1],
    height_transform[2],
    height_transform[3] - border_length,
    height_transform[4],
    height_transform[5]
  };

  int result_ew_size = ew_size - (2 * ew_border_int);
  int result_ns_size = ns_size - (2 * ns_border_int);

  printf("ew_border_int: %ld\n", ew_border_int);
  printf("result_ew_size: %d\n", result_ew_size);
  printf("ns_border_int: %ld\n", ns_border_int);
  printf("result_ns_size: %d\n", result_ns_size);

  int* result_array_buffer;

  GDALDataset* result_ds;
  GDALRasterBand* result_band;

  std::string save_tif_wo_path; 
  std::string save_tif;


  printf("calc east west...\n");

  vertex_array.CalcEW(height_nodata_value, gap_height, max_distance);


  printf("save gap size east...\n");

  result_array_buffer
    = (int*) CPLMalloc(sizeof(int) * result_ew_size * result_ns_size);

  save_tif_wo_path = save_prefix + "_east_gap_size.tif";
  save_tif = save_dir + "/" + save_tif_wo_path;

  for ( int row = 0; eErr == CE_None && row < (result_ns_size); row++ ) {
    for ( int col = 0; col < result_ew_size; col++ ) {
      Vertex* current_vertex = vertex_array.GetVertex(row + ns_border_int,
                                                      col + ns_border_int);

      result_array_buffer[(row * result_ew_size) + col]
        = current_vertex->GetGapSizeRight();
    }
  }

  result_ds = geotiff_driver->Create( save_tif.c_str(),
                                      result_ew_size, result_ns_size,
                                      default_band, GDT_Int32,
                                      options_char );

  result_ds->SetGeoTransform( result_transform );

  result_ds->SetProjection( height_projection.c_str() );

  result_band = result_ds->GetRasterBand(default_band);

  fflush(stdout);

  eErr = result_band->SetNoDataValue(save_nodata_value);  

  eErr = result_band->RasterIO( GF_Write,
                                ew_off, ns_off, result_ew_size, result_ns_size,
                                result_array_buffer,
                                result_ew_size, result_ns_size,
                                GDT_Int32,
                                nLineSpace, psExtraArg );

  GDALClose( (GDALDatasetH) result_ds );

  free (result_array_buffer);


  printf("save gap distance east...\n");

  result_array_buffer
    = (int*) CPLMalloc(sizeof(int) * result_ew_size * result_ns_size);

  save_tif_wo_path = save_prefix + "_east_gap_distance.tif";
  save_tif = save_dir + "/" + save_tif_wo_path;

  for ( int row = 0; eErr == CE_None && row < (result_ns_size); row++ ) {
    for ( int col = 0; col < result_ew_size; col++ ) {
      Vertex* current_vertex = vertex_array.GetVertex(row + ns_border_int,
                                                      col + ns_border_int);

      result_array_buffer[(row * result_ew_size) + col]
        = current_vertex->GetDistanceToGapRight();
    }
  }

  result_ds = geotiff_driver->Create( save_tif.c_str(),
                                       result_ew_size, result_ns_size,
                                       default_band, GDT_Int32,
                                       options_char );

  result_ds->SetGeoTransform( result_transform );

  result_ds->SetProjection( height_projection.c_str() );

  result_band = result_ds->GetRasterBand(default_band);

  fflush(stdout);

  eErr = result_band->SetNoDataValue(save_nodata_value);  

  eErr = result_band->RasterIO( GF_Write,
                                ew_off, ns_off, result_ew_size, result_ns_size,
                                result_array_buffer,
                                result_ew_size, result_ns_size,
                                GDT_Int32,
                                nLineSpace, psExtraArg );

  GDALClose( (GDALDatasetH) result_ds );
  free (result_array_buffer);


  printf("save gap size west...\n");

  result_array_buffer
    = (int*) CPLMalloc(sizeof(int) * result_ew_size * result_ns_size);

  save_tif_wo_path = save_prefix + "_west_gap_size.tif";
  save_tif = save_dir + "/" + save_tif_wo_path;

  for ( int row = 0; eErr == CE_None && row < (result_ns_size); row++ ) {
    for ( int col = 0; col < result_ew_size; col++ ) {
      Vertex* current_vertex = vertex_array.GetVertex(row + ns_border_int,
                                                      col + ns_border_int);

      result_array_buffer[(row * result_ew_size) + col]
        = current_vertex->GetGapSizeLeft();
    }
  }

  result_ds = geotiff_driver->Create( save_tif.c_str(),
                                      result_ew_size, result_ns_size,
                                      default_band, GDT_Int32,
                                      options_char );

  result_ds->SetGeoTransform( result_transform );

  result_ds->SetProjection( height_projection.c_str() );

  result_band = result_ds->GetRasterBand(default_band);

  fflush(stdout);

  eErr = result_band->SetNoDataValue(save_nodata_value);  

  eErr = result_band->RasterIO( GF_Write,
                                ew_off, ns_off, result_ew_size, result_ns_size,
                                result_array_buffer,
                                result_ew_size, result_ns_size,
                                GDT_Int32,
                                nLineSpace, psExtraArg );

  GDALClose( (GDALDatasetH) result_ds );
  free (result_array_buffer);


  printf("save gap distance west...\n");

  result_array_buffer
    = (int*) CPLMalloc(sizeof(int) * result_ew_size * result_ns_size);

  save_tif_wo_path = save_prefix + "_west_gap_distance.tif";
  save_tif = save_dir + "/" + save_tif_wo_path;

  for ( int row = 0; eErr == CE_None && row < (result_ns_size); row++ ) {
    for ( int col = 0; col < result_ew_size; col++ ) {
      Vertex* current_vertex = vertex_array.GetVertex(row + ns_border_int,
                                                      col + ns_border_int);

      result_array_buffer[(row * result_ew_size) + col]
        = current_vertex->GetDistanceToGapLeft();
    }
  }

  result_ds = geotiff_driver->Create( save_tif.c_str(),
                                       result_ew_size, result_ns_size,
                                       default_band, GDT_Int32,
                                       options_char );

  result_ds->SetGeoTransform( result_transform );

  result_ds->SetProjection( height_projection.c_str() );

  result_band = result_ds->GetRasterBand(default_band);

  fflush(stdout);

  eErr = result_band->SetNoDataValue(save_nodata_value);  

  eErr = result_band->RasterIO( GF_Write,
                                ew_off, ns_off, result_ew_size, result_ns_size,
                                result_array_buffer,
                                result_ew_size, result_ns_size,
                                GDT_Int32,
                                nLineSpace, psExtraArg );

  GDALClose( (GDALDatasetH) result_ds );
  free (result_array_buffer);


  printf("calc north south...\n");

  vertex_array.CalcNS(height_nodata_value, gap_height, max_distance);


  printf("save gap size south...\n");

  result_array_buffer
    = (int*) CPLMalloc(sizeof(int) * result_ew_size * result_ns_size);

  save_tif_wo_path = save_prefix + "_south_gap_size.tif";
  save_tif = save_dir + "/" + save_tif_wo_path;

  for ( int row = 0; eErr == CE_None && row < (result_ns_size); row++ ) {
    for ( int col = 0; col < result_ew_size; col++ ) {
      Vertex* current_vertex = vertex_array.GetVertex(row + ns_border_int,
                                                      col + ns_border_int);

      result_array_buffer[(row * result_ew_size) + col]
        = current_vertex->GetGapSizeRight();
    }
  }

  result_ds = geotiff_driver->Create( save_tif.c_str(),
                                       result_ew_size, result_ns_size,
                                       default_band, GDT_Int32,
                                       options_char );

  result_ds->SetGeoTransform( result_transform );

  result_ds->SetProjection( height_projection.c_str() );

  result_band = result_ds->GetRasterBand(default_band);

  fflush(stdout);

  eErr = result_band->SetNoDataValue(save_nodata_value);  

  eErr = result_band->RasterIO( GF_Write,
                                ew_off, ns_off, result_ew_size, result_ns_size,
                                result_array_buffer,
                                result_ew_size, result_ns_size,
                                GDT_Int32,
                                nLineSpace, psExtraArg );

  GDALClose( (GDALDatasetH) result_ds );
  free (result_array_buffer);


  printf("save gap distance south...\n");

  result_array_buffer
    = (int*) CPLMalloc(sizeof(int) * result_ew_size * result_ns_size);

  save_tif_wo_path = save_prefix + "_south_gap_distance.tif";
  save_tif = save_dir + "/" + save_tif_wo_path;

  for ( int row = 0; eErr == CE_None && row < (result_ns_size); row++ ) {
    for ( int col = 0; col < result_ew_size; col++ ) {
      Vertex* current_vertex = vertex_array.GetVertex(row + ns_border_int,
                                                      col + ns_border_int);

      result_array_buffer[(row * result_ew_size) + col]
        = current_vertex->GetDistanceToGapRight();
    }
  }

  result_ds = geotiff_driver->Create( save_tif.c_str(),
                                       result_ew_size, result_ns_size,
                                       default_band, GDT_Int32,
                                       options_char );

  result_ds->SetGeoTransform( result_transform );

  result_ds->SetProjection( height_projection.c_str() );

  result_band = result_ds->GetRasterBand(default_band);

  fflush(stdout);

  eErr = result_band->SetNoDataValue(save_nodata_value);  

  eErr = result_band->RasterIO( GF_Write,
                                ew_off, ns_off, result_ew_size, result_ns_size,
                                result_array_buffer,
                                result_ew_size, result_ns_size,
                                GDT_Int32,
                                nLineSpace, psExtraArg );

  GDALClose( (GDALDatasetH) result_ds );
  free (result_array_buffer);


  printf("save gap size north...\n");

  result_array_buffer
    = (int*) CPLMalloc(sizeof(int) * result_ew_size * result_ns_size);

  save_tif_wo_path = save_prefix + "_north_gap_size.tif";
  save_tif = save_dir + "/" + save_tif_wo_path;

  for ( int row = 0; eErr == CE_None && row < (result_ns_size); row++ ) {
    for ( int col = 0; col < result_ew_size; col++ ) {
      Vertex* current_vertex = vertex_array.GetVertex(row + ns_border_int,
                                                      col + ns_border_int);

      result_array_buffer[(row * result_ew_size) + col]
        = current_vertex->GetGapSizeLeft();
    }
  }

  result_ds = geotiff_driver->Create( save_tif.c_str(),
                                       result_ew_size, result_ns_size,
                                       default_band, GDT_Int32,
                                       options_char );

  result_ds->SetGeoTransform( result_transform );

  result_ds->SetProjection( height_projection.c_str() );

  result_band = result_ds->GetRasterBand(default_band);

  fflush(stdout);

  eErr = result_band->SetNoDataValue(save_nodata_value);  

  eErr = result_band->RasterIO( GF_Write,
                                ew_off, ns_off, result_ew_size, result_ns_size,
                                result_array_buffer,
                                result_ew_size, result_ns_size,
                                GDT_Int32,
                                nLineSpace, psExtraArg );

  GDALClose( (GDALDatasetH) result_ds );
  free (result_array_buffer);


  printf("save gap distance north...\n");

  result_array_buffer
    = (int*) CPLMalloc(sizeof(int) * result_ew_size * result_ns_size);

  save_tif_wo_path = save_prefix + "_north_gap_distance.tif";
  save_tif = save_dir + "/" + save_tif_wo_path;

  for ( int row = 0; eErr == CE_None && row < (result_ns_size); row++ ) {
    for ( int col = 0; col < result_ew_size; col++ ) {
      Vertex* current_vertex = vertex_array.GetVertex(row + ns_border_int,
                                                      col + ns_border_int);

      result_array_buffer[(row * result_ew_size) + col]
        = current_vertex->GetDistanceToGapLeft();
    }
  }

  result_ds = geotiff_driver->Create( save_tif.c_str(),
                                       result_ew_size, result_ns_size,
                                       default_band, GDT_Int32,
                                       options_char );

  result_ds->SetGeoTransform( result_transform );

  result_ds->SetProjection( height_projection.c_str() );

  result_band = result_ds->GetRasterBand(default_band);

  fflush(stdout);

  eErr = result_band->SetNoDataValue(save_nodata_value);  

  eErr = result_band->RasterIO( GF_Write,
                                ew_off, ns_off, result_ew_size, result_ns_size,
                                result_array_buffer,
                                result_ew_size, result_ns_size,
                                GDT_Int32,
                                nLineSpace, psExtraArg );

  GDALClose( (GDALDatasetH) result_ds );
  free (result_array_buffer);



  free (height_array_buffer);

  GDALDestroyDriverManager();

  return 0;
}

