/******************************************************************************
 *
 * Project:  FITS Driver
 * Purpose:  Implement FITS raster read/write support
 * Author:   Simon Perkins, s.perkins@lanl.gov
 *
 ******************************************************************************
 * Copyright (c) 2001, Simon Perkins
 * Copyright (c) 2008-2018, Even Rouault <even dot rouault at spatialys.com>
 * Copyright (c) 2018, Chiara Marmo <chiara dot marmo at u-psud dot fr>
 * Copyright (c) 2019, Jean-Christophe Malapert <jean-christophe dot malapert at cnes dot fr>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/

#include "crs.h"
#include "constants.h"
#include "fitsWrapper.h"
#include "wcslib.h"
#include <exception>
#include <iomanip>
#include <iostream>
#include <math.h>



Crs::Crs(fitsfile *headerFits)
    : hFITS(headerFits), naxes(getNaxes(headerFits)) {}

int Crs::getNaxes(fitsfile *headerFits) const {
  int nbAxes = fits_read_key_as_int(headerFits, GlobalConstants::WCS_NAXIS);
  return nbAxes >= 2 ? 2 : 1;
}

std::vector<int> Crs::getNaxis() const {
  std::vector<int> naxis(naxes);
  for (int i = 1; i <= naxes; i++) {
    std::string keyword = GlobalConstants::WCS_NAXIS;
    keyword += std::to_string(i);
    naxis[i - 1] = fits_read_key_as_int(hFITS, keyword);
  }
  return naxis;
}

std::vector<double> Crs::getCrval() const {
  std::vector<double> crval(naxes);
  for (int i = 1; i <= naxes; i++) {
    std::string keyword = GlobalConstants::WCS_CRVAL;
    keyword += std::to_string(i);
    crval[i - 1] = fits_read_key_as_double(hFITS, keyword);
  }
  return crval;
}

std::vector<float> Crs::getCrpix() const {
  std::vector<float> crpix(naxes);
  for (int i = 1; i <= naxes; i++) {
    std::string keyword = GlobalConstants::WCS_CRPIX;
    keyword += std::to_string(i);
    crpix[i - 1] = fits_read_key_as_float(hFITS, keyword);
  }
  return crpix;
}

std::vector<double> Crs::getCd() const {
  std::vector<double> cd(naxes * 2);
  cd[0] = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CD11);
  cd[1] = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CD12);
  cd[2] = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CD21);
  cd[3] = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CD22);
  return cd;
}

void Crs::init() {
  crs = createCrs();  
}

const OGRSpatialReference* Crs::getCrs() {  
  return crs.IsEmpty() ? nullptr : &crs;
}

Planet::Planet(fitsfile *headerFits) : Crs(headerFits) {}

bool Planet::double_equals(double a, double b, double epsilon) {
  return std::abs(a - b) < epsilon;
}

std::vector<double> Planet::getGeoTransform() {
  std::vector<double> adfGeoTransform(6);
  double crpix1 = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CRPIX1.c_str());
  double crpix2 = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CRPIX2.c_str());
  double crval1 = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CRVAL1.c_str());
  double crval2 = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CRVAL2.c_str());
  std::vector<double> cd(4);
  try {
    double cdelt1 = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CDELT1.c_str());
    double cdelt2 = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CDELT2.c_str());
    double pc11 = fits_read_key_as_double(hFITS, GlobalConstants::WCS_PC11.c_str());
    double pc12 = fits_read_key_as_double(hFITS, GlobalConstants::WCS_PC12.c_str());
    double pc21 = fits_read_key_as_double(hFITS, GlobalConstants::WCS_PC21.c_str());
    double pc22 = fits_read_key_as_double(hFITS, GlobalConstants::WCS_PC22.c_str());
    cd[0] = cdelt1 * pc11;
    cd[1] = cdelt1 * pc12;
    cd[2] = cdelt2 * pc21;
    cd[3] = cdelt2 * pc22;
  } catch (const FitsException &e) {
    cd[0] = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CD11.c_str());
    cd[1] = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CD12.c_str());
    cd[2] = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CD21.c_str());
    cd[3] = fits_read_key_as_double(hFITS, GlobalConstants::WCS_CD22.c_str());
  }
  double aRadius = fits_read_key_as_double(hFITS, GlobalConstants::FITS_A_RADIUS.c_str());  
  double radfac = GlobalConstants::DEG2RAD * aRadius;

  adfGeoTransform[1] = cd[0] * radfac;
  adfGeoTransform[2] = cd[1] * radfac;
  adfGeoTransform[4] = cd[2] * radfac;
  adfGeoTransform[5] = -cd[3] * radfac;
  if (crval1 > 180.) {
    crval1 = crval1 - 180.;
  }

  /* NOTA BENE: FITS standard define pixel integers at the center of the pixel,
     0.5 must be subtract to have UpperLeft corner */
  adfGeoTransform[0] = crval1 * radfac - adfGeoTransform[1] * (crpix1 - 0.5);
  // assuming that center latitude is also the origin of the coordinate
  // system: this is not always true.
  // More generic implementation coming soon
  adfGeoTransform[3] = -adfGeoTransform[5] * (crpix2 - 0.5);
  //+ crval2 * radfac;

  return adfGeoTransform;
}

std::string Planet::getCrsName() {
  std::string geogName;
  std::string target = fits_read_key_as_str(hFITS, GlobalConstants::FITS_OBJECT.c_str());
  geogName.assign("GCS_");
  geogName.append(target);
  return geogName;  
}

std::string Planet::getAxisName(const int axisNumber) {
  std::string axis;
  std::string ctypeName =
      GlobalConstants::WCS_CTYPE + std::to_string(axisNumber);
  std::string ctype = fits_read_key_as_str(hFITS, ctypeName.c_str());
  std::string axisValue = ctype.substr(0, 3);
  if (axisValue == "LN-") {
    axis = GlobalConstants::WKT_AXIS_LONG;
  } else if (axisValue == "LT-") {
    axis = GlobalConstants::WKT_AXIS_LAT;
  } else {
    throw CrsException("Unknow axis name");
  }
  return axis;  
}

std::string Planet::getAxisUnit(const int axisNumber) {
  return "";  
}

OGRSpatialReference Planet::createCrs() {
  OGRSpatialReference oSRS = OGRSpatialReference();
  createDatum(&oSRS);
  createCoordFrame(&oSRS);
  return oSRS;
}

void Planet::createDatum(OGRSpatialReference *oSRS) {
  double invFlattening;
  std::string geogName;
  std::string datumName;
  std::string target = fits_read_key_as_str(hFITS, GlobalConstants::FITS_OBJECT.c_str());
  geogName.assign("GCS_");
  geogName.append(target);
  datumName.assign("D_");
  datumName.append(target);
  double aRadius = fits_read_key_as_double(hFITS, GlobalConstants::FITS_A_RADIUS.c_str());
  double cRadius = fits_read_key_as_double(hFITS, GlobalConstants::FITS_C_RADIUS.c_str());
  if (!double_equals(aRadius, cRadius)) {
    invFlattening = aRadius / (aRadius - cRadius);
  } else {
    invFlattening = 0.0;
  }
  oSRS->SetGeogCS(geogName.c_str(), datumName.c_str(), target.c_str(), aRadius,
                  invFlattening, "Reference_Meridian", 0.0, SRS_UA_DEGREE,
                  GlobalConstants::DEG2RAD);
}

void Planet::createCoordFrame(OGRSpatialReference *oSRS) {
  double falseEast = 0.0, falseNorth = 0.0, scale = 1.0;
  std::string projName;

  // Define projected CRS
  const std::string ctypeName = GlobalConstants::WCS_CTYPE + std::to_string(1);
  std::string ctype = fits_read_key_as_str(hFITS, ctypeName);
  char *ctypePtr = &ctype[0];
  char *pstr = strrchr(ctypePtr, '-');
  if (!pstr) {
    throw CrsException("Unknown projection.");
  } else {
    pstr += 1;
  }

  std::string target = fits_read_key_as_str(hFITS, GlobalConstants::FITS_OBJECT.c_str());
  std::vector<double> crval = getCrval();

  /* Sinusoidal / SFL projection */
  if (strcmp(pstr, GlobalConstants::WCS_PROJ_SFL.c_str()) == 0) {
    projName.assign("Sinusoidal_");
    oSRS->SetSinusoidal(crval[0], falseEast, falseNorth);

    /* Mercator, Oblique (Hotine) Mercator, Transverse Mercator */
    /* Mercator / MER projection */
  } else if (strcmp(pstr, GlobalConstants::WCS_PROJ_MER.c_str()) == 0) {
    projName.assign("Mercator_");
    oSRS->SetMercator(crval[1], crval[0], scale, falseEast, falseNorth);

    /* Equirectangular / CAR projection */
  } else if (strcmp(pstr, GlobalConstants::WCS_PROJ_CAR.c_str()) == 0) {
    projName.assign("Equirectangular_");
    /*
    The standard_parallel_1 defines where the local radius is calculated
    not the center of Y Cartesian system (which is latitude_of_origin)
    But FITS WCS only supports projections on the sphere
    we assume here that the local radius is the one computed at the projection
    center
    */
    oSRS->SetEquirectangular2(crval[1], crval[0], crval[1], falseEast, falseNorth);
    /* Lambert Azimuthal Equal Area / ZEA projection */
  } else if (strcmp(pstr, GlobalConstants::WCS_PROJ_ZEA.c_str()) == 0) {
    projName.assign("Lambert_Azimuthal_Equal_Area_");
    oSRS->SetLAEA(crval[1], crval[0], falseEast, falseNorth);

    /* Lambert Conformal Conic 1SP / COO projection */
  } else if (strcmp(pstr, GlobalConstants::WCS_PROJ_COO.c_str()) == 0) {
    projName.assign("Lambert_Conformal_Conic_1SP_");
    oSRS->SetLCC1SP(crval[1], crval[0], scale, falseEast, falseNorth);

    /* Orthographic / SIN projection */
  } else if (strcmp(pstr, GlobalConstants::WCS_PROJ_SIN.c_str()) == 0) {
    projName.assign("Orthographic_");
    oSRS->SetOrthographic(crval[1], crval[0], falseEast, falseNorth);

    /* Point Perspective / AZP projection */
  } else if (strcmp(pstr, GlobalConstants::WCS_PROJ_AZP.c_str()) == 0) {
    projName.assign("perspective_point_height_");
    oSRS->SetProjection(SRS_PP_PERSPECTIVE_POINT_HEIGHT);
    /* # appears to need height... maybe center lon/lat */

    /* Polar Stereographic / STG projection */
  } else if (strcmp(pstr, GlobalConstants::WCS_PROJ_STG.c_str()) == 0) {
    projName.assign("Polar_Stereographic_");
    oSRS->SetStereographic(crval[1], crval[0], scale, falseEast, falseNorth);
  } else {
    throw CrsException("Unknown projection");
  }
  projName.append(target);
  //oSRS->SetProjParm(SRS_PP_FALSE_EASTING, 0.0);
  //oSRS->SetProjParm(SRS_PP_FALSE_NORTHING, 0.0);
  oSRS->SetNode(GlobalConstants::WKT_PROJCS.c_str(), projName.c_str());
  oSRS->SetAxes(GlobalConstants::WKT_PROJCS.c_str(), getAxisName(1).c_str(), OAO_East,
                getAxisName(2).c_str(), OAO_North);  
}

Sky::Sky(fitsfile *headerFits) : Crs(headerFits) {
  char *header;
  int nkeyrec, status=0, nreject;

  /* Read in the FITS header, excluding COMMENT and HISTORY keyrecords. */
  if (fits_hdr2str(hFITS, 1, NULL, 0, &header, &nkeyrec, &status)) {
    throw FitsException("Cannot convert Header Fits to string", status);
  }

  /* Interpret the WCS keywords. */
  if ((status =wcspih(header, nkeyrec, WCSHDR_all, -3, &nreject, &nwcs, &wcs))) {
    throw CrsException("wcspih ERROR");
    // fprintf(stderr, "wcspih ERROR %d: %s.\n", status, wcshdr_errmsg[status]);
    // return 1;
  }

  fits_free_memory(header, &status);
  wcsset(wcs);
}

Sky::~Sky() {
  if (wcs) {
    wcsvfree(&nwcs, &wcs);
  }
}

std::vector<double> Sky::getGeoTransform() {

  const int ncoord = 4; // Performing the conversion for five positions
  const int nelem = 2;  // Input position specified by two dimensions

  // We also define our output arrays
  double imgcrd[nelem]; // You never use this
  double phi, theta;             // You never use these
  double skycrd[ncoord * nelem];  // save result
  double world[nelem];            // Output
  double pixcrd[nelem];           // input
  int stat[ncoord];              // 0 for success, 1 for error.
  int status;                    // 0 for success, non-zero for error

  const std::vector<int> naxis = getNaxis();

  /* NOTA BENE: FITS standard define pixel integers at the center of the pixel,
          0.5 must be subtract to have UpperLeft corner */
  double fov[ncoord * nelem] = {0.0, //lower left x 0
                                   0.0,//         y 1
                                   naxis[0], //lower right  x 2
                                   0, //                    y 3                                 
                                   0, // upper left         x 4
                                   naxis[1], //             y 5
                                   naxis[0], //upper right  x 6
                                   naxis[1]}; //            y 7

  // Perform the conversion
  for (int i = 0; i < ncoord*nelem; i=i+2) { 
    pixcrd[0] = fov[i];
    pixcrd[1] = fov[i+1];
    status = wcsp2s(wcs, 1, nelem, pixcrd, imgcrd, &phi, &theta, world, stat);

    if (status > 0 || stat[0] > 0) {
      throw CrsException("Error in 'wcsp2s()' - message:"+std::string(wcs_errmsg[status]));
    } else {
      status = 0;
      stat[0] = 0;
      skycrd[i] = world[0];
      skycrd[i+1] = world[1];
    }   
  } 

  // see world file on wikipedia
  // Xgeo = C + Xpixel*A + Yline*B
  // Ygeo = F + Xpixel*D + Yline*E
  // https://www.expertgps.com/tutorials/mapping-xy-sampling-data-to-real-world-gps-coordinates.asp
  double c = skycrd[0];
  double f = skycrd[1];  
  double dx = skycrd[2] - skycrd[0];
  double dy = skycrd[5] - skycrd[1];
  double scaleX_A = dx / naxis[0];
  double scaleY_E = dy / naxis[1];
  // for x = 0 => XUpperLeft = C + NAXIS2 * B
  //double screw_B = (skycrd[4] - c) / naxis[1];
  // for y = 0 => YLowerRight = F NAXIS1 * D
  //double screw_D = (skycrd[3] - f) / naxis[0]; 
  double screw_B = (skycrd[6] - c - naxis[0] * scaleX_A) / naxis[1];
  double screw_D = (skycrd[7] - f - naxis[1] * scaleY_E) / naxis[0];


  std::vector<double> geoTransform = {c, scaleX_A, screw_B,
                                      f, screw_D,  scaleY_E};
  return geoTransform;
} 

OGRSpatialReference Sky::createCrs() {
  OGRSpatialReference oSRS = OGRSpatialReference();
  createDatum(&oSRS);
  createCoordFrame(&oSRS);
  return oSRS;
}

void Sky::createDatum(OGRSpatialReference *oSRS) {
  std::string datumName;
  datumName.assign("D_");
  datumName.append(
      fits_read_key_as_str(hFITS, GlobalConstants::WCS_RADESYS, "ICRS"));
  oSRS->SetGeogCS(getCrsName().c_str(), datumName.c_str(), "celestial_sphere",
                  1.0, 0.0, "Reference_Meridian", 0.0, SRS_UA_DEGREE,
                  GlobalConstants::DEG2RAD);
}

void Sky::createCoordFrame(OGRSpatialReference *oSRS) {
  double falseEast = 0.0, falseNorth = 0.0;
  std::string projName;

  // Define projected CRS
  const std::string ctypeName = GlobalConstants::WCS_CTYPE + std::to_string(1);
  std::string ctype = fits_read_key_as_str(hFITS, ctypeName);
  char *ctypePtr = &ctype[0];
  char *pstr = strrchr(ctypePtr, '-');
  if (!pstr) {
    throw CrsException("Unknown projection.");
  } else {
    pstr += 1;
  }

  std::vector<double> crval = getCrval();

  /* Defining projection type
  Following http://www.gdal.org/ogr__srs__api_8h.html (GDAL)
  and
  http://www.aanda.org/component/article?access=bibcode&bibcode=&bibcode=2002A%2526A...395.1077CFUL
  (FITS)
  */

  /* Point Perspective / AZP projection */
  // if( strcmp(pstr,"AZP" ) == 0 ) {
  //        projName.assign("perspective_point_height_");
  //        oSRS.SetProjection(SRS_PP_PERSPECTIVE_POINT_HEIGHT);
  /* # appears to need height... maybe center lon/lat */
  // TODO ????
  // The orthographic projection is a perspective azimuthal projection centered
  // around a given latitude and longitude.
  // PV

  //}TODO : STG
  // TODO : SZP ?
  /* Gnomonic / TAN projection */
  if (strcmp(pstr, GlobalConstants::WCS_PROJ_TAN.c_str()) == 0) {
    projName.assign("Gnomonic_");
    oSRS->SetGnomonic(crval[1], crval[0], falseEast, falseNorth);

    /* Orthographic / SIN projection */
  } else if (strcmp(pstr, GlobalConstants::WCS_PROJ_SIN.c_str()) == 0) {
    projName.assign("Orthographic_");
    oSRS->SetOrthographic(crval[1], crval[0], falseEast, falseNorth);
    // check PV

    /* 	Azimuthal_Equidistant / ARC projection */
  } else if (strcmp(pstr, GlobalConstants::WCS_PROJ_ARC.c_str()) == 0) {
    projName.assign("Azimuthal_Equidistant_");
    oSRS->SetAE(crval[1], crval[0], falseEast, falseNorth);
  }
  // TODO : ZPN ?
  /* Lambert Azimuthal Equal Area / ZEA projection */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_ZEA.c_str()) == 0) {
    projName.assign("Lambert_Azimuthal_Equal_Area_");
    oSRS->SetLAEA(crval[1], crval[0], falseEast, falseNorth);
  }
  // TODO : AIRY - not supported by GDAL ?
  // TODO : CYP ?

  /* Cylindrical Equal Area / CEA projection */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_CEA.c_str()) == 0) {
    projName.assign("Cylindrical_Equal_Area_");
    oSRS->SetCEA(crval[1], crval[0], falseEast, falseNorth);
    double pv21 = fits_read_key_as_double(hFITS, "PV2_1", 1.0);
    oSRS->SetProjParm(SRS_PP_SCALE_FACTOR, pv21);
  }

  /* Equirectangular / CAR projection */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_CAR.c_str()) == 0) {
    projName.assign("Equirectangular_");
    /*
    The standard_parallel_1 defines where the local radius is calculated
    not the center of Y Cartesian system (which is latitude_of_origin)
    But FITS WCS only supports projections on the sphere
    we assume here that the local radius is the one computed at the projection
    center
    */
    oSRS->SetEquirectangular2(crval[1], crval[0], crval[1], falseEast,
                              falseNorth);
  }

  /* Mercator, Oblique (Hotine) Mercator, Transverse Mercator */
  /* Mercator / MER projection */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_MER.c_str()) == 0) {
    double scale = 1.0;
    projName.assign("Mercator_");
    oSRS->SetMercator(crval[1], crval[0], scale, falseEast, falseNorth);
  }

  /* Sinusoidal / SFL projection */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_SFL.c_str()) == 0) {
    projName.assign("Sinusoidal_");
    oSRS->SetSinusoidal(crval[0], falseEast, falseNorth);
  }

  /* Mollweide / MOL projection */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_MOL.c_str()) == 0) {
    projName.assign("Mollweide_");
    oSRS->SetMollweide(crval[0], falseEast, falseNorth);
  }

  /* Hammer-Aitoff / AIT projection */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_AIT.c_str()) == 0) {
    projName.assign("Aitoff_");
    oSRS->SetProjection(SRS_PT_AITOFF);
    oSRS->SetProjParm(SRS_PP_LONGITUDE_OF_CENTER, crval[0]);
  }

  /* TODO COP: Conic perspective - GDAL ? */

  /* Conic Equal Area / COE projection */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_COE.c_str()) == 0) {
    projName.assign("Albers_Conic_Equal_Area_");
    double pv21 = fits_read_key_as_double(hFITS, "PV2_1");
    double pv22 = fits_read_key_as_double(hFITS, "PV2_2");
    oSRS->SetACEA(pv21, pv22, crval[1], crval[0], falseEast, falseNorth);
  }

  /* Equidistant Conic / COD */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_COD.c_str()) == 0) {
    projName.assign("Albers_Conic_Equal_Area_");
    double pv21 = fits_read_key_as_double(hFITS, "PV2_1");
    double pv22 = fits_read_key_as_double(hFITS, "PV2_2");
    oSRS->SetEC(pv21, pv22, crval[1], crval[0], falseEast, falseNorth);
  }

  /* Lambert Conformal Conic 1SP / COO projection */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_COO.c_str()) == 0) {
    projName.assign("Lambert_Conformal_Conic_");
    double pv21 = fits_read_key_as_double(hFITS, "PV2_1");
    double pv22 = fits_read_key_as_double(hFITS, "PV2_2");
    oSRS->SetLCC(pv21, pv22, crval[1], crval[0], falseEast, falseNorth);
  }

  /* Bonne / BON projection */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_BON.c_str()) == 0) {
    projName.assign("Bonn_");
    double pv21 = fits_read_key_as_double(hFITS, "PV2_1");
    oSRS->SetBonne(pv21, crval[0], falseEast, falseNorth);
  }

  /* Polyconic / PCO projection */
  else if (strcmp(pstr, GlobalConstants::WCS_PROJ_PCO.c_str()) == 0) {
    projName.assign("Polyconic_");
    oSRS->SetPolyconic(crval[1], crval[0], falseEast, falseNorth);
  } else {
    throw CrsException("Unknown projection");
  }
  projName.append(getCrsName());
  oSRS->SetNode("PROJCS", projName.c_str());
  oSRS->SetAxes("PROJCS", getAxisName(1).c_str(), OAO_East,
                getAxisName(2).c_str(), OAO_North);

  // Fix axis  
  OGR_SRSNode *poAxis1 = oSRS->GetAttrNode("AXIS");
  std::string axis1Wkt;
  axis1Wkt.assign("AXIS[\"");
  axis1Wkt.append(getAxisName(1));
  axis1Wkt.append("\",east,ORDER[1],");
  axis1Wkt.append(getAxisUnit(1));
  axis1Wkt.append("]");
  poAxis1->SetValue(axis1Wkt.c_str());

  OGR_SRSNode *poAxis2 = oSRS->GetAttrNode("AXIS");
  std::string axis2Wkt;
  axis2Wkt.assign("AXIS[\"");
  axis2Wkt.append(getAxisName(2));
  axis2Wkt.append("\",north,ORDER[2],");
  axis2Wkt.append(getAxisUnit(2));
  axis2Wkt.append("]");
  poAxis2->SetValue(axis2Wkt.c_str());
  
}

std::string Sky::getAxisUnit(const int axisNumber) {
  std::string cunitName =
      GlobalConstants::WCS_CUNIT + std::to_string(axisNumber);
  std::string cunit = fits_read_key_as_str(hFITS, cunitName,GlobalConstants::WCS_DEG);
  std::string unit = "ANGLEUNIT[";
  std::stringstream ss;
  if (cunit == GlobalConstants::WCS_DEG) {
    ss << std::fixed << std::setprecision(16) << GlobalConstants::DEG2RAD;
    unit.append("\"degree\"," + ss.str() + ",ID[\"EPSG\",9102]]");

  } else if (cunit == GlobalConstants::WCS_ARCMIN) {
    ss << std::fixed << std::setprecision(16) << GlobalConstants::WCS_ARCMIN;
    unit.append("\"arcmin\"," + ss.str() + ",ID[\"EPSG\",9103]]");

  } else if (cunit == GlobalConstants::WCS_ARCSEC) {
    ss << std::fixed << std::setprecision(16) << GlobalConstants::WCS_ARCSEC;
    unit.append("\"arcsec\"," + ss.str() + ",ID[\"EPSG\",9104]]");

  } else if (cunit == GlobalConstants::WCS_MAS) {
    ss << std::fixed << std::setprecision(16) << GlobalConstants::WCS_MAS;
    unit.append("\"mas\"," + ss.str() + ",ID[\"EPSG\",1031]]");

  } else {
    throw CrsException("Unknow unit name");
  }
  return unit;
}

Equatorial::Equatorial(fitsfile *headerFits) : Sky(headerFits) {}

std::string Equatorial::getCrsName() { return GlobalConstants::COORD_NAME_EQ; }

std::string Equatorial::getAxisName(const int axisNumber) {
  std::string axis;
  std::string ctypeName =
      GlobalConstants::WCS_CTYPE + std::to_string(axisNumber);
  std::string ctype = fits_read_key_as_str(hFITS, ctypeName.c_str());
  std::string axisValue = ctype.substr(0, 3);
  if (axisValue == "RA-") {
    axis = GlobalConstants::WKT_AXIS_RA;
  } else if (axisValue == "DEC") {
    axis = GlobalConstants::WKT_AXIS_DEC;
  } else {
    throw CrsException("Unknow axis name");
  }
  return axis;
}

Ecliptic::Ecliptic(fitsfile *headerFits) : Sky(headerFits) {}

std::string Ecliptic::getCrsName() { return GlobalConstants::COORD_NAME_EC; }

std::string Ecliptic::getAxisName(const int axisNumber) {
  std::string axis;
  std::string ctypeName =
      GlobalConstants::WCS_CTYPE + std::to_string(axisNumber);
  std::string ctype = fits_read_key_as_str(hFITS, ctypeName.c_str());
  std::string axisValue = ctype.substr(0, 3);
  if (axisValue == "ELO") {
    axis = GlobalConstants::WKT_AXIS_LONG;
  } else if (axisValue == "ELA") {
    axis = GlobalConstants::WKT_AXIS_LAT;
  } else {
    throw CrsException("Unknow axis name");
  }
  return axis;
}

HelioEcliptic::HelioEcliptic(fitsfile *headerFits) : Sky(headerFits) {}

std::string HelioEcliptic::getCrsName() {
  return GlobalConstants::COORD_NAME_HE;
}

std::string HelioEcliptic::getAxisName(const int axisNumber) {
  std::string axis;
  std::string ctypeName =
      GlobalConstants::WCS_CTYPE + std::to_string(axisNumber);
  std::string ctype = fits_read_key_as_str(hFITS, ctypeName.c_str());
  std::string axisValue = ctype.substr(0, 3);
  if (axisValue == "HLO") {
    axis = GlobalConstants::WKT_AXIS_LONG;
  } else if (axisValue == "HLA") {
    axis = GlobalConstants::WKT_AXIS_LAT;
  } else {
    throw CrsException("Unknow axis name");
  }
  return axis;
}

Galactic::Galactic(fitsfile *headerFits) : Sky(headerFits) {}

std::string Galactic::getCrsName() { return GlobalConstants::COORD_NAME_GA; }

std::string Galactic::getAxisName(const int axisNumber) {
  std::string axis;
  std::string ctypeName =
      GlobalConstants::WCS_CTYPE + std::to_string(axisNumber);
  std::string ctype = fits_read_key_as_str(hFITS, ctypeName.c_str());
  std::string axisValue = ctype.substr(0, 3);
  if (axisValue == "GLO") {
    axis = GlobalConstants::WKT_AXIS_LONG;
  } else if (axisValue == "GLA") {
    axis = GlobalConstants::WKT_AXIS_LAT;
  } else {
    throw CrsException("Unknow axis name");
  }
  return axis;
}

SuperGalactic::SuperGalactic(fitsfile *headerFits) : Sky(headerFits) {}

std::string SuperGalactic::getCrsName() {
  return GlobalConstants::COORD_NAME_SG;
}

std::string SuperGalactic::getAxisName(const int axisNumber) {
  std::string axis;
  std::string ctypeName =
      GlobalConstants::WCS_CTYPE + std::to_string(axisNumber);
  std::string ctype = fits_read_key_as_str(hFITS, ctypeName.c_str());
  std::string axisValue = ctype.substr(0, 3);
  if (axisValue == "SLO") {
    axis = GlobalConstants::WKT_AXIS_LONG;
  } else if (axisValue == "SLA") {
    axis = GlobalConstants::WKT_AXIS_LAT;
  } else {
    throw CrsException("Unknow axis name");
  }
  return axis;
}

Crs *Crs::CreateFactory(fitsfile *hFITS) {
  Crs *crs;
  const std::string ctype1Name = GlobalConstants::WCS_CTYPE + std::to_string(1);
  const std::string ctype2Name = GlobalConstants::WCS_CTYPE + std::to_string(2);
  std::string ctype1 = fits_read_key_as_str(hFITS, ctype1Name.c_str());
  std::string ctype2 = fits_read_key_as_str(hFITS, ctype2Name.c_str());
  std::string coordFrameCtype1 = ctype1.substr(0, 4);
  std::string coordFrameCtype2 = ctype2.substr(0, 4);

  /* Equatorial */
  if (coordFrameCtype1 == "RA--" || coordFrameCtype2 == "RA--") {
    crs = new Equatorial(hFITS);

    /* Galactic */
  } else if (coordFrameCtype1 == "GLON-" || coordFrameCtype2 == "GLON-") {
    crs = new Galactic(hFITS);

    /* Ecliptic */
  } else if (coordFrameCtype1 == "ELON-" || coordFrameCtype2 == "ELON-") {
    crs = new Ecliptic(hFITS);

    /* Helioecliptic */
  } else if (coordFrameCtype1 == "HLON-" || coordFrameCtype2 == "HLON-") {
    crs = new HelioEcliptic(hFITS);

    /* Super Galactic */
  } else if (coordFrameCtype1 == "SLON-" || coordFrameCtype2 == "SLON-") {
    crs = new SuperGalactic(hFITS);

    /* Planets */
  } else if (coordFrameCtype1.substr(2,3) == "LN" || coordFrameCtype2.substr(2,3) == "LN") {
    crs = new Planet(hFITS);

  } else {
    throw CrsException("Unknown coordinates frame : "+coordFrameCtype1+" / "+coordFrameCtype2);
  }

  return crs;
}

CrsException::CrsException(const std::string &message)
    : m_msg(std::string("CRS error: ") + message) {}

CrsException::CrsException(const std::string &message, const int status)
    : m_msg(std::string("CRS error: ") + message), status_msg(status) {}

const char *CrsException::what() const throw() { return m_msg.c_str(); }

const int CrsException::getStatus() const throw() { return status_msg; }