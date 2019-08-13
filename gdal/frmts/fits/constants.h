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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

class GlobalConstants {
public:
  static const double DEG2RAD;
  static const double ARCMIN2RAD;
  static const double ARCSEC2RAD;
  static const double MAS2RAD;

  static const std::string WCS_DEG;
  static const std::string WCS_ARCMIN;
  static const std::string WCS_ARCSEC;
  static const std::string WCS_MAS;

  static const std::string WKT_METER;
  static const std::string WKT_DEGREE;

  static const std::string WCS_CTYPE;
  static const std::string WCS_CTYPE1;
  static const std::string WCS_CTYPE2;
  static const std::string WCS_CUNIT;
  static const std::string WCS_NAXIS;
  static const std::string WCS_CRVAL;
  static const std::string WCS_CRVAL1;
  static const std::string WCS_CRVAL2;
  static const std::string WCS_CRPIX;
  static const std::string WCS_CRPIX1;
  static const std::string WCS_CRPIX2;
  static const std::string WCS_CDELT1;
  static const std::string WCS_CDELT2;  
  static const std::string WCS_CD11;
  static const std::string WCS_CD12;
  static const std::string WCS_CD21;
  static const std::string WCS_CD22;
  static const std::string WCS_PC11;
  static const std::string WCS_PC12;
  static const std::string WCS_PC21;
  static const std::string WCS_PC22;  
  static const std::string WCS_RADESYS;
  static const std::string WCS_PROJ_AZP;
  static const std::string WCS_PROJ_TAN;
  static const std::string WCS_PROJ_SIN;
  static const std::string WCS_PROJ_ARC;
  static const std::string WCS_PROJ_ZEA;
  static const std::string WCS_PROJ_CEA;
  static const std::string WCS_PROJ_CAR;
  static const std::string WCS_PROJ_MER;
  static const std::string WCS_PROJ_SFL;
  static const std::string WCS_PROJ_MOL;
  static const std::string WCS_PROJ_AIT;
  static const std::string WCS_PROJ_COE;
  static const std::string WCS_PROJ_COD;
  static const std::string WCS_PROJ_COO;
  static const std::string WCS_PROJ_BON;
  static const std::string WCS_PROJ_PCO;
  static const std::string WCS_PROJ_STG;

  static const std::string FITS_A_RADIUS;
  static const std::string FITS_B_RADIUS;
  static const std::string FITS_C_RADIUS;
  static const std::string FITS_OBJECT;
  static const std::string FITS_BSCALE;
  static const std::string FITS_BZERO;    
  static const std::string FITS_BLANK;    

  static const std::string WKT_AXIS_RA;
  static const std::string WKT_AXIS_DEC;
  static const std::string WKT_AXIS_LONG;
  static const std::string WKT_AXIS_LAT;
  static const std::string WKT_PROJCS;

  static const std::string COORD_NAME_EQ;
  static const std::string COORD_NAME_EC;
  static const std::string COORD_NAME_HE;
  static const std::string COORD_NAME_GA;
  static const std::string COORD_NAME_SG;
};

#endif