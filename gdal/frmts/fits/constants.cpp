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

#include "constants.h"
#include <math.h>

const double GlobalConstants::DEG2RAD = M_PI / 180.0;
const double GlobalConstants::ARCMIN2RAD = DEG2RAD / 60.0;
const double GlobalConstants::ARCSEC2RAD = DEG2RAD / 3600.0;
const double GlobalConstants::MAS2RAD = ARCSEC2RAD / 1000.0;

const std::string GlobalConstants::WCS_DEG = "deg";
const std::string GlobalConstants::WCS_ARCMIN = "arcmin";
const std::string GlobalConstants::WCS_ARCSEC = "arcsec";
const std::string GlobalConstants::WCS_MAS = "mas";

const std::string GlobalConstants::WKT_METER = "metre";
const std::string GlobalConstants::WKT_DEGREE = "degree";

const std::string GlobalConstants::WCS_CTYPE = "CTYPE";
const std::string GlobalConstants::WCS_CTYPE1 = "CTYPE1";
const std::string GlobalConstants::WCS_CTYPE2 = "CTYPE2";
const std::string GlobalConstants::WCS_CUNIT = "CUNIT";
const std::string GlobalConstants::WCS_NAXIS = "NAXIS";
const std::string GlobalConstants::WCS_CRVAL = "CRVAL";
const std::string GlobalConstants::WCS_CRVAL1 = "CRVAL1";
const std::string GlobalConstants::WCS_CRVAL2 = "CRVAL2";
const std::string GlobalConstants::WCS_CRPIX = "CRPIX";
const std::string GlobalConstants::WCS_CRPIX1 = "CRPIX1";
const std::string GlobalConstants::WCS_CRPIX2 = "CRPIX2";
const std::string GlobalConstants::WCS_CDELT1 = "CDELT1";
const std::string GlobalConstants::WCS_CDELT2 = "CDELT2";
const std::string GlobalConstants::WCS_CD11 = "CD1_1";
const std::string GlobalConstants::WCS_CD12 = "CD1_2";
const std::string GlobalConstants::WCS_CD21 = "CD2_1";
const std::string GlobalConstants::WCS_CD22 = "CD2_2";
const std::string GlobalConstants::WCS_PC11 = "PC1_1";
const std::string GlobalConstants::WCS_PC12 = "PC1_2";
const std::string GlobalConstants::WCS_PC21 = "PC2_1";
const std::string GlobalConstants::WCS_PC22 = "PC2_2";
const std::string GlobalConstants::WCS_RADESYS = "RADESYS";
const std::string GlobalConstants::WCS_PROJ_AZP = "AZP";
const std::string GlobalConstants::WCS_PROJ_TAN = "TAN";
const std::string GlobalConstants::WCS_PROJ_SIN = "SIN";
const std::string GlobalConstants::WCS_PROJ_ARC = "ARC";
const std::string GlobalConstants::WCS_PROJ_ZEA = "ZEA";
const std::string GlobalConstants::WCS_PROJ_CEA = "CEA";
const std::string GlobalConstants::WCS_PROJ_CAR = "CAR";
const std::string GlobalConstants::WCS_PROJ_MER = "MER";
const std::string GlobalConstants::WCS_PROJ_SFL = "SFL";
const std::string GlobalConstants::WCS_PROJ_MOL = "MOL";
const std::string GlobalConstants::WCS_PROJ_AIT = "AIT";
const std::string GlobalConstants::WCS_PROJ_COE = "COE";
const std::string GlobalConstants::WCS_PROJ_COD = "COD";
const std::string GlobalConstants::WCS_PROJ_COO = "COO";
const std::string GlobalConstants::WCS_PROJ_BON = "BON";
const std::string GlobalConstants::WCS_PROJ_PCO = "PCO";
const std::string GlobalConstants::WCS_PROJ_STG = "STG";

const std::string GlobalConstants::FITS_A_RADIUS = "A_RADIUS";
const std::string GlobalConstants::FITS_B_RADIUS = "B_RADIUS";
const std::string GlobalConstants::FITS_C_RADIUS = "C_RADIUS";
const std::string GlobalConstants::FITS_OBJECT = "OBJECT";
const std::string GlobalConstants::FITS_BSCALE = "BSCALE";
const std::string GlobalConstants::FITS_BZERO = "BZERO";    
const std::string GlobalConstants::FITS_BLANK = "BLANK";

const std::string GlobalConstants::WKT_AXIS_RA = "ra";
const std::string GlobalConstants::WKT_AXIS_DEC = "dec";
const std::string GlobalConstants::WKT_AXIS_LONG = "longitude";
const std::string GlobalConstants::WKT_AXIS_LAT = "latitude";
const std::string GlobalConstants::WKT_PROJCS = "PROJCS";

const std::string GlobalConstants::COORD_NAME_EQ = "Equatorial";
const std::string GlobalConstants::COORD_NAME_EC = "Ecliptic";
const std::string GlobalConstants::COORD_NAME_HE = "HelioEcliptic";
const std::string GlobalConstants::COORD_NAME_GA = "Galactic";
const std::string GlobalConstants::COORD_NAME_SG = "SuperGalactic";
