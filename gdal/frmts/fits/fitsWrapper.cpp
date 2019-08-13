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
#include "fitsWrapper.h"
#include "cpl_string.h"

void fits_read_key_longstr(fitsfile *hFITS, const char* keyword, char** longString) {
  int status = 0;
  fits_read_key_longstr(hFITS, keyword, longString, nullptr, &status);
  if (status) {
    throw FitsException("Error while reading long string for key", status); 
  }    
}

void fits_read_keyn(fitsfile *hFITS, int keyNum, char* keyword, char* value) {
  int status = 0;
  fits_read_keyn(hFITS, keyNum, keyword, value, nullptr, &status);
  if (status) {
    throw FitsException("Couldn't determine image parameters of FITS file", status); 
  }  
}

void fits_get_img_param(fitsfile *hFITS, const int maxdim, int *bitpix, int *naxis, long *naxes) {
  int status = 0;
  fits_get_img_param(hFITS, maxdim, bitpix, naxis, naxes, &status);
  if (status) {
    throw FitsException("Couldn't determine image parameters of FITS file", status); 
  }
}

void fits_close_file(fitsfile *hFITS) {
  int status = 0;
  fits_close_file(hFITS, &status);
  if (status) {
    throw FitsException("Cannot close the FITS file", status); 
  }   
}

void fits_movabs_hdu(fitsfile *hFITS, int number) {
  int status = 0;  
  fits_movabs_hdu(hFITS, number, nullptr, &status);
  if (status) {
    throw FitsException("Couldn't move to first HDU in FITS file", status);        
  } 
}

void fits_write_key_longwarn(fitsfile *hFITS) {
  int status = 0;
  fits_write_key_longwarn(hFITS, &status); 
}

void fits_update_key_longstr(fitsfile *hFITS, std::string keyword, std::string value) {
  int status = 0;
  fits_update_key_longstr(hFITS, keyword.c_str(), value.c_str(), nullptr, &status);
  if (status) {
    throw FitsException("Cannot update " + keyword, status);
  }    
}          

void fits_update_key_string(fitsfile *hFITS, std::string keyword, std::string value){
  int status = 0;
  fits_update_key(hFITS, TSTRING, keyword.c_str(), &value[0], nullptr, &status);
  if (status) {
    throw FitsException("Cannot update " + keyword, status);
  }  
}

void fits_update_key_double(fitsfile *hFITS, std::string keyword, double value){
  int status = 0;
  fits_update_key(hFITS, TDOUBLE, keyword.c_str(), &value, nullptr, &status);
  if (status) {
    throw FitsException("Cannot update " + keyword, status);
  }  
}

double fits_read_key_as_double(fitsfile *hFITS, std::string keyword) {
  double value;
  int status = 0;
  fits_read_key(hFITS, TDOUBLE, keyword.c_str(), &value, nullptr, &status);
  if (status) {
    throw FitsException("Cannot get " + keyword, status);
  }
  return value;
}

double fits_read_key_as_double(fitsfile *hFITS, std::string keyword,
                               double defaultValue) {
  double value;
  try {
    value = fits_read_key_as_double(hFITS, keyword);
  } catch (const FitsException &e) {
    value = defaultValue;
  }
  return value;
}

float fits_read_key_as_float(fitsfile *hFITS, std::string keyword) {
  float value;
  int status = 0;
  fits_read_key(hFITS, TFLOAT, keyword.c_str(), &value, nullptr, &status);
  if (status) {
    throw FitsException("Cannot get " + keyword, status);
  }
  return value;
}

int fits_read_key_as_int(fitsfile *hFITS, std::string keyword) {
  int value;
  int status = 0;
  fits_read_key(hFITS, TINT, keyword.c_str(), &value, nullptr, &status);
  if (status) {
    throw FitsException("Cannot get " + keyword, status);
  }
  return value;
}

std::string fits_read_key_as_str(fitsfile *hFITS, std::string keyword) {
  // std::string value;
  // value.reserve(FLEN_VALUE);
  char value[81];
  int status = 0;
  fits_read_key(hFITS, TSTRING, keyword.c_str(), &value, nullptr, &status);
  // fits_read_key(hFITS, TSTRING, keyword.c_str(), &value[0], nullptr,
  // &status);
  if (status) {
    throw FitsException("Cannot get " + keyword, status);
  }
  return value;
}

std::string fits_read_key_as_str(fitsfile *hFITS, std::string keyword,
                                 std::string defaultValue) {
  std::string value;
  try {
    value = fits_read_key_as_str(hFITS, keyword);
  } catch (const FitsException &e) {
    value = defaultValue;
  }
  return value;
}

FitsException::FitsException(const std::string &message)
    : m_msg(std::string("FITS error: ") + message) {}

FitsException::FitsException(const std::string &message, const int status)
    : m_msg(std::string("FITS error: ") + message), status_msg(status) {}

const char *FitsException::what() const throw() { return m_msg.c_str(); }

const int FitsException::getStatus() const throw() { return status_msg; }
