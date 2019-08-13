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
#ifndef _FITS_PLUGIN_FITSWRAPPER_
#define _FITS_PLUGIN_FITSWRAPPER_

#include <exception>
#include <fitsio.h>
#include <string>

void fits_read_key_longstr(fitsfile *hFITS, const char* keyword, char** longString);

void fits_read_keyn(fitsfile *hFITS, int keyNum, char* keyword, char* value);

void fits_get_img_param(fitsfile *hFITS, const int maxdim, int *bitpix, int *naxis, long *naxes);

void fits_close_file(fitsfile *hFITS);

void fits_movabs_hdu(fitsfile *hFITS, int number);

void fits_write_key_longwarn(fitsfile *hFITS);

void fits_update_key_longstr(fitsfile *hFITS, std::string keyword, std::string value);

void fits_update_key_string(fitsfile *hFITS, std::string keyword, std::string value);

void fits_update_key_double(fitsfile *hFITS, std::string keyword, double value);

double fits_read_key_as_double(fitsfile *hFITS, std::string keyword);

double fits_read_key_as_double(fitsfile *hFITS, std::string keyword,
                               double defaultValue);

float fits_read_key_as_float(fitsfile *hFITS, std::string keyword);

int fits_read_key_as_int(fitsfile *hFITS, std::string keyword);

std::string fits_read_key_as_str(fitsfile *hFITS, std::string keyword);

std::string fits_read_key_as_str(fitsfile *hFITS, std::string keyword,
                                 std::string defaultValue);

/*
Exception for Fits.
 */
class FitsException : public std::exception {
  std::string m_msg;
  int status_msg;

public:
  FitsException(const std::string &message);
  FitsException(const std::string &message, const int status);
  virtual const char *what() const throw() override;
  virtual const int getStatus() const throw();
};
#endif