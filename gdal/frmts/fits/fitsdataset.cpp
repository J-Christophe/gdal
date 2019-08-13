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
 * Copyright (c) 2019, Jean-Christophe Malapert <jean-christophe dot malapert at
 *cnes dot fr>
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

#include <iostream>

#include "constants.h"
#include "cpl_string.h"
#include "crs.h"
#include "fitsWrapper.h"
#include "gdal_frmts.h"
#include "gdal_pam.h"
#include "ogr_spatialref.h"

#include <cstring>
#include <fitsio.h>
#include <string.h>
#include <string>
CPL_CVSID("$Id$")

/************************************************************************/
/* ==================================================================== */
/*                              FITSDataset                             */
/* ==================================================================== */
/************************************************************************/

class FITSRasterBand;

class FITSDataset final : public GDALPamDataset {

  friend class FITSRasterBand;

public:
  ~FITSDataset();
  static int Identify(GDALOpenInfo *);
  static GDALDataset *Open(GDALOpenInfo *);
  static GDALDataset *Create(const char *pszFilename, int nXSize, int nYSize,
                             int nBands, GDALDataType eType,
                             char **papszParmList);
  const OGRSpatialReference *GetSpatialRef() const override;
  CPLErr SetSpatialRef(const OGRSpatialReference *poSRS) override;
  virtual CPLErr GetGeoTransform(double *) override;
  virtual CPLErr SetGeoTransform(double *) override;

private:
  fitsfile *hFITS;
  GDALDataType gdalDataType; // GDAL code for the image type
  int fitsDataType;          // FITS code for the image type
  bool isExistingFile;
  long highestOffsetWritten; // How much of image has been written
  bool bGeoTransformValid;
  double adfGeoTransform[6];
  OGRSpatialReference oSRS{};

  FITSDataset(); // Others should not call this constructor explicitly
  FITSDataset(GDALAccess eAccessFile);
  FITSDataset(GDALAccess eAccessFile, int xSize, int ySize);

  void Init(fitsfile *hFITS_, bool isExistingFile_);
  static void WriteFITSInfo(fitsfile *hFITS, const OGRSpatialReference *crs,
                            const double *adfGeoTransform);
  static bool is_number(const std::string& s);
  void setCrsInfos();
  void writeFitsKeywordsFromMetadata();
  static bool isIgnorableFITSHeader(const char *name);
  void setDataType(int bitpix, double offset);
  void setRasterDim(int naxis, long naxes[]);
  void fillDataItem(fitsfile *hFITS);
  static int getBitpixFrom(GDALDataType dataType);
};

/************************************************************************/
/* ==================================================================== */
/*                            FITSRasterBand                           */
/* ==================================================================== */
/************************************************************************/

class FITSRasterBand : public GDALPamRasterBand {

  friend class FITSDataset;

private:
  bool bHaveOffsetScale;
  double dfOffset;
  double dfScale;
  FITSDataset *poFDS;
  bool bNoDataSet;
  double dfNoDataValue;

public:
  FITSRasterBand(FITSDataset *, int);
  virtual ~FITSRasterBand();

  virtual CPLErr IReadBlock(int, int, void *) override;
  virtual CPLErr IWriteBlock(int, int, void *) override;

  virtual double GetNoDataValue(int *) override final;
  virtual CPLErr SetNoDataValue(double) override final;
  virtual CPLErr DeleteNoDataValue() override final;

  virtual double GetOffset(int *pbSuccess = nullptr) override final;
  virtual CPLErr SetOffset(double dfNewValue) override final;
  virtual double GetScale(int *pbSuccess = nullptr) override final;
  virtual CPLErr SetScale(double dfNewValue) override final;
};

bool FITSDataset::is_number(const std::string& s){
  char* p;
  bool result;  
  long converted = strtod(s.c_str(), &p);
  if (*p) {
    result = false;
  } else {
    result = true;
  }
  return result;
}

void FITSDataset::writeFitsKeywordsFromMetadata() {
  char **metaData = GetMetadata();
  int count = CSLCount(metaData);
  for (int i = 0; i < count; ++i) {
    const char *field = CSLGetField(metaData, i);
    if (strlen(field) == 0)
      continue;
    else {
      char *key = nullptr;
      const char *value = CPLParseNameValue(field, &key);
      // FITS keys must be less than 8 chars
      if (key != nullptr && strlen(key) <= 8 &&
          !FITSDataset::isIgnorableFITSHeader(key)) {
        // Although FITS provides support for different value
        // types, the GDAL Metadata mechanism works only with
        // string values. Prior to about 2003-05-02, this driver
        // would attempt to guess the value type from the metadata
        // value string amd then would use the appropriate
        // type-specific FITS keyword update routine. This was
        // found to be troublesome (e.g. a numeric version string
        // with leading zeros would be interpreted as a number
        // and might get those leading zeros stripped), and so now
        // the driver writes every value as a string. In practice
        // this is not a problem since most FITS reading routines
        // will convert from strings to numbers automatically, but
        // if you want finer control, use the underlying FITS
        // handle. Note: to avoid a compiler warning we copy the
        // const value string to a non const one.
        // char *valueCpy = CPLStrdup(value);
        std::string valueStr(value);
        std::string keyStr(key);
        if(FITSDataset::is_number(value)) {
          std::string::size_type sz;
          double valueNumber = std::stod (value,&sz);
          fits_update_key_double(hFITS, keyStr, valueNumber);
        } else {
          fits_update_key_longstr(hFITS, keyStr, valueStr);
        }
      }
      // Must free up key
      CPLFree(key);
    }
  }
}

/************************************************************************/
/*                          Identity()                           */
/************************************************************************/
int FITSDataset::Identify(GDALOpenInfo *poOpenInfo) {
  const char *fitsID = "SIMPLE  =                    T"; // Spaces important!
  size_t fitsIDLen = strlen(fitsID); // Should be 30 chars long
  int result;
  if ((size_t)poOpenInfo->nHeaderBytes < fitsIDLen)
    result = FALSE;
  else if (memcmp(poOpenInfo->pabyHeader, fitsID, fitsIDLen))
    result = FALSE;
  else
    result = TRUE;
  return result;
}

int FITSDataset::getBitpixFrom(GDALDataType dataType) {
  // Determine FITS type of image
  int bitpix;
  if (dataType == GDT_Byte) {
    bitpix = BYTE_IMG;
  } else if (dataType == GDT_UInt16) {
    bitpix = USHORT_IMG;
  } else if (dataType == GDT_Int16) {
    bitpix = SHORT_IMG;
  } else if (dataType == GDT_UInt32) {
    bitpix = ULONG_IMG;
  } else if (dataType == GDT_Int32) {
    bitpix = LONG_IMG;
  } else if (dataType == GDT_Float32)
    bitpix = FLOAT_IMG;
  else if (dataType == GDT_Float64)
    bitpix = DOUBLE_IMG;
  else {
    CPLError(CE_Failure, CPLE_AppDefined,
             "GDALDataType (%d) unsupported for FITS", dataType);
    throw FitsException("Unknown bitpix", CE_Failure);
  }
  return bitpix;
}

/************************************************************************/
/*                          FITSRasterBand()                           */
/************************************************************************/

FITSRasterBand::FITSRasterBand(FITSDataset *poDSIn, int nBandIn)
    : bHaveOffsetScale(false), dfOffset(0.0), dfScale(1.0), poFDS(poDSIn),
      bNoDataSet(false), dfNoDataValue(-9999.0) {
  poDS = poDSIn;
  nBand = nBandIn;
  eDataType = poDSIn->gdalDataType;
  nBlockXSize = poDSIn->nRasterXSize;
  nBlockYSize = 1;
}

/************************************************************************/
/*                          ~FITSRasterBand()                           */
/************************************************************************/

FITSRasterBand::~FITSRasterBand() { FlushCache(); }

/************************************************************************/
/*                             IReadBlock()                             */
/************************************************************************/

CPLErr FITSRasterBand::IReadBlock(CPL_UNUSED int nBlockXOff, int nBlockYOff,
                                  void *pImage) {
  // A FITS block is one row (we assume BSQ formatted data)
  FITSDataset *dataset = (FITSDataset *)poDS;
  fitsfile *hFITS = dataset->hFITS;
  int status = 0;

  // Since a FITS block is a whole row, nBlockXOff must be zero
  // and the row number equals the y block offset. Also, nBlockYOff
  // cannot be greater than the number of rows
  CPLAssert(nBlockXOff == 0);
  CPLAssert(nBlockYOff < nRasterYSize);

  // Calculate offsets and read in the data. Note that FITS array offsets
  // start at 1...
  LONGLONG offset =
      static_cast<LONGLONG>(nBand - 1) * nRasterXSize * nRasterYSize +
      static_cast<LONGLONG>(nBlockYOff) * nRasterXSize + 1;
  long nElements = nRasterXSize;

  // If we haven't written this block to the file yet, then attempting
  // to read causes an error, so in this case, just return zeros.
  if (!dataset->isExistingFile && offset > dataset->highestOffsetWritten) {
    memset(pImage, 0,
           nBlockXSize * nBlockYSize * GDALGetDataTypeSize(eDataType) / 8);
    return CE_None;
  }

  // Otherwise read in the image data
  fits_read_img(hFITS, dataset->fitsDataType, offset, nElements, nullptr,
                pImage, nullptr, &status);

  // Capture special case of non-zero status due to data range
  // overflow Standard GDAL policy is to silently truncate, which is
  // what CFITSIO does, in addition to returning NUM_OVERFLOW (412) as
  // the status.
  if (status == NUM_OVERFLOW)
    status = 0;

  CPLErr result;
  if (status) {
    CPLError(CE_Failure, CPLE_AppDefined,
             "Couldn't read image data from FITS file (%d).", status);
    result = CE_Failure;
  } else {
    result = CE_None;
  }

  return result;
}

/************************************************************************/
/*                            IWriteBlock()                             */
/*                                                                      */
/************************************************************************/

CPLErr FITSRasterBand::IWriteBlock(CPL_UNUSED int nBlockXOff, int nBlockYOff,
                                   void *pImage) {
  FITSDataset *dataset = (FITSDataset *)poDS;
  fitsfile *hFITS = dataset->hFITS;
  int status = 0;

  // Since a FITS block is a whole row, nBlockXOff must be zero
  // and the row number equals the y block offset. Also, nBlockYOff
  // cannot be greater than the number of rows

  // Calculate offsets and read in the data. Note that FITS array offsets
  // start at 1 at the bottom left...
  LONGLONG offset =
      static_cast<LONGLONG>(nBand - 1) * nRasterXSize * nRasterYSize +
      static_cast<LONGLONG>(nBlockYOff) * nRasterXSize + 1;
  long nElements = nRasterXSize;
  fits_write_img(hFITS, dataset->fitsDataType, offset, nElements, pImage,
                 &status);

  // Capture special case of non-zero status due to data range
  // overflow Standard GDAL policy is to silently truncate, which is
  // what CFITSIO does, in addition to returning NUM_OVERFLOW (412) as
  // the status.
  if (status == NUM_OVERFLOW)
    status = 0;

  // Check for other errors
  if (status) {
    CPLError(CE_Failure, CPLE_AppDefined,
             "Error writing image data to FITS file (%d).", status);
    return CE_Failure;
  }

  // When we write a block, update the offset counter that we've written
  if (offset > dataset->highestOffsetWritten)
    dataset->highestOffsetWritten = offset;

  return CE_None;
}

/************************************************************************/
/* ==================================================================== */
/*                             FITSDataset                             */
/* ==================================================================== */
/************************************************************************/

// Some useful utility functions

// Simple static function to determine if FITS header keyword should
// be saved in meta data.
constexpr int ignorableHeaderCount = 18;
static const char *const ignorableFITSHeaders[ignorableHeaderCount] = {
    "SIMPLE",  "BITPIX",   "NAXIS",    "NAXIS1", "NAXIS2", "NAXIS3",
    "END",     "XTENSION", "PCOUNT",   "GCOUNT", "EXTEND", "CONTINUE",
    "COMMENT", "",         "LONGSTRN", "BZERO",  "BSCALE", "BLANK"};

bool FITSDataset::isIgnorableFITSHeader(const char *name) {
  for (int i = 0; i < ignorableHeaderCount; ++i) {
    if (strcmp(name, ignorableFITSHeaders[i]) == 0)
      return true;
  }
  return false;
}

/************************************************************************/
/*                            FITSDataset()                            */
/************************************************************************/

FITSDataset::FITSDataset()
    : hFITS(nullptr), gdalDataType(GDT_Unknown), fitsDataType(0),
      isExistingFile(false), highestOffsetWritten(0),
      bGeoTransformValid(false) {
  oSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
}

FITSDataset::FITSDataset(GDALAccess eAccessFile) : FITSDataset() {
  eAccess = eAccessFile;
}

FITSDataset::FITSDataset(GDALAccess eAccessFile, int xSize, int ySize)
    : FITSDataset(eAccessFile) {
  nRasterXSize = xSize;
  nRasterYSize = ySize;
}

/************************************************************************/
/*                           ~FITSDataset()                            */
/************************************************************************/

FITSDataset::~FITSDataset() {

  if (!hFITS) {
    return;
  }

  if (eAccess == GA_Update) {

    try {
      fits_write_key_longwarn(hFITS);
      writeFitsKeywordsFromMetadata();

      // Access to the source file (typically when doing a gdal_translate from
      // tif to fits)
      GDALRasterBand *poSrcBand = GDALPamDataset::GetRasterBand(1);
      int pbSuccess;

      // Writing nodata value
      if (gdalDataType != GDT_Float32 && gdalDataType != GDT_Float64) {
        fits_update_key_double(hFITS, GlobalConstants::FITS_BLANK.c_str(),
                               poSrcBand->GetNoDataValue(&pbSuccess));
      }

      // Writing Scale and offset if defined
      double dfScale = poSrcBand->GetScale(&pbSuccess);
      if (pbSuccess) {
        fits_update_key_double(hFITS, GlobalConstants::FITS_BSCALE.c_str(),
                               dfScale);
      }
      double dfOffset = poSrcBand->GetOffset(&pbSuccess);
      if (pbSuccess) {
        fits_update_key_double(hFITS, GlobalConstants::FITS_BZERO.c_str(),
                               dfOffset);
      }

      // Copy georeferencing info to PAM if the profile is not FITS
      GDALPamDataset::SetSpatialRef(GDALPamDataset::GetSpatialRef());

      // Write geographic info
      try {
        double geoTransform[6];
        GDALPamDataset::GetGeoTransform(geoTransform);

        FITSDataset::WriteFITSInfo(hFITS, GetSpatialRef(), geoTransform);
      } catch (const FitsException &e) {
        fits_report_error(stderr, e.getStatus());
        CPLError(CE_Failure, CPLE_AppDefined, e.what(), e.getStatus());
      } 

    } catch (const FitsException &e) {
      fits_report_error(stderr, e.getStatus());
      CPLError(CE_Failure, CPLE_AppDefined, e.what(), e.getStatus());
    }
  }

  // Make sure we flush the raster cache before we close the file!
  FlushCache();

  // Close the FITS handle
  fits_close_file(hFITS);
}

// Set the data type according to the bitpix and offset
void FITSDataset::setDataType(int bitpix, double offset) {
  // Determine data type and nodata value if BLANK keyword is absent
  switch (bitpix) {
  case BYTE_IMG:
    gdalDataType = GDT_Byte;
    fitsDataType = TBYTE;
    break;
  case SHORT_IMG:
    if (offset == 32768.) {
      gdalDataType = GDT_UInt16;
      fitsDataType = TUSHORT;
    } else {
      gdalDataType = GDT_Int16;
      fitsDataType = TSHORT;
    }
    break;
  case LONG_IMG:
    if (offset == 2147483648.) {
      gdalDataType = GDT_UInt32;
      fitsDataType = TUINT;
    } else {
      gdalDataType = GDT_Int32;
      fitsDataType = TINT;
    }
    break;
  case FLOAT_IMG:
    gdalDataType = GDT_Float32;
    fitsDataType = TFLOAT;
    break;
  case DOUBLE_IMG:
    gdalDataType = GDT_Float64;
    fitsDataType = TDOUBLE;
    break;
  default:
    throw FitsException("Unknown datatype", CE_Failure);
  }
}

// Set image dimensions - we assume BSQ ordering
void FITSDataset::setRasterDim(int naxis, long naxes[]) {
  if (naxis == 2) {
    nRasterXSize = static_cast<int>(naxes[0]);
    nRasterYSize = static_cast<int>(naxes[1]);
    nBands = 1;
  } else if (naxis == 3) {
    nRasterXSize = static_cast<int>(naxes[0]);
    nRasterYSize = static_cast<int>(naxes[1]);
    nBands = static_cast<int>(naxes[2]);
  } else {
    throw FitsException("Cannot the image dimensions", CE_Failure);
  }
}

// Read header information from file and use it to set metadata
// This process understands the CONTINUE standard for long strings.
// We don't bother to capture header names that duplicate information
// already captured elsewhere (e.g. image dimensions and type)
void FITSDataset::fillDataItem(fitsfile *hFITS) {
  char key[100];
  char value[100];

  int nKeys = 0;
  int nMoreKeys = 0;
  int status = 0;

  fits_get_hdrspace(hFITS, &nKeys, &nMoreKeys, &status);
  for (int keyNum = 1; keyNum <= nKeys; keyNum++) {
    fits_read_keyn(hFITS, keyNum, key, value);

    if (strcmp(key, "END") == 0) {
      // We should not get here in principle since the END
      // keyword shouldn't be counted in nKeys, but who knows.
      break;
    } else if (FITSDataset::isIgnorableFITSHeader(key)) {
      // Ignore it
    } else { // Going to store something, but check for long strings etc

      // Strip off leading and trailing quote if present
      char *newValue = value;
      if (value[0] == '\'' && value[strlen(value) - 1] == '\'') {
        newValue = value + 1;
        value[strlen(value) - 1] = '\0';
      }

      // Check for long string
      if (strrchr(newValue, '&') == newValue + strlen(newValue) - 1) {
        // Value string ends in "&", so use long string conventions
        char *longString = nullptr;
        fits_read_key_longstr(hFITS, key, &longString, nullptr, &status);
        // Note that read_key_longstr already strips quotes
        SetMetadataItem(key, longString);
        free(longString);
      } else { // Normal keyword
        SetMetadataItem(key, newValue);
      }
    }
  }
}

/************************************************************************/
/*                           Init()                                     */
/************************************************************************/

void FITSDataset::Init(fitsfile *hFITS_, bool isExistingFile_) {

  hFITS = hFITS_;
  isExistingFile = isExistingFile_;
  highestOffsetWritten = 0;

  // Move to the primary HDU
  fits_movabs_hdu(hFITS, 1);

  // Get the image info for this dataset (note that all bands in a FITS dataset
  // have the same type)
  int bitpix;
  int naxis;
  const int maxdim = 3;
  long naxes[maxdim];
  double offset;

  // Get dimension and datatype of the image.
  fits_get_img_param(hFITS, maxdim, &bitpix, &naxis, naxes);

  // BZERO is not mandatory offset defaulted to 0 if BZERO is missing
  offset =
      fits_read_key_as_double(hFITS, GlobalConstants::FITS_BZERO.c_str(), 0.0);

  bool bNoDataSet;
  double dfNoDataValue;
  try {
    dfNoDataValue =
        fits_read_key_as_double(hFITS, GlobalConstants::FITS_BLANK.c_str());
    bNoDataSet = true;
  } catch (const FitsException &e) {
    bNoDataSet = false;
  };

  // Set both gdal and fits datatype
  setDataType(bitpix, offset);

  // Set raster dimension
  setRasterDim(naxis, naxes);

  // Create the bands
  double dfScale =
      fits_read_key_as_double(hFITS, GlobalConstants::FITS_BSCALE.c_str(), 1.);
  double dfOffset =
      fits_read_key_as_double(hFITS, GlobalConstants::FITS_BZERO.c_str(), 0.);
  for (int i = 0; i < nBands; ++i) {
    FITSRasterBand *fitsRasterBand = new FITSRasterBand(this, i + 1);
    fitsRasterBand->SetScale(dfScale);
    fitsRasterBand->SetOffset(dfOffset);
    if(bNoDataSet) {
      fitsRasterBand->SetNoDataValue(dfNoDataValue);
    }
    SetBand(i + 1, fitsRasterBand);
  }

  // Fill metadata in dataset
  fillDataItem(hFITS);
}

/************************************************************************/
/*                                Open()                                */
/************************************************************************/

GDALDataset *FITSDataset::Open(GDALOpenInfo *poOpenInfo) {

  if (!Identify(poOpenInfo))
    return nullptr;

  // Get access mode and attempt to open the file
  int status = 0;
  fitsfile *hFITS = nullptr;
  if (poOpenInfo->eAccess == GA_ReadOnly)
    fits_open_file(&hFITS, poOpenInfo->pszFilename, READONLY, &status);
  else
    fits_open_file(&hFITS, poOpenInfo->pszFilename, READWRITE, &status);

  GDALDataset *result;
  if (status) {
    CPLError(CE_Failure, CPLE_AppDefined,
             "Error while opening FITS file %s (%d).\n",
             poOpenInfo->pszFilename, status);
    fits_close_file(hFITS, &status);
    result = nullptr;
  } else {
    // Create a FITSDataset object and initialize it from the FITS handle
    FITSDataset *dataset = new FITSDataset(poOpenInfo->eAccess);

    // Set up the description
    dataset->SetDescription(poOpenInfo->pszFilename);

    try {
      dataset->Init(hFITS, true);
      /* -------------------------------------------------------------------- */
      /*      Initialize any information.                                 */
      /* -------------------------------------------------------------------- */
      dataset->setCrsInfos();
      dataset->TryLoadXML();

      /* -------------------------------------------------------------------- */
      /*      Check for external overviews.                                   */
      /* -------------------------------------------------------------------- */
      dataset->oOvManager.Initialize(dataset, poOpenInfo->pszFilename,
                                     poOpenInfo->GetSiblingFiles());
      result = dataset;
    } catch (const FitsException &e) {
      fits_report_error(stderr, e.getStatus());
      CPLError(CE_Failure, CPLE_AppDefined, e.what(), e.getStatus());
      delete dataset;
      result = nullptr;
    }
  }

  return result;
}

/************************************************************************/
/*                               Create()                               */
/*                                                                      */
/*      Create a new FITS file.                                         */
/************************************************************************/

GDALDataset *FITSDataset::Create(const char *pszFilename, int nXSize,
                                 int nYSize, int nBands, GDALDataType eType,
                                 CPL_UNUSED char **papszParmList) {
  int status = 0;

  // No creation options are defined. The BSCALE/BZERO options were
  // removed on 2002-07-02 by Simon Perkins because they introduced
  // excessive complications and didn't really fit into the GDAL
  // paradigm.
  // 2018 - BZERO BSCALE keywords are now set using SetScale() and
  // SetOffset() functions

  if (nXSize < 1 || nYSize < 1 || nBands < 1) {
    CPLError(CE_Failure, CPLE_AppDefined,
             "Attempt to create %dx%dx%d raster FITS file, but width, height "
             "and bands"
             " must be positive.",
             nXSize, nYSize, nBands);

    return nullptr;
  }

  // Determine FITS type of image
  int bitpix = FITSDataset::getBitpixFrom(eType);

  // Create the file - to force creation, we prepend the name with '!'
  CPLString extFilename("!");
  extFilename += pszFilename;
  fitsfile *hFITS = nullptr;
  fits_create_file(&hFITS, extFilename, &status);
  if (status) {
    CPLError(CE_Failure, CPLE_AppDefined,
             "Couldn't create FITS file %s (%d).\n", pszFilename, status);
    return nullptr;
  }

  // Now create an image of appropriate size and type
  long naxes[3] = {nXSize, nYSize, nBands};
  int naxis = (nBands == 1) ? 2 : 3;
  fits_create_img(hFITS, bitpix, naxis, naxes, &status);

  // Check the status
  if (status) {
    CPLError(CE_Failure, CPLE_AppDefined,
             "Couldn't create image within FITS file %s (%d).", pszFilename,
             status);
    fits_close_file(hFITS, &status);
    return nullptr;
  }

  FITSDataset *dataset = new FITSDataset(GA_Update, nXSize, nYSize);
  dataset->eAccess = GA_Update;
  dataset->SetDescription(pszFilename);

  GDALDataset *result;
  try {
    dataset->Init(hFITS, false);
    result = dataset;
  } catch (const FitsException &e) {
    fits_report_error(stderr, e.getStatus());
    CPLError(CE_Failure, CPLE_AppDefined, e.what(), e.getStatus());
    delete dataset;
    result = nullptr;
  }
  return result;
}

/************************************************************************/
/*                          WriteFITSInfo()                          */
/************************************************************************/

void FITSDataset::WriteFITSInfo(fitsfile *hFITS, const OGRSpatialReference *crs,
                                const double *adfGeoTransform)

{
  int status = 0;

  const double PI = std::atan(1.0) * 4;
  const double DEG2RAD = PI / 180.;

  double falseEast = 0;
  double falseNorth = 0;

  double cfactor, mres, mapres, UpperLeftCornerX, UpperLeftCornerY;
  double crpix1, crpix2;

  /* -------------------------------------------------------------------- */
  /*      Write out projection definition.                                */
  /* -------------------------------------------------------------------- */
  if (crs != nullptr) {
    // Set according to coordinate system (thanks to Trent Hare - USGS)
    std::string object, ctype1, ctype2;

    const char *target = crs->GetAttrValue("DATUM", 0);
    if (target) {    
      if (strstr(target, "Moon")) {
        object.assign("Moon");
        ctype1.assign("SE");
        ctype2.assign("SE");
      } else if (strstr(target, "Mercury")) {
        object.assign("Mercury");
        ctype1.assign("ME");
        ctype2.assign("ME");
      } else if (strstr(target, "Venus")) {
        object.assign("Venus");
        ctype1.assign("VE");
        ctype2.assign("VE");
      } else if (strstr(target, "Mars")) {
        object.assign("Mars");
        ctype1.assign("MA");
        ctype2.assign("MA");
      } else if (strstr(target, "Jupiter")) {
        object.assign("Jupiter");
        ctype1.assign("JU");
        ctype2.assign("JU");
      } else if (strstr(target, "Saturn")) {
        object.assign("Saturn");
        ctype1.assign("SA");
        ctype2.assign("SA");
      } else if (strstr(target, "Uranus")) {
        object.assign("Uranus");
        ctype1.assign("UR");
        ctype2.assign("UR");
      } else if (strstr(target, "Neptune")) {
        object.assign("Neptune");
        ctype1.assign("NE");
        ctype2.assign("NE");
      } else if (strstr(target, "Earth")) {
        object.assign("Earth");
        ctype1.assign("EA");
        ctype2.assign("EA");
      } else {
        return;
      }

      char *cstrobj = new char[object.length() + 1];
      std::strcpy(cstrobj, object.c_str());

      fits_update_key(hFITS, TSTRING, "OBJECT", cstrobj, nullptr, &status);
    }

    double aradius = crs->GetSemiMajor();
    double bradius = aradius;
    double cradius = crs->GetSemiMinor();

    cfactor = aradius * DEG2RAD;

    fits_update_key_double(hFITS, GlobalConstants::FITS_A_RADIUS.c_str(),
                           aradius);
    fits_update_key_double(hFITS, GlobalConstants::FITS_B_RADIUS.c_str(),
                           bradius);
    fits_update_key_double(hFITS, GlobalConstants::FITS_C_RADIUS.c_str(),
                           cradius);

    const char *unit = crs->GetAttrValue("UNIT", 0);

    ctype1.append("LN-");
    ctype2.append("LT-");

    // strcat(ctype1a, "PX-");
    // strcat(ctype2a, "PY-");

    std::string fitsproj;
    const char *projection = crs->GetAttrValue("PROJECTION", 0);
    double centlon = 0, centlat = 0;

    if (projection) {
      
      if (strstr(projection, "Sinusoidal")) {
        fitsproj.assign(GlobalConstants::WCS_PROJ_SFL.c_str());
        centlon = crs->GetProjParm("central_meridian", 0, nullptr);
      } else if (strstr(projection, "Equirectangular")) {
        fitsproj.assign(GlobalConstants::WCS_PROJ_CAR.c_str());
        centlat = crs->GetProjParm("standard_parallel_1", 0, nullptr);
        centlon = crs->GetProjParm("central_meridian", 0, nullptr);
      } else if (strstr(projection, "Orthographic")) {
        fitsproj.assign(GlobalConstants::WCS_PROJ_SIN.c_str());
        centlat = crs->GetProjParm("standard_parallel_1", 0, nullptr);
        centlon = crs->GetProjParm("central_meridian", 0, nullptr);
      } else if (strstr(projection, "Mercator_1SP") ||
                 strstr(projection, "Mercator")) {
        fitsproj.assign(GlobalConstants::WCS_PROJ_MER.c_str());
        centlat = crs->GetProjParm("standard_parallel_1", 0, nullptr);
        centlon = crs->GetProjParm("central_meridian", 0, nullptr);
      } else if (strstr(projection, "Polar_Stereographic") ||
                 strstr(projection, "Stereographic_South_Pole") ||
                 strstr(projection, "Stereographic_North_Pole")) {
        fitsproj.assign(GlobalConstants::WCS_PROJ_STG.c_str());
        centlat = crs->GetProjParm("latitude_of_origin", 0, nullptr);
        centlon = crs->GetProjParm("central_meridian", 0, nullptr);
      } 

      /*
                      #Transverse Mercator is supported in FITS via specific MER
         parameters. # need some more testing... #if
         EQUAL(mapProjection,"Transverse_Mercator"): #    mapProjection = "MER"
                      #    centLat = hSRS.GetProjParm('standard_parallel_1')
                      #    centLon = hSRS.GetProjParm('central_meridian')
                      #    TMscale = hSRS.GetProjParm('scale_factor')
                      #    #Need to research when TM actually applies false
         values #    #but planetary is almost always 0.0 #    falseEast =
         hSRS.GetProjParm('false_easting') #    falseNorth =
         hSRS.GetProjParm('false_northing')
      */

      ctype1.append(fitsproj);
      ctype2.append(fitsproj);

      fits_update_key_string(hFITS, GlobalConstants::WCS_CTYPE1.c_str(),
                             ctype1);
      fits_update_key_string(hFITS, GlobalConstants::WCS_CTYPE2.c_str(),
                             ctype2);
    }

    UpperLeftCornerX = adfGeoTransform[0] - falseEast;
    UpperLeftCornerY = adfGeoTransform[3] - falseNorth;

    if (centlon > 180.) {
      centlon = centlon - 180.;
    }

    if (strstr(unit, GlobalConstants::WKT_METER.c_str())) {
      // convert degrees/pixel to m/pixel
      mapres = 1. / adfGeoTransform[1];    // mapres is pixel/meters
      mres = adfGeoTransform[1] / cfactor; // mres is deg/pixel
      crpix1 = -(UpperLeftCornerX * mapres) + centlon / mres + 0.5;
      // assuming that center latitude is also the origin of the coordinate
      // system: this is not always true.
      // More generic implementation coming soon
      crpix2 = (UpperLeftCornerY * mapres) + 0.5; // - (centlat / mres);
    } else if (strstr(unit, GlobalConstants::WKT_DEGREE.c_str())) {
      // convert m/pixel to pixel/degree
      mapres = 1. / adfGeoTransform[1] / cfactor; // mapres is pixel/deg
      mres = adfGeoTransform[1];                  // mres is meters/pixel
      crpix1 = -(UpperLeftCornerX * mres) + centlon / mapres + 0.5;
      // assuming that center latitude is also the origin of the coordinate
      // system: this is not always true.
      // More generic implementation coming soon
      crpix2 = (UpperLeftCornerY * mres) + 0.5; // - (centlat / mapres);
    } else {
      throw new FitsException("Unit is not supprted", 10000);
    }

    /// Write WCS CRPIXia CRVALia CTYPEia here
    fits_update_key_double(hFITS, GlobalConstants::WCS_CRVAL1.c_str(), centlon);
    fits_update_key_double(hFITS, GlobalConstants::WCS_CRVAL2.c_str(), centlat);
    fits_update_key_double(hFITS, GlobalConstants::WCS_CRPIX1.c_str(), crpix1);
    fits_update_key_double(hFITS, GlobalConstants::WCS_CRPIX2.c_str(), crpix2);

    /* -------------------------------------------------------------------- */
    /*      Write the transform.                                            */
    /* -------------------------------------------------------------------- */

    /// Write WCS CDELTia and PCi_ja here

    double cd[4];
    cd[0] = adfGeoTransform[1] / cfactor;
    cd[1] = adfGeoTransform[2] / cfactor;
    cd[2] = adfGeoTransform[4] / cfactor;
    cd[3] = adfGeoTransform[5] / cfactor;

    double pc[4];
    pc[0] = 1.;
    pc[1] = cd[1] / cd[0];
    pc[2] = cd[2] / cd[3];
    pc[3] = -1.;

    fits_update_key_double(hFITS, GlobalConstants::WCS_CDELT1.c_str(), cd[0]);
    fits_update_key_double(hFITS, GlobalConstants::WCS_CDELT2.c_str(), cd[3]);
    fits_update_key_double(hFITS, GlobalConstants::WCS_PC11.c_str(), pc[0]);
    fits_update_key_double(hFITS, GlobalConstants::WCS_PC12.c_str(), pc[1]);
    fits_update_key_double(hFITS, GlobalConstants::WCS_PC21.c_str(), pc[2]);
    fits_update_key_double(hFITS, GlobalConstants::WCS_PC22.c_str(), pc[3]);
  }
}

/************************************************************************/
/*                          GetSpatialRef()                             */
/************************************************************************/

const OGRSpatialReference *FITSDataset::GetSpatialRef() const {
  return oSRS.IsEmpty() ? nullptr : &oSRS;
}

/************************************************************************/
/*                           SetSpatialRef()                            */
/************************************************************************/

CPLErr FITSDataset::SetSpatialRef(const OGRSpatialReference *poSRS) {
  if (poSRS == nullptr || poSRS->IsEmpty()) {
    oSRS.Clear();
  } else {
    oSRS = *poSRS;
    oSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
  }

  return CE_None;
}

/************************************************************************/
/*                          GetGeoTransform()                           */
/************************************************************************/

CPLErr FITSDataset::GetGeoTransform(double *padfTransform) {
  memcpy(padfTransform, adfGeoTransform, sizeof(double) * 6);
  CPLErr status;
  if (!bGeoTransformValid)
    status = CE_Failure;
  else
    status = CE_None;
  return status;
}

/************************************************************************/
/*                          SetGeoTransform()                           */
/************************************************************************/

CPLErr FITSDataset::SetGeoTransform(double *padfTransform) {
  memcpy(adfGeoTransform, padfTransform, sizeof(double) * 6);
  bGeoTransformValid = true;
  return CE_None;
}

/************************************************************************/
/*                             GetOffset()                              */
/************************************************************************/

double FITSRasterBand::GetOffset(int *pbSuccess) {
  if (pbSuccess)
    *pbSuccess = bHaveOffsetScale;
  return dfOffset;
}

/************************************************************************/
/*                             SetOffset()                              */
/************************************************************************/

CPLErr FITSRasterBand::SetOffset(double dfNewValue) {
  bHaveOffsetScale = true;
  dfOffset = dfNewValue;
  return CE_None;
}

/************************************************************************/
/*                              GetScale()                              */
/************************************************************************/

double FITSRasterBand::GetScale(int *pbSuccess) {
  if (pbSuccess)
    *pbSuccess = bHaveOffsetScale;
  return dfScale;
}

/************************************************************************/
/*                              SetScale()                              */
/************************************************************************/

CPLErr FITSRasterBand::SetScale(double dfNewValue) {
  bHaveOffsetScale = true;
  dfScale = dfNewValue;
  return CE_None;
}

/************************************************************************/
/*                           GetNoDataValue()                           */
/************************************************************************/

double FITSRasterBand::GetNoDataValue(int *pbSuccess) {
  double result;
  if (bNoDataSet) {
    if (pbSuccess)
      *pbSuccess = TRUE;

    result = dfNoDataValue;
  } else {
    result = GDALPamRasterBand::GetNoDataValue(pbSuccess);
  }

  return result;
}

/************************************************************************/
/*                           SetNoDataValue()                           */
/************************************************************************/

CPLErr FITSRasterBand::SetNoDataValue(double dfNoData) {
  bNoDataSet = true;
  dfNoDataValue = dfNoData;
  return CE_None;
}

/************************************************************************/
/*                        DeleteNoDataValue()                           */
/************************************************************************/

CPLErr FITSRasterBand::DeleteNoDataValue() {
  bNoDataSet = false;
  dfNoDataValue = -9999.0;
  return CE_None;
}

/************************************************************************/
/*                     setCrsInfos()                                   */
/************************************************************************/

void FITSDataset::setCrsInfos() {
  try {
    Crs *crs = Crs::CreateFactory(hFITS);
    crs->init();
    SetSpatialRef(crs->getCrs());
    std::vector<double> geoTransform = crs->getGeoTransform();
    SetGeoTransform(geoTransform.data());
  } catch (const CrsException &e) {
    CPLError(CE_Failure, CPLE_AppDefined, e.what(), 1);
  } catch (const FitsException &e) {
    fits_report_error(stderr, e.getStatus());
    CPLError(CE_Failure, CPLE_AppDefined, e.what(), e.getStatus());
  }
}

/************************************************************************/
/*                          GDALRegister_FITS()                         */
/************************************************************************/

void GDALRegister_FITS() {
  if (GDALGetDriverByName("FITS") != nullptr)
    return;

  GDALDriver *poDriver = new GDALDriver();

  poDriver->SetDescription("FITS");
  poDriver->SetMetadataItem(GDAL_DCAP_RASTER, "YES");
  poDriver->SetMetadataItem(GDAL_DMD_LONGNAME,
                            "Flexible Image Transport System");
  poDriver->SetMetadataItem(GDAL_DMD_HELPTOPIC, "frmt_various.html#FITS");
  poDriver->SetMetadataItem(GDAL_DMD_CREATIONDATATYPES,
                            "Byte UInt16 Int16 UInt32 Int32 Float32 Float64");
  poDriver->SetMetadataItem(GDAL_DMD_EXTENSION, "fits");
  poDriver->SetMetadataItem(GDAL_DMD_MIMETYPE, "application/fits");

  poDriver->pfnOpen = FITSDataset::Open;
  poDriver->pfnCreate = FITSDataset::Create;
  poDriver->pfnIdentify = FITSDataset::Identify;
  poDriver->pfnCreateCopy = nullptr;

  GetGDALDriverManager()->RegisterDriver(poDriver);
}
