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

#ifndef _FITS_PLUGIN_CRS_
#define _FITS_PLUGIN_CRS_
#include "ogr_spatialref.h"
#include <fitsio.h>
#include <string>

/**
 * A Coordinate Reference System (CRS) contains two different main elements:
 * the coordinate reference frame and the coordinate system.
 *
 * The coordinate referenece frame defines how the CRS is related to origin 
 * (position and the date of the origin - equinox, epoch, date of observation) 
 * and the coordinate system describes how the coordinates is expressed in the 
 * coordinate reference frame (e.g cartesisan coordinates, spherical coordinates,
 * coordinates of a map projection).
 */
class Crs {

/*******************************************************************
 *                        Public
 *******************************************************************/

public:

  /**
   * Init the CRS
   * 
   */
  void init();

  /**
   * Returns the coordinates reference system description.
   * @return the coordinate reference system as a OGRSpatialReference
   */
  const OGRSpatialReference* getCrs();

  /**
   * Returns the geotransform matrix.
   * Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
   * Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
   * @return the geotransform matrix   
   */
  virtual std::vector<double> getGeoTransform() = 0;

  /**
   * Default destructor.   
   */
  virtual ~Crs() = default;

  /**
   * Factory to create a CRS based on keywords in FITS.
   * @param hFITS pointer to the fits file
   * @return the coordinate reference system
   */
  static Crs *CreateFactory(fitsfile *hFITS);

/*******************************************************************
 *                        protected
 *******************************************************************/
protected:
  /*
   * Pointer to the fits file.
   */
  fitsfile *hFITS;
  /**   
   * Number of axes in the HDU
   */
  int naxes;

  /**   
   * Constructs a new coordinate reference system based
   * on the pointer to the fitsfile.
   */
  Crs(fitsfile *hFITS);

  /**
   * Returns the CRS name.
   * @return the CRS name
   */
  virtual std::string getCrsName() = 0;

  /**
   * Returns the axis name based on the axis number.
   * @param axisNumber  axis number of the HDU with axisNumber >=1
   * @return the axis name
   */
  virtual std::string getAxisName(const int axisNumber) = 0;

  /**
   * Returns the axis unit based on the axis number.
   * @param axisNumber  axis number of the HDU with axisNumber >=1
   * @return the axis unit
   */  
  virtual std::string getAxisUnit(const int axisNumber) = 0;

  /**
   * Create the coordinate reference system
   * @return the CRS description
   */
  virtual OGRSpatialReference createCrs() = 0;

  /**
   * Returns the number of axis of the current HDU of the fits file
   * @param hFITS the pointer to the fits file
   * @return the number of axis
   */
  int getNaxes(fitsfile *hFITS) const;

  /**
   * Returns the array describing the number of pixels for each axis
   * @return the array describing the number of pixels for each axis
   */
  std::vector<int> getNaxis() const;

  /**
   * Returns the array describing the coordinates of the reference point in world coordinates
   * @return the array describing the coordinates of the reference point in world coordinates
   */  
  std::vector<double> getCrval() const;

  /**
   * Returns the array describing the coordinates of the reference point in pixel coordinates
   * @return the array describing the coordinates of the reference point in pixel coordinates
   */    
  std::vector<float> getCrpix() const;

  /**
   * Returns the array describing the rotation and the scale of the field of view
   * @return the array describing the rotation and the scale of the field of view
   */   
  std::vector<double> getCd() const;

/*******************************************************************
 *                        Private
 *******************************************************************/
private:
  /**
   * CRS description   
   */
  OGRSpatialReference crs;
};





/**
 * Exception for CRS 
 */
class CrsException : public std::exception {
  std::string m_msg;
  int status_msg;

public:
  CrsException(const std::string &message);
  CrsException(const std::string &message, const int status);
  virtual const char *what() const throw() override;
  virtual const int getStatus() const throw();
};





/**
 * Planet is a concrete implementation of the CRS.
 */
class Planet: public Crs {
/*******************************************************************
 *                        Public
 *******************************************************************/  
public:
  /**
   * Constructs a planet based on the pointer to the fitsfile
   * @param headerFits pointer to the fitsfile
   */
  Planet(fitsfile *headerFits);
  
  /**
   * Default destructor   
   */
  virtual ~Planet() = default;

  /**
   * Returns the geotransform matrix
   * @return the geotransform matrix
   */
  std::vector<double> getGeoTransform() override;	


/*******************************************************************
 *                        Protected
 *******************************************************************/
protected:
  /**
   * Returns the concrete implementation of the CRS name.
   * @return the CRS name
   */
  std::string getCrsName() override;

  /**
   * Returns the concrete axis name based on the axis number.
   * @param axisNumber  axis number of the HDU with axisNumber >=1
   * @return the axis name
   */  
  std::string getAxisName(const int axisNumber) override;

  /**
   * Returns the concrete axis unit based on the axis number.
   * @param axisNumber  axis number of the HDU with axisNumber >=1
   * @return the axis unit
   */    
  std::string getAxisUnit(const int axisNumber) override;

  /**
   * Create the coordinate reference system
   * @return the CRS description
   */  
  OGRSpatialReference createCrs() override;

/*******************************************************************
 *                        Private
 *******************************************************************/
private:
  /**
   * Create the datum based on the CRS description
   * @param oSRS pointer to the CRS description
   */
  void createDatum(OGRSpatialReference *oSRS);

  /**
   * Create the coordinates frame based on the CRS description
   * @param oSRS pointer to the CRS description
   */
  void createCoordFrame(OGRSpatialReference *oSRS);

  /**
   * Test the equality between two doubles.   
   */
  bool double_equals(double a, double b, double epsilon = 0.00001);
};





/**
 * abstract CRS for celestial sphere
 */
class Sky : public Crs {

/*******************************************************************
 *                        Public
 *******************************************************************/
public:
  /**
   * Returns the geotransform matrix
   * @return the geotransform matrix
   */
  std::vector<double> getGeoTransform() override;

/*******************************************************************
 *                        Protected
 *******************************************************************/
protected:
  /**
   * Constructs a CRS for the sky
   * @param headerFits Pointer to the fits file
   */
  Sky(fitsfile *headerFits);

  /**
   * Destructor  
   */
  virtual ~Sky();

  /**
   * Returns the concrete axis unit based on the axis number.
   * @param axisNumber  axis number of the HDU with axisNumber >=1
   * @return the axis unit
   */    
  std::string getAxisUnit(const int axisNumber) override;
  
  /**
   * Create the coordinate reference system
   * @return the CRS description
   */
  OGRSpatialReference createCrs() override;

/*******************************************************************
 *                        Private
 *******************************************************************/
private:
  /**
   * wcs structure   
   */
  struct wcsprm *wcs;
  int nwcs;

  /**
   * Create the datum based on the CRS description
   * @param oSRS pointer to the CRS description
   */
  void createDatum(OGRSpatialReference *oSRS);

  /**
   * Create the coordinates frame based on the CRS description
   * @param oSRS pointer to the CRS description
   */  
  void createCoordFrame(OGRSpatialReference *oSRS);
};




/**
 * CRS based on a equatorial coordinates
 * 
 */
class Equatorial : public Sky {

/*******************************************************************
 *                        Public
 *******************************************************************/  
public:
  Equatorial(fitsfile *headerFits);

/*******************************************************************
 *                        Protected
 *******************************************************************/
protected:
  std::string getCrsName() override;
  std::string getAxisName(const int axisNumber) override;
};




/**
 * CRS based on a galactic coordinates
 * 
 */
class Galactic : public Sky {

/*******************************************************************
 *                        Public
 *******************************************************************/
public:
  Galactic(fitsfile *headerFits);

/*******************************************************************
 *                        Protected
 *******************************************************************/
protected:
  std::string getCrsName() override;
  std::string getAxisName(const int axisNumber) override;
};




/**
 * @brief CRS based on a ecliptic coordinates
 * 
 */
class Ecliptic : public Sky {

/*******************************************************************
 *                        Public
 *******************************************************************/  
public:
  Ecliptic(fitsfile *headerFits);

/*******************************************************************
 *                        Protected
 *******************************************************************/
protected:
  std::string getCrsName() override;
  std::string getAxisName(const int axisNumber) override;
};




/**
 * @brief CRS based on a helioecliptic coordinates
 * 
 */
class HelioEcliptic : public Sky {

/*******************************************************************
 *                        Public
 *******************************************************************/
public:
  HelioEcliptic(fitsfile *headerFits);

/*******************************************************************
 *                        Protected
 *******************************************************************/
protected:
  std::string getCrsName() override;
  std::string getAxisName(const int axisNumber) override;
};




/**
 * CRS based on super galactic coordinates
 */
class SuperGalactic : public Sky {

/*******************************************************************
 *                        Public
 *******************************************************************/  
public:
  SuperGalactic(fitsfile *headerFits);

/*******************************************************************
 *                        Protected
 *******************************************************************/
protected:
  std::string getCrsName() override;
  std::string getAxisName(const int axisNumber) override;
};

#endif