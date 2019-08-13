
include ../../GDALmake.opt

OBJ	=	fitsdataset.o crs.o fitsWrapper.o constants.o

FITS_OPTS	=

CPPFLAGS	:=	 $(FITS_OPTS) $(CPPFLAGS)

default:	$(OBJ:.o=.$(OBJ_EXT))

clean:
	rm -f *.o

install-obj:	$(O_OBJ:.o=.$(OBJ_EXT))
