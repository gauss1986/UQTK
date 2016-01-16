SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config3.site

TARGET = Duffing.x
TARGET2 = Lorenz.x
TARGET3 = VDP.x

SRCS = ticktock.cpp KL.cpp MCS.cpp GhanemSpanos.cpp AAPG.cpp Duffing.cpp Utils.cpp Utilsave.cpp 
SRCS2 = ticktock.cpp KL.cpp MCS.cpp GhanemSpanos.cpp Lorenz.cpp UtilsLorenz.cpp Utils.cpp Utilsave.cpp 
SRCS3 = ticktock.cpp KL.cpp MCS.cpp GhanemSpanos.cpp AAPG.cpp VDP.cpp Utils.cpp Utilsave.cpp 

OBJ = $(SRCS:.cpp=.o)
OBJ2 = $(SRCS2:.cpp=.o)
OBJ3 = $(SRCS3:.cpp=.o)

INCS2 =  UtilsLorenz.h 

LIBINCDIR = $(TOPDIR)/src_cpp/include/
LIBBINDIR = $(TOPDIR)/src_cpp/lib/
INCDIRS   = -I. -I$(LIBINCDIR)

LIBS = -L$(LIBBINDIR) -lUQTk -lquad -luqtkmcmc -luqtktools -llbfgs -lcvode-2.6.0 -ldsfmt \
       -l$(LAPACK) -l$(SLATEC) -l$(BLAS) -lxmlutils -lexpat $(FCLIB)
LDEP = $(LIBBINDIR)/libUQTk.a $(LIBBINDIR)/libquad.a $(LIBBINDIR)/liblbfgs.a $(LIBBINDIR)/libbcs.a \
       $(LIBBINDIR)/libuqtktools.a $(LIBBINDIR)/libcvode-2.6.0.a \
       $(LIBBINDIR)/libdsfmt.a $(LIBBINDIR)/lib$(LAPACK).a \
       $(LIBBINDIR)/lib$(SLATEC).a $(LIBBINDIR)/lib$(BLAS).a

DEFS = -D__$(FTNNAME)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -c $*.cpp

duffing: $(TARGET)
lorenz: $(TARGET2)
vdp: $(TARGET3)

$(TARGET): $(OBJ) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET) $(OBJ) $(LIBS)

$(TARGET2): $(OBJ2) $(INCS2) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET2) $(OBJ2) $(LIBS)

$(TARGET3): $(OBJ3) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET3) $(OBJ3) $(LIBS)

clean:
	rm -f $(OBJ) $(TARGET)*~ $(OBJ2) $(TARGET2) *~ $(OBJ3) $(TARGET3)
	rm -f *.dat



