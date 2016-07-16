SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config3.site

TARGET = Duffing.x
TARGET2 = Lorenz.x
TARGET3 = VDP.x
TARGET4 = nDuffing.x

SRCS = ticktock.cpp KL.cpp MCS.cpp GhanemSpanos.cpp AAPG.cpp Duffing.cpp Utils.cpp Utilsave.cpp KL_exp.cpp 
SRCS2 = ticktock.cpp KL.cpp MCS.cpp GhanemSpanos.cpp Lorenz.cpp UtilsLorenz.cpp Utils.cpp Utilsave.cpp 
SRCS3 = ticktock.cpp KL.cpp MCS.cpp GhanemSpanos.cpp AAPG.cpp VDP.cpp Utils.cpp Utilsave.cpp 
SRCS4 = ticktock.cpp KL.cpp nMCS.cpp nGhanemSpanos.cpp GhanemSpanos.cpp MCS.cpp nAAPG.cpp AAPG.cpp nDuffing.cpp Utils.cpp Utilsave.cpp KL_exp.cpp 

OBJ = $(SRCS:.cpp=.o)
OBJ2 = $(SRCS2:.cpp=.o)
OBJ3 = $(SRCS3:.cpp=.o)
OBJ4 = $(SRCS4:.cpp=.o)

INCS2 =  UtilsLorenz.h 

LIBINCDIR = $(TOPDIR)/src_cpp/include/
LIBBINDIR = $(TOPDIR)/src_cpp/lib/
INCDIRS   = -I. -I$(LIBINCDIR) -I$(SCINET_GSL_INC)

LIBS = -L$(LIBBINDIR) -lUQTk -lquad -luqtkmcmc -luqtktools -llbfgs -lcvode-2.6.0 -ldsfmt \
       -l$(LAPACK) -l$(SLATEC) -l$(BLAS) -lxmlutils -lexpat $(FCLIB) -L$(SCINET_GSL_LIB) -lgsl -lgslcblas
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
nduffing: $(TARGET4)

$(TARGET): $(OBJ) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET) $(OBJ) $(LIBS)

$(TARGET2): $(OBJ2) $(INCS2) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET2) $(OBJ2) $(LIBS)

$(TARGET3): $(OBJ3) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET3) $(OBJ3) $(LIBS)

$(TARGET4): $(OBJ4) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET4) $(OBJ4) $(LIBS)

test:
	g++ -Wall -O2 -D__wsu -I. -I../../src_cpp/include/ -I$(SCINET_GSL_INC) -o test.x test.cpp  -L$(SCINET_GSL_LIB) -lgsl -lgslcblas 

clean:
	rm -f $(OBJ) $(TARGET)*~ $(OBJ2) $(TARGET2) *~ $(OBJ3) $(TARGET3) *~ $(OBJ4) $(TARGET4)
	rm -f *.dat
	rm -f *.o *.x


