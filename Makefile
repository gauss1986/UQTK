SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config3.site

TARGET = Duffing.x

SRCS = ticktock.cpp KL.cpp MCS.cpp GhanemSpanos.cpp AAPG.cpp Duffing.cpp UtilsDuffing.cpp 

OBJ = $(SRCS:.cpp=.o)

INCS =  UtilsDuffing.h 

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
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -c $*.cpp -pg

all: $(TARGET)

$(TARGET): $(OBJ) $(INCS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET) $(OBJ) $(LIBS) -pg

clean:
	rm -f $(OBJ) $(TARGET)*~
	rm -f *.dat



