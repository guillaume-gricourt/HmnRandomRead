TARGET = HmnRandomRead

DIRINC = ./include
DIRSRC = ./src
DIROBJ = ./obj

SRC = $(wildcard $(DIRSRC)/*.cpp)
OBJ = $(patsubst %.cpp,$(DIROBJ)/%.o,$(notdir $(SRC)))

CXX = g++
LIBS = -I $(DIRINC)
CFLAGS = -std=c++14 -Wall -Wextra -O3 -g3

############
## HTSLIB ##
############
DIR_HTSLIB = lib/htslib-1.16
DIREXEC_HTSLIB = $(DIR_HTSLIB)/bin
DIRHEADER_HTSLIB = $(DIR_HTSLIB)/htslib

STATICLIBS_HTSLIB = $(DIR_HTSLIB)/libhts.a
LIBS_HTSLIB = -I $(DIRHEADER_HTSLIB)
CFLAGS_HTSLIB = -lpthread -lz
#-llzma -lbz2 -lcurl -lcrypto -lncurses -lm
###########
## Rules ##
###########
all: $(TARGET) $(STATICLIBS_HTSLIB)

$(STATICLIBS_HTSLIB):
	@[ -d $(DIREXEC_HTSLIB) ] || mkdir -p $(DIREXEC_HTSLIB)
	cd $(DIR_HTSLIB) && \
	./configure --disable-bz2 --disable-lzma --enable-libcurl=no --prefix=$(PWD)/$(DIREXEC_HTSLIB) && \
	$(MAKE) && \
	$(MAKE) install

$(TARGET):$(OBJ) $(STATICLIBS_HTSLIB)
	$(CXX) $(OBJ) $(STATICLIBS_HTSLIB) $(LIBS) $(LIBS_HTSLIB) $(CFLAGS) $(CFLAGS_HTSLIB) -o $@

$(DIROBJ)/%.o:$(DIRSRC)/%.cpp
	@[ -d $(DIROBJ) ] || mkdir $(DIROBJ)
	$(CXX) $(LIBS) $(LIBS_HTSLIB) $(CFLAGS) -c $< -o $@

test:
	$(MAKE)
	cp $(TARGET) tests/
	python -m pytest tests
	$(RM) tests/$(TARGET)

#.PHONY:clean
clean:
	cd $(DIR_HTSLIB) && $(MAKE) clean && $(MAKE) distclean
	$(RM) -r $(DIREXEC_HTSLIB) $(DIR_HTSLIB)/lib
	$(RM) -r $(DIROBJ)
	$(RM) $(TARGET)
