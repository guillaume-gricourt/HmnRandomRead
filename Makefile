TARGET = HmnRandomRead

DIRINC = ./include
DIRSRC = ./src
DIROBJ = ./obj
DIRLIB=./lib

SRC = $(wildcard $(DIRSRC)/*.cpp)
OBJ = $(patsubst %.cpp,$(DIROBJ)/%.o,$(notdir $(SRC)))

CXX = g++
LIBS = -I $(DIRINC) -I lib/htslib-1.16
CXXFLAGS = -std=c++14 -Wall -Wextra -O3 -g3 -lhts

# Rules
all: $(TARGET)

$(TARGET):$(OBJ)
	$(CXX) $(OBJ) $(LIBS) $(CXXFLAGS) -o $@

$(DIROBJ)/%.o:$(DIRSRC)/%.cpp
	@[ -d $(DIROBJ) ] || mkdir $(DIROBJ)
	$(CXX) $(LIBS) $(CXXFLAGS) -c $< -o $@

test:
	$(MAKE)
	cp $(TARGET) tests/
	python -m pytest tests
	$(RM) tests/$(TARGET)

#.PHONY:clean
clean:
	$(RM) -r $(DIROBJ)
	$(RM) $(TARGET)
