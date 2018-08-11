#########################
# Compiler and standard #
#########################

CMP = icpc
STD = -std=c++11

#########
# Flags #
#########

CMPFLAGS = -fopenmp
OPTFLAGS = -O2 -ansi-alias -xAVX -DASSUME_ALIGNED
#OPTFLAGS = -O2 -ansi-alias -xCORE-AVX2 -DASSUME_ALIGNED
#OPTFLAGS = -O2 -ansi-alias -xCORE-AVX512 -qopt-zmm-usage=high -DASSUME_ALIGNED
REPFLAGS = -g -qopt-report=5 -qopt-report-filter="GSimulation.cpp,251-366"

FLAGS = $(STD) $(CMPFLAGS) $(OPTFLAGS) $(REPFLAGS)

#########
# Files #
#########

SRC = GSimulation.cpp main.cpp
OBJ = $(SRC:.cpp=.o)

EXE = ./nbody.x
ADV = ./adv-v7-openmp

#################
# GENERAL RULES #
#################

all: clean $(EXE)

run: $(EXE)
	$(info )
	$(info Running the CPU executable: )
	$(EXE)

$(EXE): $(OBJ)
	$(info )
	$(info Linking the CPU executable:)
	$(CMP) $(FLAGS) -o $(EXE) $(OBJ)

%.o: %.cpp
	$(info )
	$(info Compiling $< for CPU: )
	$(CMP) $(FLAGS) -o $@ -c $<

asm: GSimulation.cpp
	$(info )
	$(info Generating assembly code: )
	$(CMP) $(FLAGS) -S $^

clean:
	$(info )
	$(info Removing temporary files: )
	rm -f $(EXE) $(OBJ) *.s *.optrpt

#################
# ADVISOR RULES #
#################

survey:
	advixe-cl -collect survey -project-dir $(ADV) -- $(EXE)

roofline:
	advixe-cl -collect survey -project-dir $(ADV) -- $(EXE)
	advixe-cl -collect tripcounts -flop -project-dir $(ADV) -- $(EXE)

map:
	advixe-cl -collect map -mark-up-list=1 -project-dir $(ADV) -- $(EXE)

dependencies:
	advixe-cl -collect dependencies -mark-up-list=1 -project-dir $(ADV) -- $(EXE) 1000 1

open-gui:
	advixe-gui $(ADV)/$(ADV).advixeproj >/dev/null 2>&1 &

clean-results:
	rm -rf $(ADV)
