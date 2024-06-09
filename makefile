# Define variables
SRC_DIR = ./src/
SRC_FILES = $(wildcard $(SRC_DIR)/*.f)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.f, %.o, $(SRC_FILES))
EXECUTABLE = lcpfct2D

# Compiler
FC = gfortran
#FCFLAGS = -O2 -Wall

# Rules
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ_FILES)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: $(SRC_DIR)/%.f
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ_FILES) $(EXECUTABLE)

.PHONY: all clean

