################################################################################
#partie définition de variable makefile

#compilateur
CC=gcc

#option
CFLAGS=-Wall -Wextra -std=c99

#dossiers
SRC_DIR=src
BUILD=build
GRAPHE=graphe


################################################################################
#partie compilation

all: driver

driver: $(BUILD)/driver.o $(BUILD)/util.o
	$(CC) $(CFLAGS) $^ -o test

$(BUILD)/driver.o: $(SRC_DIR)/driver.c
	$(CC) $(CFLAGS) -c $^ -o $@

$(BUILD)/util.o: $(SRC_DIR)/util.c
	$(CC) $(CFLAGS) -c $^ -o $@


################################################################################
#partie exécution

test: driver
	@./test

graphe: driver
	@./graphe.sh


################################################################################
#partie clean

.PHONY: clean clean_graphe clean_data clean_exec clean_all

clean:
	@rm -f $(BUILD)/*.o

clean_graphe:
	@rm -f $(GRAPHE)/*.png

clean_data:
	@rm -f donnees_graphe.txt

clean_exec: clean
	@rm -f test

clean_all: clean_graphe clean_data clean_exec
