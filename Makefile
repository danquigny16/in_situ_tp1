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
	$(CC) $(CFLAGS) $^ -o driver

$(BUILD)/driver.o: $(SRC_DIR)/driver.c
	$(CC) $(CFLAGS) -c $^ -o $@

$(BUILD)/util.o: $(SRC_DIR)/util.c
	$(CC) $(CFLAGS) -c $^ -o $@


################################################################################
#partie exécution

test: driver
	@./driver

graphe: driver
	@./driver | grep Performance | tr -s ' ' | cut -d ' ' -f 8,10 > donnees_graphe.txt
	@python3 tracer_graphe.py donnees_graphe.txt


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
