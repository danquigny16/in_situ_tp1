################################################################################
#partie définition de variable makefile

#compilateur
CC=gcc

#option
CFLAGS=-Wall -Wextra -std=c99 -O3

#dossiers
SRC_DIR=src
BUILD=build
GRAPHE=graphe


################################################################################
#partie compilation

all: $(BUILD)/driver

$(BUILD)/driver.o: $(SRC_DIR)/driver.c $(SRC_DIR)/util.h
	$(CC) -o $@ $(CFLAGS) -c $<

$(BUILD)/util.o: $(SRC_DIR)/util.c $(SRC_DIR)/util.h
	$(CC) -o $@ $(CFLAGS) -c $<

$(BUILD)/driver: $(BUILD)/driver.o $(BUILD)/util.o
	$(CC) $^ -o $@


################################################################################
#partie exécution

test: $(BUILD)/driver
	#tee permet de diriger le flux d'entrée vers plusieurs sorties, stdout + les fichiers spécifiés
	@$< | tee result.txt

result.txt: $(BUILD)/driver
	@$< | tee result.txt


$(BUILD)/result.txt: $(BUILD)/driver
	@$< | tee $@

data: $(BUILD)/result.txt
	@cat $< | grep "Performance obtenu pour des vecteurs de taille" | grep -v "unroll" | tr -s ' ' | cut -d ' ' -f 8,10 > $(GRAPHE)/produit_vect_graphe.txt
	@cat $< | grep "Performance obtenu pour des vecteurs de taille" | grep "unroll" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_vect_unroll_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "kij" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_kij_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "ijk" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_ijk_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "jik" | grep -v "unroll" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_jik_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "kji" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_kji_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "unroll" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_unroll_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "bloc" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_bloc_graphe.txt

graphe: data
	@python3 tracer_graphe.py $(GRAPHE)/produit_vect_graphe.txt $(GRAPHE)/produit_vect_unroll_graphe.txt $(GRAPHE)/produit_mat_ijk_graphe.txt $(GRAPHE)/produit_mat_jik_graphe.txt $(GRAPHE)/produit_mat_kij_graphe.txt $(GRAPHE)/produit_mat_kji_graphe.txt $(GRAPHE)/produit_mat_bloc_graphe.txt $(GRAPHE)/produit_mat_unroll_graphe.txt

################################################################################
#partie clean

.PHONY: clean clean_graphe clean_data clean_exec clean_all

clean:
	@rm -f $(BUILD)/*.o

clean_graphe:
	@rm -f $(GRAPHE)/*.png

clean_data:
	@rm -f $(GRAPHE)/*.txt

clean_exec: clean
	@rm -f $(BUILD)/driver

clean_result:
	@rm -f $(BUILD)/result.txt

clean_all: clean_graphe clean_data clean_exec clean_result
