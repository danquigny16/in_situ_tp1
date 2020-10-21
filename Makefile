################################################################################
#partie définition de variable makefile

#compilateur
CC=gcc

#bibliothèque statique
AR=ar

#option
CFLAGS=-Wall -Wextra -std=c99 -O3
ARFLAGS=-rc

#dossiers
SRC_DIR=src
BUILD=build
LIB=lib
GRAPHE=graphe


################################################################################
#partie compilation

all: $(BUILD)/driver

$(BUILD)/driver.o: $(SRC_DIR)/driver.c $(SRC_DIR)/util.h
	$(CC) -o $@ $(CFLAGS) -c $<

$(BUILD)/util.o: $(SRC_DIR)/util.c $(SRC_DIR)/util.h
	$(CC) -o $@ $(CFLAGS) -c $<

$(BUILD)/myblas.o: $(SRC_DIR)/myblas.c $(SRC_DIR)/util.h $(SRC_DIR)/myblas.h
	$(CC) -o $@ $(CFLAGS) -c $<

$(BUILD)/mylapack.o: $(SRC_DIR)/mylapack.c $(SRC_DIR)/util.h $(SRC_DIR)/mylapack.h
	$(CC) -o $@ $(CFLAGS) -c $<

$(LIB)/libmyblas.a: $(BUILD)/util.o $(BUILD)/myblas.o
	$(AR) $(ARFLAGS) $@ $^

$(LIB)/libmylapack.a: $(BUILD)/util.o $(BUILD)/mylapack.o
	$(AR) $(ARFLAGS) $@ $^

$(BUILD)/driver: $(BUILD)/driver.o $(BUILD)/util.o $(LIB)/libmylapack.a $(LIB)/libmyblas.a
	$(CC) $(BUILD)/driver.o $(BUILD)/util.o -L$(LIB)/ -lmyblas -lmylapack -o $@ 

################################################################################
#partie exécution

test: $(BUILD)/driver
	@$<

$(BUILD)/result.txt: $(BUILD)/driver
	@$< > $@

data: $(BUILD)/result.txt
	@cat $< | grep "Performance obtenu pour des vecteurs de taille" | grep -v "unroll" | tr -s ' ' | cut -d ' ' -f 8,10 > $(GRAPHE)/produit_vect_graphe.txt
	@cat $< | grep "Performance obtenu pour des vecteurs de taille" | grep "unroll" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_vect_unroll_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "kij" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_kij_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "ijk" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_ijk_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "jik" | grep -v "unroll" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_jik_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "kji" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_kji_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "unroll" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_unroll_graphe.txt
	@cat $< | grep "Performance obtenu pour des matrices de taille" | grep "bloc" | tr -s ' ' | cut -d ' ' -f 8,12 > $(GRAPHE)/produit_mat_bloc_graphe.txt

$(GRAPHE): data
	@python3 tracer_graphe.py $(GRAPHE)/produit_vect_graphe.txt $(GRAPHE)/produit_vect_unroll_graphe.txt $(GRAPHE)/produit_mat_ijk_graphe.txt $(GRAPHE)/produit_mat_jik_graphe.txt $(GRAPHE)/produit_mat_kij_graphe.txt $(GRAPHE)/produit_mat_kji_graphe.txt $(GRAPHE)/produit_mat_bloc_graphe.txt $(GRAPHE)/produit_mat_unroll_graphe.txt

################################################################################
#partie clean

.PHONY: clean clean_graphe clean_data clean_exec clean_all

clean:
	@rm -f $(BUILD)/*.o $(LIB)/*.a $(BUILD)/driver

clean_graphe:
	@rm -f $(GRAPHE)/*.png

clean_data:
	@rm -f $(GRAPHE)/*.txt

clean_exec: clean
	@rm -f $(BUILD)/driver

clean_result:
	@rm -f $(BUILD)/result.txt

clean_all: clean_graphe clean_data clean_exec clean_result
