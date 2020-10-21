# TP1 Traitement des données in-situ : HPC + traitement des données in-situ

### Auteur

Jérôme Faure  
Antoine Danquigny  
Pierre Jiang

### Utilisation

make : compile le TP  
make test : lance les tests  
make data : recupère les données du tests et les parse pour chacun de nos graphe  
make graphe : construit le graphe avec les données du test puis le met dans le dossier "graphes"

make clean : efface les fichiers objets  
make clean_graphe :  efface les graphes  
make clean_data : efface les fichiers contenant les données des graphes  
make clean_exec : efface l'exécutable et les fichiers objets  
make clean_all : efface l'exécutable, les fichiers objets, les fichiers contenant les données du graphes et les graphes

### Architecture

build : contient les fichiers objets  
graphe : contient les graphes qu'on trace  
lib : contient nos bibliothèques  
src : contient nos sources
