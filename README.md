# Devoir 3 - LINMA1170

## Description

Ce projet implémente une simulation de déformation d’un solide en 2D par la méthode des éléments finis (FEM).  
Il permet de calculer l’évolution du déplacement et de la vitesse des nœuds d’un maillage au cours du temps, selon différents modèles de comportement mécanique.

## Compilation

Pour compiler le projet, utilisez simplement :

```bash
make
```

## Utilisation

L’exécutable principal est `deformation`.  
Voici la syntaxe d’utilisation :

```bash
./deformation <model> <lc> <T> <dt> <initial.txt> <final.txt> <time.txt> <I>
```

- `<model>` : nom du modèle à utiliser (voir le dossier `models/`)
- `<lc>` : taille caractéristique du maillage
- `<T>` : temps total de la simulation (en secondes)
- `<dt>` : pas de temps (en secondes)
- `<initial.txt>` : fichier d’entrée contenant les conditions initiales (déplacements et vitesses)
- `<final.txt>` : fichier de sortie pour les déplacements et vitesses finaux
- `<time.txt>` : fichier de sortie pour l’évolution temporelle du nœud I
- `<I>` : indice du nœud à suivre dans le temps

### Exemples

```bash
./deformation fork 1 10 0.01 ./data/init/initial_fork_1.0.txt ./data/results/final.txt ./data/results/time.txt 5
./deformation fork 1 2 0.1 ./data/init/initial_fork_1.0.txt ./data/results/final.txt ./data/results/time.txt 2
```

## Fichiers importants

- `src/main.c` : point d’entrée du programme
- `models/` : différents modèles mécaniques implémentés
- `data/init/` : exemples de fichiers de conditions initiales
- `data/results/` : résultats générés par les simulations