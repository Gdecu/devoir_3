# Devoir 3 - LINMA1170 - DADI Elias & DE CUMONT Gaston

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
- `<initial.txt>` : fichier d’entrée contenant les conditions initiales (déplacements et vitesses)  (stocké dans `.\data\init\` )
- `<final.txt>` : fichier de sortie pour les déplacements et vitesses finaux
- `<time.txt>` : fichier de sortie pour l’évolution temporelle du nœud I
- `<I>` : indice du nœud à suivre dans le temps

### Exemples

```bash
./deformation fork 1 10 0.01 ./data/init/initial_fork_1.0.txt ./data/results/final.txt ./data/results/time.txt 5
./deformation fork 1 2 0.1 ./data/init/initial_fork_1.0.txt ./data/results/final.txt ./data/results/time.txt 2
```

### Test des analyse numérique

Pour obtenir les résultats de nos analyses il faut activer les fonctions `anim_enrgy(argv[5], Ksp, Msp, u, v, n, T, dt, nbr_iter)` et `convergence_complexity(Ksp, Msp, 2*n, T, argv[5])` qui sont en commentaire à la fin de la fonction main, et ensuite run chaque programmes pythons dans le dossier : `py_plot/`.

## Fichiers importants

- `src/main.c` : point d’entrée du programme
- `src/devoir_2.c` : contient le calcul de l'algorithme de Newmark et quelques autres fonctions utiles (récuperer les CI, stocké les états finaux, etc)
- `src/analyse.c` : contient toute les fonctions utilisé pour l'analyse de l'algorithme de Newmark
- `py_plot/` : contient les programmes .py pour plots les résultats des analyse de l'algorithme de Newmark
- `data/init/` : exemples de fichiers de conditions initiales
- `data/` : contients tout les autres fichiers .txt utiliser pour l'analyse numérique : animation, conservation de l'énergie, complexité temporelle, convergence, ...
- `models/` : différents modèles mécaniques implémentés
