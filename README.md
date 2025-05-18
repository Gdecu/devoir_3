# Devoir 3 - LINMA1170 
DADI Elias & DE CUMONT Gaston

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
### Exemples

```bash
./deformation fork 1 10 0.01 ./data/init/initial_fork_1.0.txt ./data/results/final.txt ./data/results/time.txt 5
./deformation fork 1 2 0.1 ./data/init/initial_fork_1.0.txt ./data/results/final.txt ./data/results/time.txt 2
```

### Test des analyse numérique

Pour obtenir les résultats de nos analyses il faut run chaque fichier python qui se trouve dans le dossier : `py_plots/`, les fichiers .txt pour ceux-ci ce trouve dans `data/`.

Sauf pour l'animation: étant donné qu'il était un peu lourd de fournir les 201 fichier.txt dans le fichier .zip, nous les avons pas inclus. Pour obtenir les données de l'animation, activez la ligne 195 dans `main.c` et ensuite run `py_plots/anim.py`.

Pour obtenir les valeurs numériques (fichiers.txt) utilisez pour nous analyse, il faut activer les lignes suivants dans la fonction `main.c` :
- `anim_enrgy(argv[5], Ksp, Msp, u, v, n, T, dt, nbr_iter)`
- `convergence_complexity(Ksp, Msp, 2*n, T, argv[5])` 
- `run_stability_vs_dt(Ksp, Msp, n);`
- `run_stability_vs_bg(Ksp, Msp, n);`
et ensuite run chaque programmes pythons dans le dossier : `py_plot/`.

## Fichiers importants

- `src/main.c` : point d’entrée du programme
- `src/devoir_2.c` : contient le calcul de l'algorithme de Newmark et quelques autres fonctions utiles (récuperer les CI, stocké les états finaux, etc)
- `src/analyse.c` : contient toute les fonctions utilisé pour l'analyse de l'algorithme de Newmark
- `py_plot/` : contient les programmes .py pour plots les résultats des analyse de l'algorithme de Newmark
- `images/` : contient toute les images du rapport
- `data/init/` : exemples de fichiers de conditions initiales
- `data/` : contients tout les autres fichiers .txt utiliser pour l'analyse numérique : animation, conservation de l'énergie, complexité temporelle, convergence & stabilité