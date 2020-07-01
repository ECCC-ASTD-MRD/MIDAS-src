Ce programme permet de sélectionner des observations qui sont dans un
domaine LAM spécifié par une grille d'un fichier standard RPN.  De
plus, il permet de splitter des observations pour les séparer en
petits domaines géographiques réguliers par rapport à une grille
définie par un champ dans un fichier standard.

Les observations peuvent être données dans un fichier BURP ou bien une
base de données SQLite.

Pour avoir une documentation complète, vous n'avez qu'à demander
l'option `-h` ou bien `--help` à l'appel du programme.

## Compilation

Pour compiler, il suffit de faire:
```bash
make splitobs_${ORDENV_PLAT}
```

## Tests unitaires

Un ensemble de tests est disponible pour vérifier que les
modifications apportées ne changent pas les résultats ou bien que ces
modifications apportent les changements attendus.  Dans ce dernier
cas, il faudrait alors mettre à jour les résultats.

Pour lancer les tests, on peut le faire interactivement avec la commande
```bash
./unittest
```

Par défaut, le programme utilisé sera celui du répertoire où le code
réside avec le nom `splitobs_${ORDENV_PLAT}`.  On peut donner un autre
programme avec la variable d'environnement `SPLITOBS`.

La commande `unittest` exécute par défaut les tests pour les BURP et
les SQLite.  On peut restreindre les tests pour faire seulement les
tests pour le BURP avec
```bash
./unittest burp
```
pour bien seulement les test pour SQLite avec
```bash
./unittest rdb
```

Un listing apparaîtra dans le répertoire de travail.

### Tests unitaires en batch

L'exécution des tests unitaires peut prendre un certain temps alors on
peut vouloir les lancer en batch.  On peut le faire avec la commande
suivante:
```bash
./unittest_batch
```
Comme pour `unittest`, on peut spécifier un programme avec la variable
d'environnement `SPLITOBS` de même qu'une famille de tests avec
l'argument `burp` ou `rdb`.  Par défaut, on lance les tests sur
`eccc-ppp4`, mais on peut modifier cela avec la variable
d'environnement `HOST_TO_SUBMIT`.
