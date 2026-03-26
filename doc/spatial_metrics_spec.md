# CoverM — Métriques spatiales de couverture

## À qui s'adresse ce document ?

Aux utilisateurs de CoverM qui veulent comprendre les nouvelles métriques spatiales
et les options associées, afin de mieux interpréter leurs résultats de détection de MAGs.

Aucune connaissance du code source n'est nécessaire.

---

## Table des matières

1. [Le problème : pourquoi les métriques existantes ne suffisent pas](#1-le-problème)
2. [Les concepts clés : îlots et gaps](#2-les-concepts-clés)
3. [Nouvelle métrique : `islands_per_mbp`](#3-islands_per_mbp)
4. [Nouvelle métrique : `max_gap`](#4-max_gap)
5. [Nouvelle métrique : `gap_fraction`](#5-gap_fraction)
6. [Nouvelle option : `--min-island-length`](#6---min-island-length)
7. [Nouveau side-output : `--coverage-profile`](#7---coverage-profile)
8. [Utilisation combinée](#8-utilisation-combinée)
9. [Grille de lecture combinée des métriques](#9-grille-de-lecture)
10. [Ce que ces métriques ne disent pas](#10-limites)
11. [Cas limites et interprétation](#11-cas-limites)
12. [Glossaire](#12-glossaire)

---

## 1. Le problème

### Pourquoi ajouter des métriques à CoverM ?

CoverM donne déjà `mean` et `covered_fraction`, qui résument la **quantité** de signal.
Pour les cas évidents (couverture très forte ou très faible), ces métriques suffisent.

Mais dans les **cas ambigus**, elles échouent. Prenons deux MAGs dans un même échantillon :

| MAG | mean | covered_fraction |
|-----|------|------------------|
| MAG_A | 5.0 | 0.30 |
| MAG_B | 5.0 | 0.30 |

Avec les métriques classiques, **ces deux MAGs sont indistinguables**. Même profondeur
moyenne, même fraction couverte.

Pourtant, la réalité est très différente :
- **MAG_A** est un **vrai positif** — l'organisme est partiellement présent, avec un signal
  réparti de manière continue sur les contigs couverts.
- **MAG_B** est un **faux positif** — le signal provient de reads d'autres organismes qui
  mappent sur des régions conservées (rRNA, gènes housekeeping, etc.). Ce n'est pas un bug
  d'alignement : c'est du cross-mapping biologiquement plausible, mais insuffisant pour
  conclure à la présence du MAG.

Ce qui les distingue, c'est la **géographie** du signal :

| MAG | islands_per_mbp | max_gap | gap_fraction |
|-----|-----------------|---------|-------------|
| MAG_A | 2.1 | 800 bp | 0.05 |
| MAG_B | 18.7 | 15 000 bp | 0.65 |

MAG_A a un signal **continu** (peu d'îlots, petits gaps, presque pas de trous internes).
MAG_B a un signal **fragmenté** (beaucoup de petits îlots, un trou de 15 kb, 65% de gaps
dans le span interne).

Les métriques spatiales séparent ce que `mean` et `covered_fraction` ne peuvent pas
distinguer.

Dans la pratique, les différences sont souvent moins tranchées que dans cet exemple,
d'où l'intérêt de combiner plusieurs métriques (voir la
[grille de lecture combinée](#9-grille-de-lecture)).

### Que proposent les nouvelles métriques ?

Trois métriques qui décrivent la **structure spatiale** du signal de couverture :

| Métrique | Question à laquelle elle répond |
|----------|-------------------------------|
| `islands_per_mbp` | Le signal est-il fragmenté en beaucoup de petits morceaux ? |
| `max_gap` | Y a-t-il un grand trou dans le signal ? |
| `gap_fraction` | Quelle proportion du span interne est constituée de trous ? |

> Les métriques existantes (`mean`, `covered_fraction`) restent utiles pour les cas évidents.
> Les métriques spatiales apportent un pouvoir de discrimination supplémentaire dans les
> cas ambigus — là où la quantité de signal ne suffit pas à conclure.

---

## 2. Les concepts clés

### Qu'est-ce qu'un îlot de couverture ?

Un **îlot** est un segment contigu de bases couvertes par au moins un read.
C'est une zone où le signal de mapping est présent, sans interruption.

```
Contig :  --------XXXXXXXXX--------XXXXX--------
                  ^^^^^^^^^        ^^^^^
                  îlot 1           îlot 2
```

Un vrai MAG tend à présenter un signal relativement continu et peu fragmenté, avec
**peu d'îlots, longs** — même si la profondeur n'est pas parfaitement uniforme partout
(régions multicopies, biais GC, profondeur variable).

Un faux positif tend à avoir **beaucoup d'îlots, courts** (un par gène conservé ou
région partagée).

### Qu'est-ce qu'un gap ?

Un **gap interne** est une zone sans couverture **entre deux îlots**, au sein d'un
même contig.

```
Contig :  --------XXXXXXXXX--------XXXXX--------
          ^^^^^^^^                       ^^^^^^^^
          pas un gap                     pas un gap
          (avant le 1er îlot)            (après le dernier)

                           ^^^^^^^^
                           GAP INTERNE
                           (entre deux îlots)
```

**Seuls les gaps internes comptent.** Les zones non couvertes aux extrémités du contig
(avant le premier îlot ou après le dernier) ne sont pas des gaps — elles peuvent refléter
un bord de contig mal assemblé ou une couverture qui s'atténue naturellement, et ne sont
pas informatives pour la détection.

### Pourquoi travailler uniquement à l'intérieur des contigs ?

Un MAG est généralement constitué de plusieurs contigs. L'ordre de ces contigs est
**arbitraire** — il ne reflète pas une continuité biologique. Les cassures entre contigs
sont des artefacts d'assemblage, pas des trous dans le signal de mapping.

Donc :
- Les gaps **entre** contigs sont **ignorés** (ils ne disent rien sur la qualité du signal)
- Les gaps **au sein** d'un contig sont **informatifs** (un trou dans un contig couvert
  signale une discontinuité réelle du signal)

### Que se passe-t-il pour les contigs sans aucune couverture ?

Ils sont **ignorés** dans les métriques spatiales. Un contig entièrement non couvert
ne contient ni îlot, ni gap interne. Son absence de signal est déjà captée par
`covered_fraction`.

### Qu'est-ce que le span interne ?

C'est la zone entre le premier et le dernier îlot d'un contig. C'est la région
où des gaps internes **peuvent** exister.

```
Contig :  --------XXX-----XXX-----XXX--------
                  |<--- span interne --->|
```

Le span interne est le dénominateur de `gap_fraction`. Il ne faut pas le confondre
avec "la zone couverte" (qui ne contient que les bases à profondeur > 0) ni avec
"la longueur du contig" (qui inclut les extrémités).

---

## 3. `islands_per_mbp`

### Que mesure cette métrique ?

La **densité de fragmentation** du signal : combien de segments couverts distincts
par mégabase (Mbp) de génome analysé.

### Comment la lire ?

| Valeur | Interprétation |
|--------|---------------|
| ~1-3 | Signal très continu. Peu de segments distincts. Compatible avec un vrai MAG. |
| ~5-15 | Fragmentation modérée. À croiser avec `max_gap` et `gap_fraction`. |
| >20 | Signal très fragmenté — beaucoup de petits segments isolés. Suspect. |

> **Important** : ces ordres de grandeur dépendent fortement de la profondeur, de la
> fragmentation des MAGs, de la technologie de séquençage et de `--min-island-length`.
> Ils devront être calibrés empiriquement sur vos données.

### Qu'est-ce qui sert de dénominateur ?

La **longueur totale des contigs qui portent du signal** (au moins un îlot détecté).

Les contigs entièrement non couverts sont exclus du dénominateur. Ainsi, `islands_per_mbp`
mesure uniquement la fragmentation du signal là où il existe, sans être diluée par
l'absence globale de couverture (déjà captée par `covered_fraction`).

### Que vaut `islands_per_mbp` si aucun contig n'a de signal ?

**0.0.**

### Pourquoi c'est utile pour la détection de MAGs ?

Un vrai MAG couvert à 10× aura typiquement peu d'îlots par contig (la couverture
est quasi continue). Un faux positif par cross-mapping aura des dizaines de petits îlots
correspondant à des régions conservées. `islands_per_mbp` sépare ces deux situations.

Cette métrique dépend aussi de la fragmentation des contigs : un MAG très fragmenté
(beaucoup de petits contigs) peut masquer certaines discontinuités du signal, car les
gaps entre contigs ne sont pas comptés. Un contig n'est pris en compte dans le
dénominateur que s'il contient au moins un îlot valide (après application de
`--min-island-length`).

---

## 4. `max_gap`

### Que mesure cette métrique ?

La **plus grande zone sans couverture à l'intérieur d'un contig**, en nombre de bases.
C'est le pire trou interne parmi tous les contigs du génome.

### Comment la lire ?

| Valeur | Interprétation |
|--------|---------------|
| 0 | Pas de gap interne. Signal continu partout. |
| < 500 bp | Petits trous compatibles avec du bruit stochastique. |
| 1 000 - 5 000 bp | Trous de taille comparable à un gène. À investiguer. |
| > 10 000 bp | Grand trou. Peut indiquer un faux positif, un problème de mapping, une structure génomique particulière (repeats, plasmides) ou un artefact d'assemblage. |

### Pourquoi est-ce souvent la métrique la plus facile à interpréter ?

- Elle ne dépend pas de la taille du génome (pas de normalisation)
- Elle n'est pas affectée par les contigs sans signal
- Elle capte directement le **pire cas** : un trou géant dans un contig par ailleurs couvert
- Un vrai MAG bien couvert n'a que de petits gaps stochastiques
- Un faux positif a des trous de plusieurs kilobases entre les régions conservées

### Que vaut `max_gap` si tout le signal est continu ?

**0.** (Pas de gap interne.)

---

## 5. `gap_fraction`

### Que mesure cette métrique ?

La **proportion de trous internes** dans le span interne — c'est-à-dire dans la zone
entre le premier et le dernier îlot de chaque contig.

```
Contig :  --------XXX-----XXX-----XXX--------
          ^^^^^^^^                   ^^^^^^^^
          ignoré                     ignoré
          (avant 1er îlot)           (après dernier)

                  |<--- span interne --->|
                  XXX-----XXX-----XXX
                      ^^^^^   ^^^^^
                      gaps internes

gap_fraction = bases de gaps internes / longueur du span interne
```

### Comment la lire ?

| Valeur | Interprétation |
|--------|---------------|
| 0.00 | Aucun trou interne. Signal continu dans le span interne. |
| 0.01 - 0.10 | Quelques petits trous. Normal pour un vrai MAG à profondeur modérée. |
| 0.10 - 0.30 | Fragmentation notable. Signal suspect. |
| > 0.50 | Plus de trous que de signal dans le span interne. Très probable faux positif. |

### Pourquoi pas simplement `1 - covered_fraction` ?

`covered_fraction` mesure la fraction du **génome entier** qui est couverte.
`gap_fraction` mesure la fragmentation **à l'intérieur du span interne**.

Ce n'est pas la même chose :
- Un MAG couvert à 30% mais de manière continue aura `covered_fraction = 0.30`
  et `gap_fraction ≈ 0.00`
- Un MAG couvert à 30% mais en 50 petits segments éparpillés aura
  `covered_fraction = 0.30` et `gap_fraction ≈ 0.60`

Les deux ont la même couverture globale, mais des géographies très différentes.

### Que vaut `gap_fraction` si le génome n'a qu'un seul îlot ?

**0.0.** Un seul îlot = pas de gap interne possible.

---

## 6. `--min-island-length`

### Que fait cette option ?

Elle définit la **longueur minimale** pour qu'un segment couvert soit compté comme un îlot.
Les segments plus courts sont reclassés comme "non couverts" et absorbés dans les gaps
adjacents.

`--min-island-length` ne change pas le mapping ni les alignements — il change uniquement
la façon dont les segments couverts sont **interprétés** pour les métriques spatiales.

### Pourquoi en a-t-on besoin ?

À faible profondeur ou avec des reads longs (ONT, PacBio), on peut observer des
micro-segments couverts de 1-2 bases par bruit de mapping. Ces artefacts :
- gonflent `islands_per_mbp` artificiellement (beaucoup de faux îlots)
- fragmentent les gaps en petits morceaux (sous-estimation de `max_gap`)

Le seuil filtre ce bruit.

### Quelle valeur choisir ?

| Technologie | Valeur recommandée | Pourquoi |
|-------------|-------------------|---------|
| Illumina short reads | 1 - 10 | Peu de bruit ponctuel |
| ONT / PacBio | 50 - 100 | Bruit de mapping local plus fréquent |
| Très faible profondeur (< 5×) | 50 - 100 | Micro-gaps stochastiques fréquents |

**Valeur par défaut : 1** (pas de filtrage).

> **Important** : ces valeurs sont indicatives. Le bon seuil dépend de votre jeu de données,
> de la technologie de séquençage et de la profondeur. Une exploration sur un sous-ensemble
> est recommandée avant de fixer cette valeur.

### Concrètement, que se passe-t-il avec le filtrage ?

```
Avant filtrage (--min-island-length 50) :

  |---gap---|XXX|---gap---|XXXXXXXXXXXXXXXXX|---gap---|XX|---gap---|
             ^^^                                       ^^
             3 bases                                   2 bases
             trop court                                trop court

Après filtrage :

  |---------gap élargi-----|XXXXXXXXXXXXXXXXX|--------gap élargi--------|
                            ^^^^^^^^^^^^^^^^^
                            seul îlot retenu
```

Les micro-segments sont absorbés dans les gaps. Les vrais îlots ne sont pas affectés.

### Est-ce que `--min-island-length` affecte les métriques existantes ?

**Non.** `mean`, `covered_fraction`, `variance` et les autres métriques CoverM classiques
ne sont pas modifiées. Seules les 3 nouvelles métriques spatiales sont concernées.

### Est-ce que ça affecte le profil BedGraph ?

**Non.** Le BedGraph contient toujours le profil brut, sans filtrage.

---

## 7. `--coverage-profile`

### Que fait cette option ?

Elle génère un **profil de couverture** pour chaque sample, sous forme de segments BedGraph,
en plus du TSV de métriques habituel.

### Quel est le format ?

**BedGraph compressé** (bgzf + index tabix) — un format standard compatible avec
`mosdepth`, `bedtools`, `IGV`, `JBrowse`, `tabix`, etc.

```
sample1.bedgraph.gz       ← profil compressé (bgzf)
sample1.bedgraph.gz.tbi   ← index pour accès par région
```

Contenu (4 colonnes : contig, début, fin, profondeur) :

```
contig_1    0       150     0
contig_1    150     480     12
contig_1    480     520     45
contig_1    520     800     12
contig_1    800     1000    0
```

Les zones de même profondeur sont fusionnées en un seul segment (run-length encoding).
La compression bgzf et l'indexation tabix sont créées automatiquement. L'index permet
l'accès par région : `tabix sample1.bedgraph.gz contig_1:1000-2000`.

### Pourquoi vouloir un profil en plus des métriques ?

Les métriques donnent un **résumé** par génome. Le profil donne le **détail complet**.

C'est utile pour :
- **Visualiser** le signal dans un genome browser (IGV, JBrowse)
- **Auditer** un résultat suspect (comprendre pourquoi un MAG a un `max_gap` élevé)
- **Explorer** de nouvelles métriques à partir des données brutes
- **Convertir** vers d'autres formats (BigWig, D4) pour des analyses spécifiques

### Comment l'utiliser ?

```bash
coverm genome -b sample.bam --genome-definition genomes.tsv \
  -m mean covered_fraction islands_per_mbp max_gap gap_fraction \
  --coverage-profile profiles/
```

L'option prend un **dossier** en argument. Deux fichiers sont créés par sample
(profil compressé + index) :

```
profiles/
  sample1.bedgraph.gz
  sample1.bedgraph.gz.tbi
  sample2.bedgraph.gz
  sample2.bedgraph.gz.tbi
```

### Le profil contient-il les mêmes données que les métriques ?

Non. Le profil est le **signal brut** :
- Pas de filtrage par `--min-island-length`
- Pas d'exclusion des extrémités de contigs (`--contig-end-exclusion`)
- Le profil complet de profondeur le long des contigs (représenté sous forme compressée
  en segments BedGraph)

Les métriques TSV sont des **statistiques dérivées** de ce signal, avec filtrage.
Les deux se complètent.

### Peut-on utiliser `--coverage-profile` sans les métriques spatiales ?

Oui. L'option est indépendante du choix de métriques. On peut écrire :

```bash
coverm genome -b sample.bam -m mean --coverage-profile profiles/
```

Le profil sera généré même si aucune métrique spatiale n'est demandée.

---

## 8. Utilisation combinée

### Comment tout utiliser ensemble ?

```bash
coverm genome \
  -b sample1.bam sample2.bam \
  --genome-definition genomes.tsv \
  -m mean covered_fraction islands_per_mbp max_gap gap_fraction \
  --min-island-length 50 \
  --coverage-profile profiles/ \
  --min-covered-fraction 0 \
  -t 8
```

Cette commande produit **en un seul passage sur les BAM** :

**1. Le TSV de métriques (stdout) :**

```
Genome   Mean   Covered Fraction   Islands per Mbp   Max Gap   Gap Fraction
MAG_001  12.3   0.92               0.7               2300      0.03
MAG_002  0.8    0.15               12.8              28000     0.72
```

**2. Les profils BedGraph (fichiers) :**

```
profiles/sample1.bedgraph
profiles/sample2.bedgraph
```

### Comment lire le TSV de résultats ?

Dans l'exemple ci-dessus :

**MAG_001** — vrai positif probable :
- `mean = 12.3` → bonne profondeur
- `covered_fraction = 0.92` → 92% du génome couvert
- `islands_per_mbp = 0.7` → très peu d'îlots (signal quasi continu)
- `max_gap = 2300` → le plus grand trou fait 2.3 kb (petit)
- `gap_fraction = 0.03` → 3% de trous internes (très faible)

**MAG_002** — faux positif probable :
- `mean = 0.8` → faible profondeur
- `covered_fraction = 0.15` → seulement 15% du génome couvert
- `islands_per_mbp = 12.8` → beaucoup d'îlots (signal fragmenté)
- `max_gap = 28000` → un trou de 28 kb (très grand)
- `gap_fraction = 0.72` → 72% de trous internes (plus de trous que de signal)

### Quelles combinaisons de métriques sont possibles ?

Les nouvelles métriques se combinent librement avec les existantes :

```bash
# Combinaison minimale pour la détection
-m mean covered_fraction max_gap

# Combinaison complète spatiale
-m mean covered_fraction islands_per_mbp max_gap gap_fraction

# Avec d'autres métriques CoverM
-m mean trimmed_mean covered_fraction variance islands_per_mbp max_gap gap_fraction rpkm
```

Les métriques spatiales ne peuvent **pas** être utilisées en mode `coverm contig`
(elles n'ont de sens qu'au niveau génome).

### `--min-covered-fraction` interagit-il avec les métriques spatiales ?

Oui. Si la couverture globale d'un génome est en dessous du seuil fixé par
`--min-covered-fraction`, **toutes** les métriques (existantes et spatiales) retournent
0.0 pour ce génome.

Dans ce cas, **0.0 signifie que le génome a été filtré par le seuil global**, pas
nécessairement qu'aucun signal brut n'existait. Le profil BedGraph, lui, contiendra
toujours les données brutes si `--coverage-profile` est activé — on peut s'y référer
pour comprendre ce qui a été filtré.

---

## 9. Grille de lecture combinée des métriques

Les métriques spatiales prennent tout leur sens **combinées entre elles et avec les
métriques existantes**. Voici une grille de lecture pour les cas les plus fréquents :

| Pattern observé | Lecture pratique |
|-----------------|-----------------|
| `covered_fraction` élevé + `gap_fraction` faible | Signal continu et étendu. Compatible avec un vrai MAG. |
| `covered_fraction` faible + `islands_per_mbp` élevé | Signal fragmenté et dispersé. Probable cross-mapping sur des régions conservées. |
| `max_gap` élevé + `gap_fraction` élevé | Forte discontinuité interne. Le signal est concentré en quelques zones isolées. |
| `covered_fraction` modéré + `gap_fraction` faible | Signal localisé mais continu. Peut indiquer une présence partielle ou un MAG incomplet. À investiguer. |
| `islands_per_mbp` faible + `max_gap` = 0 | Signal très continu (un ou deux gros blocs). Bon indicateur de présence réelle. |
| `mean` élevé + `islands_per_mbp` élevé | Profondeur élevée mais fragmentée. Possible mapping sur régions multicopies. Croiser avec taxonomie. |

> Ces patterns sont des guides d'interprétation, pas des règles absolues.
> L'interprétation finale doit toujours tenir compte du contexte biologique
> et des caractéristiques de l'échantillon.

---

## 10. Ce que ces métriques ne disent pas

Les métriques spatiales aident à **prioriser et orienter l'interprétation**, mais elles
ont des limites :

- Elles **n'identifient pas la cause** d'un signal fragmenté. Un `islands_per_mbp` élevé
  peut refléter du cross-mapping, un MAG chimère, un problème d'assemblage ou simplement
  une profondeur insuffisante.

- Elles **ne remplacent pas l'inspection taxonomique**. Un signal continu sur un MAG
  contaminé reste un signal continu. La cohérence taxonomique des reads doit être vérifiée
  par d'autres moyens.

- Elles **ne prouvent pas l'absence**. Un génome à 0.0 sur toutes les métriques peut
  être réellement absent, ou simplement filtré par `--min-covered-fraction`.

- Elles **dépendent de la qualité de l'assemblage**. Un MAG très fragmenté (beaucoup de
  petits contigs) aura mécaniquement moins de gaps internes détectables, même si le signal
  est réellement discontinu. La longueur des contigs influence ce que les métriques peuvent
  capturer.

- Elles **dépendent des choix de mapping**. Les paramètres d'alignement (identité minimale,
  filtrage des reads secondaires/supplémentaires, aligneur utilisé) influencent la structure
  du signal. Un changement de paramètres peut modifier les métriques spatiales. En
  particulier, le filtre d'identité (`--min-read-percent-identity`) a un impact majeur :
  à 0% d'identité, des reads d'organismes distants mappent sur les régions conservées et
  créent des îlots artificiels ; à 95%, ces reads disparaissent et les métriques changent
  radicalement. Il est recommandé de comparer les métriques spatiales à plusieurs niveaux
  d'identité pour évaluer la robustesse du signal.

- Elles **ne distinguent pas la nature biologique des îlots**. Un îlot sur un gène rRNA
  (probable cross-mapping) et un îlot sur une région intergénique (probable vrai signal)
  contribuent de la même manière aux métriques. Pour aller plus loin, le profil BedGraph
  (`--coverage-profile`) peut être croisé avec les annotations GFF du génome afin
  d'identifier quelles régions portent le signal.

- Elles **ne sont pas calibrées universellement**. Les seuils d'interprétation dépendent
  de la technologie de séquençage, de la profondeur, de la diversité de l'échantillon et
  de `--min-island-length`. Une calibration sur vos données est nécessaire.

---

## 11. Cas limites

### Que se passe-t-il si un génome n'a aucune couverture ?

Toutes les métriques spatiales valent **0.0**.

### Et si un génome n'a qu'un seul segment couvert sur tout son génome ?

- `islands_per_mbp` = 1 / (longueur en Mbp des contigs avec signal)
- `max_gap` = 0 (pas de gap interne)
- `gap_fraction` = 0.0 (pas de trou entre deux îlots)

### Et si `--min-island-length` est supérieur à tous les segments couverts ?

Tous les segments sont reclassés. Il n'y a plus d'îlot valide.
Toutes les métriques spatiales valent **0.0**.

### Les métriques spatiales fonctionnent-elles avec des reads courts ET des reads longs ?

Oui. Les métriques ne dépendent pas de la technologie de séquençage.
Cependant, le bruit de mapping varie selon la technologie, d'où l'intérêt de
`--min-island-length` pour filtrer les micro-segments, en particulier avec les reads longs.

---

## 12. Glossaire

| Terme | Définition |
|-------|-----------|
| **Îlot** | Segment contigu de bases couvertes par ≥1 read, de longueur ≥ `--min-island-length` |
| **Gap interne** | Zone sans couverture située entre deux îlots, au sein d'un même contig |
| **Span interne** | Zone entre le premier et le dernier îlot d'un contig — c'est la région où des gaps internes peuvent exister. Dénominateur de `gap_fraction`. |
| **Préfixe / Suffixe** | Bases non couvertes aux extrémités d'un contig, avant le premier ou après le dernier îlot. Ignorées dans les métriques spatiales. |
| **Cross-mapping** | Mapping biologiquement plausible de reads sur des régions conservées d'un génome auquel ils n'appartiennent pas réellement. Source principale de faux positifs. |
| **`islands_per_mbp`** | Nombre d'îlots par mégabase de génome couvert. Métrique de fragmentation du signal. |
| **`max_gap`** | Plus grand gap interne (en bases) parmi tous les contigs du génome. |
| **`gap_fraction`** | Fraction du span interne constituée de gaps internes. |
| **`--min-island-length`** | Seuil de longueur minimale pour qu'un segment couvert soit un îlot. Filtre le bruit de mapping. Défaut : 1. Ne modifie pas les alignements. |
| **`--coverage-profile`** | Option pour générer des fichiers BedGraph de profil de couverture. |
| **BedGraph** | Format standard pour les profils de couverture (4 colonnes : contig, début, fin, profondeur). |
| **RLE** | Run-Length Encoding — les zones de même profondeur sont fusionnées en un segment. |
| **MAG** | Metagenome-Assembled Genome — génome reconstruit à partir de données métagénomiques. |
