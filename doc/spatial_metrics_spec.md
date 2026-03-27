# CoverM — Métriques spatiales de couverture

Guide utilisateur pour les métriques spatiales ajoutées à `coverm genome`.
Aucune connaissance du code source n'est nécessaire.

---

## Le problème

`mean` et `covered_fraction` résument la **quantité** de signal, mais pas sa **répartition**.
Dans les cas ambigus, deux MAGs peuvent avoir les mêmes valeurs avec des profils très différents :

| MAG | mean | covered_fraction | Réalité |
|-----|------|------------------|---------|
| MAG_A | 5.0 | 0.30 | Signal continu — vrai positif |
| MAG_B | 5.0 | 0.30 | Pics isolés sur gènes conservés — faux positif (cross-mapping) |

Les métriques spatiales séparent ces deux cas en décrivant la **géographie** du signal.

---

## Les 3 nouvelles métriques

### `islands_per_mbp` — Fragmentation du signal

Nombre de segments couverts distincts (îlots) par mégabase de génome avec signal.

- Peu d'îlots → signal continu → probable vrai MAG
- Beaucoup d'îlots → signal fragmenté → suspect

### `max_gap` — Pire discontinuité

Plus grande zone sans couverture entre deux îlots, au sein d'un contig (en bases).

- Petit max_gap → pas de gros trou → bon signe
- Grand max_gap → trou de plusieurs kb → suspect

C'est souvent la métrique la plus facile à interpréter.

### `gap_fraction` — Proportion de trous internes

Fraction de la zone entre le premier et le dernier îlot (span interne)
qui est constituée de bases sans couverture.

- 0.0 → signal continu dans le span
- \> 0.5 → plus de trous que de signal → probable faux positif

---

## Concepts clés

**Îlot** : segment contigu de bases couvertes (depth > 0), de longueur ≥ `--min-island-length`.

**Gap interne** : zone sans couverture **entre** deux îlots, au sein d'un même contig.
Les extrémités non couvertes (préfixe/suffixe) ne sont pas des gaps.

```
Contig :  --------XXX-----XXX-----XXX--------
          ^^^^^^^^                   ^^^^^^^^
          ignoré                     ignoré

                  |<--- span interne --->|
                       ^^^   ^^^
                       gaps internes
```

**Analyse intra-contig uniquement** : les jonctions entre contigs ne sont jamais
interprétées comme des gaps (ce sont des artefacts d'assemblage).

**Contigs sans couverture** : ignorés (déjà captés par `covered_fraction`).

---

## Les 2 nouvelles options

### `--min-island-length <INT>` (défaut : 1)

Longueur minimale pour qu'un segment couvert soit un îlot.
Les segments plus courts sont reclassés comme non couverts.

```
Avant filtrage (--min-island-length 50) :
  |--gap--|XXX|--gap--|XXXXXXXXXXXXXXXXX|--gap--|XX|--gap--|
           ^^^  (3bp, trop court)        ^^  (2bp, trop court)

Après filtrage :
  |------gap élargi------|XXXXXXXXXXXXXXXXX|------gap élargi------|
```

Ne change **ni** le mapping, **ni** les métriques existantes, **ni** le profil BigWig.
Affecte uniquement les 3 métriques spatiales.

### `--coverage-profile <DIR>`

Génère un fichier **BigWig** par sample dans le dossier spécifié.

```
profiles/
  sample1.bw
  sample2.bw
```

Le BigWig contient le profil de profondeur brut (pas de filtrage, pas d'exclusion
d'extrémités). Directement visualisable dans IGV, JBrowse, UCSC genome browser.
Indépendant du choix de métriques.

---

## Utilisation

```bash
coverm genome \
  -b sample1.bam sample2.bam \
  --genome-definition genomes.tsv \
  -m mean covered_fraction islands_per_mbp max_gap gap_fraction \
  --min-island-length 50 \
  --coverage-profile profiles/ \
  -t 8
```

Produit en un seul passage :
1. Le TSV de métriques (stdout)
2. Les profils BigWig (fichiers)

Les métriques spatiales se combinent librement avec toutes les métriques existantes.
Elles ne sont disponibles qu'en mode `coverm genome` (pas `coverm contig`).

---

## Grille de lecture

| Pattern | Interprétation |
|---------|---------------|
| `covered_fraction` élevé + `gap_fraction` faible | Signal continu et étendu — probable vrai MAG |
| `covered_fraction` faible + `islands_per_mbp` élevé | Signal fragmenté et dispersé — probable cross-mapping |
| `max_gap` élevé + `gap_fraction` élevé | Forte discontinuité interne — signal concentré en zones isolées |
| `covered_fraction` modéré + `gap_fraction` faible | Signal localisé mais continu — présence partielle possible |
| `islands_per_mbp` faible + `max_gap` = 0 | Signal très continu — bon indicateur de présence |

Ces patterns sont des guides, pas des règles absolues.

---

## Limites

- **Ne distinguent pas la nature des îlots.** Un îlot sur un gène rRNA (cross-mapping)
  et un îlot sur une région intergénique (vrai signal) comptent pareil. Le profil BigWig
  peut être croisé avec les annotations GFF pour aller plus loin.

- **Dépendent des choix de mapping.** Le filtre d'identité (`--min-read-percent-identity`)
  a un impact majeur : à 0% les reads distants créent des îlots artificiels, à 95% ils
  disparaissent. Comparer les métriques à plusieurs niveaux d'identité est recommandé.

- **Dépendent de la qualité de l'assemblage.** Un MAG très fragmenté (500+ petits contigs)
  a mécaniquement moins de gaps internes détectables. Les métriques sont plus fiables
  sur des MAGs avec des contigs longs.

- **`gap_fraction` ≈ `1 - covered_fraction`** quand le génome a un seul contig long
  et des reads répartis sur toute sa longueur. C'est un cas dégénéré attendu — la
  métrique prend son sens sur des génomes multi-contigs.

- **Pas de calibration universelle.** Les seuils dépendent de la technologie, de la
  profondeur et de `--min-island-length`. Une calibration sur vos données est nécessaire.

---

## Glossaire

| Terme | Définition |
|-------|-----------|
| **Îlot** | Segment contigu couvert, de longueur ≥ `--min-island-length` |
| **Gap interne** | Zone sans couverture entre deux îlots, dans un même contig |
| **Span interne** | Zone entre premier et dernier îlot d'un contig (dénominateur de `gap_fraction`) |
| **Cross-mapping** | Mapping de reads sur des régions conservées d'un génome non présent |
| **BigWig** | Format binaire standard pour les profils de couverture (compatible IGV, JBrowse, deepTools) |
