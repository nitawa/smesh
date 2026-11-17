# `getPeriodicDiagonals` — Split Volumes avec périodicité

Fichier : `src/SMESH/SMESH_MeshEditor.cxx` (espace de noms anonyme, juste avant
`SMESH_MeshEditor::SplitVolumes`).

Fonctions liées :

| Fonction | Rôle |
|---|---|
| `SMESH_MeshEditor::SplitVolumes(...)` | Point d'entrée. Vérifie la méthode. Si `ensurePeriodicity` est vrai et que la méthode est `HEXA_TO_5` ou `HEXA_TO_6`, appelle `getPeriodicDiagonals` une seule fois, avant toute modification du maillage. |
| `getPeriodicDiagonals(...)` | Valide les entrées et calcule une diagonale imposée pour chaque quadrangle inlet/outlet. Ne modifie ni le maillage ni l'éditeur. |
| `isQuadDiagonal(...)` | Indique si deux nœuds sont les coins opposés d'un quadrangle donné par ses 4 nœuds. |
| `getTetraSplitMethod(...)` | Pour chaque volume, choisit le découpage en tétraèdres en respectant les diagonales imposées. |
| `TPeriodicDiagonals` | `std::set< SMESH_TLink >` : ensemble des diagonales imposées. |

Couches d'interface :

| Couche | Point d'entrée |
|---|---|
| IDL (`idl/SMESH_MeshEditor.idl`) | `SplitVolumesIntoTetraWithPeriodicity( elems, methodFlags, inlet, outlet, avoidOverConstrainedVolumes, tolerance )` |
| CORBA (`src/SMESH_I/SMESH_MeshEditor_i.cxx`) | Convertit `elems` en volumes (facette `-1`) et `inlet`/`outlet` en `TIDSortedElemSet` de faces, puis appelle `SplitVolumes( ..., ensurePeriodicity = true, inlet, outlet, tolerance )`. |
| Python (`src/SMESH_SWIG/smeshBuilder.py`) | `Mesh.SplitVolumesIntoTetra( elems, method, avoidOverConstrainedVolumes, ensurePeriodicity=False, inlet=None, outlet=None, tolerance=1e-7 )` |

---

## 1. Principe

La fonctionnalité découpe des volumes linéaires (hexaèdres, prismes, pyramides)
en tétraèdres. Elle garantit en plus que deux surfaces du maillage, **Inlet** et
**Outlet**, identiques à une translation près avant le découpage, le restent
après.

Seul le choix des diagonales peut casser cette périodicité. Un quadrangle se
coupe en deux triangles par la diagonale 0-2 ou par la diagonale 1-3 (les
triangles ne sont pas modifiés). Si un quadrangle de l'inlet et son image sur
l'outlet sont coupés par des diagonales différentes, les deux surfaces ne se
correspondent plus.

- **`HEXA_TO_24`** : chaque quadrangle est coupé en 4 triangles autour du centre
  de la face, sans diagonale. Le découpage normal conserve donc la périodicité.
  `getPeriodicDiagonals` n'est pas appelée, et le noyau C++ ne vérifie pas
  l'inlet et l'outlet.
- **`HEXA_TO_5` et `HEXA_TO_6`** : `getPeriodicDiagonals` valide les données,
  apparie l'inlet et l'outlet, puis choisit une diagonale pour chaque
  quadrangle de l'inlet et impose son image au quadrangle correspondant de
  l'outlet. Toutes les erreurs sont levées (`SALOME_Exception`) avant de
  modifier le maillage.

Les diagonales sont ensuite passées à `getTetraSplitMethod`, qui choisit pour
chaque volume un découpage qui les respecte (voir section 4).

---

## 2. Signature

```cpp
typedef std::set< SMESH_TLink > TPeriodicDiagonals;

TPeriodicDiagonals getPeriodicDiagonals( const SMESHDS_Mesh*                   meshDS,
                                         const SMESH_MeshEditor::TFacetOfElem& elems,
                                         const TIDSortedElemSet&               inlet,
                                         const TIDSortedElemSet&               outlet,
                                         const double                          tolerance = Precision::Confusion() );
```

- `meshDS` : le maillage. Il sert seulement à vérifier que les éléments en font partie.
- `elems` : volumes à découper (map volume → index de facette ; `-1` = découpage en tétraèdres).
- `inlet`, `outlet` : faces 2D (triangles/quadrangles linéaires) des deux surfaces périodiques.
- `tolerance` : distance absolue (unités du maillage) pour apparier un nœud de
  l'inlet et son translaté sur l'outlet. Elle doit être strictement positive.
- Valeur de retour : deux diagonales par paire de quadrangles appariés, une
  pour l'inlet et une pour l'outlet. Les triangles n'en ont pas.

La tolérance est optionnelle à tous les niveaux. Sa valeur par défaut est `Precision::Confusion()` (`1e-7`, la précision d'OCC pour les superpositions). `SplitVolumes` la transmet par son paramètre `theTolerance`, CORBA par `tolerance`, et Python par l'argument `tolerance=1e-7`.

`SMESH_TLink` est une paire de nœuds triée : le sens de la diagonale ne compte
pas. Une diagonale suffit à identifier son quadrangle, il n'y a donc pas de clé
de face.

---

## 3. Algorithme

Les étapes suivent les commentaires `STEP n` du code.

### Étape 1 — Arguments

- `elems`, `inlet` et `outlet` ne sont pas vides.
- `inlet.size() == outlet.size()`.

### Étape 2 — Volumes sélectionnés

Pour chaque volume de `elems` :

1. Il appartient au maillage et doit être découpé en tétraèdres
   (`elemFacet.second < 0`), pas en prismes.
2. Son type est `Tetra`, `Hexa`, `Penta` ou `Pyramid` (linéaire). Les éléments
   quadratiques et les polyèdres sont refusés. Les tétraèdres sont acceptés
   mais `SplitVolumes` les ignore ensuite.

### Étape 3 — Faces de l'inlet et de l'outlet

Pour chaque face de l'inlet (`iSide = 0`) puis de l'outlet (`iSide = 1`) :

1. La face appartient au maillage et est un `Triangle` ou un `Quadrangle`
   linéaire.
2. Une face de l'outlet n'est pas aussi dans l'inlet.
3. La face n'a pas de nœud répété.
4. Ses nœuds sont ajoutés à `nodes[iSide]`, et chacune de ses arrêtes
   (`SMESH_TLink`, paire de nœuds triée) est comptée dans `nbFacesOfLink[iSide]`
   (utilisé à l'étape 4).
5. **Volume propriétaire** : on cherche les volumes qui ont cette face
   (`GetElementsByNodes` + `GetFaceIndex`). Chacun doit :
   - appartenir à `elems` ;
   - n'avoir aucun voisin sur cette face (`IsFreeFace`). La face est donc au
     bord du domaine, et il ne peut pas y avoir plusieurs propriétaires ;
   - pour un quadrangle : avoir les mêmes diagonales que la face 2D
     (`isQuadDiagonal`). L'orientation et le nœud de départ peuvent différer.
     Sinon, la diagonale choisie sur la face 2D ne serait pas une diagonale de
     la face du volume, et `getTetraSplitMethod` l'ignorerait ;
   - pour un quadrangle : n'avoir aucun tétraèdre accolé à l'un des triangles
     possibles (`TTriangleFacet::hasAdjacentVol`). La diagonale serait déjà
     fixée et on ne pourrait plus l'imposer.

   Il faut ensuite `nbOwners == 1`. Comme `IsFreeFace` exclut déjà plusieurs
   propriétaires, ce test détecte en pratique une face sans volume.

Après la boucle : `nodes[0].size() == nodes[1].size()`.

Deux faces distinctes avec les mêmes nœuds sur un même côté sont détectées à
l'étape 6.

### Étape 4 — Translation

La translation est le vecteur entre les centres des nœuds du **bord** des deux
surfaces.

**Nœuds du bord.** Une arrête du bord est une arrête qui n'appartient qu'à une seule
face de la surface (`nbFacesOfLink[iSide][link] == 1`). Les nœuds de ces liens
forment `borderNodes[iSide]`.

- Il faut `borderNodes[0]` non vide et
  `borderNodes[0].size() == borderNodes[1].size()`. Sinon (surface fermée,
  bords incompatibles), une exception est levée.

Note : `SMESH_MeshAlgos::FindFreeBorders` et `SMESH::Controls::FreeEdges` ne sont pas
utilisés : ils cherchent les bords libres de tout le maillage 2D, pas ceux d'un
sous-ensemble de faces. D'autres faces accolées à l'inlet fausseraient le
résultat.

**Intérêt du bord.** Le centre ne dépend pas des nœuds intérieurs des
surfaces.

Notations : $i \in \{0, 1\}$ (0 = inlet, 1 = outlet), $\mathcal{N}_i$ les
nœuds du bord de la surface $i$, $n_i = |\mathcal{N}_i|$, $\mathbf{x}_k$ les
coordonnées du nœud $k$, $\mathbf{o}_i$ le premier nœud de $\mathcal{N}_i$
(`origin`).

$$
\mathbf{c}_i = \mathbf{o}_i + \frac{1}{n_i} \sum_{k \in \mathcal{N}_i} \left( \mathbf{x}_k - \mathbf{o}_i \right)
\;=\; \frac{1}{n_i} \sum_{k \in \mathcal{N}_i} \mathbf{x}_k
\qquad
\mathbf{t} = \mathbf{c}_1 - \mathbf{c}_0
$$

$\mathbf{c}_i$ est la moyenne simple des positions des nœuds du bord (pas le
centre géométrique pondéré par les longueurs). Le biais dû à un maillage non uniforme est
le même sur les deux surfaces et disparaît dans la différence.

**Intérêt du calcul relatif à $\mathbf{o}_i$.** Les deux écritures sont égales
en arithmétique exacte, pas en virgule flottante ($\varepsilon \approx 2.2 \times 10^{-16}$).
Soit $D_i$ le diamètre de la surface $i$ (plus grande distance entre deux de
ses nœuds). Comme $\mathbf{o}_i$ appartient à la surface,
$\|\mathbf{x}_k - \mathbf{o}_i\| \le D_i$ quel que soit le nœud choisi.

| Calcul | Erreur d'arrondi sur la somme |
|---|---|
| Somme directe des $\mathbf{x}_k$ | $\sim \varepsilon\, n_i \Vert\mathbf{x}\Vert$ |
| Somme des $\mathbf{x}_k - \mathbf{o}_i$ | $\sim \varepsilon\, n_i D_i$ |

Le calcul relatif retire le facteur $n_i$ devant les grandes coordonnées.
L'erreur totale sur $\mathbf{t}$ devient $\sim \varepsilon\,(\|\mathbf{c}\| + n_i D_i)$
au lieu de $\sim \varepsilon\, n_i \|\mathbf{c}\|$. Le terme $\varepsilon \|\mathbf{c}\|$
reste dans tous les cas : il est déjà dans les coordonnées stockées.
Le gain n'apparaît que si la surface est petite devant sa distance à l'origine
du repère ($D_i \ll \|\mathbf{x}\|$, par exemple un tube de 10 cm placé à 1 km
de l'origine). Sinon les deux calculs se valent.

### Étape 5 — Appariement des nœuds

- Un `SMESH_OctreeNode` est construit sur les nœuds de l'outlet.
- Pour chaque nœud `P` de l'inlet, `NodesAround` cherche les nœuds de l'outlet
  à moins de `tolerance` de `P + translation`.
- Il faut exactement un candidat, pas encore apparié (`matchedNodes`).
- Résultat : la bijection `nodeMap : nœud inlet → nœud outlet`.

### Étape 6 — Appariement des faces

Pour chaque face de l'inlet :

1. On prend l'image de ses nœuds par `nodeMap`, dans l'ordre de la face.
2. Parmi les faces qui contiennent ces nœuds (`SMDS_Mesh::GetElementsByNodes`),
   on garde celles qui ont le même nombre de nœuds et appartiennent à
   l'outlet. Il en faut exactement une, pas encore appariée (`matchedFaces`).
3. Pour chaque quadrangle, on vérifie que son image est également un quandrangle, et que la numérotation correspond (pour placer correctement la diagonale des deux côtés).
4. La paire `(faceInlet, faceOutlet)` est rangée dans `facePairs`.


### Étape 7 — Choix des diagonales

Pour chaque paire dont la face inlet est un quadrangle `n[0..3]` :

1. On construit les deux triangulations possibles avec des `SMDS_FaceOfNodes`
   temporaires : `(012, 230)` pour la diagonale 0-2 et `(123, 301)` pour la
   diagonale 1-3.
2. On garde celle dont la somme des `getBadRate(..., AspectRatio)` est la plus
   petite, comme `SMESH_MeshEditor::BestSplit`. En cas d'égalité, la
   diagonale 0-2.
3. On ajoute `SMESH_TLink( n[i], n[i+2] )` (inlet) et
   `SMESH_TLink( nodeMap[n[i]], nodeMap[n[i+2]] )` (outlet).

---

## 4. Utilisation dans `getTetraSplitMethod`

Pour chaque face quadrangulaire `iF` du volume :

1. Comme avant, si un tétraèdre déjà créé est accolé à la face, sa diagonale
   est ajoutée à `triaSplits` (`t012` ou `t123`).
2. Si `SMESH_TLink( fNodes[0], fNodes[2] )` (`iCommon = 0`) ou
   `SMESH_TLink( fNodes[1], fNodes[3] )` (`iCommon = 1`) est dans
   `periodicDiagonals`, le triangle
   `(nInd[iCommon], nInd[iCommon+1], nInd[iCommon+2])` est ajouté à
   `triaSplits`. Un `ASSERT( triaSplits.empty() )` (actif en debug seulement)
   vérifie qu'aucun tétraèdre n'était accolé, ce que l'étape 3 garantit.
3. Les variantes standard (`theHexTo5` ou `theHexTo6`, `thePyraTo2`,
   `thePentaTo3`) sont essayées dans l'ordre. On garde la première qui contient
   toutes les facettes de `triaSplits` (`TSplitMethod::hasFacet`) et qui, si
   demandé, n'est pas sur-contrainte.
4. Sinon, on passe au **découpage barycentrique** : un nœud est créé au centre
   du volume, chaque face est triangulée et reliée à ce nœud. Une face qui a
   une entrée dans `triaSplits` est triangulée selon cette entrée : la
   diagonale imposée est donc respectée. Un hexaèdre donne alors 12 tétraèdres.

Ensuite, pendant le découpage de chaque volume, `SplitVolumes` remplace les
faces 2D posées sur ses faces quadrangulaires par des triangles qui suivent
les tétraèdres créés (`splitMethod.hasFacet`), et met à jour les groupes
(`ReplaceElemInGroups`). Les triangles de l'inlet et de l'outlet portent donc
les diagonales imposées.

---

## 5. Schéma récapitulatif

```
SplitVolumes(elems, method, avoidOC, ensurePeriodicity=true, inlet, outlet, tolerance)
│
├─ méthode HEXA_TO_5 / HEXA_TO_6 / HEXA_TO_24 ?  sinon erreur
├─ HEXA_TO_24 → pas de diagonale, découpage normal
│
├─ getPeriodicDiagonals                      (aucune modification du maillage)
│   ├─ 1. arguments non vides, même nombre de faces, tolérance > 0
│   ├─ 2. volumes : dans le maillage, découpage en tétraèdres, types linéaires
│   ├─ 3. faces   : type, distinctes, un seul propriétaire, face libre,
│   │               mêmes diagonales que la face du volume, pas de tétra accolé
│   ├─ 4. translation = centre(bord outlet) − centre(bord inlet)
│   ├─ 5. nœud inlet  → nœud outlet   (octree, tolerance, bijection)
│   ├─ 6. face inlet  → face outlet   (nœuds translatés + diagonales compatibles)
│   └─ 7. quadrangle  : AspectRatio sur l'inlet → diagonale (n1,n2)
│                       outlet ← (nodeMap[n1], nodeMap[n2])
│
└─ pour chaque volume :
    getTetraSplitMethod(vol, …, periodicDiagonals)
      ├─ contraintes = diagonales des voisins déjà découpés ∪ diagonales périodiques
      ├─ variante standard compatible ?  oui → elle
      └─ sinon → découpage barycentrique (respecte aussi les contraintes)
    création des tétraèdres, triangulation des faces 2D, mise à jour des groupes
```

---

## 6. Points d'attention

- **`Hex_5Tet` et parité.** Dans `theHexTo5`, les faces opposées d'un hexaèdre
  ont des diagonales croisées : variante 1, 0-2 en bas et 5-7 en haut ;
  variante 2, 1-3 en bas et 4-6 en haut. Or l'image par translation de 0-2 est
  4-6. Si l'inlet et l'outlet sont aux deux bouts d'une colonne d'hexaèdres,
  la compatibilité dépend donc de la parité du nombre de couches. En cas
  d'échec, l'hexaèdre concerné passe au découpage barycentrique (12 tétraèdres
  et un nœud de plus). `Hex_6Tet` n'a pas ce problème : ses variantes ont des
  diagonales parallèles sur les faces opposées (variante 1 : 0-2 et 4-6).
- **Translation uniquement.** Les périodicités par rotation ne sont pas gérées.
  La translation est estimée par les centres des nœuds du bord : elle est
  fausse si le bord de l'outlet n'a pas exactement la même distribution de
  nœuds que celui de l'inlet.
- **Tolérance absolue.** `tolerance` n'est pas mise à l'échelle du maillage.
  Avec la valeur par défaut (`1e-7`), un très grand maillage, ou des coordonnées
  arrondies à l'import, peuvent échouer à l'étape 5 : il faut alors augmenter
  `tolerance`. Une tolérance trop grande donne plusieurs candidats, et
  l'étape 5 échoue aussi.
- **Toutes les erreurs avant la première modification.** `SplitVolumes`
  modifie le maillage volume par volume, sans retour arrière possible. Une
  exception levée pendant la boucle de découpage laisserait le maillage à
  moitié découpé. C'est pourquoi toutes les vérifications sont faites dans
  `getPeriodicDiagonals`, avant la boucle, et pourquoi `getTetraSplitMethod`
  n'utilise qu'un `ASSERT`. Toute nouvelle vérification doit suivre la même
  règle.
- **Recherche par diagonale.** `getTetraSplitMethod` cherche seulement la paire
  de nœuds. Chaque quadrange est identifié de manière unique par sa diagonale (vrai pour un maillage conforme).
- **Accès à `nodeMap` sans test** Lors des étapes 6 et 7, l'accès :`nodeMap.find(...)->second` est assuré par la contruction à l'étape 5.

---

## 7. Messages d'erreur

### `SALOME_Exception` (`SMESH_MeshEditor.cxx`, préfixe « SplitVolumes: »)

| Message | Étape | Cause |
|---|---|---|
| `periodicity requires Hex_5Tet, Hex_6Tet or Hex_24Tet` | avant l'appel | Méthode non supportée. |
| `volumes, Inlet and Outlet must not be empty` | 1 | Une des entrées est vide. |
| `the periodicity tolerance must be positive` | 1 | `tolerance` nulle, négative ou NaN. |
| `Inlet and Outlet have different numbers of faces` | 1 | Tailles différentes. |
| `invalid volume selection for periodicity` | 2 | Volume absent du maillage, ou découpage en prismes demandé. |
| `periodicity supports only linear tetrahedra, hexahedra, prisms and pyramids` | 2 | Élément quadratique, polyèdre, etc. |
| `Inlet and Outlet must contain linear triangles or quadrangles of this mesh` | 3 | Face invalide, quadratique ou d'un autre maillage. |
| `Inlet and Outlet must be distinct surfaces` | 3 | Une face est dans les deux côtés. |
| `repeated node in a periodic face` | 3 | Face dégénérée. |
| `a periodic face is not on the boundary of the selected volumes` | 3 | Propriétaire non sélectionné, ou face partagée par deux volumes. |
| `the nodes of a periodic quadrangle are not in the same order as in the face of its volume` | 3 | Les nœuds du quadrangle 2D ne sont pas dans le même ordre cyclique que ceux de la face du volume, donc leurs diagonales diffèrent. |
| `a periodic quadrangle already has adjacent tetrahedra` | 3 | Un tétraèdre existant touche déjà l'un des 4 triangles du quadrangle, donc sa diagonale est déjà fixée. |
| `each periodic face must have one selected owner` | 3 | Aucun volume ne porte la face. |
| `Inlet and Outlet have different numbers of nodes` | 3 | Nombre de nœuds différent. |
| `Inlet and Outlet borders are missing or have different numbers of nodes` | 4 | Surface fermée (aucun lien du bord) ou nombre de nœuds du bord différent. |
| `missing or ambiguous translated node on Outlet` | 5 | Pas d'image unique d'un nœud à la tolérance. |
| `periodic face connectivities do not match` | 6 | Zéro ou plusieurs faces outlet pour une face inlet, ou face outlet déjà appariée. |
| `periodic faces have incompatible node ordering` | 6 | Bons nœuds, mais diagonales différentes. |

### Couches d'interface

Ces vérifications s'appliquent aussi avec `Hex_24Tet`.

| Couche | Message | Cause |
|---|---|---|
| CORBA (`BAD_PARAM`) | `Inlet and outlet must be FACE sources of this mesh` | Source invalide ou d'un autre maillage. |
| CORBA (`BAD_PARAM`) | `Inlet and outlet must not be empty` | Source vide. |
| Python (`ValueError`) | `ensurePeriodicity requires inlet and outlet FACE sources` | `inlet` ou `outlet` vaut `None`. |
