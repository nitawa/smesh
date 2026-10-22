# `getPeriodicDiagonals` — Split Volumes avec périodicité

Fichier : `src/SMESH/SMESH_MeshEditor.cxx` (espace de noms anonyme, juste avant
`SMESH_MeshEditor::SplitVolumes`).

Fonctions liées :

| Fonction | Rôle |
|---|---|
| `SMESH_MeshEditor::SplitVolumes(...)` | Point d'entrée. Appelle `getPeriodicDiagonals` une seule fois, avant toute modification, si `ensurePeriodicity` est vrai. |
| `getPeriodicDiagonals(...)` | Valide les entrées et calcule une diagonale imposée pour chaque quadrangle inlet/outlet. |
| `getTetraSplitMethod(...)` | Pour chaque volume, choisit le découpage en tétraèdres en respectant les diagonales imposées. |
| `TVolumeFaceKey` | Clé d'identification d'une face par l'ensemble trié de ses nœuds (indépendante de l'orientation et du point de départ). |

---

## 1. Résumé

La fonctionnalité « Split Volumes avec périodicité » découpe des volumes linéaires
(hexaèdres, prismes, pyramides, tétraèdres) en tétraèdres. Elle garantit en plus
que deux surfaces du maillage, **Inlet** et **Outlet**, restent identiques à une
translation près **après** le découpage, comme elles l'étaient **avant**.

La seule chose qui peut casser cette périodicité est le choix des diagonales.
Un quadrangle se coupe en deux triangles de deux façons : par la diagonale 0-2
ou par la diagonale 1-3. Les triangles d'origine restent des triangles. Si un
quadrangle de l'inlet est coupé par une diagonale et que son image sur l'outlet
est coupée par l'autre, les deux surfaces ne se correspondent plus.

`getPeriodicDiagonals` règle ce problème **en amont** du découpage. Elle :

1. **vérifie** que les données permettent un découpage périodique : méthode,
   types d'éléments, faces au bord, nombre de faces et de nœuds identiques, etc.
   Toutes les vérifications ont lieu avant de modifier le maillage. En cas
   d'erreur, une `SALOME_Exception` est levée et le maillage reste intact ;
2. **déduit la translation** inlet → outlet à partir des barycentres des nœuds
   des deux surfaces ;
3. **apparie** chaque nœud de l'inlet à un nœud de l'outlet (recherche par
   octree à la tolérance `1e-7`), puis chaque face de l'inlet à une face de
   l'outlet ;
4. **choisit une diagonale** pour chaque quadrangle de l'inlet (critère de
   qualité `AspectRatio` via `BestSplit`) et **reporte son image** sur le
   quadrangle correspondant de l'outlet.

Le résultat est une table `TPeriodicDiagonals` : clé de face → paire d'IDs de
nœuds de la diagonale imposée. Elle est ensuite passée à `getTetraSplitMethod`,
qui l'ajoute aux contraintes déjà présentes (diagonales imposées par les
tétraèdres voisins). Le découpage standard (variantes `Hex_5Tet`/`Hex_6Tet`) est
alors choisi pour respecter toutes ces contraintes. Si aucune variante ne
convient, on passe au découpage barycentrique, qui respecte lui aussi la
diagonale imposée.

Les faces 2D inlet/outlet sont ensuite remplacées par leurs triangles à la fin
de `SplitVolumes`, avec `ReplaceElemInGroups`. Les groupes Inlet/Outlet
contiennent donc des triangles qui se correspondent par translation.

---

## 2. Signature et sorties

```cpp
void getPeriodicDiagonals( SMESH_MeshEditor&                     editor,
                           const SMESH_MeshEditor::TFacetOfElem& elems,
                           const TIDSortedElemSet&               inlet,
                           const TIDSortedElemSet&               outlet,
                           const int                             methodFlags,
                           const double                          tolerance,
                           TPeriodicDiagonals&                   diagonals );
```

- `elems` : volumes à découper (map volume → index de facette ; `-1` = découpage en tétraèdres).
- `inlet`, `outlet` : faces 2D (triangles/quadrangles linéaires) des deux surfaces périodiques.
- `methodFlags` : `HEXA_TO_5` ou `HEXA_TO_6` uniquement.
- `tolerance` : tolérance **absolue** d'appariement des nœuds (fixée à `1e-7` dans `SplitVolumes`).
- `diagonals` (sortie) :

```cpp
typedef pair< smIdType, smIdType >           TFaceDiagonal;      // IDs des 2 nœuds de la diagonale
typedef map< TVolumeFaceKey, TFaceDiagonal > TPeriodicDiagonals; // face -> diagonale imposée
```

Seuls les quadrangles ont une entrée (un triangle n'a rien à choisir). Pour
chaque paire de faces appariées, il y a **deux** entrées : une pour la face
inlet, une pour la face outlet.

`TVolumeFaceKey` construit la clé à partir des IDs triés des nœuds de coin.
On obtient donc la même clé à partir d'un `SMDS_MeshElement` face, d'un
`TIDSortedNodeSet` ou d'une face de volume (`SMDS_VolumeTool`, indice `iF`).
C'est ce qui permet à `getTetraSplitMethod` de retrouver la contrainte en
partant de la face du volume.

---

## 3. Détail de l'algorithme

### Étape 0 — Préconditions globales

- `methodFlags` doit valoir `HEXA_TO_5` ou `HEXA_TO_6`. `HEXA_TO_24` découpe
  chaque face autour de son barycentre et n'a pas de diagonale à choisir. Il
  n'est pas pris en charge ici.
- `elems`, `inlet` et `outlet` ne doivent pas être vides.
- `inlet.size() == outlet.size()`.

### Étape 1 — Validation des volumes sélectionnés

Pour chaque volume de `elems` :

1. Il doit exister dans le maillage (`FindElement(id) == elem`) et être marqué
   pour un découpage en tétraèdres (`elemFacet.second < 0`). Un découpage en
   prismes est refusé.
2. Son type doit être linéaire : `Tetra`, `Hexa`, `Penta` ou `Pyramid`. Les
   éléments quadratiques et les polyèdres sont refusés.
3. `SMDS_VolumeTool::Set` doit réussir.
4. **Conformité** : pour chaque face quadrangulaire du volume, on cherche les
   autres volumes qui partagent cette face (`GetElementsByNodes` +
   `GetFaceIndex`). Si l'un d'eux **n'est pas** dans la sélection, on lève une
   erreur. Sinon, le découpage introduirait une diagonale sur une face partagée
   avec un volume non découpé, et le maillage deviendrait non conforme.

### Étape 2 — Validation et indexation des faces inlet/outlet

Pour chaque côté (`iSide = 0` inlet, `1` outlet) et pour chaque face :

1. La face doit appartenir au maillage et être un `Triangle` ou un
   `Quadrangle` linéaire.
2. Pas de nœud répété dans la face.
3. Pas de doublon dans le même côté. La face est rangée dans
   `faces[iSide]`, une `map<TVolumeFaceKey, face>`.
4. Pour l'outlet : la face ne doit pas exister aussi dans l'inlet. Les deux
   surfaces doivent être distinctes.
5. Ses nœuds sont ajoutés à `nodes[iSide]`.
6. **Propriétaire unique au bord** : on parcourt les volumes qui contiennent
   la face (`GetElementsByNodes` + `GetFaceIndex`). Chaque volume
   « propriétaire » trouvé doit :
   - appartenir à la sélection `elems` ;
   - avoir cette face libre (`IsFaceExternal`/`IsFreeFace`). La face est donc
     au bord du domaine ;
   - avoir une connectivité de face compatible avec celle de l'élément 2D. On
     part du nœud `ownerNodes[0]` et on vérifie que les nœuds de la face du
     volume suivent ceux de la face 2D dans un sens ou dans l'autre
     (`forward`/`backward`). L'orientation peut donc être inversée ;
   - pour un quadrangle : ne pas déjà avoir de tétraèdre adjacent sur l'un des
     triangles possibles (`TTriangleFacet::hasAdjacentVol`). Une diagonale
     serait déjà fixée par ailleurs, et on ne pourrait plus l'imposer.

   Il faut exactement **un** propriétaire (`nbOwners == 1`).

Après la boucle : `nodes[0].size() == nodes[1].size()`.

### Étape 3 — Calcul de la translation

On calcule le barycentre des nœuds de chaque surface **relativement au premier
nœud** (`origin`). Cela évite d'additionner de grandes coordonnées absolues et
donc de perdre en précision :

```
centre[s] = origin_s + (1/N) * Σ (P_i - origin_s)
translation = centre[outlet] - centre[inlet]
```

Ce calcul n'est valable que si l'outlet est l'image exacte de l'inlet par une
translation (mêmes nœuds, même maillage). Une rotation n'est pas prise en
compte.

### Étape 4 — Appariement des nœuds

- Un `SMESH_OctreeNode` est construit sur les nœuds de l'outlet.
- Pour chaque nœud `P` de l'inlet, on cherche les nœuds de l'outlet à moins de
  `tolerance` de `P + translation` (`NodesAround`).
- Il faut **exactement un** candidat, pas encore apparié (`matchedNodes`).
  Sinon on lève l'erreur « missing or ambiguous translated node on Outlet ».
- On obtient la bijection `nodeMap : nœud inlet → nœud outlet`.

### Étape 5 — Appariement des faces

Pour chaque face de l'inlet :

1. On calcule l'ensemble de ses nœuds translatés (via `nodeMap`) et on cherche
   la face de l'outlet de même clé dans `faces[1]`. Elle doit exister et ne pas
   être déjà appariée (`matchedFaces`).
2. On vérifie l'ordre cyclique des nœuds : l'image de la face inlet doit suivre
   la face outlet dans un sens ou dans l'autre (même test `forward`/`backward`
   qu'à l'étape 2). Cela exclut par exemple un quadrangle dont l'image aurait
   les bons nœuds dans un ordre « croisé ».
3. La paire `(faceInlet, faceOutlet)` est rangée dans `facePairs`.

Comme `inlet.size() == outlet.size()` et que l'appariement est injectif, on
obtient une bijection entre les faces.

### Étape 6 — Choix et report des diagonales

Cette étape vient **après** toutes les validations, parce que
`SMESH_MeshEditor::BestSplit` appelle `ClearLastCreated()` et modifie donc
l'état de l'éditeur.

Pour chaque paire dont la face inlet est un quadrangle :

1. `BestSplit(face, AspectRatio)` renvoie `1` (diagonale 0-2) ou `2`
   (diagonale 1-3). On en déduit `i = résultat - 1`, puis
   `n1 = face->GetNode(i)` et `n2 = face->GetNode(i+2)`.
2. Entrée pour la face inlet : `TVolumeFaceKey(faceInlet) → (n1, n2)`.
3. Entrée pour la face outlet : `TVolumeFaceKey(faceOutlet) → (nodeMap[n1], nodeMap[n2])`.

La diagonale est choisie selon la qualité des triangles **de l'inlet**.
L'outlet reçoit l'image exacte de cette diagonale. Comme les deux faces sont
identiques à une translation près, la qualité est la même des deux côtés.

---

## 4. Utilisation des diagonales dans `getTetraSplitMethod`

Pour chaque face quadrangulaire `iF` du volume en cours :

1. Comme avant, on cherche si un tétraèdre adjacent (déjà créé) impose une
   diagonale. On remplit alors `triaSplits` avec `t012` ou `t123`.
2. Si `periodicDiagonals` n'est pas vide, on cherche une contrainte avec
   `TVolumeFaceKey(vol, iF)`. Si elle existe :
   - `iCommon = 0` si le nœud local `nInd[0]` est une extrémité de la
     diagonale, `1` sinon. La diagonale est alors `nInd[iCommon]`–`nInd[iCommon+2]` ;
   - on construit le triangle `(nInd[iCommon], nInd[iCommon+1], nInd[iCommon+2])` ;
   - si `triaSplits` contient déjà un triangle qui ne porte pas cette
     diagonale, on lève l'erreur « conflicting periodic and adjacent
     diagonals » ;
   - si `triaSplits` est vide, on y ajoute ce triangle.
3. Les variantes standard (`theHexTo5`, `theHexTo6`, `thePyraTo2`,
   `thePentaTo3`) sont essayées dans l'ordre. On garde la première qui contient
   **toutes** les facettes de `triaSplits` (`TSplitMethod::hasFacet`) et qui,
   si demandé, n'est pas sur-contrainte.
4. Si aucune variante ne convient, on passe au **découpage barycentrique** : un
   nœud est créé au centre du volume, et chaque face est triangulée puis reliée
   à ce nœud. Pour une face qui a une entrée dans `triaSplits`, `iCommon` est
   tiré du premier triangle de la liste. La diagonale imposée est donc aussi
   respectée dans ce cas. Un hexaèdre donne alors 12 tétraèdres et un nœud
   intérieur.

Ensuite, `SplitVolumes` triangule les faces 2D existantes avec les facettes des
tétraèdres créés (`splitMethod.hasFacet`). Les triangles des faces inlet/outlet
portent donc les diagonales imposées, et `ReplaceElemInGroups` met à jour les
groupes.

Remarque : l'ancienne version forçait `facetCreated = true` quand
`ensurePeriodicity` était vrai. La première variante était alors acceptée sans
vérifier la compatibilité avec les voisins. Ce n'est plus le cas.

---

## 5. Schéma récapitulatif

```
SplitVolumes(elems, method, avoidOC, ensurePeriodicity=true, inlet, outlet)
│
├─ getPeriodicDiagonals                      (aucune modification du maillage
│   ├─ 0. méthode HEXA_TO_5/6, tailles        avant la fin des validations)
│   ├─ 1. volumes : types linéaires, conformité des quadrangles partagés
│   ├─ 2. faces   : type, unicité, un seul propriétaire, face libre, ordre des nœuds
│   ├─ 3. translation = barycentre(outlet) − barycentre(inlet)
│   ├─ 4. nœud inlet  → nœud outlet   (octree, tolérance 1e-7, bijection)
│   ├─ 5. face inlet  → face outlet   (clé de nœuds + ordre cyclique)
│   └─ 6. quadrangle  : BestSplit sur l'inlet → diagonale (n1,n2)
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

## 6. Points d'attention pour les développeurs

- **Hexaèdres et `Hex_5Tet`.** Dans les deux variantes de `theHexTo5`, les
  faces opposées d'un hexaèdre ont des diagonales « croisées ». Variante 1 :
  0-2 en bas, 5-7 en haut. Variante 2 : 1-3 en bas, 4-6 en haut. L'image par
  translation de 0-2 est 4-6, pas 5-7. Si l'inlet et l'outlet sont deux faces
  opposées d'une même colonne d'hexaèdres, la compatibilité dépend donc de la
  parité du nombre de couches. En cas d'échec, on retombe sur le découpage
  barycentrique (12 tétraèdres et un nœud de plus par hexaèdre concerné).
  Avec `Hex_6Tet`, toutes les variantes ont des diagonales parallèles sur les
  faces opposées (variante 1 : 0-2 en bas, 4-6 en haut). Ce problème n'existe
  donc pas.
- **Tolérance absolue** (`1e-7`, en unités du maillage). Elle n'est pas mise à
  l'échelle de la taille du maillage. Un maillage très grand, ou importé avec
  un arrondi des coordonnées, peut échouer à l'étape 4.
- **Translation uniquement.** Les périodicités par rotation ne sont pas gérées.
  La translation est estimée par les barycentres. Elle est donc fausse si les
  deux surfaces n'ont pas exactement la même distribution de nœuds.
- **Ordre des appels.** `BestSplit` vide l'historique de création de l'éditeur.
  Il doit rester après toutes les vérifications, et aucune nouvelle
  vérification qui lève une exception ne doit être ajoutée après l'étape 6.
- **Exception dans `getTetraSplitMethod`.** L'erreur « conflicting periodic and
  adjacent diagonals » serait levée **pendant** la boucle de découpage, donc
  avec un maillage partiellement modifié. En pratique, les vérifications des
  étapes 1 et 2 rendent ce cas inatteignable : une face périodique est libre,
  n'a qu'un seul propriétaire et n'a pas de tétraèdre adjacent. Ce test sert de
  garde-fou.
- **Accès à `nodeMap` sans test** (`nodeMap.find(...)->second`). C'est sûr
  parce que tous les nœuds des faces inlet sont dans `nodes[0]` et ont tous été
  appariés à l'étape 4. Toute modification de ces étapes doit conserver cette
  propriété.
- **Côtés asymétriques.** La qualité est évaluée sur l'inlet seulement. Comme
  les deux faces sont identiques à une translation près, cela ne change rien,
  mais on pourrait aussi choisir côté outlet sans changer le résultat.

---

## 7. Messages d'erreur (`SALOME_Exception`)

| Message | Étape | Cause |
|---|---|---|
| `periodicity requires Hex_5Tet or Hex_6Tet` | 0 | Méthode non supportée (ex. `Hex_24Tet`). |
| `volumes, Inlet and Outlet must not be empty` | 0 | Une des entrées est vide. |
| `Inlet and Outlet have different numbers of faces` | 0 | Tailles différentes. |
| `invalid volume selection for periodicity` | 1 | Volume absent du maillage, ou demande de découpage en prismes. |
| `periodicity supports only linear tetrahedra, hexahedra, prisms and pyramids` | 1 | Élément quadratique, polyèdre, etc. |
| `invalid periodic volume` | 1 | `SMDS_VolumeTool::Set` a échoué. |
| `periodic selection leaves an unsplit volume across a quadrangular face` | 1 | Un voisin non sélectionné partage un quadrangle. |
| `Inlet and Outlet must contain linear triangles or quadrangles of this mesh` | 2 | Face invalide, quadratique ou d'un autre maillage. |
| `repeated node in a periodic face` | 2 | Face dégénérée. |
| `duplicate periodic face` | 2 | Même face deux fois dans un côté. |
| `Inlet and Outlet must be distinct surfaces` | 2 | Une face est dans les deux côtés. |
| `a periodic face is not on the boundary of the selected volumes` | 2 | Propriétaire non sélectionné, ou face intérieure. |
| `a periodic face has incompatible volume connectivity` | 2 | Ordre des nœuds incohérent entre la face 2D et la face du volume. |
| `a periodic quadrangle already has adjacent tetrahedra` | 2 | Diagonale déjà fixée par un tétraèdre existant. |
| `each periodic face must have one selected owner` | 2 | Zéro ou plusieurs volumes propriétaires. |
| `Inlet and Outlet have different numbers of nodes` | 2 | Nombre de nœuds différent. |
| `missing or ambiguous translated node on Outlet` | 4 | Pas d'image unique d'un nœud à la tolérance. |
| `periodic face connectivities do not match` | 5 | Pas de face outlet correspondant à une face inlet. |
| `periodic faces have incompatible node ordering` | 5 | Bons nœuds, mais ordre cyclique différent. |
| `conflicting periodic and adjacent diagonals` | `getTetraSplitMethod` | Garde-fou (voir §6). |
