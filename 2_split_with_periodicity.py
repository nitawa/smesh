"""Importez le MED, selectionnez son maillage dans SALOME, puis executez ce script.

Le maillage doit contenir deux groupes de faces nommes Inlet et Outlet.
Le decoupage modifie directement le maillage selectionne.
"""

import SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
import salome

salome.salome_init()


def get_face_group(mesh, name):
    groups = mesh.GetGroupByName(name, SMESH.FACE)
    if len(groups) != 1:
        raise RuntimeError(
            "Un unique groupe de faces '%s' est requis ; %d trouve(s)."
            % (name, len(groups))
        )
    group = groups[0]
    if group.Size() == 0:
        raise RuntimeError("Le groupe de faces '%s' est vide." % name)
    return group


print("Imoprt mesh")
# get the selected mesh from file
# ([Mesh_1], status) = smesh.CreateMeshesFromMED(r'Mesh_1_before_split.med')

# get the selected mesh from selected in app
selected_mesh_entry = salome.sg.getSelected(0)
if not selected_mesh_entry:
    raise RuntimeError("Selectionnez le maillage importe dans SALOME.")
selected_mesh = salome.IDToObject(selected_mesh_entry)
if selected_mesh is None or selected_mesh._narrow(SMESH.SMESH_Mesh) is None:
    raise RuntimeError("L'objet selectionne doit etre un maillage SMESH.")
Mesh_1 = smesh.Mesh(selected_mesh)  # Convert MeshProxy into mesh class

Inlet = get_face_group(Mesh_1, "Inlet")
Outlet = get_face_group(Mesh_1, "Outlet")

print("Split volumes avec periodicite entre Inlet et Outlet")
Mesh_1.SplitVolumesIntoTetra(
    Mesh_1,
    method=smesh.Hex_6Tet,
    avoidOverConstrainedVolumes=False,
    ensurePeriodicity=True,
    inlet=Inlet,
    outlet=Outlet,
)

print("Split reussi avec: %d tetraedres" % Mesh_1.NbTetras())
if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()
