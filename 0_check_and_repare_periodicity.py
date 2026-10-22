import math
import numpy as np

import builtins

import SalomePyQt
sgPyQt=SalomePyQt.SalomePyQt()

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()

import salome
salome.salome_init()


def checkPeriodicity(mesh, group_1, group_2, translation_coords):
    if group_1.GetNumberOfNodes() == group_2.GetNumberOfNodes():
        Face_translated = mesh.TranslateObjectMakeMesh(group_1, translation_coords, 0, group_1.GetName() + "_Translated")
        Compound_Mesh = smesh.Concatenate([Face_translated.GetMesh(), group_2], 1, 0, 1e-06, False)
        Compound_Mesh.SetName("%s_%s_%s"%(mesh.GetName(), group_1.GetName(), group_2.GetName()))
        coincident_nodes = Compound_Mesh.FindCoincidentNodesOnPart([Compound_Mesh], 1e-06, [], 0)

        if len(coincident_nodes) != Face_translated.NbNodes():
            print(f"Nodes {mesh.GetName()} : {group_1.GetName()} and {group_2.GetName()} are not periodic.")
            return False
        else:
            Compound_Mesh.MergeNodes(coincident_nodes,[],0)
            equal_elements = Compound_Mesh.FindEqualElements([Compound_Mesh], [])
            if len(equal_elements) != Face_translated.NbElements():
                print(f"Elements {mesh.GetName()} : {group_1.GetName()} and {group_2.GetName()} are not periodic.")
                return False
            else:
                print(f"{mesh.GetName()} : {group_1.GetName()} and {group_2.GetName()} are periodic.")
                Compound_Mesh.MergeElements(equal_elements,[])
                return True


#### FIXING PERIODICITY

def repare_periodicity(mesh, group_1, group_2, translation_vector, tol=1e-7):

  repared = False

  group_1_nodes = group_1.GetNodeIDs()
  group_2_nodes = group_2.GetNodeIDs()

  assert len(group_1_nodes) == len(group_2_nodes), "group_1 and group_2 do not have the same number of nodes. We can't fix periodicity"

  # for each group_1 node, compute the opposite group_2 node coordinates
  # and move the closest node to the exact coord
  for node in group_1_nodes:
    # Gets the node coodinates
    x, y, z = mesh.GetNodeXYZ(node)
    # Translates node coordinates
    x_opp = x + translation_vector[0]
    y_opp = y + translation_vector[1]
    z_opp = z + translation_vector[2]
    # Finds the closest node to the translated node
    closest_opp_node = mesh.FindNodeClosestTo(x_opp, y_opp, z_opp)
    # Checking for error
    if closest_opp_node not in group_2_nodes:
      mesh.MakeGroupByIds("node_%i"%node, SMESH.NODE, [node])
      mesh.MakeGroupByIds("opp_node_%i"%closest_opp_node, SMESH.NODE, [closest_opp_node])
      raise Exception("Opposite node is not in top_nodes... something is wrong")
    # Gets the coordinates of the closest node
    x_closest, y_closest, z_closest = mesh.GetNodeXYZ(closest_opp_node)
    # Calculates the distance
    dist = math.sqrt((x_closest-x_opp)**2+(y_closest-y)**2+(z_closest-z_opp)**2)
    
    if dist > tol:
      mesh.MoveNode(closest_opp_node, x_opp, y, z_opp)
      repared = True

  if repared:
    print(f"{mesh.GetName()} : Periodicity between {group_1.GetName()} and {group_2.GetName()} has been repared.")



print("Importing mesh")
# get the selected mesh from selected in app

selected_mesh_entry = salome.sg.getSelected(0)
selected_mesh = salome.IDToObject(selected_mesh_entry)

# Convert MeshProxy into mesh class
Mesh_1 = smesh.Mesh(selected_mesh)

d_groups = {}
mesh_groups = Mesh_1.GetGroups()
for gr in mesh_groups:
    name = gr.GetName()
    d_groups[name] = gr

#First define translation vector that puts Inlet onto Outlet
print("WARNING: Be sure that translation vector is correct.")
translation_vector = [-2.99761, 0, -3] # Automatiquement à partid du "barycentre" de Inlet/Outlet ?

p1 = checkPeriodicity(Mesh_1, d_groups["Inlet"], d_groups["Outlet"], translation_vector)
if p1:
   print(f"{Mesh_1.GetName()} is periodic.")
else:
   repare_periodicity(Mesh_1, d_groups["Inlet"], d_groups["Outlet"], translation_vector)
   p2 = checkPeriodicity(Mesh_1, d_groups["Inlet"], d_groups["Outlet"], translation_vector)
   if p2:
      print(f"After reparation, {Mesh_1.GetName()} is periodic.")
   else:
      print(f"{Mesh_1.GetName()} groups may be too different to be repaired.")

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()
