import streamlit as st

import numpy as np
import stl
import math
from stl import mesh  # pip install numpy-stl

import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy
from urllib.request import urlopen
from urllib.parse import quote

def CIRconvert(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'C'

def getH(cx,cy,cz):
  scale=1.6
  meshH=mesh.Mesh.from_file('esferauno.stl')
  meshH.x += -0.4
  meshH.y += -0.4
  meshH.z += -0.5

  meshH.x *= scale
  meshH.y *= scale
  meshH.z *= scale


  meshH.x += cx
  meshH.y += cy
  meshH.z += cz
  return meshH

def getC(cx,cy,cz):
  scale=2.4
  meshC=mesh.Mesh.from_file('esferauno.stl')
  meshC.x += -0.4
  meshC.y += -0.4
  meshC.z += -0.5
  meshC.x *= scale
  meshC.y *= scale
  meshC.z *= scale
  meshC.x += cx
  meshC.y += cy
  meshC.z += cz

  return meshC

#exmp=st.sidebar.selectbox('',['CCCCCc1cc(c2c(c1)OC([C@H]3[C@H]2C=C(CC3)C)(C)C)O','CCCCCCCC','OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2'])

#smi=st.text_input("smile", exmp)
st.markdown('# STL files from molecules')
st.markdown('## Ready for 3D Printing 🖨️')
ids=st.sidebar.text_input("identifier, name, SMILES, etc", 'Caffeine')

smi=CIRconvert(ids)
st.sidebar.write('SMILES: '+smi)
m=Chem.MolFromSmiles(smi)
m=Chem.AddHs(m)
AllChem.EmbedMolecule(m,randomSeed=0xf00d)
mol = m
atomos=[]
for atom in m.GetAtoms():
    #print( atom.GetSymbol())
    atomos.append(atom.GetSymbol())

coords=[]
for i in range(0, mol.GetNumAtoms()):
  pos = mol.GetConformer().GetAtomPosition(i)
  #print ('{0:12.4f}{1:12.4f}{2:12.4f}'.format(pos.x, pos.y, pos.z))
  coords.append([pos.x, pos.y, pos.z])

meshes=[]
for a in range(len(atomos)):
  #print(atomos[a])
  if atomos[a]=='H':
    meshes.append(getH(coords[a][0],coords[a][1],coords[a][2]))
  else:
    meshes.append(getC(coords[a][0],coords[a][1],coords[a][2]))

mesh_list=[]
for i in meshes:
  mesh_list.append(i.data)

filename=ids+'.stl'
combined = mesh.Mesh(numpy.concatenate(mesh_list))
#combined.save('combined.stl', mode=stl.Mode.ASCII)  
combined.save(filename, mode=stl.Mode.ASCII) 




with open(filename, 'rb') as my_file:
    st.sidebar.download_button(label = 'Download .stl file', data = my_file, file_name = filename, )   
