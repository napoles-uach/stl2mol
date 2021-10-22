import streamlit as st
import plotly
import numpy as np
import stl
import math
from stl import mesh  # pip install numpy-stl
import plotly.graph_objects as go
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
st.markdown('## Ready for 3D Printing or AR')
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


# Put your Python+Streamlit code here ...
# you can modify it by double cliking on the folder icon at the left
def stl2mesh3d(stl_mesh):
    # stl_mesh is read by nympy-stl from a stl file; it is  an array of faces/triangles (i.e. three 3d points) 
    # this function extracts the unique vertices and the lists I, J, K to define a Plotly mesh3d
    p, q, r = stl_mesh.vectors.shape #(p, 3, 3)
    # the array stl_mesh.vectors.reshape(p*q, r) can contain multiple copies of the same vertex;
    # extract unique vertices from all mesh triangles
    vertices, ixr = np.unique(stl_mesh.vectors.reshape(p*q, r), return_inverse=True, axis=0)
    I = np.take(ixr, [3*k for k in range(p)])
    J = np.take(ixr, [3*k+1 for k in range(p)])
    K = np.take(ixr, [3*k+2 for k in range(p)])
    return vertices, I, J, K
#my_mesh = mesh.Mesh.from_file('combined.stl')
my_mesh = mesh.Mesh.from_file(filename)
vertices, I, J, K = stl2mesh3d(my_mesh)
x, y, z = vertices.T

colorscale= [[0, '#e5dee5'], [1, '#e5dee5']]  

mesh3D = go.Mesh3d(
            x=x,
            y=y,
            z=z, 
            i=I, 
            j=J, 
            k=K, 
            flatshading=True,
            colorscale=colorscale, 
            intensity=z, 
            name='AT&T',
            showscale=False)


title = ""
layout = go.Layout(paper_bgcolor='rgb(0.999,0.999,0.999)',
            title_text=title, title_x=0.5,
                   font_color='brown',
            width=800,
            height=800,
            scene_camera=dict(eye=dict(x=1.25, y=-1.25, z=1)),
            scene_xaxis_visible=False,
            scene_yaxis_visible=False,
            scene_zaxis_visible=False)

fig = go.Figure(data=[mesh3D], layout=layout)
fig.data[0].update(lighting=dict(ambient= 0.18,
                                 diffuse= 1,
                                 fresnel=  .1,
                                 specular= 1,
                                 roughness= .1,
                                 facenormalsepsilon=0))
fig.data[0].update(lightposition=dict(x=3000,
                                      y=3000,
                                      z=10000));

fig.show()
raw_html=fig.to_html()
#raw_html = base64.b64encode(raw_html).decode()
components.html(raw_html,height=900)

with open(filename, 'rb') as my_file:
    st.sidebar.download_button(label = 'Download .stl file', data = my_file, file_name = filename, )   
