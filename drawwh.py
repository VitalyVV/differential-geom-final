from pythreejs import *
import colorsys
from math import atan
from IPython.display import display
import solution

def obj_read(filename):
    with open(filename,'r') as obj:
        lines = [ [f for f in s.split(' ') if len(f)>0] for s in obj.read().split('\n') ]
    vertices = [[float(coord) for coord in l[1:4]] for l in lines if len(l)>3 and l[0]=='v']
    faces = [[int(coord.split('/')[0])-1 for coord in l[1:4]] for l in lines if len(l)>3 and l[0]=='f']
    return faces, vertices



def heatmap(value):
    return '#{:02x}{:02x}{:02x}'.format(*[int(i*255) for i in colorsys.hls_to_rgb(0.5, .3-atan(value)/15, 1)])


def draw(faces, vertices, colorfunc=None):
    # Create the geometry:
    if colorfunc is None:
         colorfunc=lambda x:int(x)%10 - 5

    vertexcolors = 0
    if not callable(colorfunc):
        vertexcolors=[heatmap(colorfunc[i]) for  i in  range(len(vertices))]
    else:
        vertexcolors=[heatmap(colorfunc(i)) for  i in  range(len(vertices))]
    faces = [f + [None, [vertexcolors[i] for i in f], None] for f in faces]
    geometry = Geometry(faces=faces,vertices=vertices,colors=vertexcolors)
    # Calculate normals per face, for nice crisp edges:
    #geometry.exec_three_obj_method('computeFaceNormals')

    object1 = Mesh(
        geometry=geometry,
        material=MeshLambertMaterial( side="FrontSide", vertexColors='VertexColors'),
    )

    object2 = Mesh(
        geometry=geometry,
        material=MeshLambertMaterial(color= "gray", side="BackSide"),
    )

# Set up a scene and render it:
    camera = PerspectiveCamera(position=[2*max(v[0] for v in vertices), 2*max(v[1] for v in vertices), 2*max(v[2] for v in vertices)], fov=40,
                      children=[DirectionalLight(color='#cccccc', position=[-3, 5, 1], intensity=0.5)])
    scene = Scene(children=[object1, object2, camera, AmbientLight(color='#dddddd')])

    renderer = Renderer(camera=camera, background='black', background_opacity=1,
                        scene=scene, controls=[OrbitControls(controlling=camera)])
    display(renderer)


if __name__=='__main__':
    import scipy.sparse.linalg as lin
    from functools import reduce
    import numpy as np

    def kernel(x, y, t):
        integral = 0
        for v, f in zip(e_values, e_vectors.T):
            ex = np.exp(t * v)
            if np.isinf(ex):
                integral += f[x] * f[y]
            else:
                integral += ex * f[x] * f[y]
        return integral

    def heat_function(x, t, f):
        integral = 0
        heat = np.zeros(mesh.n)
        for y in range(mesh.n):
            heat[y] = kernel(x, y, t) * f[y]

        return heat


    mesh = solution.Mesh.fromobj("oloid64_tri.obj")
    mesh.draw()
    mesh.LaplaceOperator()
    L = mesh.laplace_operator
    e_values, e_vectors = lin.eigs(L)
    print(f'eigen values:\n{e_values}')
    print(f'eigen vectors:\n{e_vectors.shape}')

    f_init = np.ones(mesh.n)
    x = 0
    for t in np.arange(1.2, 2, 0.01):
        print(t)
        heat = heat_function(x, t, f_init)
        mesh.draw(heat)
