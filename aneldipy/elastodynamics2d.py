import numpy as np
from scipy.sparse.linalg import spsolve
from scipy import sparse
import aneldipy.gmsh as gmsh
import aneldipy.element2dof as element2dof
import aneldipy.assemble2dof as assemble2dof
import aneldipy.boundaryconditions2dof as boundaryconditions2dof
import aneldipy.processing as processing
import aneldipy.output as output


def solver(mesh_name, material, body_forces, traction_imposed,
           displacement_imposed,
           plot_undeformed, plot_stress, plot_deformed, print_report, **kwargs):

    mesh = gmsh.parse(mesh_name)

    s = mesh.surfaces
    mat_dic = {s[i]: material[j] for i, j in enumerate(material)}

    k_ele = element2dof.stiffness(mesh, mat_dic)

    p0q_ele = element2dof.body_forces(mesh, body_forces)

    m_ele = element2dof.mass(mesh, mat_dic)

    k = assemble2dof.globalMatrix(k_ele, mesh)

    m = assemble2dof.globalMatrix(m_ele, mesh)
    print(np.size(m), m)

    p0q = assemble2dof.globalVector(p0q_ele, mesh)

    p0t = boundaryconditions2dof.neumann(mesh, traction_imposed)

    p0 = p0q + p0t

    km, p0m = boundaryconditions2dof.dirichlet(k, p0, mesh,
                                               displacement_imposed)

    ks = sparse.csc_matrix(km)

    u = spsolve(ks, p0m)

    s_node, s_ele, e_ele = processing.stress_recovery_simple2(mesh, u,
                                                              mat_dic)
    print(mesh.num_ele)

    principal_max = processing.principal_stress_max(s_node[0], s_node[1],
                                                    s_node[2])
    principal_min = processing.principal_stress_min(s_node[0], s_node[1],
                                                    s_node[2])

    output.data(plot_stress, plot_deformed, plot_undeformed,
                print_report, mesh, s_ele, e_ele, s_node,
                principal_max, principal_min, u)
