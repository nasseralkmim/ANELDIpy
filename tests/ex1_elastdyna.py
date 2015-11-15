import numpy as np
from eldypy import elastodynamics2d

mesh_name = 'patch'

material = {'E-nu-rho': [1000.0, 0.2, 1.0]}


def body_forces(x1, x2):
    return np.array([
        0.0,
        0.0,
    ])


def traction_imposed(x1, x2):
    return {
        ('line', 3): [-1.0, 0.0],
        ('line', 1): [1.0, 0.0]
    }


def displacement_imposed(x1, x2):
    return {
        ('node', 0): [0.0, 0.0],
        ('node', 1): ['free', 0.0]
    }


period = 0.1
steps = 300
u0 = 0.0
v0 = 0.0

elastodynamics2d.solver(mesh_name,
                        material,
                        body_forces,
                        traction_imposed,
                        displacement_imposed,
                        period, steps, u0, v0,
                        plot_undeformed={'Domain': True,
                                         'Elements': True,
                                         'NodeLabel': True,
                                         'EdgesLabel': False,
                                         'ElementLabel': False,
                                         'SurfaceLabel': False},
                        plot_stress={'s11': False,
                                     's22': False,
                                     's12': False,
                                     'sPmax': False,
                                     'sPmin': False},
                        plot_deformed={'DomainUndeformed': False,
                                       'ElementsUndeformed': True,
                                       'DomainDeformed': False,
                                       'ElementsDeformed': True,
                                       'DeformationMagf': 100},
                        )
