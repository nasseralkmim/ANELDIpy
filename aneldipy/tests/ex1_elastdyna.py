__author__ = 'Nasser'

import numpy as np
import elastodynamics2d

meshName = 'patch'

material = {'E-nu-rho':[1000.0, 0.2, 1.0]}

def body_forces(x1, x2):
    return np.array([
        0.0,
        0.0,
    ])

def traction_imposed(x1, x2):
    return {
        ('line', 3):[-1.0, 0.0],
        ('line', 1):[1.0, 0.0]
    }

def displacement_imposed(x1, x2):
    return {
        ('node', 0):[0.0, 0.0],
        ('node', 1):['free', 0.0]

    }

elastodynamics2d.solver(meshName,
                    material,
                    body_forces,
                    traction_imposed,
                    displacement_imposed,
                    plotUndeformed={'Domain':True,
                                  'Elements':True,
                                  'NodeLabel':True,
                                  'EdgesLabel':False,
                                  'ElementLabel':False,
                                  'SurfaceLabel':False},
                    plotStress={'s11':False,
                              's22':False,
                              's12':False,
                              'sPmax':False,
                              'sPmin':False},
                    plotDeformed={'DomainUndeformed':False,
                                'ElementsUndeformed':False,
                                'DomainDeformed':False,
                                'ElementsDeformed':False,
                                'DeformationMagf': 100},
                    printReport={'U':False,
                                 'stress':False,
                                 'strain':False}
                    )