'''
RDX QSS Species Finder
'''
import numpy as np

from RDX_jac import rdx_jac

Jac = rdx_jac(0.5, 538)

det_Jac = np.linalg.det(Jac)

inv_Jac = np.linalg.inv(Jac)
