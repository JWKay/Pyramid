import numpy as np
import scipy
import math
import dit
from dit.pid.idep import  PID_dep
from dit.pid.iccs import PID_CCS
from dit.pid.imin import PID_WB
from dit.pid.ibroja import PID_BROJA

from dit.pid.iproj import PID_Proj


from dit import Distribution
from dit.multivariate import entropy

# Enter the non-zero frequency data from table in Figure4.R

mat = np.loadtxt('/Users/jk/Desktop/KPAGL/freqs.txt')

output = np.zeros((16, 26))

for i in range(16):

    print(i)

    freq = mat[i,]

    aa = freq / np.sum(freq)

    dat = np.reshape(aa, (4, 4, 3))

    dist1 = Distribution.from_ndarray(dat)


    p = PID_BROJA(dist1, [[0], [1]], [2])

    output[i, 0:4] =  [p.get_partial(n) for n in sorted(p._lattice, key=dit.pid.lattice.sort_key(p._lattice))]

    p = PID_WB(dist1, [[0], [1]], [2])

    output[i, 4:8] = [p.get_partial(n) for n in sorted(p._lattice, key=dit.pid.lattice.sort_key(p._lattice))]

    p = PID_Proj(dist1, [[0], [1]], [2])

    output[i, 8:12] = [p.get_partial(n) for n in sorted(p._lattice, key=dit.pid.lattice.sort_key(p._lattice))]

    p = PID_CCS(dist1, [[0], [1]], [2])

    output[i, 12:16] = [p.get_partial(n) for n in sorted(p._lattice, key=dit.pid.lattice.sort_key(p._lattice))]

    p = PID_dep(dist1, [[0], [1]], [2])

    output[i, 16:20] =  [p.get_partial(n) for n in sorted(p._lattice, key=dit.pid.lattice.sort_key(p._lattice))]






# Computation of the entropies

    H_OBA = entropy(dist1)
    H_B= entropy(dist1, [0])
    H_A= entropy(dist1, [1])
    H_O=entropy(dist1, [2])

    H_OB=entropy(dist1, [[0], [2]])
    H_OA=entropy(dist1, [[1], [2]])
    H_BA=entropy(dist1, [[0], [1]])

# Computation of the mutual informations and interaction information

    I_OB = H_B + H_O - H_OB
    I_OA = H_A + H_O - H_OA
    I_BA = H_B + H_A - H_BA
    I_OBgA = H_OA - H_A - H_OBA + H_BA
    I_OAgB = H_OB - H_B - H_OBA + H_BA
    I_OBA = H_O + H_BA - H_OBA

    II_OBA =  I_OBA -  I_OB - I_OA

    output[i, 20:26] = [I_OB, I_OA, I_OBgA, I_OAgB, I_OBA, II_OBA]


np.savetxt("/Users/jk/Desktop/KPAGL/PIDout.txt", output)

