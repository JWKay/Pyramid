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

# Enter the non-zero frequency data from table in Inf.R

dat = [42, 33,  2, 19, 15,  1, 14, 20,  1, 48, 12, 20,  8, 33,  7, 25, 15, 23, 25,  2,  6, 32, 40, 40, 16, 32, 40, 40, 40 ]



# Enter the category combinations in the order basal, apical, ap.

outcomes = [ '000', '010', '011', '020', '021', '022', '030', '031', '032', '100', '110', '111', '112','121', '122',
             '131', '132', '200', '201', '210', '211', '212', '222', '232' ,'301', '302', '312', '322',
             '332']


pmf = dat/np.sum(dat)

# Define the probability distribution

dist1 = Distribution(outcomes, pmf)

print("  ")
print("Probability distribution")
print("  ")
print(dist1)
print("  ")

# Compute five PIDs and print the output

print("Partial Information Decompositions")
print("  ")
print("{0}{1} is Shd, {0} is UnqB, {1} is UnqA, {0:1} is Syn")
print( "  ")

rcyBROJA = PID_BROJA(dist1, [[0], [1]], [2])

print("  ")
print(" The I_broja PID is plotted in Fig. 6 ")
print(rcyBROJA)
print("  ")

rcyCCS = PID_CCS(dist1, [[0], [1]], [2])

print(rcyCCS)

rcyDEP = PID_dep(dist1, [[0], [1]], [2])

print(rcyDEP)


rcyWB = PID_WB(dist1, [[0], [1]], [2])

print(rcyWB)


rcyPROJ = PID_Proj(dist1, [[0], [1]], [2])

print(rcyPROJ)
print("  ")

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

# Output for Table 1

print("Shannon measures for Table 1 ")
print(" ")

print(I_OB, I_OA, I_OBgA, I_OAgB, I_OBA, II_OBA)





