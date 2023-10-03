import math
import numpy as np

a = np.array([[1],[2],[3]])
b = np.array([[2],[8],[4]])
print(np.cross(a,b,axis=0))
print(np.cross(a,b,axis=1))

