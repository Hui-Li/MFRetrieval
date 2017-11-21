import numpy as np
from sklearn.preprocessing import normalize

q = np.random.rand(100, 50)
q = normalize(q, axis=1, norm='l2')

np.savetxt("q.txt", q, delimiter=",")

p = np.random.rand(100, 50)
p = normalize(p, axis=1, norm='l2')

np.savetxt("p.txt", p, delimiter=",")