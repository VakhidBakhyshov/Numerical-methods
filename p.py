import numpy as np
import matplotlib.pyplot as plt

N1 = 10
err1 = 0.000428483

N2 = 50
err2 = 2.11665e-05

N3 = 100
err3 = 5.23754e-06

N4 = 500
err4 = 2.06797e-07

N5 = 1000
err5 = 5.16007e-08

p1 = np.log(err1/err2)/np.log(N2/N1)
p1 = np.log(err2/err3)/np.log(N3/N2)

plt.figure(figsize=(10, 6))
plt.scatter([-np.log(1/N1), -np.log(1/N2), -np.log(1/N3), -np.log(1/N4), -np.log(1/N5)], [-np.log(err1), -np.log(err2), -np.log(err3), -np.log(err4), -np.log(err5)], color='blue')
plt.plot([-np.log(1/N1), -np.log(1/N2), -np.log(1/N3), -np.log(1/N4), -np.log(1/N5)], [-np.log(err1), -np.log(err2), -np.log(err3), -np.log(err4), -np.log(err5)], color='blue')

plt.xlabel('-np.log(1/N)')
plt.ylabel('-np.log(err)')
