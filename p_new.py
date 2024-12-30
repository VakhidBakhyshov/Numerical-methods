import numpy as np
import matplotlib.pyplot as plt

# для ф-ции x * (1 - x) * y * (1 - y)

N1 = 10
err1 = 0.000131604
N2 = 20
err2 = 3.07734e-05
N3 = 40
err3 = 7.37912e-06
N4 = 80
err4 = 1.80357e-06

p1 = np.log(err1/err2)/np.log(N2/N1)
p2 = np.log(err2/err3)/np.log(N3/N2)

plt.figure(figsize=(10, 6))
plt.scatter([-np.log(1/N1), -np.log(1/N2), -np.log(1/N3), -np.log(1/N4)], [-np.log(err1), -np.log(err2), -np.log(err3), -np.log(err4)],  color='blue')
plt.plot([-np.log(1/N1), -np.log(1/N2), -np.log(1/N3), -np.log(1/N4)], [-np.log(err1), -np.log(err2), -np.log(err3), -np.log(err4)],  color='blue')

plt.xlabel('-np.log(1/N)')
plt.ylabel('-np.log(err)')

plt.title('для ф-ции x * (1 - x) * y * (1 - y)')