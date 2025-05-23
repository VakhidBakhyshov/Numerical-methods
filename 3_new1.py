import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

scheme_1 = pd.read_csv("rate1.txt", sep=r'\s+')

scheme_1['Convergence rate'] = (np.log(scheme_1['Norma'].shift(1) / scheme_1['Norma'])) / (np.log(scheme_1['M'] / scheme_1['M'].shift(1))) 

print(scheme_1)

scheme_1_log_h = np.log(scheme_1['h'])
scheme_1_log_error = np.log(scheme_1['Norma'])
scheme_1_log_tau = np.log(scheme_1['tau'])

# Сходимость по x
plt.subplot(121)
plt.plot(scheme_1_log_h, scheme_1_log_error, 'ro-', label='Явная схема')
plt.plot(scheme_1_log_h, 2 * scheme_1_log_h, 'k--', label='Теоретическая O(h²)')

plt.xlabel('log(h)')
plt.ylabel('log(error)')
plt.title('Пространственная сходимость')
plt.legend()
plt.grid(True)

# Сходимость по t
plt.plot(scheme_1_log_tau, scheme_1_log_error, 'b^-', label='Неявная схема')
plt.plot(scheme_1_log_tau, scheme_1_log_tau, 'k--', label='Теоретическая O(τ)')
plt.xlabel('log(τ)')
plt.ylabel('log(error)')
plt.title('Временная сходимость')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

#plt.plot(scheme_2['error'][1:], scheme_2['Convergence rate'])
 
