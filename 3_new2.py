import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

scheme_2 = pd.read_csv("rate2.txt", sep=r'\s+')

scheme_2['Convergence rate'] = (np.log(scheme_2['Norma'].shift(1) / scheme_2['Norma'])) / (np.log(scheme_2['M'] / scheme_2['M'].shift(1))) 

print(scheme_2)

scheme_2_log_h = np.log(scheme_2['h'])
scheme_2_log_error = np.log(scheme_2['Norma'])
scheme_2_log_tau = np.log(scheme_2['tau'])

plt.figure(figsize=(12, 5))

# Сходимость по x
plt.subplot(121)
plt.plot(scheme_2_log_h, scheme_2_log_error, 'bo-', label='Неявная схема')
plt.plot(scheme_2_log_h, scheme_2_log_h, 'k--', label='Теоретическая O(τ)')
plt.plot(scheme_2_log_h, 2 * scheme_2_log_h, label='Теоретическая O(h²)')
plt.xlabel('log(h)')
plt.ylabel('log(error)')
plt.title('Пространственная сходимость')
plt.legend()
plt.grid(True)

# Сходимость по t
plt.plot(scheme_2_log_tau, scheme_2_log_error, 'b^-', label='Неявная схема')
plt.plot(scheme_2_log_tau, scheme_2_log_tau, 'k--', label='Теоретическая O(τ)')
plt.xlabel('log(τ)')
plt.ylabel('log(error)')
plt.title('Временная сходимость')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

#plt.plot(scheme_2['error'][1:], scheme_2['Convergence rate'])
 
