import pandas as pd
import matplotlib.pyplot as plt

file_path = 'Pn.txt'  
df = pd.read_csv(file_path, sep=r'\s+', header=0, names=['x', 'f(x)', 'Pn', 'err'])

plt.figure(figsize=(9, 4))
plt.plot(df['x'][:], df['f(x)'][:], marker='o', linestyle='-', color='b', label='f(x)')
plt.plot(df['x'][:], df['Pn'][:], marker='s', linestyle='--', color='r', label='P_n')

plt.title('Graphs Pn')
plt.xlabel('x')
plt.ylabel('Values')

plt.legend()
plt.show()

