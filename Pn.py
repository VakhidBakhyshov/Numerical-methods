import pandas as pd
import matplotlib.pyplot as plt

file_path = 'out1.txt'  
df = pd.read_csv(file_path, sep='\s+', header=0, names=['x', 'f(x)', 'Pn', 'err'])

plt.figure(figsize=(9, 4))
plt.plot(df['x'][:-1], df['f(x)'][:-1], marker='o', linestyle='-', color='b', label='f(x)')
plt.plot(df['x'][:-1], df['Pn'][:-1], marker='s', linestyle='--', color='r', label='Pn')


plt.title('Graphs Pn')
plt.xlabel('x')
plt.ylabel('Values')

plt.legend()

#plt.savefig('graphsPn.png')
plt.show()
