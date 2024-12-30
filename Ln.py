import pandas as pd
import matplotlib.pyplot as plt

file_path = 'out2.txt'  
df = pd.read_csv(file_path, sep='\s+', header=0, names=['x', 'f(x)', 'Ln', 'err'])


plt.figure(figsize=(9, 4))
plt.plot(df['x'], df['f(x)'], marker='o', linestyle='-', color='b', label='f(x)')
plt.plot(df['x'], df['Ln'], marker='s', linestyle='--', color='r', label='Ln')


plt.title('Graphs Ln')
plt.xlabel('x')
plt.ylabel('Values')

plt.legend()

#plt.savefig('graphsLn.png')
plt.show()
