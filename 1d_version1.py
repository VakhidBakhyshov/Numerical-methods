import pandas as pd
import matplotlib.pyplot as plt

file_path = 'out.txt'
df = pd.read_csv(file_path, sep='\\s+', header=0, names=['x', 'y1', 'y2'])

plt.figure(figsize=(12, 5))
plt.plot(df['x'], df['y1'], marker='o', linestyle='-', color='b', label='y1')
plt.plot(df['x'], df['y2'], marker='s', linestyle='--', color='r', label='y2')

plt.title('Graphs')
plt.xlabel('x')
plt.ylabel('Values')
plt.xlabel("(exp(x) - e) * sin(x), N = 10")

plt.legend()

plt.savefig('graphs.png')
plt.show()
