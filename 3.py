import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

spatial_nodes = np.loadtxt('x.txt')
time_steps = np.loadtxt('t.txt')
numerical_solution = np.loadtxt('matrix.txt')

# Создание сетки для визуализации
x, t = np.meshgrid(spatial_nodes, time_steps)

# Аналитическое/точное решение
# analytical_solution = np.exp(-np.pi**2 * t) * x * (x - 1)
# analytical_solution = np.exp(-np.pi**2 * t) * np.sin(np.pi * x)
analytical_solution = np.exp(-np.pi**2 * t) * np.sin(x) * (x - 1)

solution_figure = plt.figure(figsize=(18, 6))

# График численного решения
num_sol_ax = solution_figure.add_subplot(131, projection='3d')
surface_plot = num_sol_ax.plot_surface(
    x, 
    t, 
    numerical_solution, 
    cmap='summer',
    edgecolor='none'
)
num_sol_ax.set_xlabel('Пространственная координата')
num_sol_ax.set_ylabel('Временная координата')
num_sol_ax.set_title('Численное решение')
solution_figure.colorbar(surface_plot, ax=num_sol_ax, shrink=0.5)

# График аналитического/точного решения
exact_sol_ax = solution_figure.add_subplot(132, projection='3d')
exact_plot = exact_sol_ax.plot_surface(
    x, 
    t, 
    analytical_solution, 
    cmap='autumn',
    edgecolor='none'
)
exact_sol_ax.set_xlabel('Пространственная координата')
exact_sol_ax.set_title('Аналитическое решение')
solution_figure.colorbar(exact_plot, ax=exact_sol_ax, shrink=0.5)

# График ошибки
error_ax = solution_figure.add_subplot(133, projection='3d')
error_plot = error_ax.plot_surface(
    x, 
    t, 
    numerical_solution - analytical_solution, # solution_error
    cmap='winter',
    edgecolor='none'
)
error_ax.set_xlabel('Пространственная координата')
error_ax.set_title('Ошибка решения')
solution_figure.colorbar(error_plot, ax=error_ax, shrink=0.5)

plt.tight_layout()
solution_figure.suptitle('Сравнение решений уравнения теплопроводности', y=1.02)
plt.savefig("3.png")
plt.show()
