import matplotlib.pyplot as plt
import numpy as np

def load_convergence_data(input_file):
    raw_data = np.genfromtxt(input_file)
    x_vals = raw_data[:, 0]  # Логарифмы количества точек
    y_vals = raw_data[:, 1:] # Логарифмы ошибок для каждого метода
    return x_vals, y_vals

def visualize_methods_convergence(x_data, y_data):
    fig, ax = plt.subplots(figsize=(12, 8))
    
    method_names = ['Явный метод Эйлера (1ая схема)', 
                  'Неявный метод Эйлера (2ая схема)', 
                  'Метод трапеций (3ая схема)', 
                  'Двухшаговый метод (4ая схема)']
    line_styles = ['b-o', 'r--s', 'g-.^', 'm:d']
    
    for method_idx in range(y_data.shape[1]):
        ax.plot(x_data, 
               y_data[:, method_idx], 
               line_styles[method_idx],
               linewidth=2,
               markersize=8,
               label=method_names[method_idx])
    
    ax.set(xlabel='log(Количество узлов сетки)',
          ylabel='log(Погрешность)',
          title='Сравнение сходимости численных методов')
    
    ax.grid(which='both', linestyle=':', alpha=0.7)
    ax.legend(fontsize=10, framealpha=0.9)
    
    plt.savefig('2.png')
    plt.show()

def convergence_rate(errors, N_array):
    rate = []
    for i in range(1, len(errors)):
        rate_i = np.log(errors[i - 1]/errors[i])/np.log((N_array[i] - 1) / (N_array[i - 1] - 1))
        rate.append(rate_i)
    return rate


input_data_file = 'output.txt'
log_points, log_errors = load_convergence_data(input_data_file)
visualize_methods_convergence(log_points, log_errors)
rate_1_scheme = convergence_rate(log_errors[:, 0], log_points)
rate_2_scheme = convergence_rate(log_errors[:, 1], log_points)
rate_3_scheme = convergence_rate(log_errors[:, 2], log_points)
rate_4_scheme = convergence_rate(log_errors[:, 3], log_points)

for i in range(0, len(log_points) - 1):
    print(f'Для N = {log_points[i+1]}, показатель сходимости 1ой схемы = {rate_1_scheme[i]}, 2ой схемы = {rate_2_scheme[i]}, 3ей схемы = {rate_3_scheme[i]}, 4ой схемы = {rate_4_scheme[i]}')
