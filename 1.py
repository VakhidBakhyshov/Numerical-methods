import numpy as np
import matplotlib.pyplot as plt

def process_data_file(file_path):
    sizes, errors, rates = [], [], []
    
    with open(file_path, 'r') as data_source:
        for record in data_source:
            if not record.strip():
                continue
                
            parts = record.split()
        
            if len(parts) == 2:
                size = int(float(parts[0]))
                error = float(parts[1])
                sizes.append(size)
                errors.append(error)
                print(f"Размер сетки = {size}, погрешность = {error}")
            
            elif len(parts) == 3:
                size = int(float(parts[0]))
                error = float(parts[1])
                rates = float(parts[2])
                sizes.append(size)
                errors.append(error)
                print(f"Размер сетки = {size}, погрешность = {error}, показатель сходимости = {rates}")

            else: continue
        
    return sizes, errors, rates

def create_visualization(grid_dimensions, error_metrics):

    plt.figure(figsize=(10, 6))  
    plt.plot(grid_dimensions, error_metrics, marker='o', linestyle='-', color='r', label='Погрешность вычислений')  
    plt.xlabel('Размер вычислительной сетки (N)')  
    plt.ylabel('Погрешность вычислений')  
    plt.title('Анализ погрешности численного метода')  
    plt.grid(True) 
    plt.legend()
    plt.savefig('1.png')
    plt.show()

    plt.figure(figsize=(10, 6))  
    plt.plot(np.log(grid_dimensions), np.log(error_metrics), marker='o', linestyle='-', color='b', label='Погрешность вычислений')  
    plt.plot(np.log(grid_dimensions), -2*np.log(grid_dimensions))
    plt.plot(np.log(grid_dimensions), -3*np.log(grid_dimensions))
    plt.xlabel('Размер вычислительной сетки (N)')  
    plt.ylabel('Погрешность вычислений')  
    plt.title('Анализ погрешности численного метода (в логарифмических координатах)')  
    plt.grid(True) 
    plt.legend()
    plt.show()


input_file = "output.txt"
mesh_sizes, calculation_errors, convergence_rates = process_data_file(input_file)
create_visualization(mesh_sizes, calculation_errors)
