import numpy as np
import matplotlib.pyplot as plt

class Nucleus:
    # Константы для формулы Вайцзекера
    a1 = 15.75  # МэВ
    a2 = 17.8   # МэВ
    a3 = 0.711  # МэВ
    a4 = 23.7   # МэВ
    r0 = 1.2    # Фемтометр

    def __init__(self, name, A, Z):
        self.name = name
        self.A = A  # Массовое число
        self.Z = Z  # Число протонов
        self.N = A - Z  # Число нейтронов

    def binding_energy(self):
        """Вычисляет удельную энергию связи для ядра."""
        d = 0
        
        # Определяем значение a5 в зависимости от четности A и Z
        if self.A % 2 == 0 and self.Z % 2 == 0:
            d = 12.0  # четно-четное
        elif self.A % 2 == 1 and self.Z % 2 == 1:
            d = -12.0  # нечетно-нечетное
        else:
            d = 0      # четно-нечетное
        
        B = (self.a1 * self.A - 
             self.a2 * self.A**(2/3) - 
             self.a3 * self.Z**2 / self.A**(1/3) - 
             self.a4 * ((self.A - 2 * self.Z)**2) / self.A + 
             d * self.A**(-3/4))
        
        return B

    def mass(self):
        """Вычисляет массу атома."""
        m_proton = 1.007276466812  # Масса протона в атомных единицах
        m_neutron = 1.00866491588   # Масса нейтрона в атомных единицах
        
        binding_energy_value = self.binding_energy()
        
        M = (self.Z * m_proton + 
             (self.A - self.Z) * m_neutron - 
             binding_energy_value / (931.5))  # Преобразуем МэВ в атомные единицы
        return M

    def radius(self):
        """Вычисляет радиус атома."""
        R = self.r0 * self.A**(1/3)
        return R

    def is_stable(self):
        """Определяет устойчивость изотопа к бета-распаду."""
        return self.N <= self.Z + 1

    def can_split_even_even(self):
        """Проверяет возможность деления на два четно-четных осколка."""
        return self.A % 2 == 0


def main():
    # Примеры ядер для анализа
    nuclei_data = {
        'U-238': (238, 92),
        'Pu-239': (239, 94),
        'Cf-252': (252, 98),
        'Pu-238': (238, 94),
        'Te-135': (135, 52),
    }

    results = {}

    for name, (A, Z) in nuclei_data.items():
        nucleus = Nucleus(name, A, Z)
        
        results[name] = {
            'Binding Energy (MeV)': nucleus.binding_energy(),
            'Mass (atomic mass units)': nucleus.mass(),
            'Radius (fm)': nucleus.radius(),
            'Stable': nucleus.is_stable(),
            'Can Split Even-Even': nucleus.can_split_even_even()
        }

    # Вывод результатов
    for nucleus, result in results.items():
        print(f"{nucleus}:")
        for key, value in result.items():
            print(f"   {key}: {value}")
    
    # Графики
    A_values = np.array([data[0] for data in nuclei_data.values()])
    Z_values = np.array([data[1] for data in nuclei_data.values()])
    
    radii = np.array([Nucleus(name, A, Z).radius() for name, (A, Z) in nuclei_data.items()])
    binding_energies = np.array([Nucleus(name, A, Z).binding_energy() for name, (A, Z) in nuclei_data.items()])

    plt.figure(figsize=(12, 6))

    # График радиуса в зависимости от массового числа A
    plt.subplot(1, 2, 1)
    plt.plot(A_values, radii, marker='o')
    plt.title('Радиус атома в зависимости от массового числа A')
    plt.xlabel('Массовое число A')
    plt.ylabel('Радиус (фм)')
    plt.grid()

    # График удельной энергии связи в зависимости от заряда Z
    plt.subplot(1, 2, 2)
    plt.plot(Z_values, binding_energies / A_values, marker='o')
    plt.title('Удельная энергия связи в зависимости от заряда Z')
    plt.xlabel('Заряд Z')
    plt.ylabel('Удельная энергия связи (МэВ)')
    plt.grid()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
