import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from numpy import ma
import matplotlib.pyplot as plt
from tkinter import messagebox

import random
from random import choices, choice
import math
from matplotlib import pyplot as plt
import numpy as np
from sympy import lambdify, simplify, cos, symbols
from statistics import mean
import copy


class AlgoritmoGenetico:

    def __init__(self,
             precision: float,
             rango: tuple,
             limiteGeneraciones: int,
             limitePoblacion: int,
             tamanioPoblcionInicial: int,
             probabilidadMutIndiv: float,
             probabilidadMutGen: float
             ):
        # --------------
        x = symbols('x')
        expresion = ((x * 2.718 ** (x/2)) * cos(x))
        # expresion = (0.75 * sin(0.50 * x) * sin(0.25 * x) * -0.75 * sin(0.75 * x))
        expresion2 = simplify(expresion)
        self.function = lambdify((x), expresion2)
        # -----------------
        self.precision = precision
        self.rango = rango
        self.limiteGeneraciones = limiteGeneraciones
        self.limitePoblacion = limitePoblacion
        self.tamanioPoblacionInicial = tamanioPoblcionInicial

        self.Rx = self.rango[1] - self.rango[0]

        self.nPx = math.ceil(self.Rx / self.precision) + 1

        self.nBx = len(bin(self.nPx)) - 2

        self.rango_i = (0, self.nPx - 1)

        self.poblacion = []
        self.mejoresCasos = []
        self.peoresCasos = []
        self.promedioCasos = []

        self.probabilidadMutIndiv = probabilidadMutIndiv
        self.probabilidadMutGen = probabilidadMutGen

        self.first_generation = []

    def mutacion_intercambio(self, individual):
        new_individual = copy.deepcopy(individual) # copia del individuo
        punto_corte_1 = random.randint(0, self.nBx-1)
        punto_corte_2 = random.randint(0, self.nBx-1)
        if punto_corte_1 > punto_corte_2:
            punto_corte_1, punto_corte_2 = punto_corte_2, punto_corte_1
        segmento = new_individual[0][punto_corte_1:punto_corte_2]
        new_individual[0][punto_corte_1:punto_corte_2] = segmento[::-1]
        return self.generarIndividuo(new_individual[0])

    def mutacion(self, individual):
        p = random.random()
        if p <= self.probabilidadMutIndiv:
            individual = self.mutacion_intercambio(individual)
        return individual


    def generarIndividuo(self, genotipo):
        i = int("".join([str(i) for i in genotipo]), 2)
        fenotipo = self.rango[0] + (i * self.precision)

        if fenotipo > self.rango[-1]:
            fenotipo = self.rango[-1]

        aptitud = self.function(fenotipo)
        return [genotipo, i,fenotipo, aptitud]

    def poda(self):
        self.poblacion = self.poblacion[:self.limitePoblacion]

    def cruza_dos_puntos(self, a, b):
        punto_corte_1 = random.randint(0, self.nBx-1)
        punto_corte_2 = random.randint(0, self.nBx-1)
        if punto_corte_1 > punto_corte_2:
            punto_corte_1, punto_corte_2 = punto_corte_2, punto_corte_1
        genotipo_a = a[0][:punto_corte_1] + b[0][punto_corte_1:punto_corte_2] + a[0][punto_corte_2:]
        genotipo_b = b[0][:punto_corte_1] + a[0][punto_corte_1:punto_corte_2] + b[0][punto_corte_2:]
        hijo_a = self.generarIndividuo(genotipo_a)
        hijo_b = self.generarIndividuo(genotipo_b)
        return hijo_a, hijo_b

    @staticmethod
    def seleccion_ruleta(poblacion):
        aptitudes = [individuo[-1] for individuo in poblacion]
        total = sum(aptitudes)
        probabilidades = [aptitud/total for aptitud in aptitudes]
        parents = choices(poblacion, probabilidades, k=2)
        return parents

    def generarPoblacionInicial(self):
        for i in range(self.tamanioPoblacionInicial):
            while True:
                genotipo = choices([0, 1], k=self.nBx)
                individual = self.generarIndividuo(genotipo)
                if self.rango_i[0] <= individual[2] <= self.rango_i[1]:
                    self.poblacion.append(individual)
                    break

    def iniciar(self, minimize: bool):
        generation = 0

        self.generarPoblacionInicial()
        self.poblacion = sorted(
            self.poblacion,
            key=lambda y: [x[3] for x in self.poblacion],
            reverse=minimize
        )
        for i in range(self.limiteGeneraciones):
            for j in range(int(len(self.poblacion) / 2)):
                padre = self.seleccion_ruleta(self.poblacion)
                padre_a, padre_b = self.cruza_dos_puntos(padre[0], padre[1])
                padre_a = self.mutacion(padre_a)
                padre_b = self.mutacion(padre_b)
                self.poblacion.append(padre_a)
                self.poblacion.append(padre_b)

            self.poblacion = sorted(
                self.poblacion,
                key=lambda y: [y[3] for _ in self.poblacion],
                reverse=minimize
            )

            self.mejoresCasos.append(self.poblacion[0])
            self.promedioCasos.append(mean([x[3] for x in self.poblacion]))
            self.peoresCasos.append(self.poblacion[-1])

            if len(self.poblacion) > self.limitePoblacion:
                self.poda()
            generation += 1


root = tk.Tk()
w, h = 650, 400
root.geometry("%dx%d+0+0" % (w, h))
root.title('Algorítmo Genético - Maximizar y Minimizar')

root.configure(bg="white")
root.columnconfigure(0, weight=1)
root.columnconfigure(1, weight=7)

title_label = ttk.Label(root, text="Parámetros inciales del algoritmo", background='#fff',
                        font=('Lucida Sands', '12', 'bold')).place(x=250, y=10)


max_x_label = ttk.Label(root, text="Valor máximo para X:", background='#fff', font=(
    'Lucida Sands', '10')).place(x=50, y=50)
max_x_entry = ttk.Entry(root)
max_x_entry.insert(0, "20")
max_x_entry.place(x=50, y=100)
min_x_label = ttk.Label(root, text="Valor mínimo para X:", background='#fff', font=(
    'Lucida Sands', '10')).place(x=200, y=50)
min_x_entry = ttk.Entry(root)
min_x_entry.insert(0, "10")
min_x_entry.place(x=200, y=100)

resolution_x_label = ttk.Label(root, text="Precisión para X:", background='#fff', font=(
    'Lucida Sands', '10')).place(x=350, y=50)
resolution_x_entry = ttk.Entry(root)
resolution_x_entry.insert(0, "0.3")
resolution_x_entry.place(x=350, y=100)

max_generations_label = ttk.Label(root, text="Límite de generaciones:", background='#fff', font=(
    'Lucida Sands', '10')).place(x=500, y=50)
max_generations_entry = ttk.Entry(root)
max_generations_entry.insert(0, "15")
max_generations_entry.place(x=500, y=100)

initial_population_label = ttk.Label(root, text="Población inicial:", background='#fff', font=(
    'Lucida Sands', '10')).place(x=50, y=150)
initial_population_entry = ttk.Entry(root)
initial_population_entry.insert(0, "3")
initial_population_entry.place(x=50, y=200)

max_population_label = ttk.Label(root, text="Población máxima:", background='#fff', font=(
    'Lucida Sands', '10')).place(x=200, y=150)
max_population_entry = ttk.Entry(root)
max_population_entry.insert(0, "100")
max_population_entry.place(x=200, y=200)

individual_mutation_prob_label = ttk.Label(
    root, text="PMI:", background='#fff', font=('Lucida Sands', '10')).place(x=350, y=150)
individual_mutation_prob_entry = ttk.Entry(root)
individual_mutation_prob_entry.insert(0, "0.25")
individual_mutation_prob_entry.place(x=350, y=200)

gen_mutation_prob_label = ttk.Label(root, text="PMG:", background='#fff', font=(
    'Lucida Sands', '10')).place(x=500, y=150)
gen_mutation_prob_entry = ttk.Entry(root)
gen_mutation_prob_entry.insert(0, "0.40")
gen_mutation_prob_entry.place(x=500, y=200)


def run(minimize: bool):
    ga = AlgoritmoGenetico(float(resolution_x_entry.get()),
                           (float(min_x_entry.get()), float(max_x_entry.get())),
                           int(max_generations_entry.get()),
                           int(max_population_entry.get()),
                           int(initial_population_entry.get()),
                           float(individual_mutation_prob_entry.get()),
                           float(gen_mutation_prob_entry.get()))
    ga.iniciar(minimize)

    figure2 = plt.figure()

    plt.plot(np.arange(0, ga.limiteGeneraciones), [
             x[3] for x in ga.mejoresCasos], label="Mejores casos")
    plt.plot(np.arange(0, ga.limiteGeneraciones), [
             x[3] for x in ga.peoresCasos], label="Peores casos")
    plt.plot(np.arange(0, ga.limiteGeneraciones),
             ga.promedioCasos, label="Promedio de casos")
    plt.legend()
    plt.title("Evolución de la población")
    plt.xlabel("Generaciones/Iteraciones")
    plt.ylabel("Valor de aptitud")
    plt.show()

    messagebox.showinfo(
        message=f"Genotipo: {ga.poblacion[0][0]}\ni:{ga.poblacion[0][1]}, Fenotipo i: {ga.poblacion[0][2]},  Aptitud: {ga.poblacion[0][3]}",
        title="Mejor individuo")


login_button = ttk.Button(root, text="Maximizar",
                          command=lambda: run(True)).place(x=250, y=250)
login_button = ttk.Button(root, text="Minimizar",
                          command=lambda: run(False)).place(x=350, y=250)

root.mainloop()
