import random
from random import choices, choice
import math
from matplotlib import pyplot as plt
import numpy as np
from sympy import lambdify, simplify, var, sin, cos, symbols
from statistics import mean

class GeneticAlgorithm:

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
        expresion = (0.75 * sin(0.50 * x) * sin(0.25 * x) * -0.75 * sin(0.75 * x))
        expresion2 = simplify(expresion)
        self.f = lambdify((x), expresion2)
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

    def mutacion(self, individual):

        p = random.random()
        if p < self.probabilidadMutIndiv:
            for _ in range(self.nBx):
                index = random.randrange(self.nBx)
                individual[0][index] = individual[0][index] if random.random() > self.probabilidadMutGen else \
                    abs(individual[0][index] - 1)
            individual = self.generarIndividuo(individual[0])
            return individual
        else:
            return individual

    def generarIndividuo(self, genotipo):
        i = int("".join([str(i) for i in genotipo]), 2)
        fenotipo = self.rango[0] + (i * self.precision)

        if fenotipo > self.rango[-1]:
            fenotipo = self.rango[-1]

        aptitud = self.f(fenotipo)
        
        return [genotipo, i,fenotipo, aptitud]

    def poda(self):
        self.poblacion = self.poblacion[:self.limitePoblacion]

    def cruza(self, a, b):
        limite = random.randint(1, self.nBx)
        genotipoa = a[0][0:limite] + b[0][limite:]
        genotipob= b[0][0:limite] + a[0][limite:]
        padre_a = self.generarIndividuo(genotipoa)
        padre_b = self.generarIndividuo(genotipob)
        return padre_a, padre_b

    @staticmethod
    def seleccionarPadre(poblacion):
        parents = []
        for _ in range(2):
            parents.append(choice(poblacion))
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
                padre = self.seleccionarPadre(self.poblacion)
                padre_a, padre_b = self.cruza(padre[0], padre[1])
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
            