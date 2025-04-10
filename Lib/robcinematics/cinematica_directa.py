import sympy as sp

class AlgoritmoDH():
    def __init__(self):

        self.grados_libertad = 0
        self.tabladh = []
        self.js = []
        self.matriceshomogeneas = []
        self.matrizT = sp.eye(4)
      
    def definir_simbolos(self, grados_libertad):

        self.grados_libertad = grados_libertad

        # Define la cantidad de js
        self.js = sp.symbols(f'J1:{self.grados_libertad+1}')

        return self.js

    def matrices_homogeneas(self, tabla, guardar=False):

        self.tabladh = tabla

        for i in range(self.grados_libertad):

            theta, d, a, alpha = self.tabladh[i]

            theta = theta * sp.pi / 180
            alpha = alpha * sp.pi / 180

            matriz_auxiliar = sp.Matrix([[sp.cos(theta), -sp.sin(theta)*sp.cos(alpha),  sp.sin(theta)*sp.sin(alpha), a*sp.cos(theta)],
                                         [sp.sin(theta),  sp.cos(theta)*sp.cos(alpha), -sp.cos(theta)*sp.sin(alpha), a*sp.sin(theta)],
                                         [0,              sp.sin(alpha),                sp.cos(alpha),               d],
                                         [0,              0,                            0,                           1]
                                        ])

            self.matriceshomogeneas.append(matriz_auxiliar)

        if guardar:
            return self.matriceshomogeneas 
    def matriz_T(self, guardar=False):

        for i in range(self.grados_libertad):
            self.matrizT = self.matrizT * self.matriceshomogeneas[i]

        if guardar:
            return self.matrizT

    def sustitucion_matriz_T(self,valores_j, guardar=False):

        dict_valores_j = {self.js[i]: valores_j[i] for i in range(self.grados_libertad)}
        resultados_matriz_T = self.matrizT.subs(dict_valores_j).evalf()

        if guardar:
            return resultados_matriz_T
    
    def ver_matriz(self,matriz):
        matriz = sp.Matrix(matriz)
        matriz = sp.simplify(matriz)
        sp.pprint(matriz)

    def resultados_xyz(self, matriz):
        x = matriz[0, 3]
        y = matriz[1, 3]
        z = matriz[2, 3]
        xyz = (x, y, z)

        return xyz