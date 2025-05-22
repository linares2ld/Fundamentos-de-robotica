import sympy as sp

# Algoritmo computacional para el modelado dinámico por Lagrange-Euler

class AlgoritmoLE():
    def __init__(self, grados_libertad):

        self.gl = grados_libertad

        self.ms = []
        self.ls = []
        self.js = []
        self.jds = []
        self.jdds = []

        self.m_homogeneas = []        
      
    def definir_simbolos(self): 
    # Define todos los simbolos necesarios para el algoritmo

        derv = '\u1d48'

        g = sp.symbols(f'g')
        self.ms = sp.symbols(f'm1:{self.gl+1}') 
        self.ls = sp.symbols(f'l1:{self.gl+1}') 
        self.js = sp.symbols(f'j1:{self.gl+1}') 
        self.jds = sp.symbols(f'j{derv}1:{self.gl+1}')
        self.jdds = sp.symbols(f'j{derv}{derv}1:{self.gl+1}')

        return g, self.ms, self.ls, self.js, self.jds, self.jdds

    def matrices_homogeneas(self, tabla_DH):

        for i in range(self.gl):
            theta, d, a, alpha = tabla_DH[i]

            m_aux = sp.Matrix([[sp.cos(theta), -sp.sin(theta)*sp.cos(alpha),  sp.sin(theta)*sp.sin(alpha), a*sp.cos(theta)],
                               [sp.sin(theta),  sp.cos(theta)*sp.cos(alpha), -sp.cos(theta)*sp.sin(alpha), a*sp.sin(theta)],
                               [0,              sp.sin(alpha),                sp.cos(alpha),               d              ],
                               [0,              0,                            0,                           1              ]])

            self.m_homogeneas.append(m_aux)        

    def L_E2(self, tabla_DH):
    # L-E 2: obtiene de las matrices de transformación 0Ai para cada elemento i

        self.matrices_homogeneas(tabla_DH)
        m_Ai = []   
        aux = self.m_homogeneas[0]

        for i in range(self.gl):
            for j in range(1, i+1):
                aux *=  self.m_homogeneas[j]

            m_Ai.append(aux)
        
        return m_Ai

    def L_E3(self, m_Ai):
    # L-E 3: obtiene las matrices Uij 

        m_Uij = []

        for A_i in m_Ai:
            for q_j in self.js:
                derivada = A_i.diff(q_j)
                m_Uij.append(derivada)

        return m_Uij

    def L_E4(self, m_Uij):
    # L-E 4: obtiene las matrices Uijk

        m_Uijk = []

        for Uij in m_Uij:
            for q_k in self.js:
                derivada = Uij.diff(q_k)
                m_Uijk.append(derivada)

        return m_Uijk
    
    def L_E5(self, tabla):
    # L-E 5: obtiene las matrices de pseudinercia Jpi para cada elemento

        m_Jpi = []

        for i in range(self.gl):
            x, y, z, trash = tabla[i]

            m_aux = sp.Matrix([[x**2, x*y, x*z, x],
                               [x*y, y**2, y*z, y],
                               [x*z, y*z, z**2, z],
                               [x,   y,   z,   1]])
            
            m_Jpi_i = m_aux.applyfunc(lambda expr: sp.integrate(expr, self.ms[i]))

            m_Jpi.append(m_Jpi_i)
        
        return m_Jpi

    def L_E6(self, m_Uij, m_Jpi):
    # L-E 6: obtiene la matriz de inercias D = [Dij]

        D = sp.zeros(self.gl, self.gl)

        for i in range(self.gl):
            for j in range(self.gl):
                suma = 0
                for k in range(self.gl):
                    idx_ki = k * self.gl + i  # indexa como si m_Uij fuera una matriz aplanada
                    idx_kj = k * self.gl + j
                    Ukj = m_Uij[idx_kj]
                    Uki = m_Uij[idx_ki]
                    Jpk = m_Jpi[k]

                    traza = (Ukj * Jpk * Uki.T).trace()
                    suma += traza
                D[i, j] = sp.simplify(suma)

        return D

    def L_E7(self, m_Uij, m_Uijk, m_Jpi):
    # L-E 7: obtiene los términos hikm definidos por una sumatoria desde j=max(i, k, m)

        hikm = []

        for i in range(self.gl):
            for k in range(self.gl):
                for m in range(self.gl):
                    suma = 0
                    max_j = max(i, k, m)
                    for j in range(max_j, self.gl):
                        # Índices para acceder a U_{jkm} y U_{ji}
                        idx_jkm = j * self.gl**2 + k * self.gl + m
                        idx_ji = j * self.gl + i
                        
                        Ujkm = m_Uijk[idx_jkm]
                        Uji = m_Uij[idx_ji]
                        Jj = m_Jpi[j]

                        traza = (Ujkm * Jj * Uji.T).trace()
                        suma += traza

                    # Guarda h_{ikm} en un diccionario con tupla como clave
                    suma = sp.simplify(suma)
                    hikm.append(suma)

        return hikm

    def L_E8(self, hikm):
    # L-E 8: obtiene la matriz columna de fuerzas de Coriolis y centrípeta H = [hi]T

        H = []

        for i in range(self.gl):
            suma = 0
            for k in range(self.gl):
                for m in range(self.gl):
                    idx = i * self.gl**2 + k * self.gl + m
                    hikm_val = hikm[idx]
                    suma += hikm_val * self.jds[k] * self.jds[m]
            H.append(sp.simplify(suma))

        return sp.Matrix(H)
    
    def L_E9(self, v_g, m_Uij, v_r):

        C = []

        for i in range(self.gl):
            suma = 0
            for j in range(self.gl):

                idx = j * self.gl + i  # índice de Uij
                Uij = m_Uij[idx]

                rj = sp.Matrix(v_r[j])

                val =  self.ms[j] * v_g * Uij * rj
                suma += val[0]

            C.append(sp.simplify(suma))

        return sp.Matrix(C)