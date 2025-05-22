import modulos.algLE as le
import sympy as sp

Robot2G = le.AlgoritmoLE(2)

Js = Robot2G.definir_simbolos()

tabla_DH = [[Js[0], 0, 0, sp.rad(-90)],
            [0, Js[1], 0,  0         ]]

m_Ai = Robot2G.L_E2(tabla_DH)
m_Uij = Robot2G.L_E3(m_Ai)
m_Uijk = Robot2G.L_E4(m_Uij)

L1 = sp.symbols('L1')
Tabla_L_E5 = [[0, 0, L1],
              [0, 0, 0 ]]

m_Jpi = Robot2G.L_E5(Tabla_L_E5)
m_d = Robot2G.L_E6(m_Uij,m_Jpi)

# ---------------- Sustitución de valores para la matriz D del paso 6 (inicio) ----------------
J2, L1, m1, m2 = sp.symbols('J2 L1 m1 m2')

sustituciones = {
    J2: 1.3,
    L1: 1.1,
    m1: 1,
    m2: 1.2
}

d_numerica = m_d.subs(sustituciones)

# Mostrar matriz resultante
print("\nMatriz simbólica D:")
sp.pprint(m_d)

print("\nMatriz D con valores sustituidos:")
sp.pprint(d_numerica)

# ---------------- Sustitución de valores para la matriz D del paso 6 (fin) -------------------

m_h = Robot2G.L_E7(m_Uij,m_Uijk,m_Jpi)