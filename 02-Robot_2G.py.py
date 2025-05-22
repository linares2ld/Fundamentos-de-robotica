import modulos.algLE as le
import sympy as sp

Robot2G = le.AlgoritmoLE(2)

# L-E 1.
ms,ls,js,jds,jdds = Robot2G.definir_simbolos()
tabla_DH = [[js[0], 0, 0, sp.rad(-90)],
            [0, js[1], 0,  0         ]]

# -------------------- Asignación de valores (inicio) --------------------
sustituciones = {ms[0]: 1, ms[1]: 2,
                 ls[0]: 1, ls[1]: 0,
                 js[0]: 0, js[1]: 1,
                 jds[0]: 2, jds[1]: 2, 
                 jdds[0]: 1, jdds[1]: 1}
# -------------------- Asignación de valores (fin) -----------------------

# L-E 2.
m_Ai = Robot2G.L_E2(tabla_DH)

# L-E 3.
m_Uij = Robot2G.L_E3(m_Ai)

# L-E 4.
m_Uijk = Robot2G.L_E4(m_Uij)

# L-E 5.
v_r = [[0, 0, ls[0]],
       [0, 0, 0 ]]
m_Jpi = Robot2G.L_E5(v_r)

# L-E 6.
D = Robot2G.L_E6(m_Uij,m_Jpi)

# ---------------- Impresión de la matriz D (inicio) ----------------
D_numerica = D.subs(sustituciones)

# Mostrar matriz resultante
print("\nMatriz simbólica D:\n")
sp.pprint(D)

print("\nMatriz D con valores sustituidos:\n")
sp.pprint(D_numerica)

# ---------------- Impresión de la matriz D (fin) -------------------

# L-E 7.
m_h = Robot2G.L_E7(m_Uij,m_Uijk,m_Jpi)

# L-E 8.
H = Robot2G.L_E8(m_h)

# ---------------- Impresión de la matriz H (inicio) ----------------
H_numerica = H.subs(sustituciones)

# Mostrar matriz resultante
print("\nMatriz simbólica H:\n")
sp.pprint(H)

print("\nMatriz H con valores sustituidos:\n")
sp.pprint(H_numerica)

# ---------------- Impresión de la matriz H (fin) -------------------
