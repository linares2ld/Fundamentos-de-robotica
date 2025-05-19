import Lib.robcinematics.algDH as rc
import pandas as pd

tabla_resultados = []
valores_de_las_J = [[33.67, 68.603, -70.13, 0, 43.83, -28.992],
                    [-28.849, 70.002, -75.066, 0, 53.427, -21.854],
                    [25.048, -12.395, 38.412, -25.919, -33.122, 41.675],
                    [-31.187, 19.02, -28.766, -19.17, 26.37, -33.561]] 

fanuc = rc.AlgoritmoDH()
Js = fanuc.definir_simbolos(6)

tabla_DH = [[Js[0], 0, 75, -90],
            [Js[1]-90, 0, 300, 180],
            [Js[1]+Js[2], 0, 75, -90],
            [Js[3], -320, 0, 90],
            [Js[4], 0, 0, -90],
            [Js[5], -80, 0, 0]]

fanuc.matrices_homogeneas(tabla_DH)
fanuc.matriz_T()

for i in range(len(valores_de_las_J)):
    matriz_resultados = fanuc.sustitucion_matriz_T(valores_de_las_J[i], True)
    resultados_xyz = fanuc.resultados_xyz(matriz_resultados)

    tabla_resultados.append(resultados_xyz)


df = pd.DataFrame(tabla_resultados, columns=["X", "Y", "Z"], index=[f"Posición {i+1}" for i in range(len(valores_de_las_J))])
df = df.map(lambda x: round(x, 3))
print(df)