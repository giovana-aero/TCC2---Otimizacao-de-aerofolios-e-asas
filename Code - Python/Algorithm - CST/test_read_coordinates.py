import matplotlib.pyplot as plt
# import pandas as pd
import numpy as np

# https://stackoverflow.com/questions/8282553/removing-character-in-list-of-strings
# f = open('coordenadas.dat','r')
# lines = f.readlines()
# result1 = []
# result2 = []
# for x in lines:
#     result1.append(x.split(' ')[0])
#     result2.append(x.split(' ')[1])
# f.close()

# Ficamos com duas listas de strings.
# No caso da segunda lista, é necessário tirar os caracteres de quebra de linha:
# result2 = 



# https://www.statology.org/pandas-read-text-file/
# coo = pd.read_csv('coordenadas.dat',sep=' ',header=None)




# https://stackoverflow.com/questions/3518778/how-do-i-read-csv-data-into-a-record-array-in-numpy
coo = np.genfromtxt('coordenadas.dat',delimiter=' ')