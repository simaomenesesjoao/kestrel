import numpy as np
import sys 
import os
sys.stdout.flush()


N = 100000

lista_data= []
count = 0
for filename in os.listdir(os.getcwd()):
    count += 1

    f = open(filename, 'r')
    A = f.readlines()
    # print(A)

    B = np.zeros(N)


    for i,a in enumerate(A):
        num = float(a.split(" ")[1])
        B[i] = num

    lista_data.append(B)


    print(count, filename, B[0])
    # if(count > 10):
        # break


# print(lista_data)
print("finished reading")


avg = np.zeros(N)
avg1    = np.zeros(N)
avg10   = np.zeros(N)
avg100  = np.zeros(N)
avg500  = np.zeros(N)
avg_all = np.zeros(N)
for i in range(count):
    print(i, lista_data[i][0])
    avg += (lista_data[i] - avg)/(i+1)
    if(i == 0):
        avg1 = avg.copy()
        np.savetxt("avg1.txt", avg1, delimiter="\n")
    if(i == 10):
        avg10 = avg.copy()
        np.savetxt("avg10.txt", avg10, delimiter="\n")
    if(i == 100):
        avg100 = avg.copy()
        np.savetxt("avg100.txt", avg100, delimiter="\n")
    if(i == 500):
        avg500 = avg.copy()
        np.savetxt("avg500.txt", avg500, delimiter="\n")
    if(i == count - 1):
        avg_all = avg.copy()
        np.savetxt("avg_all.txt", avg_all, delimiter="\n")

print("finished averaging")

