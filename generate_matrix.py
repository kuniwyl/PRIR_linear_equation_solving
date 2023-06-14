import numpy as np
import random

N = 1000

def gen():
    return [random.randint(1, 100) for _ in range(N)]

A = []
for _ in range(N):
    row = gen()
    A.append(row)

for i in range(N):
    A[i][i] = int(sum(A[i]) / len(A[i]) * 4)

A = np.array(A)
B = np.array(gen())

X = np.linalg.solve(A, B)

with open('data/matrix.txt', 'w') as file:
    file.write(str(B.size) + "\n")
    file.write(str(B.size) + "\n")
    for row in A:
        for elem in row:
            file.write(str(elem) + " ")
        file.write("\n")

with open('data/vector.txt', 'w') as file:
    file.write(str(B.size) + '\n')
    for elem in B:
        file.write(str(elem) + "\n")

with open('data/result.txt', 'w') as file:
    for elem in X:
        file.write(str(elem) + "\n")

