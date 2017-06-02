# -*- encoding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

test = np.zeros((int(w_n_points), int(h_n_points)))
test_w = test.shape[0]
test_h = test.shape[1]
size = [test_w, test_h]

T = test[:]

def is_border(i, j, side=""):
    if side == "":
        if(i == 0 or i == size[0]-1 or j == 0 or j == size[1]-1):
            return True
    if side == "T":
        if(i==0 and j !=0 and j != size[1]-1):
            return True
    if side == "B":
        if(i==size[0]-1 and j != 0 and j != size[1]-1):
            return True
    if side == "L":
        if(i != 0 and i != size[0]-1 and j == 0):
            return True
    if side == "R":
        if(i != 0 and i != size[0]-1 and j == size[1]-1):
            return True

def read_props_file(file_name):
    with open(file_name, "r") as f:
        props = {k: float(v) for (k, v) in [e[:-1].strip().split(" ") for e in f.readlines()]}

    for p in props:
        if "ISOLATED" in p:
            props[p] = False if props[p] == 0 else True

    return props

def solve(props):

    # Material Properties
    X_SIZE                 = props["X_SIZE"]
    Y_SIZE                 = props["Y_SIZE"]
    THERMAL_CONDUCTIVITY   = props["THERMAL_CONDUCTIVITY"]
    DENSITY                = props["DENSITY"]
    SPECIFIC_HEAT_CAPACITY = props["SPECIFIC_HEAT_CAPACITY"]

    # Simulation Properties
    TOTAL_TIME    = props["TOTAL_TIME"]
    TIME_STEP     = props["TIME_STEP"]
    TIME_N_POINTS = TOTAL_TIME / TIME_STEP
    XY_STEP       = props["XY_STEP"]
    X_N_POINTS    = int(X_SIZE / XY_STEP)
    Y_N_POINTS    = int(Y_SIZE / XY_STEP)

    # Boundary Conditions
    INIT_TEMPERATURE_INNER         = props["INIT_TEMPERATURE_INNER"]
    INIT_TEMPERATURE_TOP_BORDER    = props["INIT_TEMPERATURE_TOP_BORDER"]
    INIT_TEMPERATURE_BOTTOM_BORDER = props["INIT_TEMPERATURE_BOTTOM_BORDER"]
    INIT_TEMPERATURE_LEFT_BORDER   = props["INIT_TEMPERATURE_LEFT_BORDER"]
    INIT_TEMPERATURE_RIGHT_BORDER  = props["INIT_TEMPERATURE_RIGHT_BORDER"]
    ISOLATED_TOP_BORDER            = props["ISOLATED_TOP_BORDER"]
    ISOLATED_BOTTOM_BORDER         = props["ISOLATED_BOTTOM_BORDER"]
    ISOLATED_LEFT_BORDER           = props["ISOLATED_LEFT_BORDER"]
    ISOLATED_RIGHT_BORDER          = props["ISOLATED_RIGHT_BORDER"]

    # Calculate material thermal diffusivity
    MATERIAL_THERMAL_DIFFUSIVITY = THERMAL_CONDUCTIVITY / (DENSITY * SPECIFIC_HEAT_CAPACITY)

    # Calculate F0 constant
    F0 = MATERIAL_THERMAL_DIFFUSIVITY * (TIME_STEP / (XY_STEP ** 2))
    if F0 > 0.25:
        print("ERROR: F0 must be <= 0.25")
        quit()

    time = np.linspace(0.0, TOTAL_TIME, TIME_N_POINTS)

def DB():

    for i in range(test_w):
        for j in range(test_h):
            if i == 0:
                test[i][j] = INIT_TEMPERATURE_TOP_BORDER
            if j == 0 and (i != 0):
                test[i][j] = INIT_TEMPERATURE_LEFT_BORDER
            if j == (test_h-1) and (i != 0):
                test[i][j] = INIT_TEMPERATURE_RIGHT_BORDER 
            if i == test_w-1 and j != (test_h-1) and j != 0:
                test[i][j] = INIT_TEMPERATURE_BOTTOM_BORDER
    
    if T_ISOLATED:
        test[0][0] = INIT_TEMPERATURE_LEFT_BORDER
        test[0][test_w-1] = INIT_TEMPERATURE_RIGHT_BORDER
        test[test_h-1][0] = INIT_TEMPERATURE_BOTTOM_BORDER
        test[test_h-1][test_w-1] = INIT_TEMPERATURE_BOTTOM_BORDER

    for t in range(len(time)-1):
        if t != -1:
            prev_T = np.copy(T)
            for i in range(test_w):
                for j in range(test_h):
                    if not is_border(i,j):
                        T[i][j] = F0*(prev_T[i+1][j] + prev_T[i-1][j] + prev_T[i][j+1] + prev_T[i][j-1]) + prev_T[i][j]*(1 - 4*F0)
                    if ISOLATED_BOTTOM_BORDER and is_border(i,j,"B"):
                        T[i][j] = F0*(2*prev_T[test_w-2][j] + prev_T[test_w-1][j+1] + prev_T[test_w-1][j-1]) + (1-4*F0)*prev_T[test_w-1][j]
 
    return T

def DL():

    for i in range(test_w):
        for j in range(test_h):
            if i == 0:
                test[i][j] = INIT_TEMPERATURE_TOP_BORDER
            if j == 0 and (i != 0):
                test[i][j] = 0
            if j == (test_h-1) and (i != 0):
                test[i][j] = INIT_TEMPERATURE_RIGHT_BORDER 
            if i == test_w-1 and j != (test_h-1) and j != 0:
                test[i][j] = INIT_TEMPERATURE_BOTTOM_BORDER 

    for t in range(len(time)-1):
        if t != -1:
            prev_T = np.copy(T)
            for i in range(test_w):
                for j in range(test_h):
                    if not is_border(i,j):
                        T[i][j] = F0*(prev_T[i+1][j] + prev_T[i-1][j] + prev_T[i][j+1] + prev_T[i][j-1]) + prev_T[i][j]*(1 - 4*F0)
                    if ISOLATED_LEFT_BORDER and is_border(i,j,"L"):
                        T[i][j] = F0*(prev_T[i+1][0] + prev_T[i-1][0] + 2*prev_T[i][1]) + (1-4*F0)*prev_T[i][0]
 
    return T

def DT():

    superior = 0
    inferior = 100
    esquerda = 50
    direita = 75

    for i in range(test_w):
        for j in range(test_h):
            if i == 0:
                test[i][j] = superior
            if j == 0 and (i != 0):
                test[i][j] = esquerda
            if j == (test_h-1) and (i != 0):
                test[i][j] = direita 
            if i == test_w-1 and j != (test_h-1) and j != 0:
                test[i][j] = inferior 

    test[0][0]=esquerda
    test[0][test_w-1]=direita
    test[test_h-1][0]=inferior
    test[test_h-1][test_w-1]=inferior

    for t in range(len(time)-1):
        if t != -1:
            prev_T = np.copy(T)
            for i in range(test_w):
                for j in range(test_h):
                    if not is_border(i,j):
                        T[i][j] = F0*(prev_T[i+1][j] + prev_T[i-1][j] + prev_T[i][j+1] + prev_T[i][j-1]) + prev_T[i][j]*(1 - 4*F0)
                    if ISOLATED_TOP_BORDER and is_border(i,j,"T"):
                        T[i][j] = F0*(2*prev_T[1][j] + prev_T[0][j+1] + prev_T[0][j-1]) + (1-4*F0)*prev_T[0][j]
    
    return T

def DR():

    superior = 50
    inferior = 75
    esquerda = 100
    direita = 0

    for i in range(test_w):
        for j in range(test_h):
            if i == 0:
                test[i][j] = superior
            if j == 0 and (i != 0):
                test[i][j] = esquerda
            if j == (test_h-1) and (i != 0):
                test[i][j] = direita 
            if i == test_w-1 and j != (test_h-1) and j != 0:
                test[i][j] = inferior 

    test[0][0]=esquerda
    test[0][test_w-1]=superior
    test[test_h-1][0]=esquerda
    test[test_h-1][test_w-1]=inferior

    for t in range(len(time)-1):
        if t != -1:
            prev_T = np.copy(T)
            for i in range(test_w):
                for j in range(test_h):
                    if not is_border(i,j):
                        T[i][j] = F0*(prev_T[i+1][j] + prev_T[i-1][j] + prev_T[i][j+1] + prev_T[i][j-1]) + prev_T[i][j]*(1 - 4*F0)
                    if ISOLATED_RIGHT_BORDER and is_border(i,j,"R"):
                        T[i][j] = F0*(prev_T[i+1][test_h-1] + prev_T[i-1][test_h-1] + 2*prev_T[i][test_h-2]) + (1-4*F0)*prev_T[i][test_h-1]
    return T


if __name__ == '__main__':
    result = DR()
    #plt.imshow(result, cmap='bwr', interpolation='nearest')
    #plt.show()
    np.savetxt('test.out', result, fmt='%1.4e') 