# -*- encoding: utf-8 -*-

import numpy as np

def is_superior( i, j, size):
    if(i==0 and j !=0 and j != size[1]-1):
        return True

def is_inferior(i, j, size):
    if(i==size[0]-1 and j != 0 and j != size[1]-1):
        return True

def is_left(i, j, size):
    if(i != 0 and i != size[0]-1 and j == 0):
        return True

def is_right(i, j, size):
    if(i != 0 and i != size[0]-1 and j == size[1]-1):
        return True

def is_border(i, j, size):
    if(i == 0 or i == size[0]-1 or j == 0 or j == size[1]-1):
        return True



def D1():
    alfa = 1                                # difusibilidade da barra
    x = 50                                  # tamanho da barra
    nx = 10                                 # numero de discretizacoes em x
    dx = x/nx                               # passo no espaco
    tfinal = 500
    nt = 100                                # numero de discretizacoes em t
    dt = tfinal/nt                          # passo no tempo
    lmb = alfa*(dt/(dx**2))                 # lambda
    lmb = 0.2
    X = [i for i in range(0, x + dx, dx)]   #lista de posicoes
    T = [0]*len(X)                          #lista de temperaturas

    for t in range(0 , tfinal, dt):
        prev_T = T[:]
        for i in range(0, len(T)):
            if t == 0:                        # Condicao inicial
                if i ==0:
                    T[i]=0
                elif i == nx:
                    T[i]=0
                else:
                    T[i]=20
            else:
                if((i != 0) and (i != nx)):
                    #calcula a temperatura interna da barra
                    T[i] = prev_T[i] + lmb*(prev_T[i+1] - 2*prev_T[i] + prev_T[i-1])

        print(T)

def printm(m):
    for line in m:
        print(line)

def DI():
    #file = open("test.in", "r") 
    #print(file.read())

    print("Computing 2D temperature matrix")

############################################
    #top_edge = int(file.readline(5))
    #left_edge = int(file.readline(5))
    #right_edge = int(file.readline(5))
    #bottom_edge = int(file.readline(5))
    #internal_temperature = int(file.readline(5))

    #tfinal = int(file.readline(5))
    #t_n_points = int(file.readline(5))

##################################################
    top_edge = 100
    left_edge = 75
    right_edge = 50
    bottom_edge = 0
    internal_temperature = 0

    tfinal = 200.0
    t_n_points = 160000.0
    dt = tfinal/t_n_points
    #print(dt)

    time = np.linspace(0.0, tfinal, 100)
#####################################################
    #width = int(file.readline(5))
    #w_n_points = int(file.readline(5))
    #dwidth = width/w_n_points

    #height = int(file.readline(5))
    #h_n_points = int(file.readline(5))
    #dheight = height/h_n_points
#######################################################
    width = 0.5
    w_n_points = 5
    dwidth = width/w_n_points

    height = 0.5
    h_n_points = 5
    dheight = height/h_n_points

    step = np.linspace(0.0, width, 5)

    test = np.zeros((int(w_n_points), int(h_n_points)))
    test_w = test.shape[0]
    test_h = test.shape[1]
    

    #alpha = int(file.readline(5))
    alpha = 1
    F0 = alpha*(float(dt)/float((dwidth**2)))

    print(F0)

    T = test[:]

    superior = 100
    inferior = 0
    esquerda = 75
    direita = 50

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



            


    for t in range(len(time)-1):
        if t != -1:
            prev_T = np.copy(T)
            for i in range(test_w):
                for j in range(test_h):
                    if not is_border(i,j,[test_w, test_h]):
                        T[i][j] = F0*(prev_T[i+1][j] + prev_T[i-1][j] + prev_T[i][j+1] + prev_T[i][j-1]) + prev_T[i][j]*(1 - 4*F0)
                    if is_inferior(i, j, [test_w, test_h]):
                        T[i][j] = F0*(2*prev_T[test_w-2][j] + prev_T[test_w-1][j+1] + prev_T[test_w-1][j-1]) + (1-4*F0)*prev_T[test_w-1][j]
                    #if is_superior(i, j, [test_w, test_h]):
                    #    T[i][j] = F0*(2*prev_T[1][j] + prev_T[0][j+1] + prev_T[0][j-1]) + (1-4*F0)*prev_T[0][j]
                    #if is_left(i,j,[test_w, test_h]):
                    #    T[i][j] = F0*(prev_T[i+1][0] + prev_T[i-1][0] + 2*prev_T[j][1]) + (1-4*F0)*prev_T[i][0]


                    
    #file.close() 
    return T

def DL():
    #file = open("test.in", "r") 
    #print(file.read())

    print("Computing 2D temperature matrix")

############################################
    #top_edge = int(file.readline(5))
    #left_edge = int(file.readline(5))
    #right_edge = int(file.readline(5))
    #bottom_edge = int(file.readline(5))
    #internal_temperature = int(file.readline(5))

    #tfinal = int(file.readline(5))
    #t_n_points = int(file.readline(5))

##################################################
    top_edge = 100
    left_edge = 75
    right_edge = 50
    bottom_edge = 0
    internal_temperature = 0

    tfinal = 200.0
    t_n_points = 160000.0
    dt = tfinal/t_n_points
    #print(dt)

    time = np.linspace(0.0, tfinal, 100)
#####################################################
    #width = int(file.readline(5))
    #w_n_points = int(file.readline(5))
    #dwidth = width/w_n_points

    #height = int(file.readline(5))
    #h_n_points = int(file.readline(5))
    #dheight = height/h_n_points
#######################################################
    width = 0.5
    w_n_points = 5
    dwidth = width/w_n_points

    height = 0.5
    h_n_points = 5
    dheight = height/h_n_points

    step = np.linspace(0.0, width, 5)

    test = np.zeros((int(w_n_points), int(h_n_points)))
    test_w = test.shape[0]
    test_h = test.shape[1]
    

    #alpha = int(file.readline(5))
    alpha = 1
    F0 = alpha*(float(dt)/float((dwidth**2)))

    print(F0)

    T = test[:]

    superior = 75
    inferior = 50
    esquerda = 0
    direita = 100

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

    test[0][0]=superior
    test[0][test_w-1]=direita
    test[test_h-1][0]=inferior
    test[test_h-1][test_w-1]=direita

    for t in range(len(time)-1):
        if t != -1:
            prev_T = np.copy(T)
            for i in range(test_w):
                for j in range(test_h):
                    if not is_border(i,j,[test_w, test_h]):
                        T[i][j] = F0*(prev_T[i+1][j] + prev_T[i-1][j] + prev_T[i][j+1] + prev_T[i][j-1]) + prev_T[i][j]*(1 - 4*F0)
                    #if is_inferior(i, j, [test_w, test_h]):
                    #    T[i][j] = F0*(2*prev_T[test_w-2][j] + prev_T[test_w-1][j+1] + prev_T[test_w-1][j-1]) + (1-4*F0)*prev_T[test_w-1][j]
                    #if is_superior(i, j, [test_w, test_h]):
                    #    T[i][j] = F0*(2*prev_T[1][j] + prev_T[0][j+1] + prev_T[0][j-1]) + (1-4*F0)*prev_T[0][j]
                    if is_left(i,j,[test_w, test_h]):
                        T[i][j] = F0*(prev_T[i+1][0] + prev_T[i-1][0] + 2*prev_T[i][1]) + (1-4*F0)*prev_T[i][0]


                    
    #file.close() 
    return T

def DS():
    #file = open("test.in", "r") 
    #print(file.read())

    print("Computing 2D temperature matrix")

############################################
    #top_edge = int(file.readline(5))
    #left_edge = int(file.readline(5))
    #right_edge = int(file.readline(5))
    #bottom_edge = int(file.readline(5))
    #internal_temperature = int(file.readline(5))

    #tfinal = int(file.readline(5))
    #t_n_points = int(file.readline(5))

##################################################
    top_edge = 100
    left_edge = 75
    right_edge = 50
    bottom_edge = 0
    internal_temperature = 0

    tfinal = 200.0
    t_n_points = 160000.0
    dt = tfinal/t_n_points
    #print(dt)

    time = np.linspace(0.0, tfinal, 100)
#####################################################
    #width = int(file.readline(5))
    #w_n_points = int(file.readline(5))
    #dwidth = width/w_n_points

    #height = int(file.readline(5))
    #h_n_points = int(file.readline(5))
    #dheight = height/h_n_points
#######################################################
    width = 0.5
    w_n_points = 5
    dwidth = width/w_n_points

    height = 0.5
    h_n_points = 5
    dheight = height/h_n_points

    step = np.linspace(0.0, width, 5)

    test = np.zeros((int(w_n_points), int(h_n_points)))
    test_w = test.shape[0]
    test_h = test.shape[1]
    

    #alpha = int(file.readline(5))
    alpha = 1
    F0 = alpha*(float(dt)/float((dwidth**2)))

    print(F0)

    T = test[:]

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
                    if not is_border(i,j,[test_w, test_h]):
                        T[i][j] = F0*(prev_T[i+1][j] + prev_T[i-1][j] + prev_T[i][j+1] + prev_T[i][j-1]) + prev_T[i][j]*(1 - 4*F0)
                    #if is_inferior(i, j, [test_w, test_h]):
                    #    T[i][j] = F0*(2*prev_T[test_w-2][j] + prev_T[test_w-1][j+1] + prev_T[test_w-1][j-1]) + (1-4*F0)*prev_T[test_w-1][j]
                    if is_superior(i, j, [test_w, test_h]):
                        T[i][j] = F0*(2*prev_T[1][j] + prev_T[0][j+1] + prev_T[0][j-1]) + (1-4*F0)*prev_T[0][j]
                    #if is_left(i,j,[test_w, test_h]):
                    #    T[i][j] = F0*(prev_T[i+1][0] + prev_T[i-1][0] + 2*prev_T[i][1]) + (1-4*F0)*prev_T[i][0]
    
    #file.close() 
    return T

def DD():

    print("Computing 2D temperature matrix")

    tfinal = 200.0
    t_n_points = 160000.0
    dt = tfinal/t_n_points

    time = np.linspace(0.0, tfinal, 100)

    width = 0.5
    w_n_points = 5
    dwidth = width/w_n_points

    height = 0.5
    h_n_points = 5
    dheight = height/h_n_points

    step = np.linspace(0.0, width, 5)

    test = np.zeros((int(w_n_points), int(h_n_points)))
    test_w = test.shape[0]
    test_h = test.shape[1]
    
    alpha = 1
    F0 = alpha*(float(dt)/float((dwidth**2)))

    print(F0)

    T = test[:]

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
                    if not is_border(i,j,[test_w, test_h]):
                        T[i][j] = F0*(prev_T[i+1][j] + prev_T[i-1][j] + prev_T[i][j+1] + prev_T[i][j-1]) + prev_T[i][j]*(1 - 4*F0)
                    if is_right(i, j, [test_w, test_h]):
                        T[i][j] = F0*(prev_T[i+1][test_h-1] + prev_T[i-1][test_h-1] + 2*prev_T[i][test_h-2]) + (1-4*F0)*prev_T[i][test_h-1]
    return T


if __name__ == '__main__':
    # D1()
    # printm(D2())
    np.savetxt('test.out', DS(), fmt='%1.4e') 
    