# -*- encoding: utf-8 -*-

import sys
import matplotlib.pyplot as plt
import numpy as np

def is_border(i, j, side=""):
    """Check if position (i,j) is in matrix border"""
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
    """Read input file"""
    with open(file_name, "r") as f:
        props = {k: float(v) for (k, v) in [e[:-1].strip().split(" ") for e in f.readlines()]}

    for p in props:
        if "ISOLATED" in p:
            props[p] = False if props[p] == 0 else True

    return props

def solve(props):
    """Solve 2D matrix for temperature"""
    global size

    # Material Properties
    X_SIZE                 = props["X_SIZE"] # Material width
    Y_SIZE                 = props["Y_SIZE"] # Material height
    THERMAL_CONDUCTIVITY   = props["THERMAL_CONDUCTIVITY"]
    DENSITY                = props["DENSITY"]
    SPECIFIC_HEAT_CAPACITY = props["SPECIFIC_HEAT_CAPACITY"]

    # Simulation Properties
    TOTAL_TIME    = props["TOTAL_TIME"]    # Total simulation time
    TIME_STEP     = props["TIME_STEP"]     # Interval between two time points
    TIME_N_POINTS = TOTAL_TIME / TIME_STEP # Number of time points
    XY_STEP       = props["XY_STEP"]       # Interval between two size points
    X_N_POINTS    = int(X_SIZE / XY_STEP)  # Number of points in X
    Y_N_POINTS    = int(Y_SIZE / XY_STEP)  # Number of points in Y

    # Boundary Conditions
    INIT_TEMPERATURE_INNER         = props["INIT_TEMPERATURE_INNER"] # Initial temperature (C) for non-border points
    INIT_TEMPERATURE_TOP_BORDER    = props["INIT_TEMPERATURE_TOP_BORDER"]
    INIT_TEMPERATURE_BOTTOM_BORDER = props["INIT_TEMPERATURE_BOTTOM_BORDER"]
    INIT_TEMPERATURE_LEFT_BORDER   = props["INIT_TEMPERATURE_LEFT_BORDER"]
    INIT_TEMPERATURE_RIGHT_BORDER  = props["INIT_TEMPERATURE_RIGHT_BORDER"]
    ISOLATED_TOP_BORDER            = props["ISOLATED_TOP_BORDER"]    # 1 if border is isolated, 0 otherwise
    ISOLATED_BOTTOM_BORDER         = props["ISOLATED_BOTTOM_BORDER"] # 1 if border is isolated, 0 otherwise
    ISOLATED_LEFT_BORDER           = props["ISOLATED_LEFT_BORDER"]   # 1 if border is isolated, 0 otherwise
    ISOLATED_RIGHT_BORDER          = props["ISOLATED_RIGHT_BORDER"]  # 1 if border is isolated, 0 otherwise

    # Calculate material thermal diffusivity
    MATERIAL_THERMAL_DIFFUSIVITY = THERMAL_CONDUCTIVITY / (DENSITY * SPECIFIC_HEAT_CAPACITY)

    # Calculate F0 constant
    F0 = MATERIAL_THERMAL_DIFFUSIVITY * (TIME_STEP / (XY_STEP ** 2))
    if F0 > 0.25:
        print("ERROR: F0 must be <= 0.25")
        quit()

    # Configure simulation steps
    time = np.linspace(0.0, TOTAL_TIME, TIME_N_POINTS)

    # Init old and current matrix
    oldT = np.zeros((X_N_POINTS, Y_N_POINTS))
    oldT_w = oldT.shape[0]
    oldT_h = oldT.shape[1]
    size = [oldT_w, oldT_h]
    T = oldT[:]

    # Set initial temperatures
    for i in range(oldT_w):
        for j in range(oldT_h):
            if i == 0:
                oldT[i][j] = INIT_TEMPERATURE_TOP_BORDER
            if j == 0 and (i != 0):
                oldT[i][j] = INIT_TEMPERATURE_LEFT_BORDER
            if j == (oldT_h-1) and (i != 0):
                oldT[i][j] = INIT_TEMPERATURE_RIGHT_BORDER
            if i == oldT_w-1 and j != (oldT_h-1) and j != 0:
                oldT[i][j] = INIT_TEMPERATURE_BOTTOM_BORDER

    # Set vertices specifically depending on isolated border
    if ISOLATED_TOP_BORDER:
        oldT[0][0] = INIT_TEMPERATURE_LEFT_BORDER
        oldT[0][oldT_w-1] = INIT_TEMPERATURE_RIGHT_BORDER
        oldT[oldT_h-1][0] = INIT_TEMPERATURE_BOTTOM_BORDER
        oldT[oldT_h-1][oldT_w-1] = INIT_TEMPERATURE_BOTTOM_BORDER

    if ISOLATED_BOTTOM_BORDER:
        oldT[0][0] = INIT_TEMPERATURE_TOP_BORDER
        oldT[0][oldT_w-1] = INIT_TEMPERATURE_TOP_BORDER
        oldT[oldT_h-1][0] = INIT_TEMPERATURE_LEFT_BORDER
        oldT[oldT_h-1][oldT_w-1] = INIT_TEMPERATURE_RIGHT_BORDER

    if ISOLATED_LEFT_BORDER:
        oldT[0][0] = INIT_TEMPERATURE_TOP_BORDER
        oldT[0][oldT_w-1] = INIT_TEMPERATURE_RIGHT_BORDER
        oldT[oldT_h-1][0] = INIT_TEMPERATURE_BOTTOM_BORDER
        oldT[oldT_h-1][oldT_w-1] = INIT_TEMPERATURE_RIGHT_BORDER

    if ISOLATED_RIGHT_BORDER:
        oldT[0][0] = INIT_TEMPERATURE_LEFT_BORDER
        oldT[0][oldT_w-1] = INIT_TEMPERATURE_TOP_BORDER
        oldT[oldT_h-1][0] = INIT_TEMPERATURE_LEFT_BORDER
        oldT[oldT_h-1][oldT_w-1] = INIT_TEMPERATURE_BOTTOM_BORDER

    # Main FDM loop
    for t in range(len(time)-1):
        prev_T = np.copy(T)
        for i in range(oldT_w):
            for j in range(oldT_h):
                if not is_border(i,j):
                    T[i][j] = F0*(prev_T[i+1][j] + prev_T[i-1][j] + prev_T[i][j+1] + prev_T[i][j-1]) + prev_T[i][j]*(1 - 4*F0)
                if ISOLATED_BOTTOM_BORDER and is_border(i,j,"B"):
                    T[i][j] = F0*(2*prev_T[oldT_w-2][j] + prev_T[oldT_w-1][j+1] + prev_T[oldT_w-1][j-1]) + (1-4*F0)*prev_T[oldT_w-1][j]
                if ISOLATED_LEFT_BORDER and is_border(i,j,"L"):
                    T[i][j] = F0*(prev_T[i+1][0] + prev_T[i-1][0] + 2*prev_T[i][1]) + (1-4*F0)*prev_T[i][0]
                if ISOLATED_TOP_BORDER and is_border(i,j,"T"):
                    T[i][j] = F0*(2*prev_T[1][j] + prev_T[0][j+1] + prev_T[0][j-1]) + (1-4*F0)*prev_T[0][j]
                if ISOLATED_RIGHT_BORDER and is_border(i,j,"R"):
                    T[i][j] = F0*(prev_T[i+1][oldT_h-1] + prev_T[i-1][oldT_h-1] + 2*prev_T[i][oldT_h-2]) + (1-4*F0)*prev_T[i][oldT_h-1]

    return T

def main():
    props = read_props_file(sys.argv[1])
    result = solve(props)
    
    with open("Temp2D.out", "w") as outf:
      for i in result:
        for j in i:
          outf.write(str(j) + " ")
        outf.write("\n")

    # Plot heat map
    plt.imshow(result, cmap='bwr', interpolation='nearest')
    plt.show()

if __name__ == '__main__':
    main()
