# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})

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
    TIME_N_POINTS = TOTAL_TIME / TIME_STEP + 1
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

    # Init matrixes used in FDM
    oldT = np.full((int(X_N_POINTS), int(Y_N_POINTS)), INIT_TEMPERATURE_INNER)
    curT = oldT

    # Create simulation steps
    time_steps = np.linspace(0, TOTAL_TIME, TOTAL_TIME / TIME_STEP + 1)

    # Initialize border temperature and inner temperature
    for i in range(X_N_POINTS):
        for j in range(Y_N_POINTS):
            if i == Y_N_POINTS - 1:
                oldT[i][j] = INIT_TEMPERATURE_BOTTOM_BORDER
            elif i == 0:
                oldT[i][j] = INIT_TEMPERATURE_TOP_BORDER
            elif j == 0:
                oldT[i][j] = INIT_TEMPERATURE_LEFT_BORDER
            elif j == X_N_POINTS - 1:
                oldT[i][j] = INIT_TEMPERATURE_RIGHT_BORDER

    # FDM
    for t in range(len(time_steps)):
        for i in range(X_N_POINTS):
            for j in range(Y_N_POINTS):
                # Check if is border isolated; make sure we are in border; ignore corners
                if ISOLATED_TOP_BORDER and (i == 0 and j > 0 and j < X_N_POINTS - 1):
                    curT[i][j] = F0*(2*oldT[i+1][j] + oldT[i][j+1] + oldT[i][j-1]) + (1 - 4*F0)*oldT[i][j]
                if ISOLATED_LEFT_BORDER and (j == 0 and i > 0 and i < Y_N_POINTS - 1):
                    curT[i][j] = F0*(oldT[i+1][j] + oldT[i-1][j] + 2*oldT[i][j+1]) + (1 - 4*F0)*oldT[i][j]
                if ISOLATED_RIGHT_BORDER and (j == X_N_POINTS - 1 and i > 0 and i < Y_N_POINTS - 1):
                    curT[i][j] = F0*(oldT[i+1][j] + 2*oldT[i-1][j] + oldT[i][j-1]) + (1 - 4*F0)*oldT[i][j]
                if ISOLATED_BOTTOM_BORDER and (i == Y_N_POINTS - 1 and j > 0 and j < X_N_POINTS - 1):
                    curT[i][j] = F0*(2*oldT[i-1][j] + oldT[i][j+1] + oldT[i][j-1]) + (1 - 4*F0)*oldT[i][j]
                # If not in border, use another equation
                if i > 0 and j > 0 and i < Y_N_POINTS - 1 and j < X_N_POINTS - 1:
                    curT[i][j] = F0*(oldT[i+1][j] + oldT[i-1][j] + oldT[i][j+1] + oldT[i][j-1]) + (1 - 4*F0)*oldT[i][j]
        oldT = curT

    return curT

def main():
    print("STATUS: Solving 2D...")

    # Solve using file passed by user
    props = read_props_file(sys.argv[1])
    T = solve(props)

    print("STATUS: Done!")
    print(T)

    np.savetxt("tms.out", T)
    print("\nSTATUS: Saved data to tms.out")

    # Plot heat map
    plt.imshow(T, cmap='bwr', interpolation='nearest')
    plt.show()

if __name__ == '__main__':
    main()
