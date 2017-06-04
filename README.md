# Temp2D

### Input file
#### File with extension [.in] containing following properties:
```sh
X_SIZE 0.5
Y_SIZE 0.5
THERMAL_CONDUCTIVITY 1
DENSITY 1
SPECIFIC_HEAT_CAPACITY 1
XY_STEP 0.1
TOTAL_TIME 200
TIME_STEP 0.001
INIT_TEMPERATURE_INNER 0
INIT_TEMPERATURE_TOP_BORDER 100
INIT_TEMPERATURE_BOTTOM_BORDER 0
INIT_TEMPERATURE_LEFT_BORDER 75
INIT_TEMPERATURE_RIGHT_BORDER 50
ISOLATED_TOP_BORDER 0
ISOLATED_BOTTOM_BORDER 1
ISOLATED_LEFT_BORDER 0
ISOLATED_RIGHT_BORDER 0
```

| Property | Value |
| --- | --- |
| X_SIZE | Width of 2D surface (Unit: m) |
| Y_SIZE | Height of 2D surface (Unit: m) |
| THERMAL_CONDUCTIVITY | Material thermal conductivity (Unit: W/(m·K)) |
| DENSITY | Material density (Unit: kg/m<sup>3</sup>) |
| SPECIFIC_HEAT_CAPACITY | Material specific heat capacity (Unit: J/(kg·K)) |
| XY_STEP | Interval between discrete points in X and Y (Unit: m)|
| TOTAL_TIME | Total simulation time (Unit: s)|
| TIME_STEP | Interval between time points (Unit: s)|
| INIT_TEMPERATURE_INNER | Initial temperature for non-border points (Unit: ºC)|
| INIT_TEMPERATURE_*_BORDER | Initial temperature for border points (Unit: ºC)|
| ISOLATED_*_BORDER | Is border isolated? (0 or 1) |

### Output file
#### File with extension [.out] containing final temperature matrix:
```sh
100.0 100.0 100.0 100.0 100.0 
75.0 83.410925947 82.6286022085 74.2614414109 50.0 
75.0 76.0151015796 72.8420414759 64.4171634352 50.0 
75.0 72.8074388954 68.3072986803 60.5651708541 50.0 
75.0 71.9073553216 67.0145434958 59.536221301 50.0 
```
### Output plot
#### Heat map showing temperature distribution in 2D surface

![Temperature Plot](http://i.imgur.com/qWLWa7U.png)
