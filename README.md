---
description: >-
  Gammapy is a python library to calculate spatial continuity functions commonly
  used in geostatistical approach.
---

# gammapy-

## Installation

Spatial continuity functions are the representation of spatial dependencies among statistical variables and the space. They could be considered as the backbone of geostatistical analysis. Gammapy is a simple library to help spatial continuity calculations based on GamV from GSLib 90 package. 

To install the package, clone the folder from repository and run the command in the gammapy\dist folder. 

```text
pip install gammapy-0.0.0tar.gz
```

## Dependencies 

Gammapy needs several libraries to run properly. To better run gammapy functions, verify if you have installed these following dependencies.

* numpy 1.15.4
* pandas 0.23.4
* matplotlib 3.0.2
* itertools 
* scipy 1.1.0
* sklearn 0.20.1
* mpl\_toolkits 

## Experimental Continuity functions 

The 3D experimental functions are related to funcs\_3D object. To create experimental spatial continuity functions you should inform the following attributes:

* **dataset** = pandas Dataframe containing \(X,Y,Z coordinates and Head and tail properties\)
* **x\_label** = A string containing the name of X coordinates
* **y\_label** = A string containing the name of Y coordinates
* **z\_label** = A string containing the name of Z coordinates
* **head\_property** = string containing the name of the first propertie
* **tail property** = string containing the name of the second propertie
* **nlags** = integer containing the number of lags in experimental continuity functions
* **lagdistance** = float containing the lag size
* **lineartolerance** = float containing the linear tolerance
* **htolerance** = float containing the horizontal tolerance value in degrees
* **vtolerance** = float containing the vertical tolerance value in degrees
* **hband** = float containing the horizontal band width
* **vband** = float containing the vertical band width
* **azimuth** = float containing the experimental function azimuth in degrees
* **dip** = float containing the experimental function dip in degrees

#### If you have any problem with the class parameters, please use the help\(\) command 

The following methods are used to calculate spatial continuity functions:

## funcs\_3D.distances\(\)

```text

Returns a pandas dataframe with the distances of each pair of 
samples containing in dataset

>> Output 
>> ............................................................................. 
DX  = diference of x cartesian values from the head and tails of the vector 
DY  = diference of y cartesian values from the head and tails of the vector  
DZ  = diference of z cartesian values from the head and tails of the vector 
XY  = Distance projection on XY plane of the vector  
H   = Distance value from head and tail of vector  
Var 1 (head) = Value from variable 1 on the head of vector 
Var 2 (head) = Value from variable 2 on the head of vector  
Var 1 (tail) = Value from variable 1 on the tail of vector 
Var 2 (tail) = Value form variable 2 on the tail of vector 
INDEX HEAD   = Index of propertie 1 sample 
INDEX TAIL   = Index of propertie 2 sample

```



## funcs\_3D.permissible\_pairs\_omni \(lag\_multiply\)

```text
    '''
    Returns the set of permissible pairs for a lag multiplication of lag_multiply
    for a omnidirecional search strategy


    >> Input 
    >>............................................................................
    >> lag_multiply = integer containing theh multiple of lag distance 
    '''
```



## funcs\_3D.permissible\_pairs \(lag\_multiply\)

```text
    Returns the set of permissible pairs for a lag multiplication of lag_multiply
    for a directional search strategy

    >> lag_multiply = integer containing theh multiple of lag distance 
```

