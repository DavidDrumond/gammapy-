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



```text
   """"funcs_3D
   Instances: 
   dataset (pandas.Dataframe): Input Dataframe containing the dataset 
   x_label (string): Label of x coordinates contained in dataset 
   y_label (string): Label of y coordinates contained in dataset 
   z_label (string): Label of z coordinates contained in dataset 
   head_property  (string): label of head property contained in dataset 
   tail property (string): label of tail property contained in dataset  
   lagdistance (double): lag size of experimental spatial functions
   lineartolerance (double): lag linear tolerance of experimental spatial functions 
   htolerance (double): angular horizontal tolerance of experimental spatial functions in degrees   
   vtolerance  (double): angular vertical tolerance of experimental spatial functions in degrees 
   hband (double): horizontal bandwidth of spatial functions  
   vband  (double): vertical bandwidth of spatial functions  
   azimuth  (double): Azimuth value for experimental continuity function in degrees 
   dip (double): Dip value for experimental continuity functions in degrees 

   Methods:
   distances(self) : Calculate the matrix distance of all pairs 
   permissible_pairs_omni (self, lag_multiply) : Calculate the permissible sample pairs for omnidirecional functions for irregular grids 
   permissible_pairs(self , lag_multiply) : Calculate the permissible sample pairs for directional functions for irregular grids
   hscatter(self, lag_multiply) : Calculate the hscatterplot for a distance multiple of a lag
   calculate_experimental(self , lag_multiply,  type_var) : Calculate the experimental continuity function value for a distance multiple of a lag and a type of variogram 
   calculate_experimental_omini(self , lag_multiply,  type_var) : Calculate the experimental omnidirecional continuity function value for a distance multiple of a lag and a type of variogram
   calculate_experimental_function(self, type_var) : Calculate the experimental continuity function for all lag values 
   calculate_experimental_function_omni(self, type_var) : Calculate the omnidirecional experimental continuity function for all lag values

   """
```

## funcs\_3D.distances\(\)

Returns the 3Dimensional distances between all samples 

```text
		'''distances
		Returns:	
		 distance_dataframe (pandas.DataFrame): Pandas Dataframe containing all the distance metrics
		 DX (pandas.DataFrame.Series) : Difference of x cartesian coordinates 
		 DY (pandas.DataFrame.Series) : = diference of y cartesian values from the head and tails of the vector  
		 DZ (pandas.DataFrame.Series) : = diference of z cartesian values from the head and tails of the vector 
		 XY (pandas.DataFrame.Series) : = Distance projection on XY plane of the vector  
		 H  (pandas.DataFrame.Series) : = Distance value from head and tail of vector  
		 Var 1 (head) (pandas.DataFrame.Series) : Value from variable 1 on the head of vector 
		 Var 2 (head) (pandas.DataFrame.Series) : Value from variable 2 on the head of vector  
		 Var 1 (tail) (pandas.DataFrame.Series) : Value from variable 1 on the tail of vector 
		 Var 2 (tail) (pandas.DataFrame.Series) : Value form variable 2 on the tail of vector 
		 INDEX HEAD   (pandas.DataFrame.Series) : Index of propertie 1 sample 
		 INDEX TAIL   (pandas.DataFrame.Series) : Index of propertie 2 sample
		'''
```

## funcs\_3D.permissible\_pairs\_omni \(lag\_multiply\)

```text
		'''permissible_pairs_omni
		Args:
		 lag_multiply (double): Mutliple of lag distance
		Returns:	
		 distances (pandas.DataFrame): Returns the permissible sample pairs for omnidirecional functions
		'''
```



## funcs\_3D.permissible\_pairs \(lag\_multiply\)

```text
    Returns the set of permissible pairs for a lag multiplication of lag_multiply
    for a directional search strategy

    >> lag_multiply = integer containing theh multiple of lag distance 
```

