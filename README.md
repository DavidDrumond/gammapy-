---
description: >-
  Gammapy is a python library to calculate spatial continuity functions commonly
  used in geostatistical approach.
---

# gammapy-

## Installation

Spatial continuity functions are the representation of spatial dependencies among statistical variables and the space. They could be considered as the backbone of geostatistical analysis. Gammapy is a simple library to help spatial continuity calculations based on GamV from GSLib 90 package. 

To install the package, clone the folder from repository and run the command in the gammapy\dist folder. You can also just compy the gammapy file to your directory and import to your project. 

```text
pip install gammapy-0.1-py3-none-any.whl
```

## Example

This is a simple example....

```text
import numpy as np 
import pandas as pd 
from gammapy import gammapy 


data = np.loadtxt('Ferro 3d.dat', skiprows= 1, unpack=False)
df = pd.DataFrame(data, columns=['X','Y','Z','fe'])

nlags = 10
lagdistance= 80
lineartolerance = 40
htolerance = 45.0
vtolerance = 45.0
hband = 40
vband = 40 
vband = 40
azimuth = 157
dip = 30

gamma = gammapy.funcs_3D(df,'X','Y','Z','fe','fe',nlags,lagdistance,
                         lineartolerance, htolerance, vtolerance, 
                         hband, vband, azimuth,dip)

variogram = gamma.calculate_experimental_function("Variogram")


gamma.modelling(experimental_dataframe = variogram, 
                rotation_reference = [157,30,0],
                model_func = ["Spherical", "Spherical"],
                ranges =[[100,50,50],[250,50,50]],
                contribution = [10,30],
                nugget = 5,
                inverted = False)

print(help(gamma.covariogram_map_3d))
gamma.covariogram_map_3d(df, 'X', 'Y', 'Z', 'fe',plot=True, cuty=[600,650])
```

![Covariance map ](.gitbook/assets/cov_map.png)

![Variogram modelling](.gitbook/assets/modelling%20%281%29.png)

You might also check the jupyter notebook example at the directory...

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

The following methods and attributes are used to calculate spatial statistics. Because spatial continuity functions demands high computer performance, it could be not allowed, for most machineries, calculate spatial statistics for a high quantity sample size. If dataset have more than 20.000 samples, it could be possible randomly select subsets from the dataset, obtained the same results. To better achieve high sample size spatial continuity functions,  change the choice attribute, and select a random sample according your wishes. 

Obs: Dip values must be inserted with negative values

```text

	"""
       funcs_3D
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
	   choice (double): Random size of sampling dataset if data is greater than 20.000
	   self.choice (int) : Select a number of samples to perform random method for high number of data
   Methods:
       calculate_experimental_omini(self , lag_multiply,  type_var) : Calculate the experimental omnidirecional continuity function value for a distance multiple of a lag and a type of variogram
       modelling(self, experimental_dataframe, rotation_reference, model_func, ranges, contribution, nugget, inverted= False, plot_graph = True ): Modelling spatial continuity functions 
       covariogram_map_3d(self,property_value, neighbors, division = 20, alpha= 0.7,  cutx =[-np.inf, np.inf],cuty =[-np.inf,np.inf],cutz =[-np.inf,np.inf], size =20 ): Estimate covariance maps in three dimensions

	"""
```



## funcs\_3D.calculate\_experimental\_function\(self, type\_c, omni = False, plot\_graph=False, show\_pairs=False\)

This function computes the experimental continuity functions according a select model. 5 kinds of model are reliable, \(Variogram, Covariogram, Correlogram, PairWise and RelativeVariogram\).  If omnidirecional variogram will be performed, select the omni property equals True. This function could show the number of sample pairs on graph selecting the plot\_graph equals True. 

```text
'''calculate_experimental_function
Args:	
 omni (bool)       = Boolean for selecting omnidirecional experimental function, if True select omnidirecional
 plot_graph (bool) = Boolean for selecting plotting experimental values 
 show_pairs (bool) = Boolean for selecting plotting experimental number of pairs 
 type_var (string): String containing the type of spatial continuity function to calculate
 					5 admissible functions are possible:
					"Variogram"
	 				"Covariogram"
	 				"Correlogram"
	 				"PairWise"
	 				"RelativeVariogram" 
Returns:
 df (pandas.DataFrame): Pandas Dataframe containing the experimental continuity functions of all lags
'''
```



## funcs\_3D.modelling\(self, experimental\_dataframe, rotation\_reference, model\_func, ranges, contribution, nugget, inverted= False, plot\_graph = True \)

This function could be used for modelling the experimental functions. Rotation axis criteria is the same of GSLib program. Azimuth could be measure as clockwise the Y coordinate, dip values are negative to downward and Rake values are clockwise to the maximum direction.  For more information, see the geostats reference:

[http://geostatisticslessons.com/lessons/anglespecification](http://geostatisticslessons.com/lessons/anglespecification)

To plot Covariogram values, the inverted parameter could be change to True. If user want to plot the modelled graph, please change the plot\_graph to True. 

```text
'''plot_experimental_function)_omni
Args:
 experimental dataframe (pandas.DataFrame): Pandas DataFrame containing experimental continuity functions     
 rotation reference (list(azimuth, dip, rake)): List containing the reference of principal directions angles in degrees 
 model_func(list(string)) : List containing the models for all structures. size of the list must be the same of the number of structures
                             3 admissible functions are possible:
                             "Spherical"
                            "Gaussian"
                            "Exponential"
 ranges(list(list(maximum range, medium range, minimum range))) : list of lists containing the maximum, medium and minimum range for each number of structures
 contribution (list): list of contributions for each strucutre 
 nugget (double): Nugget effect value 
 inverted (bool): If true plot model according covariogram form, otherwise plot model according the variogram form
 plot_graph (bool): If true plot the experimental variogram and the spatial continuity model

Returns:
 plot (matplotlib.pyplot): Plot of omnidirecional experimental continuity function and the spatial continuity model for one direction 
'''
```

## funcs\_3D.covariogram\_map\_3d\(self, df, x\_label, y\_label, z\_label, property\_value, plot= False, division = 20, alpha= 0.7, cutx =\[-np.inf, np.inf\],cuty =\[-np.inf,np.inf\],cutz =\[-np.inf,np.inf\] \)

An easy way to see the principal components of spatial anisotropy is to performed covariogram maps in three dimensions. This function interpolate the spatial data using Knearest neighbors, and calculating the Fourier Transform from the data. To change the number of neighbors please change the neighbors parameter.

```text
		'''covariogram_map_3d
		Args:
		property_value(string): String containing the property to create the covariogram map
		neighbors (int) : Number of neighbors using in KNearest neighbors 
		division(int, optional): discretize number of covariogram map
		size(int, optional): size of bullet
		alpha(float, optional): the level of transparency (0- transparent, 1-solid)
		cutx (list, optional): list containing the minimum cutsize and the maximum cutsize for x coordinates
		cuty (list, optional): list containing the minimum cutsize and the maximum cutsize for y coordinates
		cutz (list, optional): list containing the minimum cutsize and the maximum cutsize for z coordinates
						
		Returns:
		plot (matplotlib.pyplot): Plot of Covariance map in three dimensional scale 
		'''
```



