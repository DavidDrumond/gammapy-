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

{% tabs %}
{% tab title="First Tab" %}
asdfadfadf
{% endtab %}

{% tab title="Second Tab" %}

{% endtab %}
{% endtabs %}



{% api-method method="get" host="" path="" %}
{% api-method-summary %}
distances
{% endapi-method-summary %}

{% api-method-description %}
Return the distance matrix for all samples containined in the dataset 
{% endapi-method-description %}

{% api-method-spec %}
{% api-method-request %}
{% api-method-path-parameters %}
{% api-method-parameter name="None" type="string" required=false %}

{% endapi-method-parameter %}
{% endapi-method-path-parameters %}
{% endapi-method-request %}

{% api-method-response %}
{% api-method-response-example httpCode=200 %}
{% api-method-response-example-description %}

{% endapi-method-response-example-description %}

```
All the distances between samples. 
H = Euclidian distance 
XY = projection in XY plane 
DX = Difference of X coordinates 
DY = Difference of Y coordinates 
DZ= Difference of Z coordinates 
Var 1 head = First property at the head of vector 
Var 2 head = Second property at the head of vector
Var 1 tail = First property in the tail of vector 
Var 2 tail = Second property in the tail of vector
Index Head = index of the sample  considered in the head 
Index Tail = index of the sample considered in the tail 
```
{% endapi-method-response-example %}
{% endapi-method-response %}
{% endapi-method-spec %}
{% endapi-method %}



{% api-method method="get" host="" path="" %}
{% api-method-summary %}
permissible\_pairs\_omni / permissible\_pairs\_
{% endapi-method-summary %}

{% api-method-description %}
Return the distance matrix for all pairs of samples allow to be in the direction \(azimuth/dip\) considered. This function is necessary to calculate directional experimental functions.
{% endapi-method-description %}

{% api-method-spec %}
{% api-method-request %}
{% api-method-path-parameters %}
{% api-method-parameter name="lag\_multipy" type="number" required=true %}
Number of lags to calculate the permissible pairs 
{% endapi-method-parameter %}
{% endapi-method-path-parameters %}
{% endapi-method-request %}

{% api-method-response %}
{% api-method-response-example httpCode=200 %}
{% api-method-response-example-description %}

{% endapi-method-response-example-description %}

```

```
{% endapi-method-response-example %}
{% endapi-method-response %}
{% endapi-method-spec %}
{% endapi-method %}







