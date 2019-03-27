

# IMPORT PACKAGES # 

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import itertools as it 
from scipy import signal 
from sklearn.neighbors import KNeighborsRegressor
from mpl_toolkits.mplot3d import Axes3D
import math 

plt.style.use('seaborn-muted')

import pathlib
 
script_dir = pathlib.Path(__file__).parent.resolve()


class funcs_3D:


	def __init__(self, dataset, x_label, y_label, 
				 z_label, head_property, tail_property,
				 nlags, lagdistance, lineartolerance,
				 htolerance, vtolerance, hband, vband,
				 azimuth, dip):

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
		self.dataset = dataset
		self.x_label = str(x_label)
		self.y_label = str(y_label)
		self.z_label = str(z_label)
		self.head_property = str(head_property)
		self.tail_property = str(tail_property)
		self.nlags = int(nlags) 
		self.lagdistance  = float(lagdistance)
		self.lineartolerance = lineartolerance
		self.htolerance = float(htolerance)
		self.vtolerance = float(vtolerance)
		self.hband = float(hband) 
		self.vband = float(vband) 
		self.azimuth = float(azimuth)
		self.dip = float(dip) 
		self.dist = pd.DataFrame()

	def distances(self):

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

		X = self.dataset[self.x_label]
		Y = self.dataset[self.y_label]
		Z = self.dataset[self.z_label]
		HEAD = self.dataset[self.head_property]
		TAIL = self.dataset[self.tail_property]
		points = zip(X,Y,Z,HEAD,TAIL,HEAD.index, TAIL.index)

		pairs = list(it.combinations(points, 2))
		distanceh  =[]
		distancexy = []

	
		distancex = np.array([(pair[0][0] - pair[1][0]) for pair in pairs])
		distancey = np.array([(pair[0][1] - pair[1][1]) for pair in pairs])
		distancez = np.array([(pair[0][2] - pair[1][2]) for pair in pairs])
		distancexy =np.array([np.sqrt(distancex[i]**2 + distancey[i]**2) for i in range(len(distancex))]) + 0.000001
		distanceh =np.array([np.sqrt(distancex[i]**2 + distancey[i]**2 + distancez[i]**2) for i in range(len(distancex))])+ 0.000001
		head_1 = np.array([pair[0][3] for pair in pairs])
		head_2 = np.array([pair[0][4] for pair in pairs])
		tail_1 = np.array([pair[1][3] for pair in pairs])
		tail_2 = np.array([pair[1][4] for pair in pairs])
		index_h = np.array([pair[0][5] for pair in pairs])
		index_t = np.array([pair[0][6] for pair in pairs])

		distance_dataframe  =pd.DataFrame(np.array([distancex, 
						  distancey, 
						  distancez, 
						  distancexy, 
						  distanceh, 
						  head_1,
						  head_2,
						  tail_1,
						  tail_2,
						  index_h,
						  index_t]).T, 
						  columns=['DX', 'DY', 'DZ', 'XY', 'H', 'Var 1 (head)', 'Var 2 (head)', 'Var 1 (tail)', 'Var 2 (tail)', 'INDEX HEAD', 'INDEX TAIL'])
		return distance_dataframe


	def permissible_pairs(self , lag_multiply):

		'''permissible_pairs
		Args:
		 lag_multiply (double): Mutliple of lag distance
		Returns:	
		 distances (pandas.DataFrame): Returns the permissible sample pairs for omnidirecional functions
		'''

		cos_Azimuth = np.cos(np.radians(90-self.azimuth))
		sin_Azimuth = np.sin(np.radians(90-self.azimuth))
		cos_Dip     = np.cos(np.radians(90-self.dip))
		sin_Dip     = np.sin(np.radians(90-self.dip))

		minimum_range = lag_multiply*self.lagdistance - self.lineartolerance
		maximum_range = lag_multiply*self.lagdistance + self.lineartolerance

		htol = np.abs(np.cos(np.radians(self.htolerance)))
		vtol=np.abs(np.cos(np.radians(self.vtolerance)))

		
		check_azimuth = np.abs((self.dist['DX']*cos_Azimuth + self.dist['DY']*sin_Azimuth)/self.dist['XY'])
		check_dip     = np.abs((self.dist['XY']*sin_Dip + self.dist['DZ']*cos_Dip)/self.dist['H'])
		check_bandh   = np.abs(cos_Azimuth*self.dist['DY']- sin_Azimuth*self.dist['DX'])
		check_bandv	  = np.abs(sin_Dip*self.dist['DZ'] - cos_Dip*self.dist['XY'])

	
		filter_dist = self.dist[(self.dist['H'] >= minimum_range) & 
							  (self.dist['H'] <= maximum_range) & 
							  (check_azimuth.values >= htol) &
							  (check_dip.values >= vtol) &
							  (check_bandh.values < self.hband)&
							  (check_bandv.values < self.vband)]

		filter_dist = filter_dist.dropna()

		return filter_dist


	def calculate_experimental(self , lag_multiply,  type_var):

		'''calculate_experimental
		Args:
		 lag_multiply (double): Mutliple of lag distance	
		 type_var (string): String containing the type of spatial continuity function to calculate
		 					5 admissible functions are possible:
							"Variogram"
			 				"Covariogram"
			 				"Correlogram"
			 				"PairWise"
			 				"RelativeVariogram"
		Returns:
		 value (double): Experimental continuity function value 
		 number_of_pairs (int) : number of pairs used to calculate the experimental function 
		 average_distace (double): average distance of experimental continuity function value  
		'''

		points = self.permissible_pairs(lag_multiply)
		
		if type_var not in['Variogram','Covariogram','Correlogram','Pair_Wise','Relative_Variogram']:
			raise Exception("Experimental continuity function not in admissible functions")

		if points.empty == False:
			number_of_pairs = points.shape[0]
			average_distance = points['H'].mean()
			value = 0
			if type_var == 'Variogram': 
				value = ((points['Var 1 (head)'] - points['Var 1 (tail)'])*(points['Var 2 (head)'] - points['Var 2 (tail)']))/(2*number_of_pairs)
				value = value.sum()
			if type_var == 'Covariogram': 
				value = ((points['Var 1 (head)'] - points['Var 1 (head)'].mean())*(points['Var 2 (tail)']-points['Var 2 (tail)'].mean()))/number_of_pairs
				value = value.sum()
			if type_var == 'Correlogram':
				value = ((points['Var 1 (head)'] - points['Var 1 (head)'].mean())*(points['Var 2 (tail)']-points['Var 2 (tail)'].mean()))/(number_of_pairs*points['Var 1 (head)'].var()*points['Var 2 (tail)'].var())
				value = value.sum()
			if type_var == 'Pair_Wise':
				value = 2*((points['Var 1 (head)'] - points['Var 1 (tail)'])/(points['Var 2 (head)'] + points['Var 2 (tail)']))**2/number_of_pairs
				value = value.sum()
			if type_var == 'Relative_Variogram':
				average_tail = (points['Var 1 (tail)'] +  points['Var 2 (tail)'])/2
				average_head = (points['Var 1 (head)'] +  points['Var 2 (head)'])/2
				value = 4*((points['Var 1 (head)'] - points['Var 1 (tail)'])*(points['Var 2 (head)'] - points['Var 2 (tail)']))/(number_of_pairs*(average_head + average_tail)**2)
				value = value.sum()
			return value, number_of_pairs, average_distance
		return None , None, None



	def calculate_experimental_function(self, type_var, plot=True, show_pairs=True):

		'''calculate_experimental_function
		Args:	
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

		self.dist = self.distances()
		values =[]
		number_of_pairs = []
		average_distance = []
		p = 0
		for i in range(0,self.nlags):
			value, pair, distance = self.calculate_experimental(i, type_var)
			if value != None:
				values.append(value)
				number_of_pairs.append(pair)
				average_distance.append(distance)
			p += 1
			

		df = pd.DataFrame(np.array([values,number_of_pairs,average_distance]).T, 
						  columns = ['Spatial continuity', 'Number of pairs', 'Average distance'])
		
		if plot == True:
			fig = plt.figure()
			ax = fig.add_subplot(111)

			ax.plot(df['Average distance'].values, df['Spatial continuity'].values)
			ax.set_xlabel('Lag distance (h)')
			ax.set_ylabel(type_var)
			ax.set_title('Experimental continuity function : {} , azimuth {}  and dip {} '.format(str(type_var), str(self.azimuth), str(self.dip)))
			if show_pairs == True:
				x, y = df['Average distance'].values, df['Spatial continuity'].values
				for i, j  in enumerate(df['Number of pairs'].values):
					ax.annotate(str(j), xy =(x[i], y[i]), xytext =(x[i], (y[i]+0.05*y[i])))
				ax.set_ylim((min(y),1.10*max(y)))
			plt.grid()
			plt.show()

		return df 


	def modelling(self, experimental_dataframe, rotation_reference, model_func, ranges, contribution, nugget, inverted= False, plot_graph = True ):

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


		y = math.cos(math.radians(self.dip))*math.sin(math.radians(self.azimuth))
		x = math.cos(math.radians(self.dip))*math.cos(math.radians(self.azimuth))  
		z = math.sin(math.radians(self.dip))

		angle_azimuth = math.radians(90-rotation_reference[0])
		angle_dip = -math.radians(90-rotation_reference[1])
		angle_rake = math.radians(90-rotation_reference[2])


		rotation1 = [[math.cos(angle_azimuth), - math.sin(angle_azimuth), 0],
					 [math.sin(angle_azimuth), math.cos(angle_azimuth), 0],
					 [0,0,1]]

		rotation2 = [[1, 0, 0],
					 [0, math.cos(angle_dip), math.sin(angle_dip)],
					 [0,-math.sin(angle_dip),math.cos(angle_dip)]]

		rotation3 = [[math.cos(angle_rake), 0, -math.sin(angle_rake)],
					 [0, 1, 0],
					 [math.sin(angle_rake),0,math.cos(angle_rake)]]

		rotation_transform = np.dot(rotation3,np.dot(rotation2,np.dot(rotation1, np.array([x,y,z]).T)))

		rotated_range =[]

		for i in ranges:
			rangex = float(i[0])
			rangey = float(i[1])
			rangez = float(i[2])


			rotated = (np.multiply(rotation_transform, [rangex, rangey, rangez]))
			rotated_range.append(math.sqrt(rotated[0]**2+rotated[1]**2+rotated[2]**2))

		distancemax = experimental_dataframe['Average distance'].max()
		distances = np.linspace(0, distancemax, 200)

		model = []

		if inverted == False:
			for i in distances:
				soma = 0
				for j, o, l  in zip(contribution, model_func, rotated_range):
					if o == 'Exponential':
						soma += j*(1-math.exp(-3*i/float(l)))
					if o == 'Gaussian':
						soma += j*(1-math.exp(-3*(i/float(l)**2)))
					if o == 'Spherical':
						if i <= l:
							soma += j*(1.5*i/float(l)-0.5*(i/float(l))**3)
						else:
							soma += j
				soma += nugget
				model.append(soma)
		else:
			for i in distances:
				soma = 0
				for j, o, l  in zip(contribution, model_func, rotated_range):
					if o == 'Exponential':
						soma += (j+nugget)*(math.exp(-3*i/float(l)))
					if o == 'Gaussian':
						soma += (j+nugget)*(math.exp(-3*(i/float(l)**2)))
					if o == 'Spherical':
						if i <= l:
							soma += (j+nugget)*(1 - (1.5*i/float(l)-0.5*(i/float(l))**3))
						else:
							soma += 0
				model.append(soma)

		df = pd.DataFrame(np.array([distances, model]).T, columns= ['distances', 'model']) 


		if plot_graph == True:
			x, y = experimental_dataframe['Average distance'].values, experimental_dataframe['Spatial continuity'].values
			plt.plot(distances, model, label= 'Model')
			plt.plot(x,y, label='Experimental Variogram')
			plt.legend()
			plt.xlabel("Lag distance (h)")
			plt.ylabel("Experimental continuity function value")
			plt.title("Experimental continuity function modelling")
			plt.grid()
			plt.show()


		return df 


	def covariogram_map_3d(self,property_value, neighbors, division = 20, alpha= 0.7,  cutx =[-np.inf, np.inf],cuty =[-np.inf,np.inf],cutz =[-np.inf,np.inf], size =20 ):

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
		X = self.dataset[self.x_label].values
		Y = self.dataset[self.y_label].values
		Z = self.dataset[self.z_label].values
		R = self.dataset[property_value].values 

		max_x, min_x = max(X), min(X)
		max_y, min_y = max(Y), min(Y)
		max_z, min_z = max(Z), min(Z)

		cordinatesx = np.linspace(min_x,max_x,division)
		cordinatesy = np.linspace(min_y,max_y, division)
		cordinatesz = np.linspace(min_z,max_z,division)

		cordinates = np.array([np.array([i,j,k]).T for i in cordinatesx for j in cordinatesy for k in cordinatesz])
		estimates = np.zeros(len(cordinates))

		nb = KNeighborsRegressor(n_neighbors=neighbors).fit(np.array([X,Y,Z]).T, R)

		for i, j in zip(range(len(estimates)), cordinates): 
			estimates[i] = nb.predict(j.reshape(1, -1))[0]
			
		estimates = estimates.reshape((division,division,division))
		
		fft_u  = np.conjugate(np.fft.fftn(estimates,axes=(0,1,2), norm ='ortho'))
		fft_u_roll = np.fft.fftn(estimates,axes=(0,1,2), norm ='ortho')
		product = np.multiply(fft_u_roll,fft_u)
		Covariance = np.fft.fftshift(np.fft.ifftn(product,axes=(0,1,2), norm ='ortho'))
		Covariance = np.real(Covariance)	
		Covariance = Covariance.reshape(-1,1)[:,0]/(division*division*division)

		filter_cord = (cordinates[:,0]>cutx[0]) & (cordinates[:,0]<cutx[1]) &  (cordinates[:,1]>cuty[0]) & (cordinates[:,1]<cuty[1]) &  (cordinates[:,2]>cutz[0]) & (cordinates[:,2]<cutz[1])

		cx = cordinates[:,0][filter_cord]
		cy = cordinates[:,1][filter_cord]
		cz = cordinates[:,2][filter_cord]

		Covariance = Covariance[filter_cord]
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		ax.set_xlabel("x")
		ax.set_ylabel("y")
		ax.set_zlabel("z")
		scat = ax.scatter(cx, cy, cz,  cmap='Spectral', s= size, marker="D", c=Covariance, alpha=alpha)
		fig.colorbar(scat)
		plt.show()



data = np.loadtxt("Ferro 3d.dat", skiprows= 1, unpack=False)
df = pd.DataFrame(data, columns=['X','Y','Z','fe'])



nlags = 10
lagdistance= 80
lineartolerance = 40
htolerance = 45.0
vtolerance = 45.0
hband = 40
vband = 40 
vband = 40
azimuth = 157.0
dip = 20.0

gamma = funcs_3D(df,'X','Y','Z','fe','fe',nlags,lagdistance,
                         lineartolerance, htolerance, vtolerance, 
                         hband, vband, azimuth,dip)



# Calculate variogram values 
variogram = gamma.calculate_experimental_function("Variogram")
print(variogram)



		


		