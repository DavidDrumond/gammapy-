

# IMPORT PACKAGES # 

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from sklearn.neighbors import KNeighborsRegressor
from mpl_toolkits.mplot3d import Axes3D
import math 
import numba
import warnings



plt.style.use('seaborn-muted')

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
	   choice (double): Random size of sampling dataset if data is greater than 20.000
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
		if type(dataset) is pd.DataFrame:
			self.dataset = dataset
		else:
			raise ValueError("Dataset is not a pandas DataFrame")
		try:
			self.x_label = str(x_label)
			self.y_label = str(y_label)
			self.z_label = str(z_label)
			self.head_property = str(head_property)
			self.tail_property = str(tail_property)
		except:
			raise ValueError("Labels must be a string columns of dataset Pandas DataFrame")
		try:
			self.nlags = int(nlags) 
			self.lagdistance  = float(lagdistance)
			self.lineartolerance = float(lineartolerance) 
			self.htolerance = float(htolerance)
			self.vtolerance = float(vtolerance)
			self.hband = float(hband) 
			self.vband = float(vband) 
			self.azimuth = float(azimuth)
		except:
			raise ValueError("Experimental parameters must be numbers")
		self.dip = float(dip) 
		self.__dist = []
		self.__type_var = "Variogram"
		self.__omni = False
		self.__htol = 0
		self.__vtol=0
		self.__check_azimuth = 0
		self.__check_dip     = 0
		self.__check_bandh   = 0
		self.__check_bandv	  = 0
		self.choice	  = 1500


	@numba.jit(fastmath=True)
	def hdist(self , distancex, distancey, distancez):
		dist =np.zeros(distancex.shape[0])
		for i in range(distancex.shape[0]):
			dist[i] = np.sqrt(distancex[i]**2 + distancey[i]**2 + distancez[i]**2) + 0.0000000001
		return dist

	@numba.jit(fastmath=True)
	def xydist(self , distanceh,):
		cos = np.cos(np.radians(self.dip))
		xy = np.zeros(distanceh.shape[0])
		for i in range(distanceh.shape[0]):
			xy[i] = np.abs(distanceh[i]*cos)
		return xy

	@numba.jit(fastmath=True)
	def xdist(self , pairs):
		dist =np.zeros(pairs.shape[0])
		for i in range(pairs.shape[0]):
			dist[i] = (pairs[i][0][0] - pairs[i][1][0])
		return dist

	@numba.jit(fastmath=True)
	def ydist(self , pairs):
		dist =np.zeros(pairs.shape[0])
		for i in range(pairs.shape[0]):
			dist[i] = (pairs[i][0][1] - pairs[i][1][1])
		return dist

	@numba.jit(fastmath=True)
	def zdist(self , pairs):
		dist =np.zeros(pairs.shape[0])
		for i in range(pairs.shape[0]):
			dist[i] = (pairs[i][0][2] - pairs[i][1][2])
		return dist

	@numba.jit(fastmath=True)
	def combin(self , points,n, max_dist):
		dist =[]
		p = 0
		for i in range(0,n):
			for j in range((i+1),n):
				if (points[i][0] - points[j][0]) < max_dist:
					if (points[i][1] - points[j][1]) < max_dist:
						if (points[i][2] - points[j][2]) < max_dist:
							dist.append((points[i] , points[j]))
				p += 1
		return np.array(dist)

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

		max_dist = (self.nlags + 1)*self.lagdistance
		X = self.dataset[self.x_label].values
		Y = self.dataset[self.y_label].values
		Z = self.dataset[self.z_label].values
		HEAD = self.dataset[self.head_property].values
		TAIL = self.dataset[self.tail_property].values
		if (X.shape[0] > 20000):
			warnings.warn("Warning .... dataset is too big (> 20000)! Applying random method ")
			choice = int(self.choice)
			if choice < 20000:
				points = np.array(list(zip(X,Y,Z,HEAD,TAIL)))
				points = np.array([points[np.random.randint(0,len(points))] for i in range(choice)])
			else:
				raise ValueError("Subset of random selected data is greater than datasize")
		else:
			points = np.array(list(zip(X,Y,Z,HEAD,TAIL)))
		pairs = self.combin(points,points.shape[0], max_dist)
		distancex = self.xdist(pairs) 
		distancey =  self.ydist(pairs) 
		distancez =  self.zdist(pairs)
		distanceh = self.hdist(distancex, distancey, distancez)
		distancexy = self.xydist(distanceh)
		head_1 = np.array([pair[0][3] for pair in pairs])
		head_2 = np.array([pair[0][4] for pair in pairs])
		tail_1 = np.array([pair[1][3] for pair in pairs])
		tail_2 = np.array([pair[1][4] for pair in pairs])
		distance_dataframe =  np.array([distancex, 
						distancey, 
						distancez, 
						distancexy, 
						distanceh, 
						head_1,
						head_2,
						tail_1,
						tail_2,]).T
		return distance_dataframe[distanceh[:,] < max_dist ]

	def __permissible_pairs (self , lag_multiply):

		'''permissible_pairs
		Args:
		 lag_multiply (double): Mutliple of lag distance
		Returns:	
		 distances (numpy array): Returns the permissible sample pairs for omnidirecional functions
		'''
		minimum_range = lag_multiply*self.lagdistance - self.lineartolerance
		maximum_range = lag_multiply*self.lagdistance + self.lineartolerance
		if self.omni == False:
			filter_dist = self.dist[(self.dist[:,4] >= minimum_range) & 
							  (self.dist[:,4] <= maximum_range) & 
							  (self.check_azimuth >= self.htol) &
							  (self.check_dip >= self.vtol) &
							  (self.check_bandh < self.hband)&
							  (self.check_bandv < self.vband)]
		else:
			filter_dist = self.dist[(self.dist[:,4] >= self.minimum_range) & 
							  (self.dist[:,4] <= self.maximum_range)]
		return filter_dist

	def __calculate_experimental(self , lag_multiply):

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

		points = self.__permissible_pairs(lag_multiply)
		
		if self.type_var not in['Variogram','Covariogram','Correlogram','PairWise','Relative_Variogram']:
			raise Exception("Experimental continuity function not in admissible functions")

		if points.size  != 0:
			number_of_pairs = float(points.shape[0])
			average_distance = points[:,4].mean()
			value_exp = 0
			if self.type_var == 'Variogram': 
				value_exp = ((points[:,5] - points[:,7])*(points[:,6] - points[:,8]))/(2*number_of_pairs)
				value_exp = value_exp.sum()
			if self.type_var == 'Covariogram': 
				value_exp = ((points[:,5] - points[:,5].mean())*(points[:,8]-points[:,8].mean()))/number_of_pairs
				value_exp = value_exp.sum()
			if self.type_var == 'Correlogram':
				value_exp = ((points[:,5] - points[:,5].mean())*(points[:,8]-points[:,8].mean()))/(number_of_pairs*points[:,5].var()*points[:,8].var())
				value_exp = value_exp.sum()
			if self.type_var == 'PairWise':
				value_exp = 2*((points[:,5] - points[:,7])/(points[:,6] + points[:,8]))**2/number_of_pairs
				value_exp = value_exp.sum()
			if self.type_var == 'Relative_Variogram':
				average_tail = (points[:,7] +  points[:,8])/2
				average_head = (points[:,5] +  points[:,6])/2
				value_exp = 4*((points[:,5] - points[:,7])*(points[:,6] - points[:,8]))/(number_of_pairs*(average_head + average_tail)**2)
				value_exp = value_exp.sum()
			return [value_exp, number_of_pairs, average_distance]
		return [np.nan , np.nan, np.nan]

	def calculate_experimental_function(self, type_c, omni = False, plot_graph=False, show_pairs=False):

		'''calculate_experimental_function
		Args:	
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

		self.dist = self.distances()
		self.type_var = type_c
		self.omni = omni
		cos_Azimuth = np.cos(np.radians(90-self.azimuth))
		sin_Azimuth = np.sin(np.radians(90-self.azimuth))
		cos_Dip     = np.cos(np.radians(90-self.dip))
		sin_Dip     = np.sin(np.radians(90-self.dip))
		self.htol = np.abs(np.cos(np.radians(self.htolerance)))
		self.vtol= np.abs(np.cos(np.radians(self.vtolerance)))
		self.check_azimuth = np.abs((self.dist[:,0]*cos_Azimuth + self.dist[:,1]*sin_Azimuth)/self.dist[:,3])
		self.check_dip     = np.abs((self.dist[:,3]*sin_Dip + self.dist[:,2]*cos_Dip)/self.dist[:,4])
		self.check_bandh   = np.abs(cos_Azimuth*self.dist[:,1]- sin_Azimuth*self.dist[:,0])
		self.check_bandv	  = np.abs(sin_Dip*self.dist[:,2] - cos_Dip*self.dist[:,3])

		number_of_int = range(1, (self.nlags +1))
		value_exp = np.array(list(map(self.__calculate_experimental, number_of_int)))
		df = pd.DataFrame(value_exp, 
						  columns = ['Spatial continuity', 'Number of pairs', 'Average distance'])
		df = df.dropna()
		if plot_graph == True:
			fig = plt.figure()
			ax = fig.add_subplot(111)
			ax.plot(df['Average distance'].values, df['Spatial continuity'].values)
			ax.set_xlabel('Lag distance (h)')
			ax.set_ylabel(type_c)
			ax.set_title('Experimental continuity function : {} , azimuth {}  and dip {} '.format(str(type_c), str(self.azimuth), str(self.dip)))
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

		if len(ranges[0]) != 3:
			raise ValueError("Variogram ranges must range 3 principal directions")
		if len(ranges) == len(contribution) and len(ranges) == len(model_func):
			pass
		else:
			raise ValueError("Variogram structures must be the same size")


		y = math.cos(math.radians(self.dip))*math.cos(math.radians(self.azimuth))
		x = math.cos(math.radians(self.dip))*math.sin(math.radians(self.azimuth))  
		z = math.sin(math.radians(self.dip))

		angle_azimuth = math.radians(rotation_reference[0])
		angle_dip = math.radians(rotation_reference[1])
		angle_rake = math.radians(rotation_reference[2])


		rotation1 = np.array([[math.cos(angle_azimuth), -math.sin(angle_azimuth), 0],
					 [math.sin(angle_azimuth), math.cos(angle_azimuth), 0],
					 [0,0,1]])

		rotation2 = np.array([[1, 0, 0],
					 [0, math.cos(angle_dip), math.sin(angle_dip)],
					 [0,-math.sin(angle_dip),math.cos(angle_dip)]])

		rotation3 = np.array([[math.cos(angle_rake), 0, -math.sin(angle_rake)],
					 [0, 1, 0],
					 [math.sin(angle_rake),0,math.cos(angle_rake)]])

		rot1 = np.dot(rotation1.T, np.array([x,y,z]))
		rot2 = np.dot(rotation2.T, rot1)
		rot3= np.dot(rotation3.T,rot2)


		print(rot3)
		rotated_range =[]

		for i in ranges:
			rangex = float(i[0])
			rangey = float(i[1])
			rangez = float(i[2])


			rotated = (np.multiply(rot3, np.array([rangex, rangey, rangez]).T))
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



nlags = 10
lagdistance= 10
lineartolerance = 20
htolerance = 45.0
vtolerance = 45.0
hband = 20
vband = 20
azimuth = 45
dip = 0


data = np.loadtxt("Níquel 3D.txt", skiprows= 1, unpack=False)
df = pd.DataFrame(data, columns=['X','Y','Z','fe'])
display(df.describe())
display(df['fe'].var())

gamma3 = funcs_3D(df,'X','Y','Z','fe','fe',nlags,lagdistance,
                         lineartolerance, htolerance, vtolerance, 
                         hband, vband, azimuth,dip)





gamma3.choice = 500
variogram3 = gamma3.calculate_experimental_function("Covariogram", plot_graph = False, show_pairs = True)



# Calculate average covariogram maps in three dimensions 

rotation_reference = [45,0,0]
ranges = [[45,0,0]]
model_func = ["Spherical"]
contribution = [0.08]
nugget = 0.10


model_value3 = gamma3.modelling(experimental_dataframe = variogram3, 
                rotation_reference = rotation_reference ,
                model_func = model_func,
                ranges =ranges ,
                contribution = contribution,
                nugget = nugget,
                inverted = True)
