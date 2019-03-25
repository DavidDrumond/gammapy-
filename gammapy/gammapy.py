

# IMPORT PACKAGES # 

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import itertools as it 
from scipy import signal 
from sklearn.neighbors import KNeighborsRegressor
from mpl_toolkits.mplot3d import Axes3D

plt.style.use('seaborn-muted')



class funcs_3D:

	def __help__():

		'''
		Continuity functions 3D is a Python object to perform experimental continuity functions in 3d datasets:\n
		See the below attributes:
		-----------------------------------------------------------------------------\n
			>> dataset = pandas dataframe containing  (X,Y,Z cordinates and Head and tail properties)\n 
	        >> x_label = string containing the name of X cordinante\n  
	        >> y_label = string containing the name of Y cordinate\n  
	        >> z_label = string containing the name of Z cordinante\n  
	        >> head_property = string containing the name of the first propertie\n  
	        >> tail property = string containing the name of the second propertie\n  
	        >> nlags = integer containing the number of lags in experimental continuity functions\n  
	        >> lagdistance = float containing the lag size\n  
	        >> lineartolerance = float containing the linear tolerance\n 
	        >> htolerance = float containing the horizontal tolerance value in degrees\n  
	        >> vtolerance = float containing the vertical tolerance value in degrees\n  
	        >> hband = float containing the horizontal band width\n  
	        >> vband = float containing the vertical band width\n  
	        >> azimuth = float containing the experimental function azimuth in degrees\n  
	        >> dip = float containing the experimental function dip in degrees\n  
		---------------------------------------------------------------------------------------------------------\n 
		'''

	def __init__(self, dataset, x_label, y_label, 
				 z_label, head_property, tail_property,
				 nlags, lagdistance, lineartolerance,
				 htolerance, vtolerance, hband, vband,
				 azimuth, dip):


		
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

	def distances(self):

		'''
		Returns a pandas dataframe with the distances of each pair of samples containing in dataset

		>> Output 
		>> ........................................................................................
		>> Pandas DataFrame containing 
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

	
		distancex = [(pair[0][0] - pair[1][0]) for pair in pairs]
		distancey = [(pair[0][1] - pair[1][1]) for pair in pairs]
		distancez = [(pair[0][2] - pair[1][2]) for pair in pairs]
		distancexy =[np.sqrt(distancex[i]**2 + distancey[i]**2) for i in range(len(distancex))]
		distanceh =[np.sqrt(distancex[i]**2 + distancey[i]**2 + distancez[i]**2) for i in range(len(distancex))]
		head_1 = [pair[0][3] for pair in pairs]
		head_2 = [pair[0][4] for pair in pairs]
		tail_1 = [pair[1][3] for pair in pairs]
		tail_2 = [pair[1][4] for pair in pairs]
		index_h = [pair[0][5] for pair in pairs]
		index_t = [pair[0][6] for pair in pairs]

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

	def permissible_pairs_omni (self, lag_multiply):

		'''
		Returns the set of permissible pairs for a lag multiplication of lag_multiply
		for a omnidirecional search strategy


		>> Input 
		>>............................................................................
		>> lag_multiply = integer containing theh multiple of lag distance 
		'''

		distances = self.distances()
		minimum_range = lag_multiply*self.lagdistance - self.lineartolerance
		maximum_range = lag_multiply*self.lagdistance + self.lineartolerance

	
		distances = distances[(distances['H'] >= minimum_range) & 
							  (distances['H'] <= maximum_range)]

		distances = distances.dropna()

		return distances




	def permissible_pairs(self , lag_multiply):

		'''
		Returns the set of permissible pairs for a lag multiplication of lag_multiply
		for a directional search strategy

		>> lag_multiply = integer containing theh multiple of lag distance 
		'''

		distances = self.distances()
		cos_Azimuth = np.cos(np.radians(90-self.azimuth))
		sin_Azimuth = np.sin(np.radians(90-self.azimuth))
		cos_Dip     = np.cos(np.radians(90-self.dip))
		sin_Dip     = np.sin(np.radians(90-self.dip))

		minimum_range = lag_multiply*self.lagdistance - self.lineartolerance
		maximum_range = lag_multiply*self.lagdistance + self.lineartolerance

		htol = np.abs(np.cos(np.radians(self.htolerance)))
		vtol=np.abs(np.cos(np.radians(self.vtolerance)))

		
		check_azimuth = np.abs((distances['DX']*cos_Azimuth + distances['DY']*sin_Azimuth)/distances['XY'])
		check_dip     = np.abs((distances['XY']*sin_Dip + distances['DZ']*cos_Dip)/distances['H'])
		check_bandh   = np.abs(cos_Azimuth*distances['DY']- sin_Azimuth*distances['DX'])
		check_bandv	  = np.abs(sin_Dip*distances['DZ'] - cos_Dip*distances['XY'])

	
		distances = distances[(distances['H'] >= minimum_range) & 
							  (distances['H'] <= maximum_range) & 
							  (check_azimuth.values >= htol) &
							  (check_dip.values >= vtol) &
							  (check_bandh.values < self.hband)&
							  (check_bandv.values < self.vband)]

		distances = distances.dropna()

		return distances

	def hscatter(self, lag_multiply):

		''' 
			Create a hscatterplot for a lag distance 
			>> lag_multiply = integer containing theh multiple of lag distance 
		'''


		df = self.permissible_pairs(lag_multiply)
		figure, ax = plt.subplots(nrows= 2, ncols= 2, figsize=(10,10))


		correlations = df[['Var 1 (head)','Var 2 (head)','Var 1 (tail)','Var 2 (tail)']].corr().values


		print (ax)
		ax[0][0].scatter(df['Var 1 (head)'].values, df['Var 2 (head)'].values )
		ax[0][0].set_xlabel('Var 1 (head)')
		ax[0][0].set_ylabel('Var 2 (head)')
		ax[0][0].set_title('Correlation :  {}'.format(str(correlations[0][1])))

		ax[0][1].scatter(df['Var 1 (head)'], df['Var 2 (tail)'] )
		ax[0][1].set_xlabel('Var 1 (head)')
		ax[0][1].set_ylabel('Var 2 (tail)')
		ax[0][1].set_title('Correlation :  {}'.format(str(correlations[0][3])))

		ax[1][0].scatter(df['Var 1 (tail)'], df['Var 2 (head)'] )
		ax[1][0].set_xlabel('Var 1 (tail)')
		ax[1][0].set_ylabel('Var 2 (head)')
		ax[1][0].set_title('Correlation :  {}'.format(str(correlations[2][1])))

		ax[1][1].scatter(df['Var 1 (tail)'], df['Var 2 (tail)'] )
		ax[1][1].set_xlabel('Var 1 (tail)')
		ax[1][1].set_ylabel('Var 2 (tail)')
		ax[1][1].set_title('Correlation :  {}'.format(str(correlations[2][3])))
		plt.show()

	def calculate_experimental(self , lag_multiply,  type_var):

		''''
		Calculate experimental continuity functions for a lag muntiplication and type of function type 
		for an omnidirecional variogram 

		>> lag_multiply = integer containing theh multiple of lag distance 
		>> type_var = String containing the type of spatial continuity function to calculate 
				 			5 admissible functions are possible:
				 				"Variogram"
				 				"Covariogram"
				 				"Correlogram"
				 				"PairWise"
				 				"RelativeVariogram" 
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


	def calculate_experimental_omini(self , lag_multiply,  type_var):

		'''
		Calculate omnidirecional experimental continuity functions for a lag muntiplication and type of function type 
		for an omnidirecional variogram 

		>> lag_multiply = integer containing theh multiple of lag distance 
		>> type_var = String containing the type of spatial continuity function to calculate 
				 			5 admissible functions are possible:
				 				"Variogram"
				 				"Covariogram"
				 				"Correlogram"
				 				"PairWise"
				 				"RelativeVariogram" 
		'''

		points = self.permissible_pairs_omni(lag_multiply)
		
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
			if type_var == 'PairWise':
				value = 2*((points['Var 1 (head)'] - points['Var 1 (tail)'])/(points['Var 2 (head)'] + points['Var 2 (tail)']))**2/number_of_pairs
				value = value.sum()
			if type_var == 'RelativeVariogram':
				average_tail = (points['Var 1 (tail)'] +  points['Var 2 (tail)'])/2
				average_head = (points['Var 1 (head)'] +  points['Var 2 (head)'])/2
				value = 4*((points['Var 1 (head)'] - points['Var 1 (tail)'])*(points['Var 2 (head)'] - points['Var 2 (tail)']))/(number_of_pairs*(average_head + average_tail)**2)
				value = value.sum()
			return value, number_of_pairs, average_distance
		return None , None, None

	def calculate_experimental_function(self, type_var):

		'''
		Calculate directional experimental functions for a type of continuity function type
		>> type_var = String containing the type of spatial continuity function to calculate 
				 			5 admissible functions are possible:
				 				"Variogram"
				 				"Covariogram"
				 				"Correlogram"
				 				"PairWise"
				 				"RelativeVariogram"
		'''


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
		return df 

	def calculate_experimental_function_omni(self, type_var):

		'''
		Calculate ominidrecional experimental functions for a type of continuity function type
		> type_var = String containing the type of spatial continuity function to calculate 
				 			5 admissible functions are possible:
				 				"Variogram"
				 				"Covariogram"
				 				"Correlogram"
				 				"PairWise"
				 				"RelativeVariogram"
		'''


		values =[]
		number_of_pairs = []
		average_distance = []
		p = 0
		for i in tqdm(range(0,self.nlags), desc ="Calculating lag : "):
			value, pair, distance = self.calculate_experimental_omini(i, type_var)
			if value != None:
				values.append(value)
				number_of_pairs.append(pair)
				average_distance.append(distance)
			p += 1
			

		df = pd.DataFrame(np.array([values,number_of_pairs,average_distance]).T, 
						  columns = ['Spatial continuity', 'Number of pairs', 'Average distance'])
		return df 

	def Database():
		df =pickle.load('database_example', 'wb')
		return df 

	def plot_experimental_function(self, type_var, show_pairs = False):

		'''
			Function to plot directional experimental continuity functions
			>> type_var = String containing the type of spatial continuity function to calculate 
				 			5 admissible functions are possible:
				 				"Variogram"
				 				"Covariogram"
				 				"Correlogram"
				 				"PairWise"
				 				"RelativeVariogram"
			>> show_pairs = Boolean to show the number of spatial continuity pairs for each spatial continuity value on the graph 
							False = Doesn`t show 
							True = Show the pairs on the graph  

		'''

		df = self.calculate_experimental_function(type_var)
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


	def plot_experimental_function_omni(self, type_var, show_pairs = False):

		'''
			Function to plot omnidirecional experimental continuity functions
			>> type_var = String containing the type of spatial continuity function to calculate 
				 			5 admissible functions are possible:
				 				"Variogram"
				 				"Covariogram"
				 				"Correlogram"
				 				"PairWise"
				 				"RelativeVariogram"
			>> show_pairs = Boolean to show the number of spatial continuity pairs for each spatial continuity value on the graph 
							False = Doesn`t show 
							True = Show the pairs on the graph  

		'''

		df = self.calculate_experimental_function_omni(type_var)
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

	def modelling(self, experimental_dataframe, rotation_reference, model_func, ranges, contribution, nugget, inverted= False, plot_graph = True ):

		'''
			Define the model for an experimental function
			>> experimental dataframe = pandas dataframe containing the experimental values from a spatial continuity function 
			>> rotation reference = list containing the rotation reference from principal components [Azimuth, dip, rake]
			>> model_func = list containing the model functions for each spatial structure 
							three admissible functions are possible:
								"Spherical"
								"Gaussian"
								"Exponential"
			>> ranges = list of list containing the ranges from the principal components for each structure'
						example = [[max range structure 1, med range structure 1, min range structure 1], [max range structure 2, med range structure 2, min range structure 2]]
			>> contribution = list containing the constribution of each structure 
			>> nugget = float containing the nugget effect value 
			>> inverted = Bool to choice the way of plot : False -> for variograms ,  True -> for covariograms
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


	def covariogram_map_3d(self, df, x_label, y_label, z_label, property_value, plot= False, division = 20, alpha= 0.7, cutx =[-np.inf, np.inf],cuty =[-np.inf,np.inf],cutz =[-np.inf,np.inf] ):

		''''	
		Create 3D covariogram maps and returns the value in a np.array 
		Only calculate direct covariogram values 

		>> df = pandas dataframe containing the dataset to calculate covariogram maps 
		>> xlabel = string containing the label for X cordinate 
		>> ylabel = string containing the label for the Y cordinate 
		>> zlabel = string contanining the label for the Z cordinate 
		>> propertie_value =  string containing the propertie the label of the propertie to calculate 
		>> plot = if True plot a 3D grid image for covariogram maps 
		>> division =  Integer number of discretization points of Covariogram map, standart is 20 points in x, y and z 
		>> alpha = Float number containing the trasnparency of map 
		>> cutx = list containing the x cut size of covariogram map [minx, maxx]
		>> cuty = list containing the y cut size of covariogram map [miny, maxy]
		>> cutz = list containing the z cut labels of covariogram map [minz, maxz]


		'''

		X = df[x_label].values
		Y = df[y_label].values
		Z = df[z_label].values
		R = df[property_value].values 

		max_x, min_x = max(X), min(X)
		max_y, min_y = max(Y), min(Y)
		max_z, min_z = max(Z), min(Z)

		cordinatesx = np.linspace(min_x,max_x,division)
		cordinatesy = np.linspace(min_y,max_y, division)
		cordinatesz = np.linspace(min_z,max_z,division)

		cordinates = np.array([np.array([i,j,k]).T for i in cordinatesx for j in cordinatesy for k in cordinatesz])
		estimates = np.zeros(len(cordinates))

		nb = KNeighborsRegressor().fit(np.array([X,Y,Z]).T, R)

		for i, j in zip(range(len(estimates)), cordinates): 
			estimates[i] = nb.predict(j.reshape(1, -1))[0]
			
		estimates = estimates.reshape((division,division,division))
		
		fft_u  = np.conjugate(np.fft.fftn(estimates))
		fft_u_roll = np.fft.fftn(estimates)
		product = np.multiply(fft_u_roll,fft_u)
		Covariance = np.fft.fftshift(np.fft.ifftn(product))
		Covariance = np.real(Covariance)

		if plot == True:
			Covariance = Covariance.reshape(-1,1)[:,0]
			
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
			scat = ax.scatter(cx, cy, cz,  cmap='viridis', s= 40, marker="D", c=Covariance, alpha=alpha)
			cbar = fig.colorbar(scat)
			plt.show()
		
		return Covariance

			



				










		


		