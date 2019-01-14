import netCDF4
import numpy as np

class Model():
	def __init__(self, config):
		mesh = netCDF4.Dataset(config.mesh_file)
		variables = mesh.variables
		elem_map = list(variables['elem_map'])

		self.X = variables['coordx']
		self.Y = variables['coordy']
		self.Z = variables['coordz']
		self.connect = np.array(variables['connect1']) - 1    
		self.distance_to_boundary = dict()
		self.boundary_surface = dict()
		self.boundary_internal = dict()
		
		npt = len(self.X)
		point_type = dict()

		for node in self.connect:
			for point in node:
				if not point in point_type:
					point_type[point] = 1

				else:
					point_type[point] += 1

		for i in range(len(variables['ss_status'])):
			bname = variables['ss_names'][i]
			bid = variables['ss_prop1'][i]
			name = bname.tostring().decode('ascii').strip().rstrip('\x00')
			assert name in config.boundary
			btype = config.boundary[name]

			side = []
			for node in variables['elem_ss%d' % bid]:
				side.append(int(node))

			side = np.array(side)
			if btype == 0:
				self.boundary_surface[name] = side

			elif btype == 1:
				self.boundary_internal[name] = side

			self.distance_to_boundary[name] = np.zeros(npt)
			
			point_type2 = dict()
			for elem in side:
				node_id = elem_map.index(elem)
				node = self.connect[node_id]
				
				for point in node:
					if not point in point_type2:
						point_type2[point] = 1

					else:
						point_type2[point] += 1

			for j in range(npt):
				x0 = self.X[j]
				z0 = self.Z[j]

				min_dist = float('inf')
				min_point = None
				min_point2 = None
				min_node = []

				for elem in side:
					node_id = elem_map.index(elem)
					node = self.connect[node_id]
					
					for point in node:
						is_boundary = point_type[point] < 4
						if point_type2[point] == 1 and point_type[point] == 2:
							is_boundary = False

						if is_boundary:
							x1 = self.X[point]
							z1 = self.Z[point]

							if point == min_point:
								min_node.append(node_id)

							else:
								dist = (x1 - x0) ** 2 + (z1 - z0) ** 2
								if dist < min_dist:
									min_dist = dist
									min_point = point
									min_node = [node_id]
						
				
				assert min_point != None

				x2 = self.X[min_point]
				z2 = self.Z[min_point]
				self.distance_to_boundary[name][j] = np.sqrt((x2 - x0) ** 2 + (z2 - z0) ** 2)
				
					


				

		
	