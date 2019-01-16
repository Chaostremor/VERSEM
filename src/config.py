"""This reads the input filue and source info.

Author: Congyue Cui


"""

import argparse
import re
import os.path as path
from imp import load_source
import numpy as np
import json

class Source():
	def __init__(self, src, config_dir):
		"""
		Class for source information
		stf_file: python file that defines source time function
		stf_args: arguments that are passed to stf_file
		location: location of the source
		term: force direction in 2D
		"""
		stf_file = src['stf_file']
		stf_name = re.sub('.py$', '', path.basename(stf_file))
		stf = load_source(stf_name, path.join(config_dir, stf_file))
		stf = getattr(stf, stf_name)
		stf_args = src['stf_args']
		self.stf = stf
		self.stf_args = stf_args
		self.location = src['location']
		self.term = src['term']

class Config():
	def __init__(self):
		"""
		Class for configurations
		"""
		parser = argparse.ArgumentParser()
		parser.add_argument('--config-file', nargs='?', default='./input/config.json')
		config_file = parser.parse_args().config_file
		
		assert path.exists(config_file)

		output = dict()
		boundary = dict()
		config_dir = path.dirname(config_file)

		with open(config_file) as f:
			data = json.load(f)

			self.mesh_file = path.join(config_dir, data['mesh_file'])
			self.material_file = path.join(config_dir, data['material_file'])
			self.output_dir = path.join(config_dir, data['output_dir'])
			self.boundary = data['boundary']

			self.sources = []
			self.stations = data['stations']
			for src in data['sources']:
				self.sources.append(Source(src, config_dir))
			
			self.nt = data['nt']
			self.dt = data['dt']

			self.ngll_x = data['ngll_x']
			self.ngll_y = data['ngll_y']
			self.ngll_z = data['ngll_z']
			self.dim = data['dim']
			self.solver = data['solver']
