# Makefile for the VERSEM Python suite
#
#

all: init test

# Initializing installation
init:
	pip install -r requirements.txt

# Testing installation
test:
	py.test tests

# Make Documentation
docs:
	cd docs
	make latexpdf
	make html

# Cleaning up
clean:
	rm -rf results/timesteps
	rm -f results/gll_coordinates.npy
	#cd docs 
	#make clean
	

.PHONY: init test
