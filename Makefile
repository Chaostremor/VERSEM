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
docs: *.rst
	sphinx-apidoc -o ./docs/source/ ./ [*.npy]
	make latexpdf -C docs
	make html -C docs

# Cleaning up
clean:
	rm -rf results/timesteps
	rm -f results/gll_coordinates.npy
	#cd docs 
	#make clean
	

.PHONY: init test
