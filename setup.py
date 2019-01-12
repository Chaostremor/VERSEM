""" :noindex:
Setup.py file with generic info
"""
import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "VERSEM",
    version = "0.2",
    author = "Congyue Cui, Srijan Bharati Das, Lucas Sawade, Chao Song, Fan Wu",
    author_email = "lsawade@princeton.edu",
    description =("A versatile Spectral Element Method"),
    license = "GNU",
    keywords = "SEM,Spectral Element Method, Wave Equation",
    url = "http://packages.python.org/rom2dectesttesttest",
    packages=['input','src', 'tests'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)
