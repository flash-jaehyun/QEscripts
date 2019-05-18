from setuptools import setup

setup(
   name="QEstruct"
   version="0.1"
   author_email="wwwennei@gmail.com"
   description="Mini module for analyzing structures from QE"
   packages=['struct']
   install_requires=[
         "ase",
         "numpy",
         "pymatgen"
   ],
)
