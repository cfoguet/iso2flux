#!/usr/bin/env python


from setuptools import setup, find_packages

setup(name='iso2flux',
      version='0.1001',
      description='iso2flux',
      author='Carles Foguet',
      author_email='cfoguet@outlook.com',
      url='',
      scripts=["run_iso2flux.py"],
      install_requires=["cython","python-libsbml","openpyxl","numpy","scipy","matplotlib","pandas","nose","lxml","Tkinter","tkFileDialog","cobra"],
      packages=find_packages(),
     )
