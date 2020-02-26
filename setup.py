# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

setup(name='iso2flux',
      version='0.7.1',
      description='iso2flux',
      author='Carles Foguet',
      author_email='cfoguet@ub.edu',
      url='',
      scripts=["create_iso2flux_model.py","get_intervals.py","integrate_gene_expression.py","solve_iso2flux_label.py","p13cmfa.py","run_iso2flux_analysis.py"],
      install_requires=["pygmo==2.5","cython==0.23.4","python-libsbml==5.16.0","openpyxl==2.3.3","numpy==1.13.3","scipy==0.19.0","sympy==1.0","lxml","cobra==0.9.0"],
      packages=find_packages(),
     )

#["create_and_solve_iso2flux_model.py"]
#scripts=["create_iso2flux_model.py","get_intervals.py","integrate_gene_expression.py","solve_iso2flux_label.py","p13cmfa.py","create_and_solve_iso2flux_model.py"],
