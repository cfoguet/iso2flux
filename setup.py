from setuptools import setup, find_packages

setup(name='iso2flux',
      version='0.6.1',
      description='iso2flux',
      author='Carles Foguet',
      author_email='cfoguet@ub.edu',
      url='',
      scripts=["create_iso2flux_model.py","get_intervals.py","integrate_gene_expression.py","solve_iso2flux_label.py","p13cmfa.py"],
      install_requires=["cloudpickle==0.4","pygmo==2.3","cython","python-libsbml","openpyxl","numpy","scipy","lxml","cobra==0.7"],
      packages=find_packages(),
     )
