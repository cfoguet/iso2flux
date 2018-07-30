from setuptools import setup, find_packages

setup(name='iso2flux',
      version='0.7.0',
      description='iso2flux',
      author='Carles Foguet',
      author_email='cfoguet@ub.edu',
      url='',
      scripts=["create_and_solve_iso2flux_model.py"],
      install_requires=["pygmo==2.3","cython==0.27.3","python-libsbml==5.16.0","openpyxl==2.5.0","numpy==1.14.0","scipy==0.19.0","sympy==1.0","lxml","cobra==0.7.0"],
      packages=find_packages(),
     )

#["create_and_solve_iso2flux_model.py"]
#scripts=["create_iso2flux_model.py","get_intervals.py","integrate_gene_expression.py","solve_iso2flux_label.py","p13cmfa.py","create_and_solve_iso2flux_model.py"],
