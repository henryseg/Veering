from setuptools import setup

from pathlib import Path
long_description = (Path(__file__).parent / "README.md").read_text()

setup(name='veering',
      version='0.1',
      description='Taut and veering triangulations',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Anna Parlak, Henry Segerman, Saul Schleimer',
      author_email='segerman@math.okstate.edu',
      url='https://github.com/henryseg/Veering',
      packages=['veering'],
      package_data={'veering': ['data/*']},
      install_requires=['regina', 'snappy'],
      keywords='surfaces, manifolds, geometry, taut triangulation, veering triangulation'
     )
