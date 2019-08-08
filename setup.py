from setuptools import setup

setup(name='simplicial_operators',
      version='0.1',
      description='A package to define and use simplicial operators',
      author='Anibal Medina',
      author_email='anibal.medinamardones@epfl.ch',
      url = 'https://github.com/Medina-Mardones/simplicial_operators',
      download_url = 'https://github.com/Medina-Mardones/simplicial_operators/archive/v_01.tar.gz',
      keywords = ['topology', 'simplicial complexes', 'topological data analysis', 'homology'],
      install_require=[],
      license='MIT',
      packages=['simplicial_operators'],
      calssifier=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Mathematicians',
          'Topic :: Software Development :: Build Tools',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6'
      ],
      zip_safe=False)
