from setuptools import setup


setup(name='rough_surfaces',
      version='0.1',
      description='Analysis, contact and flow - fractures and rough surfaces',
      author='Philipp S. Lang',
      author_email='plang@slb.com',
      download_url='https://github.com/plang85/rough_surfaces.git',
      install_requires=['numpy>=1.9.1',
                        'scipy>=0.14',
                        'matplotlib'], # hate this here TODO somehow get rid of plotting stuff
      extras_require={
          'test': ['pytest>=3.6.0',
                   'pytest-pep8',
                   'pytest-xdist',
                   'pytest-cov',
                   'codecov'],
      },
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Intended Audience :: Developers',
          'Intended Audience :: Engineers',
          'Programming Language :: Python :: 3.6',
          'Topic :: Software Development :: Libraries',
          'Topic :: Software Development :: Libraries :: Python Modules'
      ],
      packages=['rough_surfaces'])