import subprocess
import sys

try:
    from conans import client
except ImportError:
    subprocess.call([sys.executable, '-m', 'pip', 'install', 'conan'])
    from conans import client

try:
    from skbuild import setup
except ImportError:
    subprocess.call([sys.executable, '-m', 'pip', 'install', 'scikit-build'])
    from skbuild import setup
    
with open('README.md', 'r') as fh:
    long_description = fh.read()

setup_requirements = ['scikit-build>=0.11.1',
                      'cmake>=3.15',
                      'conan>=1.25'
                     ]

setup (name = 'QSystem',
       version='1.2.0r1',
       cmake_source_dir='.',
       author='Evandro Chagas Ribeiro da Rosa, Bruno Gouvêa Taketani',
       author_email='ev.crr97@gmail.com',
       description='A quantum computing simulator for Python',
       long_description=long_description,
       long_description_content_type='text/markdown',
       url='https://gitlab.com/evandro-crr/qsystem',
       packages=['qsystem'],
       setup_requires=setup_requirements,
       classifiers=[
            'Programming Language :: Python :: 3',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
        ]
       )

