from setuptools import setup
from setuptools.extension import Extension

with open('README.md', 'r') as fh:
    long_description = fh.read()

ext_module = Extension('_qsystem',
        sources=['src/qsystem.cpp',
                 'src/gate.cpp',
                 'src/microtar.c', 
                 'src/qs_ancillas.cpp',
                 'src/qs_errors.cpp',
                 'src/qs_evol.cpp',
                 'src/qs_utility.cpp'],
        include_dirs=['armadillo-code/include'],
        extra_compile_args=['-std=c++17']
        )

setup (name = 'QSystem',
       version='1.0.0rc3',
       author='Evandro Chagas Ribeiro da Rosa',
       author_email='ev.crr97@gmail.com',
       description='A Python quantum computer simulator',
       long_description=long_description,
       long_description_content_type='text/markdown',
       url='https://gitlab.com/evandro-crr/qsystem',
       ext_modules = [ext_module],
       packages=['qsystem'],
       classifiers=[
            'Programming Language :: Python :: 3',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
        ]
       )

