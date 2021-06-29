#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('docs/history.rst') as history_file:
    history = history_file.read()

requirements = [
    'numpy >= 1.19.1',
    'scipy >= 1.5.0',
    'matplotlib >= 3.2.2',
    'pandas >= 1.1.1',
    'mstools >= 0.1.0']

setup_requirements = []

test_requirements = []

setup(
    author="E. A. Montoya-Araque\\\A. J. Aparicio-Ortube\\\D. G. Zapata-Medina\\\L. G. Arboleda-Monsalve",
    author_email='eamontoyaa@gmail.com',
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Application software for calculating the preconsolidation pressure from the incremental loading oedometer testing",
    entry_points={
        'console_scripts': [
            'pysigmap=pysigmap.cli:main',
        ],
    },
    install_requires=requirements,
    license="BSD license",
    long_description_content_type='text/x-rst',
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords=['Preconsolidation pressure', 'Yield stress',
              'incremental loading oedometer testing', 'Python', 'application software'],
    name='pysigmap',
    packages=find_packages(include=['pysigmap', 'pysigmap.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/eamontoyaa/pysigmap',
    version='0.1.8',
    zip_safe=False,
)
