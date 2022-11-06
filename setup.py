import pathlib
from setuptools import setup, find_packages
import os


package_name='Binny'
def get_package_version():

    dir_name=os.path.dirname(os.path.abspath(__file__))
    init_path = os.path.join(dir_name, package_name.lower(), '__init__.py')
    package_version=None
    with open(init_path) as file:
        for line in file:
            if '__version__' in line:
                package_version=line.replace('__version__','')
                package_version=package_version.strip('\n')
                package_version=package_version.strip()
                package_version=package_version.strip('=')
                package_version=package_version.strip()
                package_version=package_version.strip('"')
                package_version=package_version.strip('"')
    return package_version



# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text(encoding='utf-8')
LICENSE = (HERE / "LICENSE").read_text(encoding='utf-8')

long_description='Binny - Automated binning algorithm to recover high-quality genomes from complex metagenomic datasets.'

setup(
    name=package_name,
    version=get_package_version(),
    author="Oskar Hickl, Anna Heintz-Buschart",
    author_email="oskar.hickl@uni.lu",
    description="Automated binning tool",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/a-h-b/binny",
    project_urls={
        "Bug Tracker": "https://github.com/a-h-b/binny/issues",
    },
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        'Operating System :: MacOS',
        'Operating System :: POSIX :: Linux',
    ],
    license=LICENSE,
    include_package_data=True,
    install_requires=['networkx==2.6.3',
                      'snakemake==7.16.0',
                      'pyarrow==7.0.0',
                      'hdbscan==0.8.29',
                      'scikit-bio==0.5.6',
                      'opentsne==0.6.1'
                      ],
    
    entry_points={
        "console_scripts": [
            "binny=binny.__main__:main",
        ],
    },
)