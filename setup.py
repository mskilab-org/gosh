import os
import re
from setuptools import setup, find_packages

def get_version():
    version_file = os.path.join(os.path.dirname(__file__), "gosh_cli", "__init__.py")
    with open(version_file, "r") as f:
        contents = f.read()
    version_match = re.search(r'^__version__ = [\'"]([^\'"]+)[\'"]', contents, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

setup(
    name='gosh-cli',
    version=get_version(),
    description='gOSh - gOS sHell: A CLI tool for the nf-gOS pipeline',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Shihab Dider',
    packages=find_packages(where='.'),
    package_dir={'': '.'},
    include_package_data=True,
    package_data={'gosh_cli': ['utils/*.txt']},
    install_requires=[
        'Click',
        'openai',
    ],
    entry_points={
        'console_scripts': [
            'gosh=gosh_cli.main:cli',
        ],
    },
    python_requires='>=3.6',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
