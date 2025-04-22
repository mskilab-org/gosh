from setuptools import setup, find_packages

setup(
    name='gosh-cli',
    version='0.1.4',
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
