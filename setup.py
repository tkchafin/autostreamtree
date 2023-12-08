from setuptools import setup, find_packages

setup(
    name='autostreamtree',
    version='1.0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'autostreamtree=autostreamtree.cli:main',
            'networkDimensions=autostreamtree_scripts.networkDimensions:main',
            'streeToDendrogram=autostreamtree_scripts.streeToDendrogram:main',
        ],
    },
    install_requires=[
        "numpy",
        "pandas>=2.0",
        "networkx>=3.0",
        "seaborn",
        "matplotlib",
        "geopandas",
        "pyogrio",
        "momepy",
        "scipy",
        "scikit-learn",
        "mantel",
        "pysam",
        "sortedcontainers"
    ],
    extras_require={
        'dev': ['pytest']
    },
    package_data={
        'autostreamtree': ['data/*'],
    },
    include_package_data=True,
    author='Tyler Chafin',
    author_email='tyler.chafin@bioss.ac.uk',
    description='A package for fitting genetic distances to spatial networks',
    keywords='population genetics, genetic distance, river networks',
    url='https://github.com/tkchafin/autostreamtree',
)