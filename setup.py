# Update from https://packaging.python.org/en/latest/tutorials/packaging-projects/
#
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sympkf",
    version="1.0.3",
    author="Olivier Pannekoucke ",
    author_email="olivier.pannekoucke@meteo.fr",
    description="Symbolic Calculation for the Parametric Dynamics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/opannekoucke/sympkf",
    project_urls={
        "Bug Tracker": "https://github.com/opannekoucke/sympkf/issues",
    },                
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: CeCILL-B Free Software License Agreement (CECILL-B) ",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},        
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6"        
)
