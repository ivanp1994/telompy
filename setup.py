
from setuptools import setup, find_packages

# with open("README.md", "r") as _f:
#    long_description = _f.read()

with open("requirements_dev.txt", "r") as _f:
    extra_reqs = _f.read().strip().split("\n")

with open("requirements.txt", "r") as _f:
    install_reqs = _f.read().strip().split("\n")

setup(name="TelOMpy",
      version="0.0.1",
      description="Python based calculator for length of telomers from BNGO de novo assembly",
      packages=find_packages(),
      #      long_description=long_description,
      #      long_description_content_type='text/markdown',
      #      url="",
      author="Ivan Pokrovac",
      author_email="ivan.pokrovac.fbf@gmail.com",
      license="MIT",
      classifiers=["Development Status :: 3 - Alpha",
                   "License :: OSI Approved :: MIT License",
                   "Programming Language :: Python :: 3.8",
                   "Programming Language :: Python :: 3 :: Only",
                   "Intended Audience :: Science/Research"],
      install_requires=install_reqs,
      python_requires=">=3.8",
      extras_require={"extra": extra_reqs,
                      },
      entry_points={
          'console_scripts': [
              'telompy = telompy.cli:command_line_target'
          ]}
      )
