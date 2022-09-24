from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as rdme:
    readme = rdme.read()

setup(name = "vector",
      version = "1.0.0",
      author = "prb1509",
      description = "A small package for vector operations in various coordinate systems.",
      long_description = readme,
      long_description_content_type = "text/markdown",
      url = "https://github.com/prb1509/vector",
      packages = find_packages(),
      license = "GNU General Public License v3 (GPLv3)",
      )
      

