import os
import os.path
from setuptools import setup, find_packages


def get_package_data():
    package_data = ['resources/*.txt', 'resources/testing/*']
    for root, dirnames, filenames in os.walk('./eva_cttv_pipeline/resources/json_schema'):
        root = root.replace("./eva_cttv_pipeline/", "")
        for filename in filenames:
            new_fn = os.path.join(root, filename)
            package_data.append(new_fn)
    return package_data


def get_requires():
    requires = []
    with open("requirements.txt", "rt") as f:
        for line in f:
            requires.append(line.rstrip())
    return requires


setup(name='eva_cttv_pipeline',
      version='0.1',
      packages=find_packages(),
      install_requires=get_requires(),
      package_data={
          'eva_cttv_pipeline': get_package_data()
      },
      tests_require=get_requires(),
      setup_requires=get_requires()
      )
