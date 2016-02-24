import os
import os.path
from setuptools import setup, find_packages


def get_package_data():
    package_data = ['resources/*.xls', 'resources/variant_summary.txt']
    for root, dirnames, filenames in os.walk('./eva_cttv_pipeline/resources/json_schema'):
        root = root.replace("./eva_cttv_pipeline/", "")
        for filename in filenames:
            new_fn = os.path.join(root, filename)
            package_data.append(new_fn)
    return package_data


setup(name='eva_cttv_pipeline',
      version='0.1',
      packages=find_packages(),
      install_requires=[
          'datetime',
          'jsonschema>=v2.5.0',
          'setuptools-git',
          'xlrd'
      ],
      package_data={
          'eva_cttv_pipeline': get_package_data()
      }
      )
