from setuptools import setup

setup(name='eva_cttv_pipeline',
      version='0.1',
      packages=['eva_cttv_pipeline'],
      install_requires=[
          'datetime',
          'jsonschema>=v2.5.0',
          'setuptools-git',
          'xlrd'
      ],
      package_data={
        'eva_cttv_pipeline': [
            'resources/*.xls',
            'resources/variant_summary.txt'
        ]
      },
      include_package_data=True)
