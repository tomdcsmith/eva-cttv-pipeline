from setuptools import setup

setup(name='eva_cttv_pipeline',
      version='0.1',
      packages=['eva_cttv_pipeline'],
      install_requires=[
          'codecs',
          'datetime',
          'http.client',
          'json',
          'jsonschema>=v2.5.0',
          'optparse',
          'os',
          'setuptools-git',
          'sys',
          'time',
          'urllib.error',
          'urllib.parse',
          'urllib.request',
          'xlrd',
      ],
      package_data={
        'eva_cttv_pipeline': [
            'resources/*.xls',
            'resources/variant_summary.txt'
        ]
      },
      include_package_data=True)
