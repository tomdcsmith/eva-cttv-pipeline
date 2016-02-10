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
          'sys',
          'time',
          'urllib.error',
          'urllib.parse',
          'urllib.request',
          'xlrd'
      ])
