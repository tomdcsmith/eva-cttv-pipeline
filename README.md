## README ##

[![Build Status](https://travis-ci.org/EBIvariation/eva-cttv-pipeline.svg?branch=master)](https://travis-ci.org/EBIvariation/eva-cttv-pipeline)
[![Coverage Status](https://coveralls.io/repos/github/EBIvariation/eva-cttv-pipeline/badge.svg?branch=master)](https://coveralls.io/github/EBIvariation/eva-cttv-pipeline?branch=master)


Minimum Python version needed: 3.5


Building and (optional) Setting up virtual environment
-------

1. "git clone --recursive git@github.com:EBIvariation/eva-cttv-pipeline.git"
2. "cd eva-cttv-pipeline"
3. [OPTIONAL] "virtualenv -p python3.5 venv"
4. [OPTIONAL] "source venv/bin/activate" ("venv/bin/deactivate" to deactivate virtualenv)
5. pip install -r requirements.txt
6. and then one of:
    7. to install: "python3 setup.py install"
    8. to install to develop: "python3 setup.py develop"
    9. to build a source distribution: "python3 setup.py sdist"


Usage
-------

Please see the [GitHub wiki](https://github.com/EBIvariation/eva-cttv-pipeline/wiki/How-to-submit-an-OpenTargets-batch) for usage