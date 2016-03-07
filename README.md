## README ##

Building
-------
1. "git clone git@github.com:EBIvariation/eva-cttv-pipeline.git"
2. "cd eva-cttv-pipeline"
3. and then:
	4. to install: "python3 setup.py install"
	5. to install to develop: "python3 setup.py develop"
	6. to build a source distribution: "python3 setup.py sdist"

Setting up virtual environment
-------
For a Python virtual environment to work with the pipeline:

1. "cd eva-cttv-pipeline"
2. "virtualenv -p python3.5 venv"
3. "source venv/bin/activate" ("venv/bin/deactivate" to deactivate virtualenv)
4. pip install -r /path/to/requirements.txt
