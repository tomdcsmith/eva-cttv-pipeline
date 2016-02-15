import importlib
import importlib._bootstrap
import importlib.util
import os
import sys


def get_resource_file(package, resource):
    spec = importlib.util.find_spec(package)
    if spec is None:
        return None
    mod = (sys.modules.get(package) or importlib._bootstrap._load(spec))
    if mod is None or not hasattr(mod, '__file__'):
        return None

    parts = resource.split('/')
    parts.insert(0, os.path.dirname(mod.__file__))
    resource_name = os.path.join(*parts)
    return resource_name
