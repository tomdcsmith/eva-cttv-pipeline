import errno
import subprocess
import importlib
import importlib._bootstrap
import importlib.util
import os
import sys
import shutil
import eva_cttv_pipeline.config as config


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


def copy_and_overwrite(from_path, to_path):
    if os.path.exists(to_path):
        shutil.rmtree(to_path)
    shutil.copytree(from_path, to_path)


def copy_dir(src, dest):
    try:
        copy_and_overwrite(src, dest)
    except OSError as e:
        # If the error was caused because the source wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dest)
        else:
            print('Directory not copied. Error: %s' % e)


def change_json_refs(local_schema_dir):

    command = "find " + local_schema_dir + " -type f -exec sed -i -e \"s/https:\/\/raw.githubusercontent.com\/CTTV\/json_schema\/master/file:\/\/" + local_schema_dir.replace("/", "\/") + "/g\" {} \;"
    subprocess.check_output(command, shell=True)

    evidence_base_json = os.path.join(local_schema_dir, "evidence/base.json")
    evidence_base_json_temp = os.path.join(local_schema_dir, "evidence/base_temp.json")
    command = "grep -v '\"id\": \"base_evidence\"' " + evidence_base_json + " > " + evidence_base_json_temp + "; mv " + evidence_base_json_temp + " " + evidence_base_json
    subprocess.check_output(command, shell=True)

    command = "grep -v '\"id\": \"#single_lit_reference\"' " + evidence_base_json + " > " + evidence_base_json_temp + "; mv " + evidence_base_json_temp + " " + evidence_base_json
    subprocess.check_output(command, shell=True)

    command = "sed -i -e \"s/evidence\/base.json#base_evidence\/definitions\/single_lit_reference/evidence\/base.json#definitions\/single_lit_reference/g\""
    subprocess.check_output(command, shell=True)


def create_local_schema():
    json_schema_dir = get_resource_file(__package__, "resources/json_schema")
    # local_schema_dir = str(os.path.join(str(Path(json_schema_dir).parent), "schema_copy"))
    local_schema_dir = str(os.path.join(os.path.dirname(json_schema_dir), config.local_schema))

    os.makedirs(local_schema_dir, exist_ok=True)
    print("Copying json schema")
    print(json_schema_dir + " to " + local_schema_dir)
    copy_dir(json_schema_dir, local_schema_dir)

    change_json_refs(local_schema_dir)


def check_for_local_schema():
    local_schema = get_resource_file(__package__, "resources/" + config.local_schema)
    if not os.path.exists(local_schema):
        create_local_schema()
