import os
import os.path
from setuptools import setup, find_packages


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


def process_schema():
    json_schema_dir = utilities.get_resource_file("eva_cttv_pipeline", "resources/json_schema")
    # local_schema_dir = str(os.path.join(str(Path(json_schema_dir).parent), "schema_copy"))
    local_schema_dir = str(os.path.join(os.path.dirname(json_schema_dir), "schema_copy"))

    os.makedirs(local_schema_dir, exist_ok=True)
    print("Copying json schema")
    print(json_schema_dir + " to " + local_schema_dir)
    copy_dir(json_schema_dir, local_schema_dir)

    change_json_refs(local_schema_dir)


class Install(_install):
    def run(self):
        _install.run(self)
        print("############################## MAIN INSTALL FINISHED ##############################")
        process_schema()


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
