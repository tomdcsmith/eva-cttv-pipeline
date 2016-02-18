import os
from setuptools import setup, find_packages
from setuptools.command.install import install as _install


# def copy_and_overwrite(from_path, to_path):
#     if os.path.exists(to_path):
#         shutil.rmtree(to_path)
#     shutil.copytree(from_path, to_path)
#
#
# def copy_dir(src, dest):
#     try:
#         copy_and_overwrite(src, dest)
#     except OSError as e:
#         # If the error was caused because the source wasn't a directory
#         if e.errno == errno.ENOTDIR:
#             shutil.copy(src, dest)
#         else:
#             print('Directory not copied. Error: %s' % e)


# def process_schema():
#     copy_dir("./eva_cttv_pipeline/resources/json_schema/src", "./eva_cttv_pipeline/resources/schema_copy")


# class CustomInstallCommand(install):
#
#     def process_schema(self):
#         json_schema_dir = utilities.get_resource_file("eva_cttv_pipeline", "resources/json_schema")
#         local_schema_dir = str(os.path.join(str(Path(json_schema_dir).parent), "schema_copy/"))
#
#         if not os.path.exists(local_schema_dir):
#             os.makedirs(local_schema_dir)
#         copy_dir(json_schema_dir, local_schema_dir)
#
#     def run(self):
#         self.process_schema()
#         install.run(self)


# class Install(_install):
#     def run(self):
#         # do stuff here
#         _install.run(self)


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
