""" Specification and functions for operating on the YAML files:

* tests_config.yml
* defaults_config.yml
"""
from collections import namedtuple
import enum
import os
from pathlib import Path
from typing import Dict, Type, Union, Optional, Callable
import yaml

from ..runner.profile import CompilerBuild


# Properties/attributes of a test case.
ConfigurationDefaults = namedtuple('ConfigurationDefaults',
                                   ['files_under_test',
                                    'inputs',
                                    'repeat',
                                    'failing_builds',
                                    'group',
                                    'comments',
                                    'depends_on',
                                    'cmd_line_args']
                                   )

# PyYaml tags, for use with yaml constructors (Must start with `!`)
group_tag = u'!Group'
build_tag = u'!Build'


def dynamically_generate_test_group_enum_class() -> Type[enum.Enum]:
    """ Dynamically generate the Group enum class from YAML file attributes.

    NOTE. defaults_config.yml gets parsed twice by the test suite, which
    could be avoided. However, the file is small so the overhead is negligible.

    :return Group: Enum class containing the enums defined by `group_execution`
    in 'defaults_config.yml'.
    """
    # Get the full path of defaults_config.yml.
    this_file_path = os.path.dirname(os.path.realpath(__file__))
    defaults_config_file = Path(*Path(this_file_path).parts[:-2]) / 'defaults_config.yml'

    try:
        with open(defaults_config_file, 'r') as stream:
            data = yaml.safe_load(stream)
    except FileNotFoundError:
        raise FileNotFoundError(f'Cannot find YAML defaults file {defaults_config_file}')

    group_names = [name for name in data['group_execution'].keys()]
    return enum.Enum(value='Group', names=group_names)


# Group enum class
Group = dynamically_generate_test_group_enum_class()


# Default attributes for each test case:
# * files_under_test. Defaults are method specific, so define later
#   - Making this a mandatory field in the yaml would simplify this.
# * inputs. Default are scraped from ref/ directory of each test
#   - Test framework should definitely be modified to make this a mandatory field in the yaml
# * depends_on: default assumes a single dependency or a chain of dependencies, not multiple dependencies
#   - One would need to modify it to a list.
default_attributes = ConfigurationDefaults(files_under_test=[],
                                           inputs=[],
                                           repeat=False,
                                           failing_builds=[],
                                           depends_on=None,
                                           group=Group.NONE,  # Cannot resolve until runtime
                                           comments='',
                                           cmd_line_args=''
                                           )


def enum_group_constructor(loader: yaml.SafeLoader, node: yaml.nodes.ScalarNode) -> Group:
    """ Custom YAML Constructor: Parse group string entry to enum.

    :return Group enum.
    """
    enum_str = loader.construct_scalar(node)
    try:
        return Group[enum_str]
    except KeyError:
        raise KeyError(f"Invalid Group enum in yaml file: {enum_str}")


def enum_compiler_build_constructor(loader: yaml.SafeLoader, node: yaml.nodes.ScalarNode) -> CompilerBuild:
    """ Custom YAML Constructor: Parse compiler-build string to CompilerBuild object.

    For example:

    A YAML entry string 'intel_mpiandsmp' is returned to CompilerBuild('intel', 'mpiandsmp')
    such that an instance of CompilerBuild is present in the parsed dictionary.

    :return CompilerBuild instance.
    """
    compiler_build_str = loader.construct_scalar(node)
    try:
        compiler, build = compiler_build_str.split('_')[0:2]
    except ValueError:
        raise ValueError(f'String representation of CompilerBuild enum is invalid: {compiler_build_str}. \n'
                         f'See src/runner/profile.py or config file documentation for valid strings.')
    return CompilerBuild(compiler, build)


def robust_yaml_arg(load_func: Callable) -> Callable:
    """Decorator for YAML loader to handle file or YAML-formatted string as an argument.
    Note, there may be cleaner ways to implement this.

    :param load_func: Yaml loader.
    :return modified_function: Modified Yaml loader.
    """
    def modified_function(file_or_yaml_str: Union[str, Path], custom_constructor: Optional[dict] = None):

        # 1. Must be a file
        if isinstance(file_or_yaml_str, Path):

            if not file_or_yaml_str.is_file():
                raise FileNotFoundError(f'File not found: {file_or_yaml_str}')

            with open(file_or_yaml_str, "r") as fid:
                yaml_str = fid.read()
            return load_func(yaml_str, custom_constructor)

        # 2. Can be a file name or yaml string
        elif isinstance(file_or_yaml_str, str):

            # File
            if os.path.isfile(file_or_yaml_str):
                with open(file_or_yaml_str, "r") as fid:
                    yaml_str = fid.read()
                return load_func(yaml_str, custom_constructor)

            # `Assume` string contains yaml and pass on
            return load_func(file_or_yaml_str, custom_constructor)

        # Argument is wrong type
        else:
            raise ValueError(f'First argument is not a file or string: {file_or_yaml_str}')

    return modified_function


# Signature of custom_constructor dict
CustomConType = Dict[str, Callable]


@robust_yaml_arg
def yaml_load(yaml_str, custom_constructor: Optional[CustomConType] = None) -> dict:
    """ Parse YAML file.

    Expect custom_constructor of the form:
    {u'!Group': enum_group_constructor}

    :param yaml_str: File name or yaml-formatted string
    :param custom_constructor: Optional custom constructors applied to parsed data.
    :return config: Dict of parsed data.
    """
    # Extend loader to parse the data through custom constructors
    if custom_constructor is None:
        custom_constructor = {}

    custom_loader = yaml.SafeLoader
    for tag, constructor in custom_constructor.items():
        custom_loader.add_constructor(tag, constructor)

    # Parse yaml
    try:
        config = yaml.load(yaml_str, Loader=custom_loader)
    except yaml.YAMLError:
        raise yaml.YAMLError(f'Invalid formatting in YAML file: {yaml_str}')

    return config
