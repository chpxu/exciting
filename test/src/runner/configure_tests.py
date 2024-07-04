"""Parse test definitions from the yaml configuration file.
Where fields are missing, insert default settings from the defaults yaml configuration file.
"""
import copy
from typing import List, Dict, Union
import numpy as np
import os

from .set_inputs import input_files_for_tests
from ..exciting_settings.constants import settings

from ..tolerance.tol_classes import methods as valid_methods
from ..io.file_system import get_method_from_tolerance_file
from ..io.yaml_configuration import Group, ConfigurationDefaults, default_attributes


def get_method(test_name: str) -> str:
    """Get calculation method.

    method is defined as 'groundstate', 'gw', etc.

    One needs to know the method to choose the correct default `files_under_test`,
    when they are not specified in the config file.

    The problem:
      In most test cases, the file path defines the method. For example, `BSE/LiF`.
      However, properties are all grouped under one subdirectory of the same name.
      In this case, one can infer the method from the tolerance file.
      This is annoying to unit test though, as it requires one to interact with the file system. So every routine
      that uses `get_method` (i.e. all higher-level routines in configure_tests.py) will then require
      some file mocking when unit-tested.

    Hence, the use of this solution:
       a) check the test name string and
       b) fall back to the tolerance file iff it's a `properties` test.

    Better solutions one could implement, but require a design change are:
     a) Don't interact with the file system and instead specify every method in config input
        - This is the cleanest solution with minimal extra typing
        - May also restrict developers attempting to test multiple properties in one test.
     b) Split up the properties subdirectory according to the different tolerance files.

    :param str test_name: Test name, prepended by path.
    :return str method: Method name.
    """
    method = get_method_from_test_name(test_name)
    # This will interact with the file system
    if isinstance(method, KeyError):
        method = get_method_from_tolerance_file(test_name)
    return method


def get_method_from_test_name(test_name: str) -> Union[str, KeyError]:
    """Get method string from test name.

    Specific to the config file test naming convention.
    Expect test_name = 'method/test_case' or 'test_farm/method/test_case'
    Note, relying on strings and dealing with case is not the most robust approach.

    :param str test_name: Test name, prepended by path.
    :return Union[str, KeyError] method: Method name or KeyError
    """
    method = test_name.split('/')[-2].lower()
    if method not in valid_methods:
        # Properties are not differentiated in the directory structure of the test_farm, hence this requirement
        return KeyError()
    return method


def initialise_tests(test_names: List[str], default_inputs: dict, default_files_under_test: dict) -> dict:
    """Create tests dictionary, with all tests initialised using default properties.

    :param test_names: List of test names.
    :param default_inputs: Default inputs, defined per test_name.
    :param default_files_under_test: Default files under test, defined for all methods.
    :return dict tests: Initialised dictionary of test cases.
    """
    # defaults.inputs *currently* defined from searching test_farm.
    # Expect default.inputs = {'PBE-Al': ['input.xml', 'Ar.xml'], ..., 'test_case': ['input', 'files']}
    inconsistent_test_names = set(default_inputs.keys()) - set(test_names)
    if inconsistent_test_names:
        raise KeyError(f"Default input keys (from inspecting test_farm) inconsistent with test_names (from config "
                       f"file): "
                       f"{inconsistent_test_names}")

    tests = {}
    for name in test_names:
        method = get_method(name)
        default_dict = default_attributes._asdict()
        # Default files under test are method-specific, so must define here
        default_dict['files_under_test'] = default_files_under_test[method]
        # Test framework should be modified to make this a mandatory field in the yaml, which would remove this line
        default_dict['inputs'] = default_inputs[name]
        tests[name] = default_dict

    # Check all keys are assigned default values.
    random_test_defaults = next(iter(tests.values()))
    assert set(random_test_defaults.keys()) == set(ConfigurationDefaults._fields), \
        "Field present in namedtuple is missing from this initialisation."

    return tests


def check_tests_have_all_configuration_fields(tests: dict):
    """Check that each test case has all properties defined by the Configuration tuple.

    :param dict tests: Test dictionary, with initialised values.
    """
    expected_keys = set(ConfigurationDefaults._fields)

    for name, attributes in tests.items():
        test_keys = set(attributes.keys())
        if expected_keys != test_keys:
            raise KeyError(f'Test dictionary {name}, not initialised with all ConfigurationDefaults keys')


def assign_config_keys(initialised_tests: dict, config_data: dict) -> dict:
    """ Add fields specified in config data to the initialised test dictionary.

    This could be reduced to:

    ```python
    for name in tests.keys():
        tests[name].update(config_data[name])
    ```
    if parsing:

    ```yaml
    groundstate/LDA_PW-noncollinear-Fe:

    ```
    returned {}, rather than None.

    :param initialised_tests: All test cases, with fields/attributes initialised to defaults.
    :param config_data: All test cases defined in yaml configuration file.
    :return tests: Test cases defined in yaml configuration file, with missing fields assigned default values.
    """
    tests = copy.deepcopy(initialised_tests)
    for name in tests.keys():
        specified_attributes = config_data[name]
        # For entries in tests_config.yml with no explicit options, their value parses as None, not {}
        if isinstance(specified_attributes, dict):
            tests[name].update(config_data[name])
    return tests


def configure_all_tests(config_data: dict, default_files_under_test: dict) -> dict:
    """Configure all tests in the suite.

    Initialise tests with default settings, then replace any properties that have been explicitly specified in the
    test config file.

    :param config_data: Test config file data.
    :param default_files_under_test: Default files under test (each key is a method and values are files under test).
    :return dict tests: Tests dict, containing all test cases (even those that shouldn't run) and their properties.
    """
    # Must prepend with 'test_farm' directory, as this is not specified in the naming convention of the config file
    config_data = {os.path.join(settings.test_farm, name): value for name, value in config_data.items()}

    test_names = [name for name in config_data.keys()]
    default_inputs = input_files_for_tests(test_names, subdirectory='ref')
    initialised_tests = initialise_tests(test_names, default_inputs, default_files_under_test)
    tests = assign_config_keys(initialised_tests, config_data)

    return tests


def find_tests_with_attribute(tests: dict, attribute_key: str, default_value) -> dict:
    """
    Create a dictionary of tests with a specific attribute, as specified by attribute_key.

    For example, attribute_key = failing tests
    return failing_tests = {'test1': failing_builds_i, 'testi': failing_builds_1}

    :param dict tests: Tests.
    :param str attribute_key: Key of attribute.
    :param default_value: Default value of attribute, as specified in 'default'.
    :return dict: subset of tests.
    """
    test_names = []
    specific_attributes = []

    for name, attributes in tests.items():
        test_names.append(name)
        specific_attributes.append(attributes[attribute_key])

    not_equal = [entry != default_value for entry in specific_attributes]
    indices = np.where(not_equal)[0]
    test_names_with_attribute = [test_names[i] for i in indices]

    return {name: tests[name] for name in test_names_with_attribute}


def sort_tests_by_group(tests: dict) -> Dict[Group, List[str]]:
    """
    Sort tests according to their `group` tag/attribute.

    For example:

    {SLOW_TESTS: ['test_farm/XANES/TiO2'],
     GROUP_NAME: ['test_farm/method/test_name', 'test_farm/method/another_test']
    }

    :param dict tests: Tests.
    :return Dict[Group, List[str]] tests_by_group: Tests sorted by group.
    """
    # Initialise entries for all possible groups
    tests_by_group = {group: [] for group in Group}

    # Log tests by 'group' attribute
    for name, attributes in tests.items():
        group = attributes['group']
        tests_by_group[group].append(name)

    return tests_by_group


def find_tests_to_skip_by_group(all_tests: dict, group_execution: dict) -> set:
    """
    List tests with groups that are not assigned to run.

    :param dict all_tests: All test cases.
    :param dict group_execution: All groups, and whether they should run.
     For example {SLOW_TESTS: False, GROUP_NAME: True, ...}

    :return dict tests_to_skip: Tests to skip.
    """
    tests_by_group = sort_tests_by_group(all_tests)

    tests_to_run = set()
    for group_str, to_run in group_execution.items():
        if to_run:
            group = Group[group_str]
            tests_to_run.update(tests_by_group[group])

    tests_to_skip = set(all_tests.keys()) - tests_to_run
    return tests_to_skip
