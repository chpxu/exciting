"""Reading, writing and interacting with the file system.
"""
import os
import shutil
from typing import List, Set, Generator

from ..exciting_settings.constants import settings
from ..io.tolerances import list_tolerance_files_in_directory
from ..tolerance.tol_classes import tol_file_to_method


def create_run_dir(path_to_test_case: str, run_dir: str):
    """
    Create run directory for a test case with the name test_name in test_farm.
    If a run directory already exists, it will be deleted.

    :param str path_to_test_case:  path to the test case
    :param str run_dir:    name of the run directory
    """

    # check if the test case exists
    try:
        os.path.isdir(path_to_test_case)
    except NotADirectoryError:
        raise NotADirectoryError('%s does not exist.', path_to_test_case)

    if os.path.isdir(os.path.join(path_to_test_case, run_dir)):
        try:
            shutil.rmtree(os.path.join(path_to_test_case, run_dir))
        except OSError:
            raise OSError('%s could not be removed.' % os.path.join(path_to_test_case, run_dir))

    os.mkdir(os.path.join(path_to_test_case, run_dir))


def flatten_directory(path: str):
    """
    Flatten the file structure by factor 1 for the directory at path.
    :param path: path to directory which gets flattened
    """
    try:
        os.path.isdir(path)
    except Exception:
        raise NotADirectoryError('%s does not exist.' % path)

    dirs = next(os.walk(path))[1]

    for dir in dirs:
        files = next(os.walk(os.path.join(path, dir)))[2]
        for file in files:
            try:
                shutil.move(os.path.join(path, dir, file), path)
            except shutil.Error:
                os.remove(os.path.join(path, file))
                shutil.move(os.path.join(path, dir, file), path)

        shutil.rmtree(os.path.join(path, dir))


def files_to_remove(files: List[str], remove_list: List[str]) -> List[str]:
    """
    Return the intersection between a list of files and a list of files that shall be removed.

    :param List[str] files: List of files
    :param List[str] remove_list: List of files to be removed
    :return List[str] to_remove: List of files to be removed
    """
    to_remove = []
    for remove in remove_list:
        for file in files:
            if remove in file:
                to_remove.append(file)

    return to_remove


def remove_ignored_files(path: str, ignored_output: list):
    """
    Remove files that are ignored for the tests in a directory.

    :param str path:              Path to the directory
    :param list ignored_output:    Files that are ignored by the tests.
    """
    try:
        files = next(os.walk(path))[2]
    except NotADirectoryError:
        raise NotADirectoryError('%s is not a directory.' % path)

    for file in files_to_remove(files, ignored_output):
        os.remove(os.path.join(path, file))


def copy_calculation_inputs(source: str, destination: str, input_files: List[str]):
    """
    Copy input files for a calculation from source to destination.
    Performs checking to ensure all inputs exist.

    :param str source: Path to the exciting calculation where to copy the input files from
    :param str destination: Path to the directory where the input files shall be copied to.
    :param List[str] input_files: Input files required to perform a calculation
    """
    inputs_present = input_files_in_directory(source, input_files)
    missing_files = set(input_files) - set(inputs_present)

    if missing_files:
        raise FileNotFoundError(f'Required input file(s) are missing from the test farm directory: {missing_files}')

    if not os.path.isdir(destination):
        raise NotADirectoryError('Target directory does not exist:', destination)

    for file in input_files:
        shutil.copyfile(os.path.join(source, file), os.path.join(destination, file))


def input_files_in_directory(directory: str, files: List[str]) -> List[str]:
    """
    Return files which are present in files argument and directory.

    :param str directory: Path to directory to check
    :param List[str] files: List of possible files
    :return List[str]:  Files present in directory
    """
    files_in_directory = next(os.walk(directory))[2]
    if not files_in_directory:
        raise FileNotFoundError(f'Directory "{directory}" contains no files')

    species_files = set(files_in_directory) & set(files)
    return list(species_files)


def get_test_directories(test_farm_root: str, basename=False) -> Generator[str, None, None]:
    """Get test directories from the test farm root.

    Test farm has the directory structure:
     test_farm/method/test_directory

    :param str test_farm_root: Test farm root directory name
    :param bool basename: Only return test directory base names
    :return List[str] test_dirs: List of test directories, given relative to test_farm_root
    if basename is false.
    """
    method_dirs = next(os.walk(test_farm_root))[1]

    if basename:
        for method_dir in method_dirs:
            method_dir = os.path.join(test_farm_root, method_dir)
            yield from next(os.walk(method_dir))[1]

    else:
        for method_dir in method_dirs:
            method_dir = os.path.join(test_farm_root, method_dir)
            test_dirs = next(os.walk(method_dir))[1]
            yield from (os.path.join(method_dir, t) for t in test_dirs)


def get_all_test_cases(test_farm: str) -> Set[str]:
    """Get all test cases present in the test_farm directory by file system inspection.

    This *assumes* a specific test farm subdirectory structure.

    :param str test_farm: Test farm directory.
    :return Set[str] test_cases: All test cases present in test_farm, with names prepended by `test_farm/method`.
    """
    methods = next(os.walk(test_farm))[1]

    test_cases = []
    for method in methods:
        directory = os.path.join(test_farm, method)
        test_names = next(os.walk(directory))[1]
        test_cases += [os.path.join(test_farm, method, name) for name in test_names]

    return set(test_cases)


def get_method_from_tolerance_file(test_name: str, subdirectory=settings.ref_dir) -> str:
    """Get method from tolerance file inspection.

    Uses the method specified in the {tolerance file : method map}

    :param str test_name: Test name, prepended by path.
    :param optional, str subdirectory: Subdirectory of test case, in which to look for tolerance file.
    :return str method: Method name.
    """
    directory = os.path.join(test_name, subdirectory)
    tolerance_files = list_tolerance_files_in_directory(directory)

    if len(tolerance_files) != 1:
        raise ValueError(f"Found {len(tolerance_files)} tolerance files in directory: {directory}")

    tol_file = tolerance_files[0]
    try:
        method = tol_file_to_method[tol_file]
    except KeyError:
        method_and_testname = "/".join(s for s in test_name.split('/')[-2:])
        valid_tol_files_str = "\n".join(m for m in tol_file_to_method.keys())
        raise KeyError(f'Invalid tolerance file, {tol_file}, for {method_and_testname}. \n'
                       f'Method must correspond to a tolerance:\n{valid_tol_files_str}'
                       )

    return method
