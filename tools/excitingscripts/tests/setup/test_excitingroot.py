import os, pytest, subprocess
from unittest import mock

from excitingtools.utils.test_utils import MockFile
from excitingscripts.setup.excitingroot import set_exciting_root

@pytest.fixture
def input_mock(tmp_path) -> MockFile:
    """ Mock file containing $EXCITINGROOT
    """
    string = "Testing if $EXCITINGROOT is replaced by the correct string.\n"
    input_file = tmp_path / "input.xml"
    input_file.write_text(string)
    return MockFile(input_file, string)

def test_set_exciting_root(input_mock, tmp_path):
   # May need to check I've not missed anything... treat this as pseudo-code

   exciting_root = "some/path"
   reference = input_mock.string.replace("$EXCITINGROOT", exciting_root)

   output = tmp_path / "output.xml"
   output.write_text("")

   set_exciting_root(input_mock.full_path, output, excitingroot=exciting_root)
   assert output.read_text() == reference

@mock.patch.dict(os.environ, {"EXCITINGROOT": "some/path"})
def test_cli_using_flags(input_mock):
    
    reference = input_mock.string.replace("$EXCITINGROOT", "some/path")

    # setting filepath explicitly

    subprocess.run(["python3", "-m", "excitingscripts.setup.excitingroot", "-i", input_mock.full_path, "-o", input_mock.full_path])

    with open(input_mock.full_path, "r") as f:
        altered_data = f.read()

    assert altered_data == reference, "Did not replace `$EXCTIINGROOT` correctly using options."

@mock.patch.dict(os.environ, {"EXCITINGROOT": "some/path"})
def test_cli_using_defaults(input_mock, tmp_path):

    reference = input_mock.string.replace("$EXCITINGROOT", "some/path")

    cwd_ = os.getcwd()

    os.chdir(tmp_path.as_posix())

    try:
        subprocess.run(["python3", "-m", "excitingscripts.setup.excitingroot"])
    finally:
        os.chdir(cwd_)

    with open(input_mock.full_path, "r") as f:
        altered_data = f.read()

    assert altered_data == reference, "Did not replace `$EXCTIINGROOT` correctly using defaults."
