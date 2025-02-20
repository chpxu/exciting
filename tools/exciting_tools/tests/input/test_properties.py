""" Test properties.

NOTE:
All attribute tests should assert on the XML tree content's as the attribute
order is not preserved by the ElementTree.tostring method. Elements appear to
be fine.
"""
import re

import pytest

from excitingtools.exciting_dict_parsers.input_parser import parse_element_xml
from excitingtools.input.bandstructure import band_structure_input_from_ase_atoms_obj, \
    get_bandstructure_input_from_exciting_structure
from excitingtools.input.input_classes import ExcitingPropertiesInput, ExcitingGroundStateInput, ExcitingInputXML
from excitingtools.input.structure import ExcitingStructure


@pytest.fixture
def ase_ag():
    """If we cannot import ase we skip tests with this fixture.

    :returns: bulk Silver crystal
    """
    ase_build = pytest.importorskip('ase.build')
    return ase_build.bulk('Ag')


def test_bandstructure_properties_input():
    """Test giving the bandstructure input as a nested dict."""
    properties = {
        "bandstructure": {"plot1d": {"path": {
            "steps": 100,
            "point":
                [
                    {"coord": [1, 0, 0], "label": "Gamma"},
                    {"coord": [0.625, 0.375, 0], "label": "K"},
                    {"coord": [0.5, 0.5, 0], "label": "X", "breakafter": True},
                    {"coord": [0, 0, 0], "label": "Gamma"},
                    {"coord": [0.5, 0, 0], "label": "L"}]}}}}

    properties_input = ExcitingPropertiesInput(**properties)
    xml_string = properties_input.to_xml_str()

    assert properties == parse_element_xml(xml_string)


def test_bs_from_ase(ase_ag):
    """Test getting XML bandstructure input for elemental silver.

    We test first that the bandstructure path defined by the high symmetry
    points is correct. We also check the high symmetry point names are
    correct.
    """
    ase_ag.set_cell(ase_ag.cell * 3)
    bs = band_structure_input_from_ase_atoms_obj(ase_ag)
    properties_input = ExcitingPropertiesInput(bandstructure=bs)
    bs_xml_string = properties_input.to_xml_str()
    # Get the high symmetry path.
    all_labels = re.findall(r'label=\"([A-Z])', bs_xml_string)
    assert all_labels == ['G', 'X', 'W', 'K', 'G', 'L', 'U', 'W', 'L', 'K', 'U', 'X']
    # Ensure there is only one breakafter statement.
    assert len(re.findall(r'breakafter', bs_xml_string)) == 1

    bs_xml = properties_input.to_xml()
    points = bs_xml.find("bandstructure").find("plot1d").find("path").findall("point")
    assert len(points) == 12
    assert points[0].get("coord") == "0.0 0.0 0.0"
    assert points[1].get("coord") == "0.5 0.0 0.5"
    assert points[2].get("coord") == "0.5 0.25 0.75"
    assert points[3].get("coord") == "0.375 0.375 0.75"
    assert points[4].get("coord") == "0.0 0.0 0.0"
    assert points[5].get("coord") == "0.5 0.5 0.5"
    assert points[6].get("coord") == "0.625 0.25 0.625"
    assert points[7].get("coord") == "0.5 0.25 0.75"
    assert points[8].get("coord") == "0.5 0.5 0.5"
    assert points[9].get("coord") == "0.375 0.375 0.75"
    assert points[9].get("breakafter") == "true"
    assert points[10].get("coord") == "0.625 0.25 0.625"
    assert points[11].get("coord") == "0.5 0.0 0.5"


def test_bs_input(ase_ag):
    """Test writing exciting full input xml file with bandstructure property."""
    gs = ExcitingGroundStateInput(rgkmax=5.0)
    struct = ExcitingStructure(ase_ag)
    bs = band_structure_input_from_ase_atoms_obj(ase_ag)
    properties_input = ExcitingPropertiesInput(bandstructure=bs)
    input_xml = ExcitingInputXML(
        structure=struct, groundstate=gs, title="BS exciting",
        properties=properties_input).to_xml()

    assert input_xml.tag == 'input'
    assert input_xml.keys() == []

    subelements = list(input_xml)
    assert len(subelements) == 4
    title_xml = subelements[0]
    assert title_xml.tag == 'title'
    assert title_xml.keys() == []
    assert title_xml.text == 'BS exciting'

    properties_xml = subelements[3]
    assert properties_xml.tag == 'properties'
    assert properties_xml.keys() == []
    assert properties_xml.text == ' '
    assert properties_xml[0].tag == 'bandstructure'
    assert properties_xml[0][0].tag == 'plot1d'
    assert properties_xml[0][0][0].tag == 'path'

    first_coordinates = properties_xml[0][0][0][0].get('coord')
    first_label = properties_xml[0][0][0][0].get('label')
    assert first_coordinates == '0.0 0.0 0.0'
    assert first_label == 'G'


def test_get_bandstructure_input_from_exciting_structure(ase_ag):
    structure = ExcitingStructure(ase_ag)
    bs_xml = get_bandstructure_input_from_exciting_structure(structure).to_xml()
    points = bs_xml.find("plot1d").find("path").findall("point")
    assert len(points) == 12
    assert points[0].get("coord") == "0.0 0.0 0.0"
    assert points[1].get("coord") == "0.5 0.0 0.5"


def test_get_bandstructure_input_from_exciting_structure_stretch(ase_ag):
    structure = ExcitingStructure(ase_ag)
    structure.crystal_properties.stretch = [2, 1, 1]
    bs_xml = get_bandstructure_input_from_exciting_structure(structure).to_xml()
    points = bs_xml.find("plot1d").find("path").findall("point")
    assert len(points) == 13
    assert points[0].get("coord") == "0.0 0.0 0.0"
    assert points[1].get("coord") == "0.5 0.0 0.5"
    assert points[2].get("coord") == "0.125 -0.375 0.3125"
