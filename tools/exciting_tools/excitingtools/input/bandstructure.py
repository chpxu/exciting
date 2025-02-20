""" Generate bandstructure input from atoms. """
from __future__ import annotations

import re
from typing import List

import numpy as np
from ase import Atoms
from ase.cell import Cell

from excitingtools import ExcitingStructure
from excitingtools.input.input_classes import ExcitingBandStructureInput

_default_steps = 100


def band_structure_input_from_cell(cell: Cell,
                                   steps: int = _default_steps) -> ExcitingBandStructureInput:
    """ Get band path from lattice cell using ASE.

    :param cell: ase.cell.Cell object
    :param steps: number of steps for bandstructure calculation
    :return: the exciting band structure input
    """
    bandpath = cell.bandpath()
    points = []
    pattern = r'[A-Za-z,][0-9]*'
    for point in re.findall(pattern, bandpath.path):
        if point == ",":
            points[-1]["breakafter"] = True
        else:
            points.append({"coord": list(bandpath.special_points[point]), "label": point})

    return ExcitingBandStructureInput(plot1d={"path": {"steps": steps, "point": points}})


def band_structure_input_from_lattice(lattice_vectors: List[List[float]] | np.ndarray,
                                      steps: int = _default_steps) -> ExcitingBandStructureInput:
    """ Get band path from lattice vectors as array or list of lists using ASE. Lattice vectors
    in array needs to be stored row-wise.

    :param lattice_vectors: lattice
    :param steps: number of steps for bandstructure calculation
    :return: the exciting band structure input
    """
    return band_structure_input_from_cell(Cell(lattice_vectors), steps=steps)


def band_structure_input_from_ase_atoms_obj(atoms_obj: Atoms,
                                            steps: int = _default_steps) -> ExcitingBandStructureInput:
    """ Get band path from ase object.

    :param atoms_obj: ase atoms object
    :param steps: number of steps for bandstructure calculation
    :return: the exciting band structure input
    """
    return band_structure_input_from_cell(atoms_obj.cell, steps)


def get_bandstructure_input_from_exciting_structure(structure: ExcitingStructure,
                                                    steps: int = _default_steps) -> ExcitingBandStructureInput:
    """ Use ase bandpath to get bandstructure input for exciting using ASE.

    :param structure: exciting structure input object
    :param steps: number of steps for bandstructure calculation
    :return: bandstructure input for exciting
    """
    return band_structure_input_from_lattice(structure.get_lattice(), steps=steps)
