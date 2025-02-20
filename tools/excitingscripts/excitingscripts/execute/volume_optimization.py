import os
import pathlib
from argparse import ArgumentParser
from os.path import join
from typing import Union

import numpy as np
from excitingscripts.execute.single import run_exciting
from excitingtools import parser_chooser


def execute_volume_optimization(input_file: Union[str, pathlib.Path], number_volume_values: int,
                                root_directory=os.getcwd(), excitingroot=os.getenv("EXCITINGROOT")) -> np.ndarray:
    """Execute a series of exciting calculations with different volumes obtained by varying the lattice constant.

    :param input_file: Input file.
    :param number_volume_values: Number of volume values for which structures are generated by varying the lattice
    constant.
    :param root_directory: Root directory.
    :param excitingroot: Environment variable string.
    :returns: NumPy array containing total energy values and corresponding volume values.
    """
    total_energy = []
    volume_values = []

    for i_v in range(0, number_volume_values):

        run_exciting(f"{root_directory}/volume-{i_v + 1}", excitingroot)

        results = parser_chooser(join(os.getcwd(), f"{root_directory}/volume-{i_v + 1}/INFO.OUT"))
        max_scf = max([int(i) for i in results["scl"].keys()])
        converged_results = results["scl"][str(max_scf)]
        total_energy.append(converged_results["Total energy"])
        volume_values.append(results["initialization"]["Unit cell volume"])

    energy_volume_data = np.vstack((volume_values, total_energy)).T

    return energy_volume_data


def main() -> None:
    parser = ArgumentParser(description="""Execute a series of exciting calculations with different volumes obtained by
                                        varying the lattice constant.""")

    parser.add_argument("--input-file", "-i",
                        type=Union[str, pathlib.Path],
                        default=["input.xml"],
                        nargs=1,
                        dest="infile",
                        help="name of the input file")

    parser.add_argument("number_volume_values",
                        type=int,
                        nargs=1,
                        help="number of volume values for which structures are generated")

    parser.add_argument("--root-directory", "-r",
                        type=Union[str, pathlib.Path],
                        default=[os.getcwd()],
                        nargs=1,
                        dest="root_directory",
                        help="root path for files that are created by this script")

    args = parser.parse_args()

    energy_volume_data = execute_volume_optimization(args.infile[0], args.number_volume_values[0],
                                                     args.root_directory[0])

    with open(f"{args.root_directory[0]}/energy-vs-volume", "w") as f:
        np.savetxt(f, energy_volume_data)

if __name__ == "__main__":
    main()
