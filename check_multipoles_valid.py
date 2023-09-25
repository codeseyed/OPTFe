import sys
from cspy import Molecule
from cspy.apps.dma import generate_combined_name
from cspy.chem.multipole import DistributedMultipoles
from cspy.crystal import Crystal
from cspy.potentials import available_potentials
from cspy.util.path import Path
from cspy.configuration import CONFIG



def check_multipoles_valid(xyz_files, axis, charges, multipoles):
    potential = "w99rev_6311"
    c = DistributedMultipoles.from_dma_file(charges)
    # if len(c.molecules) > 0:
        # print("Succesfully read charges file: %s mols", len(c.molecules))
    m = DistributedMultipoles.from_dma_file(multipoles)
    # if len(m.molecules) > 0:
        # print(
        #     "Succesfully read higher order multipoles file: %s mols", len(m.molecules)
        # )
    axis_contents = Path(axis).read_text()
    success = True
    for rank, multipoles in (("charges", c), ("multipoles", m)):
        crys = Crystal.from_xyz_files(xyz_files, titl="tmp")
        potential_type = available_potentials[potential][0]
        crys.neighcrys_setup(potential_type=potential_type, file_content=axis_contents)
        #print("Mapping: %s", rank)
        reduction = 0.1 if potential_type == "W" else None
        mapping = crys.map_multipoles(multipoles, foreshorten_hydrogens=reduction)
        for i, (j, rmsd, inv) in enumerate(mapping):
            print(
                f"{xyz_files[i]}: {rank} idx={j} {rmsd}: {'ok' if rmsd < 1e-4 else 'failed'}")
            if rmsd > 1e-4:
                success = False
    return success

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "xyz_files",
        type=str,
        nargs="+",
        help="Xyz files containing molecules for generation",
    )
    parser.add_argument(
        "-c",
        "--charges",
        type=str,
        default=CONFIG.get("csp.charges_file"),
        help="Rank0 multipole file",
    )
    parser.add_argument(
        "-m",
        "--multipoles",
        type=str,
        default=CONFIG.get("csp.multipoles_file"),
        help="RankN multipole file",
    )
    parser.add_argument(
        "-a",
        "--axis",
        type=str,
        default=None,
        help="Axis filename for structure minimization",
    )
    args = parser.parse_args()
    default_name = generate_combined_name(args.xyz_files)
    if args.axis is None:
        args.axis = default_name + ".mols"
        print("No axis file provided, trying %s", args.axis)
    if args.charges is None:
        args.charges = default_name + "_rank0.dma"
        print("No charges file provided, trying %s", args.charges)
    if args.multipoles is None:
        args.multipoles = default_name + ".dma"
        print("No multipole file provided, trying %s", args.multipoles)
    if not check_multipoles_valid(args.xyz_files, args.axis, args.charges, args.multipoles):
        print("No good matches for given multipoles! Exiting...")
        sys.exit(1)
    else:
        print("Multipoles looking good")

if __name__ == "__main__":
    main()
