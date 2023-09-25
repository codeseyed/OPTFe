import argparse
import logging
import zipfile
import pandas as pd
from cspy.apps.dma import generate_combined_name
from cspy.util.logging_config import FORMATS, DATEFMT
from cspy.chem import Molecule
from cspy.crystal.generate_crystal import CrystalGenerator
from cspy.deprecated.core.io.res import ResCrystalMapper, ResFileWriter
from cspy.deprecated.core.io.xyz import XYZFileReader
from cspy.deprecated.core.models import Crystal, Lattice
from cspy.deprecated.core.resources.crystallographic_space_groups import CrystallographicSpaceGroups
from cspy.deprecated.core.tasks.internal.landscape_generator.crystal import Convergence
from cspy.deprecated.core.tasks.internal.landscape_generator.crystal import CrystalLandscapeGenerator


LOG = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="App for running landscape generation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input-filenames",
        type=str,
        nargs="+",
        required=True,
        help="File containing input molecular structure (XYZ format)",
    )
    parser.add_argument(
        "-g",
        "--space-group",
        type=int,
        required=True,
        help="International tables space group number",
    )
    parser.add_argument(
        "-s",
        "--initial-seed",
        default=1,
        type=int,
        help="Initial seed for sobol sampling",
    )
    parser.add_argument(
        "-n",
        "--number-valid",
        type=int,
        default=1,
        help="Number of valid structures to generate",
    )
    parser.add_argument(
        "--clg-old",
        action="store_true",
        help="Use the old version of the CLG",
    )
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "WARNING", "ERROR", "INFO"),
        default="INFO",
        help="Set the level of detail output in logs",
    )
    args = parser.parse_args()
    logging.basicConfig(
        level=args.log_level, format=FORMATS[args.log_level], datefmt=DATEFMT
    )

    n_mols = len(args.input_filenames)
    LOG.info(
        "Generated crystals will have %d molecule%s in their asymmetric unit",
        n_mols, "" if n_mols < 2 else "s",
    )

    sg = args.space_group
    seed = args.initial_seed
    name = generate_combined_name(args.input_filenames)
    zip_name = f"{name}_{args.space_group}.zip"
    if args.clg_old:
        asymmetric_unit = [
            XYZFileReader(x).read() for x in args.input_filenames]
        model = Crystal(
            asymmetric_unit=asymmetric_unit,
            lattice=Lattice(1.0, 1.0, 1.0, 90, 90, 90),
            space_group=CrystallographicSpaceGroups.get(args.space_group)
        )
        convergence = Convergence(
            number_valid_structures=args.number_valid,
            max_trials=1000 * args.number_valid
        )

        clg = CrystalLandscapeGenerator(model=model, convergence=convergence)
        clg.run(start_index=args.initial_seed)
        results = clg.landscape.summary_dict()
        LOG.info("Summary:\n%s", pd.DataFrame(results).T)

        LOG.info("Writing structures to %s", zip_name)
        with zipfile.ZipFile(zip_name, "a") as zip_file:
            for obj in clg.landscape.accepted_objs:
                seed = obj.search_point.seed
                filename = f"{name}_{sg}_{seed}.res"
                titl = f"{name} cspy generated sg={sg} seed={seed}"
                file_content = ResCrystalMapper().map_to_content(obj.crystal)
                file_content = file_content._replace(TITL=titl)
                writer = ResFileWriter(
                    file_content=file_content,
                    output_location=""
                ).res_string_form
                zip_file.writestr(filename, writer)
    else:
        mols = [Molecule.load(x) for x in args.input_filenames]
        clg = CrystalGenerator(mols, sg)

        success = 0
        crystals = []
        while success != args.number_valid:
            crystal = clg.generate(seed)
            if crystal:
                success += 1
                crystals.append((crystal, seed))
            seed += 1
        rate = 100 * success / (seed - args.initial_seed)
        results = {"Results": {
            "Accepted": success,
            "Rejected": seed - success - args.initial_seed,
            "Success Rate": f"{rate:.06f}%",
        }}
        LOG.info("Summary:\n%s", pd.DataFrame(results).T)

        LOG.info("Writing structures to %s", zip_name)
        with zipfile.ZipFile(zip_name, "a") as zip_file:
            for crystal, seed in crystals:
                filename = f"{name}_{sg}_{seed}.res"
                titl = f"{name} cspy generated sg={sg} seed={seed}"
                shelx_string = crystal.to_shelx_string(titl=titl)
                zip_file.writestr(filename, shelx_string)


if __name__ == "__main__":
    main()
