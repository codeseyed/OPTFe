from cspy.crystal import Crystal
import argparse
from multiprocessing import Pool
from functools import partial

def niggli(input_file, p1):
    structure = Crystal.load(input_file)
    file_ext = input_file[-4:]
    name = input_file.split(file_ext)[0]
    niggli_structure = structure.standardized()

    if p1:
        niggli_P1 = niggli_structure.as_P1()
        niggli_P1.save('{}_p1.cif'.format(name))
    else:
        niggli_structure.save('{}.cif'.format(name))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "structures",
        type=str,
        nargs="+",
        help="CIF or RES files containing input structures."
             )

    parser.add_argument(
        "-j",
        "--num_threads",
        type=int,
        help="Number of parallel threads to use when doing the comparisons.",
        default=1
        )

    parser.add_argument(
        "--p1",
        action='store_true',
        help="Saves in P1 space group.",
        default=False,
        )

    args = parser.parse_args()
    structures = [item for item in args.structures]
    
    niggli_run = partial(niggli, p1=args.p1)

    with Pool(int(args.num_threads)) as pool:
        pool.map(niggli_run, structures)

if __name__ == "__main__":
    main()
