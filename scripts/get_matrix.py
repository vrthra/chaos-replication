import scipy.sparse
import pandas as pd
import argparse
from pathlib import Path
default_catch_types = ['classes', 'methods']
testsuites = ['random', 'organic', 'dynamosa']
subjects = ['commons-csv', 'commons-collections', 'commons-math', 'commons-lang', 'commons-dbcp', 'commons-configuration', 'commons-imaging',
            'commons-net', 'commons-io', 'commons-compress']


def get_killmatrix(input_dir, output_dir, subject, testsuite, sample_type):
    name = f"{subject}_{testsuite}_{sample_type}-killmatrix"
    cols_path = input_dir / subject / (name+".cols")
    rows_path = input_dir / subject / (name+".rows")
    matrix_path = input_dir / subject / (name+".npz")
    with cols_path.open() as fd:
        cols = [c.strip() for c in fd.readlines()]
    with rows_path.open() as fd:
        rows = [r.strip() for r in fd.readlines()]
    matrix = scipy.sparse.load_npz(str(matrix_path)).todense()
    df = pd.DataFrame(data=matrix, index=rows, columns=cols)
    out_path = output_dir / f"{subject}_{testsuite}_{sample_type}-killmatrix.csv"
    df.to_csv(str(out_path))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-dir", default="../data/killmatrix", help="input dir with npz killmatrix files")
    parser.add_argument("-o", "--output-dir", help="output dir")
    parser.add_argument("-st", "--sample-type", choices=default_catch_types, help="output type")
    parser.add_argument("-s", "--subject", help="subject")
    parser.add_argument("-t", "--testsuite", choices=testsuites, help="testsuite")
    args = parser.parse_args()
    get_killmatrix(Path(args.input_dir), Path(args.output_dir), args.subject, args.testsuite, args.sample_type)
