"""
a script to parse Juge logs
"""
import argparse
import csv
import re
import sys
from pathlib import Path

from tqdm import tqdm

re_line = re.compile(r'^(?P<test_method>\[?.*\]?)\((?P<test_class>[a-zA-Z0-9._]+)\) (?P<kill>(PASS|FAIL|IGNORED|INVALID)) ?(?P<err_tail>\[?.*)$')
re_err = re.compile(r'^\[(?P<err>.+?)(, (?P<line>\d+))?\]$')


# (, )?(?P<line>\d*)\]?


def parse_line(line):
    res = re_line.search(line)
    if not res:
        return None
    return res['test_class'], res['test_method'], res['kill'], res['err_tail']


def parse_err(err_tail):
    # TODO: sanitize non-unicode characters
    res = re_err.search(err_tail)
    if not res:
        return None
    return res['err'], res['line']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(epilog='Example of use: python3 collect_data.py -p commons-lang -g organic -r 01')
    parser.add_argument("-p", "--proj", help="project", required=True)
    parser.add_argument("-g", "--generator", help="generator", required=True)
    parser.add_argument("-r", "--run", help="run", required=True)
    parser.add_argument("-d", "--dir", help="root folder, not mandatory, defaults to current dir ", default=".")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    proj = args.proj
    generator = args.generator
    run = args.run
    root = args.dir
    # e.g., root/csv/random-tests/03/
    wd = Path(root).absolute() / proj / f"{generator}-tests" / run
    mutant_coverage_file = Path(root).absolute() / proj / f"{generator}-{run}-mutant-coverage.csv"  # no use now, assuming we already have it
    # "GENERATOR,RUN,CUT,ALL_MUTANTS,COVERED_MUTANTS" > ${MUTANT_COVERAGE}
    kill_map_file = Path(root).absolute() / proj / f"{generator}-{run}-killmatrix.csv"
    fieldnames = ["GENERATOR", "RUN", "CUT", "MUTANT", "TEST_CLASS", "TEST_METHOD", "KILL", "EXCEPTION", "LINE_NUMBER"]
    with kill_map_file.open('w') as csv_fd:
        writer = csv.DictWriter(csv_fd, fieldnames=fieldnames, delimiter='|', quotechar='"', quoting=csv.QUOTE_MINIMAL, doublequote=True)
        writer.writeheader()
        for f_name in tqdm(wd.rglob("test-ese*.txt")):
            mutant = f_name.parents[0].name
            cut_dir = f_name.parents[1]
            cut = f_name.parents[1].name
            with f_name.open() as fd:
                try:
                    while True:
                        line = fd.readline()
                        if not line:
                            break
                        line = line.strip()
                        parsed_line = parse_line(line)
                        if not parsed_line:
                            continue
                        test_class, test_method, kill, err_tail = parsed_line
                        if kill != 'PASS':
                            while not err_tail.endswith(']'):
                                err_tail += " " + fd.readline().strip()
                            exception, ln = parse_err(err_tail)
                        else:
                            exception = ""
                            ln = ""
                        writer.writerow({'GENERATOR': generator,
                                         'RUN': run,
                                         'CUT': cut,
                                         'MUTANT': f"{cut}_{mutant}",
                                         'TEST_CLASS': test_class,
                                         'TEST_METHOD': test_method,
                                         'KILL': kill,
                                         'EXCEPTION': exception,
                                         'LINE_NUMBER': ln
                                         })
                except:
                    print("Smth went wrong")
