import argparse
import io
import re
import tarfile
from pathlib import Path

import pandas

mutation_regex = re.compile(
    "MutationDetails \[id=MutationIdentifier \[location=Location \[clazz=(?P<clazz>.+), method=(?P<method>.+), methodDesc=(?P<desc>.+)\], indexes=\[(?P<idx>.+)\], mutator=(?P<mutator>.+)\], filename=(?P<filename>.+), block=(?P<block>.+), lineNumber=(?P<lineNumber>.+), description=(?P<description>.+), testsInOrder=\[(.*)\], isInFinallyBlock=(.+), poison=(.+)\]")
subjects = ['commons-csv', 'commons-collections', 'commons-math', 'commons-lang', 'commons-dbcp', 'commons-configuration', 'commons-imaging',
            'commons-net', 'commons-io', 'commons-compress']
testsuites = ['random', 'organic', 'evosuite', 'dynamosa']
re_subject = re.compile(fr"({'|'.join(subjects)})")
re_testsuite = re.compile(fr"({'|'.join(testsuites)})")


def get_subject_from_name(name):
    s_search = re_subject.search(name)
    return s_search.group(0)


def get_testsuite_from_name(name):
    ts_search = re_testsuite.search(str(name))
    return ts_search.group(0)


def get_mutants_from_tar(tar_path):
    tarf = tarfile.open(tar_path, "r")
    lines = []

    for file in tarf.getmembers():
        fd = io.TextIOWrapper(tarf.extractfile(file), "utf-8")
        with fd:
            m_data = fd.readline().strip()
        m_parsed = mutation_regex.search(m_data)
        fpath = Path(file.name)
        cut = f"{fpath.parent.parent.name}_{fpath.parent.name}"
        if m_parsed is None:
            print(m_data)
        lines.append({"mutatedClass": m_parsed.group("clazz"), "mutatedMethod": m_parsed.group("method"),
                      "methodDescription": m_parsed.group("desc"), "index": int(m_parsed.group("idx")),
                      "mutator": m_parsed.group("mutator"), "cut": cut, "lineNumber": int(m_parsed.group("lineNumber")),
                      "block": m_parsed.group("block"), "description": m_parsed.group("description")})
    return pandas.DataFrame(lines)


def get_total_mutants(path, out_path):
    res = []
    for file in path.iterdir():
        print(file.name)
        if not file.name.endswith("-mutations.tar.gz"):
            continue
        mutants = get_mutants_from_tar(file)
        mutants_count = len(mutants)
        subject = get_subject_from_name(file.name)
        testsuite = get_testsuite_from_name(file.name)
        res.append({"subject": subject, "testsuite": testsuite, "total_mutants": mutants_count})
    df = pandas.DataFrame(res)
    df = df.sort_values(by="subject")
    df.to_csv(out_path, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--input-dir", help="input dir with mutants.tar.gz files")
    parser.add_argument("-o", "--output-file", help="output csv file")
    args = parser.parse_args()
    get_total_mutants(Path(args.input_dir), Path(args.output_file))
