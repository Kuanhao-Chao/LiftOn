from glob import glob
import lifton.liftoff 
from lifton.liftoff import liftoff_main
import os
import pytest

def test_yeast(tmp_path):
    print("tmp_path: ", tmp_path)
    os.makedirs(tmp_path, exist_ok=True)
    # cleanup
    for tfile in glob('lifton/liftoff/tests/*.mmi') + glob('lifton/liftoff/tests/*.fai') + glob('lifton/liftoff/tests/*_db'):
        print("tfile: ", tfile)
        try:
            os.unlink(tfile)
        except OSError:
            pass

    # inputs
    asmbl = 'lifton/liftoff/tests/GCA_000146045.2_R64_genomic.fna.gz'
    annot = 'lifton/liftoff/tests/GCA_000146045.2_R64_genomic.gff.gz'
    target = 'lifton/liftoff/tests/GCA_000146045.2_R64_genomic.fna.gz'

    # outputs


    output = f'{tmp_path}/GCA_000146045.2_R64_to_GCA_000146045.2_basic.gff3'
    unmapped = f'{tmp_path}/unmapped_features_GCA_000146045.2_R64_basic.txt'
    tempdir = f'{tmp_path}/sandbox'
    expout = 'lifton/liftoff/tests/GCA_000146045.2_R64_to_GCA_000146045.2_R64_expected_basic.gff'
    expunmapped = 'lifton/liftoff/tests/unmapped_features_GCA_000146045.2_R64_expected_basic.txt'

    # run the program
    args = ['-g', annot, '-o', output, '-u', unmapped, '-dir', tempdir, target, asmbl]
    liftoff_main.main(args)

    # verify the output
    with open(output, 'r') as fh1, open(expout, 'r') as fh2:
        fh1_lines = fh1.readlines()
        fh2_lines = fh2.readlines()
        for i in range(len(fh1_lines)):
                if fh1_lines[0][0] != "#":
                    observed_output = fh1_lines[i].strip().split("\t")[0:8]
                    expected_output = fh2_lines[i].strip().split("\t")[0:8]
                    assert observed_output == expected_output
    with open(unmapped, 'r') as fh1, open(expunmapped, 'r') as fh2:
        observed_output = fh1.read().strip()
        expected_output = fh2.read().strip()
        assert observed_output == expected_output
