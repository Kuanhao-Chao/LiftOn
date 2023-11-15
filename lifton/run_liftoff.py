import subprocess
import os, copy
from lifton.liftoff import liftoff_main
from lifton.liftoff.tests import test_basic, test_advanced

# def check_liftoff_installed():
#     installed = False
#     try:
#         installed = liftoff.__version__ >= "1.6.3"
#     except:
#         pass
#     return installed


def run_liftoff(args):
    print(">> run_liftoff")
    liftoff_args = copy.deepcopy(args)
    liftoff_outdir = os.path.dirname(args.output) + "/liftoff/"
    os.makedirs(liftoff_outdir, exist_ok=True)
    liftoff_args.output = liftoff_outdir + "/liftoff.gff3"
    liftoff_main.run_all_liftoff_steps(liftoff_args)

    # test_basic.test_yeast(liftoff_outdir + "test_basic/")
    # test_advanced.test_yeast(liftoff_outdir + "test_advance/")
