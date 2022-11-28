import argparse

import subprocess as sp
from pathlib import Path

from copy import copy

FPM = "fpm"

def make_icaf_config(args):
    if (args.cores is None) or (args.program is None):
        return
    Path("icaf_config").write_text(
        f"-genvall -n {args.cores} ./{args.program}\n"
    )

PROFILES = {
    "gfortran_release": {
        "compiler": "gfortran",
        "flags": "-Ofast -Wall -fno-stack-arrays -fcoarray=single",
    },
    "gfortran_debug": {
        "compiler": "gfortran",
        "flags": "-O0 -g -Wall -fcheck=all -fimplicit-none "
                 "-fcoarray=single",
    },
    "caf_release": {
        "compiler": "caf",
        "flags": "-Ofast -Wall -fno-stack-arrays",
    },
    "caf_debug": {
        "compiler": "caf",
        "flags": "-O0 -g -Wall -fbounds-check -fimplicit-none ",
    },
    "ifort_release": {
        "compiler": "ifort",
        "flags": "-traceback -Ofast -warn all,nounused "
                 "-coarray=distributed -heap-arrays "
                f"-coarray-config-file={f'{str(Path.cwd())}/icaf_config'}",
        "extra_func": make_icaf_config,
    },
}

PROFILES["gfortran_gprof"] = copy(PROFILES["gfortran_release"])
PROFILES["gfortran_gprof"]["flags"] += " -pg"


def clear():
    for directory in ['.', "src", "test", "example"]:
        [path.unlink() for path in Path(directory).glob("*.mod")]
    sp.run(["rm", "-rf", "build"])


def compile(profile):
    sp.run([FPM, "build", "--tests",
            "--compiler", profile['compiler'],
            "--flag", profile['flags']],
           check=True)


parser = argparse.ArgumentParser(description="Compile program.")
parser.add_argument("--profile", type=str, default=None,
                    help="compile profile, specified in global dictionary")
parser.add_argument("--clear", action="store_true",
                    help="whether to clear build directory and mod files")
parser.add_argument("--cores", action="store", default=None,
                    help="USED WITH INTEL PROFILES:"
                         "create 'icaf_config' file with the specified "
                         "number of cores")
parser.add_argument("--program", action="store", default='',
                    help="USED WITH INTEL PROFILES:"
                         "create 'icaf_config' file with the specified "
                         "program to run")

if __name__ == "__main__":

    args = parser.parse_args()

    if args.clear: clear()

    if args.profile is not None:
        if args.profile not in PROFILES:
            raise ValueError(f"Invalid profile {args.profile}, either specify "
                              "an existing profile or create one in the "
                              "PROFILES dictionary")
        profile = PROFILES[args.profile]
        compile(profile)
        if "extra_func" in profile: profile["extra_func"](args)
