import sys 

def log(*argv, debug=False):
    if debug:
        print(*argv, file=sys.stderr)
