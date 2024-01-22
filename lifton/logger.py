import sys 

def log(*argv, debug=False):
    """
        This function logs the message to stderr.

        Parameters:
        - argv: message
        - debug: debug mode
    """
    if debug:
        print(*argv, file=sys.stderr)
