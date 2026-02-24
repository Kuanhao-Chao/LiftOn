import sys

# ANSI colour codes (disabled on non-TTY so log files stay clean)
_USE_COLOUR = sys.stderr.isatty()
_RESET  = "\033[0m"  if _USE_COLOUR else ""
_YELLOW = "\033[33m" if _USE_COLOUR else ""
_RED    = "\033[31m" if _USE_COLOUR else ""
_GREEN  = "\033[32m" if _USE_COLOUR else ""
_BOLD   = "\033[1m"  if _USE_COLOUR else ""


def log(*argv, debug=False):
    """
    Log a message to stderr.

    Parameters
    ----------
    *argv : any
        Message parts (joined by space, same as print).
    debug : bool
        Only print when debug=True.
    """
    if debug:
        print(*argv, file=sys.stderr)


def log_info(*argv):
    """Always print an informational message to stderr (no debug flag needed)."""
    print(*argv, file=sys.stderr)


def log_warning(*argv):
    """Print a yellow WARNING message to stderr."""
    print(f"{_YELLOW}{_BOLD}[WARNING]{_RESET}", *argv, file=sys.stderr)


def log_error(*argv):
    """Print a red ERROR message to stderr."""
    print(f"{_RED}{_BOLD}[ERROR]{_RESET}", *argv, file=sys.stderr)


def log_success(*argv):
    """Print a green SUCCESS message to stderr."""
    print(f"{_GREEN}✓{_RESET}", *argv, file=sys.stderr)


def log_section(title: str, body_lines: list, kind: str = "info") -> None:
    """
    Print a boxed section header to stderr.

    Parameters
    ----------
    title : str
        Section title shown in the header bar.
    body_lines : list[str]
        Lines to print inside the box.
    kind : str
        One of 'info', 'warning', 'error', 'success' — controls the colour.
    """
    width = 70
    colour_map = {
        "info":    "",
        "warning": _YELLOW,
        "error":   _RED,
        "success": _GREEN,
    }
    colour = colour_map.get(kind, "")
    reset = _RESET if colour else ""

    top    = f"╔{'═' * (width - 2)}╗"
    sep    = f"╠{'═' * (width - 2)}╣"
    bottom = f"╚{'═' * (width - 2)}╝"

    def _row(text=""):
        pad = max(0, width - 4 - len(text))
        return f"║  {text}{' ' * pad}  ║"

    print(colour + top, file=sys.stderr)
    print(_row(title), file=sys.stderr)
    if body_lines:
        print(sep + reset, file=sys.stderr)
        for line in body_lines:
            # Wrap long lines
            for chunk in _wrap_line(line, width - 6):
                print(colour + _row(chunk) + reset, file=sys.stderr)
    print(colour + bottom + reset, file=sys.stderr)


def _wrap_line(text: str, max_width: int) -> list:
    """Wrap a single line to max_width characters."""
    if len(text) <= max_width:
        return [text]
    chunks = []
    while len(text) > max_width:
        chunks.append(text[:max_width])
        text = text[max_width:]
    if text:
        chunks.append(text)
    return chunks
