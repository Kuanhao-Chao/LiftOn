# ---------------------------------------------------------------------------
# Author: Kuan-Hao Chao <kuanhao.chao@gmail.com>
# ---------------------------------------------------------------------------
"""Pure-Python fallback parser. Used as the correctness oracle for the Rust
implementation and as a fallback when the native extension is unavailable.
"""

from .parser import parse_file, parse_bytes, detect_dialect

__all__ = ["parse_file", "parse_bytes", "detect_dialect"]
