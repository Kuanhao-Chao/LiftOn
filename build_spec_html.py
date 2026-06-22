#!/usr/bin/env python3
"""
build_spec_html.py — render lifton_technical_spec.md into a single self-contained,
offline HTML file (embedded CSS, sticky collapsible TOC sidebar, Pygments-highlighted
code, and inlined Mermaid diagrams). No CDN / external resources.

Usage:
    python build_spec_html.py \
        --md  lifton_technical_spec.md \
        --out lifton_technical_spec.html \
        --mermaid spec_assets/mermaid.min.js
"""
import argparse
import html
import re
import sys
from pathlib import Path

import markdown
from pygments.formatters import HtmlFormatter

MERMAID_PLACEHOLDER = "@@LIFTON_MERMAID_BLOCK_{}@@"
_FENCE_MERMAID = re.compile(r"^```mermaid[ \t]*\n(.*?)\n```[ \t]*$", re.DOTALL | re.MULTILINE)


def extract_mermaid(md_text):
    """Pull ```mermaid fences out before Markdown processing so the diagram
    text is never HTML-escaped. Returns (text_with_placeholders, [raw_diagrams])."""
    blocks = []

    def _sub(m):
        idx = len(blocks)
        blocks.append(m.group(1))
        # Surround with blank lines so Markdown treats the token as its own paragraph.
        return f"\n\n{MERMAID_PLACEHOLDER.format(idx)}\n\n"

    return _FENCE_MERMAID.sub(_sub, md_text), blocks


def reinsert_mermaid(html_body, blocks):
    """Replace each placeholder paragraph with a <pre class="mermaid"> element
    holding the raw (un-escaped) diagram source mermaid.js expects."""
    for idx, raw in enumerate(blocks):
        token = MERMAID_PLACEHOLDER.format(idx)
        pre = f'<pre class="mermaid">\n{raw}\n</pre>'
        # Markdown will have wrapped the bare token in <p>...</p>.
        html_body = html_body.replace(f"<p>{token}</p>", pre)
        html_body = html_body.replace(token, pre)  # belt-and-suspenders
    return html_body


CSS = r"""
:root{
  --bg:#ffffff; --fg:#1f2328; --muted:#656d76; --accent:#0969da; --accent2:#0a7f5b;
  --border:#d0d7de; --code-bg:#f6f8fa; --sidebar-bg:#f6f8fa; --sidebar-w:320px;
  --table-stripe:#f6f8fa; --warn:#9a6700; --warn-bg:#fff8c5;
}
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{margin:0;color:var(--fg);background:var(--bg);
  font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;}
a{color:var(--accent);text-decoration:none}
a:hover{text-decoration:underline}
.layout{display:flex;align-items:flex-start}
/* Sidebar */
#sidebar{position:sticky;top:0;height:100vh;width:var(--sidebar-w);min-width:var(--sidebar-w);
  overflow-y:auto;background:var(--sidebar-bg);border-right:1px solid var(--border);padding:18px 14px 60px;}
#sidebar .toc-title{font-weight:700;font-size:15px;letter-spacing:.02em;text-transform:uppercase;
  color:var(--muted);margin:4px 6px 12px}
#sidebar .toc ul{list-style:none;margin:0;padding-left:14px}
#sidebar .toc>ul{padding-left:0}
#sidebar .toc li{margin:1px 0}
#sidebar .toc a{display:block;padding:3px 8px;border-radius:6px;color:var(--fg);font-size:13.5px;line-height:1.35}
#sidebar .toc a:hover{background:#eaeef2;text-decoration:none}
#sidebar .toc a.active{background:#ddf4ff;color:var(--accent);font-weight:600}
.toc-toggle{cursor:pointer;user-select:none;color:var(--muted);display:inline-block;width:14px;text-align:center}
.collapsed > ul{display:none}
.searchbox{width:100%;padding:7px 10px;margin:0 0 12px;border:1px solid var(--border);border-radius:8px;font-size:13px}
/* Content */
#content{flex:1;max-width:980px;margin:0 auto;padding:36px 48px 120px;min-width:0}
#content h1{font-size:30px;border-bottom:2px solid var(--border);padding-bottom:.3em;margin-top:0}
#content h2{font-size:24px;border-bottom:1px solid var(--border);padding-bottom:.25em;margin-top:2.2em}
#content h3{font-size:19px;margin-top:1.8em}
#content h4{font-size:16px;margin-top:1.4em;color:#374151}
#content h2,#content h3,#content h4{scroll-margin-top:14px}
#content p,#content li{font-size:15px}
code,kbd{font-family:ui-monospace,SFMono-Regular,"SF Mono",Menlo,Consolas,monospace;font-size:13px}
:not(pre)>code{background:var(--code-bg);padding:.15em .4em;border-radius:6px;border:1px solid var(--border)}
pre{background:var(--code-bg);border:1px solid var(--border);border-radius:8px;padding:14px 16px;overflow:auto}
pre code{background:none;border:none;padding:0}
.codehilite{background:var(--code-bg);border:1px solid var(--border);border-radius:8px;overflow:auto;margin:14px 0}
.codehilite pre{border:none;margin:0;background:none}
table{border-collapse:collapse;width:100%;margin:16px 0;display:block;overflow-x:auto}
th,td{border:1px solid var(--border);padding:7px 11px;text-align:left;font-size:13.5px;vertical-align:top}
th{background:var(--code-bg);font-weight:600}
tr:nth-child(2n) td{background:var(--table-stripe)}
blockquote{margin:14px 0;padding:.4em 1em;color:var(--muted);border-left:4px solid var(--border);background:var(--code-bg)}
hr{border:none;border-top:1px solid var(--border);margin:2em 0}
.mermaid{background:#fff;border:1px solid var(--border);border-radius:8px;padding:14px;margin:18px 0;text-align:center;overflow:auto}
.banner{background:var(--warn-bg);border:1px solid #d4a72c;color:var(--warn);border-radius:8px;padding:8px 14px;margin:0 0 20px;font-size:13px}
/* Responsive */
@media (max-width:980px){
  #sidebar{position:fixed;left:0;z-index:30;transform:translateX(-100%);transition:transform .2s;box-shadow:2px 0 12px rgba(0,0,0,.15)}
  #sidebar.open{transform:translateX(0)}
  #content{padding:24px 18px 100px}
  #menubtn{display:block}
}
#menubtn{display:none;position:fixed;top:12px;left:12px;z-index:40;background:var(--accent);color:#fff;border:none;
  border-radius:8px;padding:8px 12px;font-size:14px;cursor:pointer}
/* Print */
@media print{
  #sidebar,#menubtn{display:none}
  #content{max-width:none;padding:0}
  pre,table,.mermaid{page-break-inside:avoid}
  a{color:inherit}
}
"""

PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
{css}
{pygments_css}
</style>
</head>
<body>
<button id="menubtn" onclick="document.getElementById('sidebar').classList.toggle('open')">&#9776; Contents</button>
<div class="layout">
  <nav id="sidebar">
    <div class="toc-title">LiftOn Spec</div>
    <input class="searchbox" id="tocsearch" type="text" placeholder="Filter sections…" oninput="filterToc(this.value)">
    <div class="toc">{toc}</div>
  </nav>
  <main id="content">
{body}
  </main>
</div>
<script>
{mermaid_js}
</script>
<script>
try {{ mermaid.initialize({{ startOnLoad: true, securityLevel: 'loose', theme: 'neutral', flowchart: {{ htmlLabels: true, useMaxWidth: true }} }}); }}
catch (e) {{ console.error('mermaid init failed', e); }}

// Collapsible TOC: add toggles to every li that has a nested ul.
(function(){{
  var toc = document.querySelector('#sidebar .toc');
  if(!toc) return;
  toc.querySelectorAll('li').forEach(function(li){{
    var sub = li.querySelector(':scope > ul');
    if(sub){{
      var t = document.createElement('span');
      t.className = 'toc-toggle'; t.textContent = '▾';
      t.onclick = function(e){{ e.stopPropagation(); li.classList.toggle('collapsed');
        t.textContent = li.classList.contains('collapsed') ? '▸' : '▾'; }};
      li.insertBefore(t, li.firstChild);
    }}
  }});
}})();

// Scroll-spy: highlight the TOC entry for the heading currently in view.
(function(){{
  var links = Array.prototype.slice.call(document.querySelectorAll('#sidebar .toc a'));
  var map = {{}};
  links.forEach(function(a){{ var id = decodeURIComponent((a.getAttribute('href')||'').slice(1)); if(id) map[id]=a; }});
  var heads = Object.keys(map).map(function(id){{ return document.getElementById(id); }}).filter(Boolean);
  function onScroll(){{
    var cur=null, y=window.scrollY+90;
    for(var i=0;i<heads.length;i++){{ if(heads[i].offsetTop<=y) cur=heads[i]; else break; }}
    links.forEach(function(a){{ a.classList.remove('active'); }});
    if(cur && map[cur.id]){{ map[cur.id].classList.add('active'); }}
  }}
  document.addEventListener('scroll', onScroll, {{passive:true}}); onScroll();
}})();

function filterToc(q){{
  q=(q||'').toLowerCase();
  document.querySelectorAll('#sidebar .toc li').forEach(function(li){{
    var a=li.querySelector(':scope > a');
    var txt=a?a.textContent.toLowerCase():'';
    var self=txt.indexOf(q)>=0;
    var kid=!!li.querySelector('ul a');  // keep parents of matches visible (cheap)
    li.style.display=(q===''||self||kid)?'':'none';
  }});
}}
</script>
</body>
</html>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", default="lifton_technical_spec.md")
    ap.add_argument("--out", default="lifton_technical_spec.html")
    ap.add_argument("--mermaid", default="spec_assets/mermaid.min.js")
    ap.add_argument("--title", default="LiftOn — Technical Specification")
    args = ap.parse_args()

    md_path = Path(args.md)
    if not md_path.is_file():
        sys.exit(f"ERROR: markdown source not found: {md_path}")
    raw = md_path.read_text(encoding="utf-8")

    text, mermaid_blocks = extract_mermaid(raw)

    md = markdown.Markdown(
        extensions=["extra", "toc", "codehilite", "sane_lists", "admonition"],
        extension_configs={
            "toc": {"permalink": False, "toc_depth": "2-4"},
            "codehilite": {"guess_lang": False, "css_class": "codehilite"},
        },
    )
    body = md.convert(text)
    body = reinsert_mermaid(body, mermaid_blocks)
    toc = getattr(md, "toc", "") or ""

    pyg_css = HtmlFormatter(style="default").get_style_defs(".codehilite")

    mermaid_path = Path(args.mermaid)
    if not mermaid_path.is_file():
        sys.exit(f"ERROR: mermaid bundle not found: {mermaid_path} (run the download step)")
    mermaid_js = mermaid_path.read_text(encoding="utf-8")
    # Defensive: never let a literal </script> in the bundle close our inline tag.
    mermaid_js = mermaid_js.replace("</script>", "<\\/script>")

    out = PAGE.format(
        title=html.escape(args.title),
        css=CSS,
        pygments_css=pyg_css,
        toc=toc,
        body=body,
        mermaid_js=mermaid_js,
    )
    Path(args.out).write_text(out, encoding="utf-8")
    n_mer = out.count('<pre class="mermaid">')
    print(f"Wrote {args.out}  ({len(out):,} bytes, {len(mermaid_blocks)} mermaid blocks placed, "
          f"{n_mer} rendered in DOM, {body.count('<table')} tables)")
    if n_mer != len(mermaid_blocks):
        print(f"WARNING: placed {len(mermaid_blocks)} mermaid blocks but DOM has {n_mer} <pre class=mermaid>")


if __name__ == "__main__":
    main()
