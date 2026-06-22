## 6. External Tool Integration

LiftOn drives two external aligners and one in-process protein aligner:

1. **Liftoff** (a complete vendored fork, `lifton/liftoff/`) — DNA-level lift-over, invoked **as a Python library** (not a subprocess). Liftoff itself shells out to **minimap2** (or drives `mappy` in-process under `--native`).
2. **miniprot** — protein-to-genome aligner, invoked as a **subprocess** (`miniprot --gff-only ...`), or through a `pyminiprot`-shaped facade under `--native`.
3. **parasail** — Smith/Waterman-Gotoh-style alignment library, used **in-process** (PyO3-style C extension) for both protein-vs-protein and DNA-vs-DNA pairwise alignment inside the per-locus loop.

This part of the spec covers the LiftOn-side wrappers that build these invocations (§6.1 Liftoff wrapper, §6.2 miniprot subprocess + parsing, §6.3 native-binding facades, §6.4 consolidated command/parameter table). The internals of the vendored Liftoff fork are covered separately in §6.5.

Two key entry points sit *above* the wrappers and decide whether to run the tool at all:

- `lifton_utils.exec_liftoff(outdir, ref_db, args)` (`lifton/lifton_utils.py:97`) — short-circuits when the user supplied a pre-built Liftoff GFF3 via `-L`/`--liftoff`.
- `lifton_utils.exec_miniprot(outdir, args, tgt_genome, ref_proteins_file)` (`lifton/lifton_utils.py:117`) — short-circuits when the user supplied a pre-built miniprot GFF3 via `-M`/`--miniprot`.

Both follow the same idiom:

```python
annotation = args.liftoff            # or args.miniprot
if annotation is None or not os.path.exists(annotation):
    annotation = run_liftoff.run_liftoff(...)   # or run_miniprot.run_miniprot(...)
return annotation
```

**Gotcha (short-circuit semantics):** the guard is `args.X is None OR not os.path.exists(args.X)`. A user-supplied path that does **not** exist on disk is silently ignored and the tool is run anyway — there is no error for "you pointed me at a missing `-L` file". The path string the user passed is also the path the runner writes to only by coincidence; `run_liftoff`/`run_miniprot` compute their *own* output paths under `<outdir>/liftoff/` and `<outdir>/miniprot/` and ignore `args.liftoff`/`args.miniprot` entirely. The short-circuit is the *only* place the user path is consulted.

`exec_miniprot` additionally calls `check_miniprot_installed()` (`lifton/lifton_utils.py:130`) unconditionally — even on the short-circuit path — but discards its return value (a leftover; it does not gate anything).

---

### 6.1 `run_liftoff.run_liftoff` — the Liftoff library wrapper

`run_liftoff.run_liftoff(output_dir, ref_db, args)` (`lifton/run_liftoff.py:10`) runs the vendored Liftoff and returns **either** a filesystem path (legacy/default) **or** an in-memory GFF3 `bytes` blob (Phase 8, `--inmemory-liftoff`). The return type is therefore `str | bytes`; downstream Step 5 ingest (`lifton/lifton.py:453`) accepts both shapes transparently because `annotation.Annotation` can build a FeatureDB from a path or from a string blob.

#### 6.1.1 Inputs

| Param | Type | Meaning |
|---|---|---|
| `output_dir` | `str` | Base intermediate-files dir; Liftoff artefacts go under `<output_dir>liftoff/`. Caller passes a trailing-slash-terminated dir, so string concatenation (`output_dir + "liftoff/"`) yields a valid path. |
| `ref_db` | `annotation.Annotation` (wraps `gffutils.FeatureDB` or gffbase) | Reference annotation Liftoff lifts from. Passed straight into the vendored Liftoff entry points. |
| `args` | `argparse.Namespace` | Must carry `polish` (bool); may carry `inmemory_liftoff` (bool, default `False` via `getattr`). The whole namespace is **deep-copied** before mutation (see below). |

#### 6.1.2 Algorithm (numbered)

1. **Deep-copy args** (`run_liftoff.py:38`): `liftoff_args = copy.deepcopy(args)`. All Liftoff-specific path mutations happen on the copy so the caller's `args` is never corrupted (Liftoff reads many fields like `args.output`, `args.u`, `args.directory`).
2. **Compute output dir & make it** (`:39-40`): `liftoff_outdir = output_dir + "liftoff/"`; `os.makedirs(liftoff_outdir, exist_ok=True)`.
3. **Set Liftoff output paths on the copy** (`:41-43`):
   - `liftoff_annotation = liftoff_outdir + "liftoff.gff3"`
   - `liftoff_args.output = liftoff_annotation`
   - `liftoff_args.u = liftoff_outdir + "unmapped_features.txt"` (unmapped-features report; **always written to disk** in both modes).
4. **Decide mode** (`:45`): `inmemory = bool(getattr(args, "inmemory_liftoff", False))`. Note this reads the *original* `args`, not `liftoff_args` — equivalent here since the field is not mutated.
5. **Save & raise recursion limit** (`:56-57`):
   - `_orig_recursion_limit = sys.getrecursionlimit()`
   - `sys.setrecursionlimit(max(_orig_recursion_limit, 10000))`
   - Rationale (Phase 16 Tier 3): real reference annotations (NCBI RefSeq, GENCODE) drive Liftoff's recursive feature-hierarchy traversal past Python's default 1000-frame limit. The `max(orig, 10000)` form **never lowers** an already-higher limit. 10000 frames stays well inside an 8 MB Linux thread stack (~13K frames).
6. **`try:` — run Liftoff** (`:58-67`):
   - **In-memory branch** (`inmemory == True`, `:59-65`):
     1. `from lifton.liftoff import inmemory_emitter` (lazy import).
     2. `lifted_feature_list, feature_db, _ref_parent_order, _unmapped = liftoff_main.run_all_liftoff_steps_inmemory(liftoff_args, ref_db)` — runs all Liftoff steps but **returns the lifted-feature list in RAM** rather than serialising it to `liftoff.gff3`.
     3. `gff_bytes = inmemory_emitter.lifted_features_to_gff3_bytes(lifted_feature_list, liftoff_args, feature_db)` — serialises to a `bytes` blob using **the same serialisation helpers** as the on-disk `write_new_gff`, guaranteeing byte-identical output to the legacy file.
   - **On-disk branch** (`inmemory == False`, `:66-67`): `liftoff_main.run_all_liftoff_steps(liftoff_args, ref_db)` — runs all Liftoff steps and writes `liftoff.gff3` to disk. No return value is captured; the path is reconstructed later.
7. **`except Exception as e:` — fatal-error handling** (`:68-81`):
   1. `import traceback` (lazy).
   2. `logger.log_error(f"Liftoff encountered a fatal error during native execution: {e}")` — the short message.
   3. `logger.log_error("Full Python traceback ...\n" + traceback.format_exc().rstrip("\n"))` (Phase 16 Tier 2) — the **complete formatted traceback**, so the deepest LiftOn / vendored-Liftoff frame is visible (critical for diagnosing `RecursionError`, whose `str()` is only "maximum recursion depth exceeded …").
   4. `logger.log_error("LiftOn cannot proceed without a valid Liftoff baseline annotation.")`.
   5. `sys.exit(1)` — Liftoff failure is **fatal** for the whole run (unlike miniprot, which degrades gracefully).
8. **`finally:` — restore recursion limit** (`:82-83`): `sys.setrecursionlimit(_orig_recursion_limit)`. This restores the original even on the `sys.exit(1)` path (the `SystemExit` raised by `sys.exit` propagates *through* the `finally`), so the process-global recursion limit is never left perturbed.
9. **Polish suffix** (`:85-86`): `if args.polish: liftoff_annotation += "_polished"`. When Liftoff polishing is enabled it writes to `liftoff.gff3_polished`; the returned path is adjusted to match. **Gotcha:** this suffix logic runs *after* the run block and *only* mutates the path string — it has no effect in the in-memory branch (which returns `gff_bytes`, computed before any polish step).
10. **Return** (`:90-92`):
    - `if inmemory: return gff_bytes` (a `bytes` blob).
    - else `return liftoff_annotation` (a `str` path, possibly `..._polished`).

**Gotcha (byte-identity contract):** the in-memory branch is required to be byte-for-byte identical to the on-disk file. This is enforced by the `--inmemory-liftoff` 4-cell matrix (`tests/test_liftoff_inmemory.py`) and the 24-cell native matrix. The identity holds because `inmemory_emitter.lifted_features_to_gff3_bytes` *shares* the GFF3 line-formatting code with `lifton.liftoff.write_new_gff`. Any change to Liftoff's serialisation must touch both paths or break the gate.

**Gotcha (recursion-limit save/restore must wrap the whole call):** the `try/finally` is the *only* thing guaranteeing the global limit is restored. If a future refactor moves the `run_all_liftoff_steps` call outside the `try`, an exception would leave the process at 10000 frames permanently, perturbing every later component (including parasail callbacks and the writer).

#### 6.1.3 Liftoff-side consumers (context, not part of the wrapper)

After `run_liftoff` returns, Step 5 (`lifton/lifton.py:453`) re-ingests the result into a FeatureDB `l_feature_db` via `annotation.Annotation(liftoff_annotation, infer_genes=False, infer_transcripts=False, merge_strategy=args.merge_strategy, id_spec=None, force=args.force, verbose=args.verbose, auto_convert_gtf=False)`. The same module (`run_liftoff.py`) also hosts the per-locus Liftoff→protein processing helpers (`initialize_lifton_gene` `:95`, `lifton_add_trans_exon_cds` `:116`, `process_liftoff_with_protein` `:142`, `process_liftoff` `:186`); those belong to the Step-7 pipeline, not the external-tool boundary, and are specified elsewhere. The only external-boundary facts worth pinning here:

- `lifton_add_trans_exon_cds` (`:135`) reads CDS children with `featuretype=('CDS', 'stop_codon')` — Liftoff emits both, and LiftOn treats `stop_codon` rows as part of the CDS chain.
- `process_liftoff` carries a `_visited` cycle-detection set (`:214-224`) that raises `LiftOnInputError` on a circular `Parent=` chain *before* it can trip a `RecursionError` — a defence layered under the §6.1.2-step-5 recursion-limit bump.

---

### 6.2 `run_miniprot.run_miniprot` — the miniprot subprocess wrapper

`run_miniprot.run_miniprot(outdir, args, tgt_genome, ref_proteins_file)` (`lifton/run_miniprot.py:80`) runs miniprot and returns `str | bytes | None`:

- **path** (`str`) — default file-write branch, the legacy Phase 5 contract.
- **bytes** — `--stream` branch, stdout piped to RAM, no disk write.
- **`None`** — any failure mode; the pipeline then continues **Liftoff-only** (miniprot is *not* fatal, unlike Liftoff).

#### 6.2.1 Inputs

| Param | Type | Meaning |
|---|---|---|
| `outdir` | `str` | Base intermediate dir (trailing slash); miniprot artefacts under `<outdir>miniprot/`. |
| `args` | `argparse.Namespace` | Must carry `mp_options` (str); may carry `stream` (bool, default `False`) and `native` (bool, default `False`), both read via `getattr`. |
| `tgt_genome` | `str` | Path to target-genome FASTA (the database miniprot maps proteins **into**). |
| `ref_proteins_file` | `str` | Path to reference-proteins FASTA (the queries miniprot maps **from**). |

#### 6.2.2 Command construction

```python
miniprot_path = "miniprot"            # bare name; resolved via PATH
command = (
    [miniprot_path, "--gff-only", tgt_genome, ref_proteins_file]
    + [opt for opt in args.mp_options.split(" ") if opt]
)
```
(`run_miniprot.py:120-125`)

| Token | Source | Meaning |
|---|---|---|
| `miniprot` | hard-coded `:120` | binary name; **no configurable path** — must be on `$PATH`. |
| `--gff-only` | hard-coded `:123` | emit **GFF3 only** to stdout (no PAF/alignment block, no `##` headers except what miniprot itself emits). |
| `tgt_genome` | arg `:123` | positional 1: the genome being indexed/mapped into. |
| `ref_proteins_file` | arg `:123` | positional 2: the protein FASTA queries. |
| `*args.mp_options.split(" ")` (non-empty tokens) | `:124` | user extra flags. **Default `mp_options` is the empty string** (`lifton/lifton.py:70`), so by default **no extra flags** are appended. The `if opt` filter drops empty tokens so a leading/trailing/double space cannot inject `''` argv entries. |

**Gotcha (argv splitting is naive `.split(" ")`):** `mp_options` is split on single spaces only. A value containing tabs, or a quoted multi-word argument, will be mis-tokenised. This is the same naive splitting Liftoff uses for `mm2_options`. Threading is **not** added here — miniprot runs single-threaded unless the user passes `-t N` inside `mp_options`.

Both `args.mp_options` and the full command are `print()`-ed to stdout for visibility (`:121`, `:126`).

#### 6.2.3 Three execution branches

The branch is selected by `(stream_mode, native_mode) = (bool(getattr(args,"stream",False)), bool(getattr(args,"native",False)))` (`run_miniprot.py:128-129`):

**Branch A — `stream_mode AND native_mode`** (`:131-158`): route through the native facade.
1. `from lifton.native_bindings import MiniprotIndex`.
2. Construct `idx = MiniprotIndex(tgt_genome, mp_options=args.mp_options, miniprot_path="miniprot", ref_proteins_path=ref_proteins_file)`.
3. `bundle = idx.align_all()` (runs miniprot **once**, caches result).
4. Set `stdout_bytes = bundle.raw_bytes`, `stderr_text = ""`, `return_code = 0`, `output_size = len(stdout_bytes)`.
5. On `RuntimeError` from the facade (its way of re-raising miniprot's non-zero exit / `ERROR`): print `"[LiftOn] miniprot (native facade) failed: {exc}"` to stderr and `return None`.
   - **Gotcha:** today the facade *still* shells out to the subprocess underneath, so this path is **byte-identical** to Branch B (asserted by `tests/test_native_bindings.py::TestMiniprotFacadeStreamingParity`). It exists so a future real `pyminiprot` PyO3 binding can be swapped in with no caller change.

**Branch B — `stream_mode` only** (`:159-174`): stdout→RAM with bounded chunked drain.
1. `proc = subprocess.Popen(command, stdout=PIPE, stderr=PIPE, bufsize=1 << 20)` (1 MiB pipe buffer).
2. `stdout_bytes, stderr_bytes, return_code = _drain_stream_chunks(proc, chunk_size=1 << 16)` (64 KiB chunks).
3. `stderr_text = stderr_bytes.decode("utf-8", errors="replace") if stderr_bytes else ""`.
4. `output_size = len(stdout_bytes) if stdout_bytes else 0`.

**Branch C — default file-write** (`:175-187`, legacy Phase 5 baseline):
1. `with open(miniprot_output, "w") as fw:` where `miniprot_output = outdir + "miniprot/" + "miniprot.gff3"` (`:117-119`).
2. `proc = subprocess.run(command, stdout=fw, stderr=subprocess.PIPE, text=True)` — stdout streamed straight to the file handle; stderr captured for scanning.
3. `stderr_text = proc.stderr or ""`; `return_code = proc.returncode`.
4. `output_size = os.path.getsize(miniprot_output) if os.path.exists(miniprot_output) else 0`.

#### 6.2.4 `_drain_stream_chunks` (bounded, double-allocation-avoiding drain)

`_drain_stream_chunks(proc, *, chunk_size=65536)` (`run_miniprot.py:6`) replaces `proc.communicate()` (Phase 15c / V3.10). `communicate()` materialises `bytes(stdout) + bytes(stderr)` and holds both simultaneously — an 8-figure peak-RSS pathology for multi-GB miniprot runs.

Algorithm:
1. `import threading` (lazy).
2. `stdout_chunks: list[bytes] = []`, `stderr_chunks: list[bytes] = []`.
3. Define inner `_drain(stream, sink)`:
   - if `stream is None`: return.
   - loop: `buf = stream.read(chunk_size)`; `if not buf: return`; else `sink.append(buf)`.
   - `finally:` `stream.close()` (best-effort, swallows exceptions).
4. Start two **daemon** threads, one per stream: `t_out` on `(proc.stdout, stdout_chunks)`, `t_err` on `(proc.stderr, stderr_chunks)`.
5. `rc = proc.wait()`; then `t_out.join()`; `t_err.join()`.
6. `return b"".join(stdout_chunks), b"".join(stderr_chunks), rc`.

**Gotcha (peak-RAM contract):** during collection, peak buffer is bounded by `chunk_size` per stream (64 KiB each as called from Branch B). The **final** `b"".join(...)` is unavoidable — `gffbase.create_db(from_string=True)` needs the full blob — but the doubled-allocation interim of `communicate()` is eliminated. Two threads are needed to avoid pipe-buffer deadlock (if you read stdout to EOF first while stderr's 64 KiB pipe fills, miniprot blocks writing stderr).

**Gotcha (chunk_size mismatch):** Branch B passes `chunk_size=1 << 16` (65536) explicitly; this equals the function default but is stated independently — do not assume they track each other.

#### 6.2.5 Post-run: stderr echo and the three failure modes

After whichever branch ran (`:189-222`):

1. **Echo stderr** (`:190-191`): `if stderr_text: print(stderr_text, end="", file=sys.stderr)` — surface miniprot's own log lines verbatim.
2. **Failure mode 1 — non-zero exit** (`:194-201`): `if return_code != 0:` print a "miniprot exited with code N … continue Liftoff-only" message and `return None`.
3. **Failure mode 2 — `ERROR` on stderr with exit 0** (`:204-212`): `if stderr_text and "ERROR" in stderr_text.upper():` print a message and `return None`. This guards a **known miniprot quirk**: exit code 0 *and* an error logged. The check is **case-insensitive** (`.upper()`), so it matches `Error`, `error`, `ERROR`.
4. **Failure mode 3 — empty/absent output** (`:215-222`): `if output_size == 0:` print a message and `return None`.

All three failure messages share the contract "Miniprot output will be skipped — LiftOn will continue using Liftoff results only."

5. **Outer `except Exception as exc:`** (`:224-231`): any unanticipated exception (e.g. `FileNotFoundError` from a missing `miniprot` binary in Branch C's `subprocess.run`) prints `"[LiftOn] miniprot failed unexpectedly: {exc}"` and returns `None`. So miniprot being uninstalled degrades to Liftoff-only rather than crashing.

6. **Return** (`:233-236`): `if stream_mode: return stdout_bytes` (covers Branches A & B). Else `return miniprot_output` (the path; Branch C).

**Gotcha (`ERROR` substring is greedy):** failure mode 2 fires on *any* occurrence of the substring `ERROR` (uppercased) anywhere in stderr — including inside a sequence name, a warning, or a benign log line. A reference protein record literally named `…ERROR…` whose name leaks into stderr would abort the entire miniprot stage. This is intentional conservatism (favour Liftoff-only over a possibly-corrupt miniprot GFF) but is byte-identity-relevant: whether miniprot output is used at all hinges on this string scan.

#### 6.2.6 `check_miniprot_installed`

`check_miniprot_installed()` (`run_miniprot.py:54`) returns `True`/`False`:
1. `command = ["miniprot", "--version"]`.
2. `try: subprocess.run(command); installed = True`.
3. `except (FileNotFoundError, PermissionError, NotADirectoryError, subprocess.SubprocessError): pass` (leaves `installed = False`).

**Gotcha (narrowed except, V1.1b):** the except is deliberately narrow — `MemoryError`, `OSError` (other than the listed subtypes), `KeyboardInterrupt`, `SystemExit` propagate instead of being mis-reported as "not installed". Also note `subprocess.run` here does **not** check the return code, so a `miniprot` that exists but errors on `--version` still reports installed. The function's return value is currently **unused** by `exec_miniprot` (it is called for its side effects only).

#### 6.2.7 Re-ingest and the miniprot ID-mapping (downstream consumers)

After `run_miniprot` returns, Step 5 (`lifton/lifton.py:466`) builds `m_feature_db` from the path **or** bytes blob (same `annotation.Annotation(...)` call shape as Liftoff, `infer_genes=False`, `infer_transcripts=False`). If miniprot returned `None`, `m_feature_db = None` and the whole miniprot contribution is skipped with a stderr note (`:476-481`).

The miniprot GFF3 record model that LiftOn consumes (built by miniprot's `--gff-only`):

| Featuretype | How LiftOn reads it | Where |
|---|---|---|
| `mRNA` | top-level per-alignment record; carries `ID=<miniprot id>` and `Target=<ref_protein_id> <qstart> <qend>` | `miniprot_id_mapping` `lifton_utils.py:447`; Step 8 loop `lifton.py:551` |
| `CDS` | children of an `mRNA`; the protein-coding blocks | `lifton_miniprot_with_ref_protein` `run_miniprot.py:272`; `LiftOn_miniprot_alignment` `lifton_utils.py:322` |
| `stop_codon` | counted **with** CDS via `featuretype=('CDS','stop_codon')` when assembling the miniprot transcript's CDS chain | `LiftOn_miniprot_alignment` `lifton_utils.py:319,322` |

**ID-mapping — `miniprot_id_mapping(m_feature_db)`** (`lifton_utils.py:431`): builds two dicts by walking `m_feature_db.features_of_type("mRNA")`:
1. `miniprot_id = feature["ID"][0]` — the miniprot-assigned mRNA id.
2. `aa_trans_id = str(feature.attributes["Target"][0]).split(" ")[0]` — the **reference protein id is the first whitespace-delimited token of the `Target=` attribute** (miniprot writes `Target=<id> <start> <end> [<strand>]`).
3. `ref_id_2_m_id_trans_dict[aa_trans_id]` → **list** of miniprot ids (a reference protein can map to many target loci — appended if the key already exists, `:450-453`).
4. `m_id_2_ref_id_trans_dict[miniprot_id]` → the single `aa_trans_id` (`:454`).

When `m_feature_db is None`, both dicts are empty (`:444-445`).

**Gotcha (Target token parsing):** if a reference protein id contains a space, `.split(" ")[0]` truncates it. Reference protein FASTA headers must therefore have whitespace-free record ids for the round-trip `aa_trans_id` to match the reference annotation's transcript ids. This is the load-bearing join key between miniprot output and the reference annotation.

**Miniprot-only emission (Step 8) — `process_miniprot`** (`run_miniprot.py:286`) walks every miniprot `mRNA` not already covered by Liftoff and applies these filters before emitting it as an additional gene copy:
- **Overlap filter** (`:291`): `check_ovps_ratio(mtrans, Interval(mtrans.start, mtrans.end, mtrans_id), args.overlap, tree_dict)`. If the miniprot mRNA overlaps an existing Liftoff locus by a fraction `> args.overlap` (**default 0.1**, `lifton.py:136`), it is skipped (`is_overlapped == True` → no emission).
- **Reference resolution** (`:295-298`): `get_ref_ids_miniprot(...)`; skip if `ref_gene_id is None`.
- **Processed-pseudogene filter** (`:304-305`): `if len(CDS children)==1 and ref_trans_exon_num_dict[ref_trans_id] > 1: return None` — a single-CDS miniprot hit against a multi-CDS reference is treated as a retro-copied processed pseudogene and dropped.
- **Length-ratio filter** (`:306-307`): `miniprot_trans_ratio = (mtrans.end - mtrans.start + 1) / ref_features_len_dict[ref_gene_id]`; keep only if `args.min_miniprot < ratio < args.max_miniprot` (**defaults 0.9 and 1.5**, `lifton.py:82,86`). The kept ratio is recorded as the `miniprot_annotation_ratio` attribute formatted to 3 decimals (`:309`).

**Gotcha (strict-inequality length filter):** the comparison is `ratio > min_miniprot AND ratio < max_miniprot` — both bounds are **exclusive**. A ratio of exactly 0.9 or exactly 1.5 is rejected.

**Per-gene miniprot alignment inside Step 7 — `LiftOn_miniprot_alignment`** (`lifton_utils.py:263`) is invoked when a Liftoff transcript's protein identity is `< 1` and chooses the best miniprot candidate to chain with. Its two acceptance checks per candidate miniprot id:
1. **Overlap check** (`:291-294`): the miniprot mRNA interval must overlap the Liftoff transcript interval (`segments_overlap_length`) **and** share the same `seqid`; else `continue`.
2. **Cross-gene check** (`:303-315`): using the per-chromosome `IntervalTree`, if the miniprot interval overlaps any gene locus that the Liftoff transcript does **not** overlap (`miniprot_cross_gene_loci`), the candidate is rejected. This prevents a miniprot hit that bridges two genes from polluting one of them.
For each surviving candidate it builds a temporary `Lifton_TRANS`, aligns with `align.lifton_parasail_align`, and keeps the highest-identity one (`:328-330`), setting `has_valid_miniprot = True`.

---

### 6.3 Native bindings facade (`lifton/native_bindings/`)

Phase 10 introduced a `native_bindings` package presenting one uniform shape over two backends:

- **minimap2 → `mappy`** — a real PyO3 binding, already on PyPI; wrapped by `MinimapAligner`.
- **miniprot → `pyminiprot`** — a real native binding **does not yet exist**; `MiniprotIndex` presents the API a real binding *would* and falls back to the subprocess underneath.

Public surface (`native_bindings/__init__.py:13-31`): `MinimapAligner`, `MiniprotIndex`, `GFF3Bundle`, `GFF3Hit`, `MinimapHit`, `is_mappy_available`, `is_pyminiprot_native_available`.

#### 6.3.1 `types.py` — common data structures

**`MinimapHit`** (`types.py:11`, `@dataclass(frozen=True)`) — one minimap2 hit, projected from `mappy.Alignment` to only the fields LiftOn (or a SAM round-trip) consumes:

| Field | Type | Meaning |
|---|---|---|
| `query_name` | `str` | the query/feature name LiftOn supplied |
| `ctg` | `str` | target contig name |
| `r_st` | `int` | **0-based** reference start (mappy convention) |
| `r_en` | `int` | **exclusive** reference end |
| `q_st` | `int` | 0-based query start |
| `q_en` | `int` | exclusive query end |
| `strand` | `int` | `+1` or `-1` |
| `mapq` | `int` | mapping quality |
| `NM` | `int` | edit distance |
| `cigar_str` | `str` | CIGAR string |
| `is_primary` | `bool` | primary-alignment flag |

**`GFF3Hit`** (`types.py:33`, `@dataclass(frozen=True)`) — one miniprot GFF3 row decoded into the 9 GFF3 columns:

| Field | Type | Source column |
|---|---|---|
| `seqid` | `str` | col 0 |
| `source` | `str` | col 1 |
| `featuretype` | `str` | col 2 |
| `start` | `int` | col 3 (parsed via `int`) |
| `end` | `int` | col 4 (parsed via `int`) |
| `score` | `str` | col 5 (kept as string — may be `.`) |
| `strand` | `str` | col 6 |
| `phase` | `str` | col 7 (kept as string — may be `.`) |
| `attributes` | `str` | col 8 (**raw** attribute string, e.g. `"ID=MP1;Target=tx1 1 66"`) |

`GFF3Hit.from_gff_line(cls, line)` (`types.py:51`) — parse one tab-separated row:
1. Return `None` if `not line` **or** `line.startswith("#")` (blank/comment/directive).
2. `cols = line.rstrip("\n").split("\t")`; return `None` if `len(cols) != 9`.
3. `try: start = int(cols[3]); end = int(cols[4]) except ValueError: return None`.
4. Construct and return the `GFF3Hit` with columns mapped as in the table.

**Gotcha (silent drops):** any line that isn't exactly 9 tab-separated columns, or whose start/end aren't integers, is silently skipped (`None`). This mirrors how the subprocess path's GFF3 was line-parsed, so the facade is parse-equivalent. The `attributes` column is **not** parsed into a dict — it is kept verbatim, so downstream consumers re-parse it (preserving exact byte layout).

**`GFF3Bundle`** (`types.py:78`, mutable `@dataclass`):

| Field | Type | Default | Meaning |
|---|---|---|---|
| `hits` | `List[GFF3Hit]` | `[]` (`default_factory`) | parsed records |
| `raw_bytes` | `bytes` | `b""` | the **full** GFF3 stdout blob, carried so the gffbase ingest path can use it without reconstructing lines |

Implements `__iter__` (over `hits`) and `__len__` (`len(hits)`).

**Gotcha (raw_bytes is the byte-identity anchor):** `raw_bytes` is the exact stdout from miniprot. The streaming-native path (§6.2.3 Branch A) returns `bundle.raw_bytes`, so the *bytes* fed into gffbase are byte-identical to the subprocess streaming branch — the parsed `hits` list is essentially a convenience view that the current LiftOn path does not actually consume (it ingests `raw_bytes`).

#### 6.3.2 `minimap_facade.py` — `MinimapAligner`

`is_mappy_available()` (`minimap_facade.py:22`): `try: import mappy; return True except ImportError: return False`. Used so hermetic CI (no `mappy` installed) falls back to the subprocess path. The facade itself **never spawns a subprocess** — subprocess minimap2 is the legacy `align_features.py`.

`MinimapAligner.__init__(self, target_fa, *, preset=None, threads=1, mm2_options="", best_n=50)` (`:38`):
1. If `not is_mappy_available()`: raise `RuntimeError("mappy is not installed; install it via pip install mappy or conda install -c bioconda mappy.")`.
2. `import mappy`.
3. `self._aligner = mappy.Aligner(fn_idx_in=target_fa, preset=preset, n_threads=threads, best_n=best_n)`.
4. If `not self._aligner`: raise `RuntimeError(f"mappy.Aligner failed to open target FASTA: {target_fa}")`.
5. `self._lock = threading.Lock()` — only used in tests/debug; `mappy.Aligner.map` is itself thread-safe.

| Construction param | Default | CLI equivalence |
|---|---|---|
| `preset` | `None` | matches `minimap2 -x <preset>` (e.g. `splice`) |
| `threads` | `1` | `minimap2 -t N` (passed as `n_threads`) |
| `mm2_options` | `""` | the same string LiftOn's `args.mm2_options` held (documented as parameter-equivalent; **note** the current implementation accepts it but does not parse it into individual `mappy` knobs — it is a forward-compat surface) |
| `best_n` | `50` | matches `-N 50` from the default `mm2_options` |

`map(self, query_name, query_seq)` (`:74`) — generator yielding one `MinimapHit` per `mappy` `Alignment`. For each `hit` from `self._aligner.map(query_seq)` it projects the fields, coercing with `int(...)`/`str(...)`/`bool(...)` and defaulting `NM` to `0` and `is_primary` to `True` via `getattr` (so a future `mappy` lacking those attributes degrades gracefully). The docstring notes `mappy.Aligner.map` **releases the GIL** during the C kernel, enabling `ThreadPoolExecutor` scaling — this is *why* `--native` unlocks real threading on gffbase backends.

`map_one_best(self, query_name, query_seq)` (`:97`) — iterate `map(...)`, skip non-primary hits, return the highest-`mapq` primary hit, or `None`.

Read-only helpers: `n_seq` (property → `self._aligner.n_seq`), `seq_names()` → `list(self._aligner.seq_names)`.

#### 6.3.3 `miniprot_facade.py` — `MiniprotIndex`

`is_pyminiprot_native_available()` (`miniprot_facade.py:39`): `try: import pyminiprot; return hasattr(pyminiprot, "Index") except ImportError: return False`. Today always `False` (no such package).

`MiniprotIndex.__init__(self, target_fa, *, mp_options="", miniprot_path="miniprot", ref_proteins_path=None)` (`:74`):
1. Store `target_fa`, `mp_options`, `miniprot_path`, `_ref_proteins_path`, `_cached_bundle=None`.
2. `self._native = None`; if `is_pyminiprot_native_available()`: try `import pyminiprot; self._native = pyminiprot.Index(target_fa, mp_options=mp_options)`, swallowing any exception back to `None`. (Forward-compat; never taken today.)

`align_all(self, ref_proteins_path=None) -> GFF3Bundle` (`:101`) — run miniprot **once** over the whole protein FASTA, cache, return:
1. If `self._cached_bundle is not None: return self._cached_bundle` (O(1) on repeat).
2. `proteins = ref_proteins_path or self._ref_proteins_path`; if both `None` raise `ValueError`.
3. **Native branch** (`self._native is not None`, never taken today): iterate `self._native.align_file(proteins)`, build `raw` from `native_hit.to_gff_line()+"\n"`, parse each via `GFF3Hit.from_gff_line`, cache and return.
4. **Subprocess fallback** (`:134-166`):
   - `cmd = [miniprot_path, "--gff-only", target_fa, proteins] + [opt for opt in mp_options.split(" ") if opt]` — **same argv shape as `run_miniprot`** (§6.2.2).
   - `proc = subprocess.Popen(cmd, stdout=PIPE, stderr=PIPE, bufsize=1<<20)`; `stdout_bytes, stderr_bytes = proc.communicate()`.
   - If `proc.returncode != 0`: raise `RuntimeError(f"miniprot exited with code {rc}: {err}")`.
   - If `stderr_bytes and b"ERROR" in stderr_bytes.upper()`: raise `RuntimeError("miniprot reported an ERROR …")`.
   - Parse: for each `line` in `stdout_bytes.decode("utf-8","replace").splitlines()`, `GFF3Hit.from_gff_line(line)`, collect non-`None`.
   - `self._cached_bundle = GFF3Bundle(hits=hits, raw_bytes=stdout_bytes)`; return it.

**Gotcha (facade uses `communicate()`, not the chunked drain):** `align_all`'s subprocess fallback uses `proc.communicate()` — it does **not** use `_drain_stream_chunks`. So the native+streaming path (§6.2.3 Branch A) re-acquires the very double-allocation that Branch B's chunked drain was built to avoid. This is acceptable because Branch A is opt-in and the facade is forward-compat scaffolding; the bytes are still identical, only peak RSS differs.

**Gotcha (facade re-raises rather than returns `None`):** unlike `run_miniprot`, the facade *raises* `RuntimeError` on miniprot failure; `run_miniprot` Branch A catches it and converts to `return None` (`run_miniprot.py:152-158`), preserving the "skip miniprot, continue Liftoff-only" contract.

`align(self, protein_seq) -> Iterator[GFF3Hit]` (`:168`) — per-protein query shape. Today it requires `align_all()` to have run first (else raises `RuntimeError`), and simply returns `iter(self._cached_bundle.hits)` — i.e. **all** hits, unfiltered, because the subprocess GFF carries the protein record id (not its sequence) in `Target=`. Downstream filtering by `Target=` is the caller's job. This matches the legacy contract where `miniprot.gff3` held all hits regardless of which protein triggered them.

Properties: `is_native` (`:196`) → `self._native is not None`; `raw_bytes` (`:201`) → `self._cached_bundle.raw_bytes` (raises if `align_all()` not yet called).

#### 6.3.4 The drop-in contract for a future real PyO3 binding

The facade is engineered so a real `pyminiprot` binding swaps in with **no call-site change**:
- `is_pyminiprot_native_available()` flips to `True` when `import pyminiprot` succeeds and exposes an `Index`.
- `MiniprotIndex.__init__` then sets `self._native` to a real `pyminiprot.Index`.
- `align_all` / `align` route through the native iterator, building `GFF3Hit` from `native_hit.to_gff_line()` — the **same `GFF3Hit.from_gff_line` parser**, so the structured shape callers see is invariant across subprocess and native backends.
- `MinimapAligner` already wraps the real `mappy` binding directly; the only fallback is "mappy not installed → use subprocess minimap2 via `align_features.py`".

The unifying invariant: every backend ultimately yields **`raw_bytes` byte-identical to the subprocess stdout** (miniprot) or **`MinimapHit` field-identical to the SAM round-trip** (minimap2), which is exactly what keeps the 24-cell `--native` matrix green.

---

### 6.4 Tool command contracts & default-parameter table

Consolidated invocation contracts for every external tool, with exact argv/flags and defaults.

#### 6.4.1 minimap2 (via vendored Liftoff)

Built in `lifton/liftoff/align_features.py`. Liftoff resolves the binary via `get_minimap_path(args)` (`:106`): `"minimap2"` unless `args.m` is set (the `-m PATH` CLI flag). Two argv shapes depending on genome size (threshold `max_single_index_size = 4_000_000_000`, `:64`):

**Standard (genome ≤ 4 Gbp)** (`:76-78`):
```
minimap2 -o <output.sam> <target.mmi> <features.fa> <*mm2_options> -t <threads>
```
preceded by index build (`build_minimap2_index`, `:122`) **only if** `<target>.mmi` does not already exist:
```
minimap2 -d <target>.mmi <target> <*mm2_options> -t <threads>
```

**Split-index (genome > 4 Gbp)** (`:71-74`):
```
minimap2 -o <output.sam> <target> <features.fa> <*mm2_options> --split-prefix <prefix> -t <threads>
```

| Token | Value / default | Notes |
|---|---|---|
| binary | `minimap2` (or `-m PATH`) | resolved via PATH if bare |
| `-o` | `<features>_to_<target>.sam` under `args.directory` | SAM output |
| `<target.mmi>` / `<target>` | prebuilt index, else raw FASTA | index reused if present |
| `<features.fa>` | `<features_name>_genes.fa` | extracted reference features |
| `args.mm2_options` | **`-a --end-bonus 5 --eqx -N 50 -p 0.5`** (`lifton/lifton.py:66-68`) | split on single spaces (`.split(" ")`) |
| `-t` | `str(threads)` from `args.threads` (**default 1**, `lifton.py:107`) | appended last |
| `--split-prefix` | only on the >4 Gbp path | enables minimap2 split-index mode |

Decoding the default `mm2_options`: `-a` = SAM output; `--end-bonus 5` = bonus for reaching read ends (favours full-length feature alignment); `--eqx` = use `=`/`X` CIGAR ops (match/mismatch distinction); `-N 50` = retain up to 50 secondary alignments; `-p 0.5` = report secondaries scoring ≥ 0.5× the best.

**Gotcha (preset-equivalence under `--native`):** `MinimapAligner` exposes `preset`, `threads`, `best_n=50` (matching `-N 50`) and `mm2_options` — but the in-process `mappy` path does **not** re-parse `mm2_options` into individual knobs. The byte-identity matrix passes because the native path is validated against the subprocess SAM round-trip on the synthetic fixtures, not because the option strings are mechanically equivalent.

#### 6.4.2 miniprot (subprocess)

Built in `run_miniprot.run_miniprot` (`:122-125`) and identically in `MiniprotIndex.align_all` (`:135-140`):
```
miniprot --gff-only <target_genome.fa> <ref_proteins.fa> <*mp_options>
```

| Token | Value / default | Notes |
|---|---|---|
| binary | `miniprot` | **no path override**; must be on PATH |
| `--gff-only` | always | GFF3-only to stdout |
| `<target_genome.fa>` | positional 1 | genome to map into |
| `<ref_proteins.fa>` | positional 2 | reference protein queries |
| `args.mp_options` | **`""`** (empty; `lifton.py:70`) | split on single spaces, empty tokens dropped |
| threads | not added | single-threaded unless `-t N` placed inside `mp_options` |

Output destinations: file `<outdir>miniprot/miniprot.gff3` (default), or stdout→RAM (`--stream`). Failure → `None` → Liftoff-only.

#### 6.4.3 mappy (in-process, `--native`)

`mappy.Aligner(fn_idx_in=target_fa, preset=preset, n_threads=threads, best_n=best_n)` (`minimap_facade.py:55`).

| Param | Default | Equivalent |
|---|---|---|
| `fn_idx_in` | target FASTA path | — |
| `preset` | `None` | `-x <preset>` |
| `n_threads` | `1` | `-t` |
| `best_n` | `50` | `-N 50` |

`Aligner.map(seq)` releases the GIL during the C kernel (enables `ThreadPoolExecutor`).

#### 6.4.4 parasail (in-process, `lifton/align.py`)

Two distinct matrices/penalty sets, both run with the global Needleman-Wunsch traceback kernel `parasail.nw_trace_scan_sat(query, target, gap_open, gap_extend, matrix)` (i.e. **global** alignment with full traceback, "scan" vectorisation, saturation-checked).

**Protein alignment** — `parasail_align_protein_base(protein_seq, ref_protein_seq)` (`align.py:84`):

| Param | Value | Source |
|---|---|---|
| matrix | `parasail.Matrix("blosum62")` | `align.py:95` |
| `gap_open` | **11** | `align.py:96` |
| `gap_extend` | **1** | `align.py:97` |
| empty-input guard | raises `LiftOnAlignmentError` if either seq is `""` | `align.py:102-107` |

**DNA / transcript alignment** — `parasail_align_DNA_base(trans_seq, ref_trans_seq)` (`align.py:134`):

| Param | Value | Source |
|---|---|---|
| matrix | `parasail.matrix_create("ACGTN*", 1, -3)` (match +1, mismatch −3) | `align.py:148` |
| `gap_open` | **5** | `align.py:149` |
| `gap_extend` | **2** | `align.py:150` |
| input sanitisation | both seqs passed through `_sanitise_for_parasail_dna` | `align.py:151-152` |

`_sanitise_for_parasail_dna(seq)` (`align.py:16`): uppercases; if every char is in `{A,C,G,T,N,*}` (`_PARASAIL_DNA_ALPHABET = frozenset("ACGTN*")`, `:13`) returns as-is; otherwise replaces every out-of-alphabet char (IUPAC ambiguity codes R/Y/S/W/K/M/B/D/H/V, gaps, lowercase residue) with `N`. This prevents crashing parasail's C kernel on an unknown byte (V4.2 fix).

**Gotcha (matrix alphabet vs sanitiser must agree):** the DNA matrix is created over exactly `"ACGTN*"`. The sanitiser coerces to exactly that alphabet. If one is widened without the other, parasail either scores garbage (unknown char in matrix) or crashes (char outside matrix). The `*` is included to score the stop-codon translation token / terminal sentinel.

**Identity computation:** both protein and DNA paths compute `identity = matches / length` from the traceback via `get_id_fraction.get_AA_id_fraction` / `get_DNA_id_fraction` (`align.py:128-129`, `:176-177`). This identity drives the Step-7 branch logic: a Liftoff transcript with protein `identity == 1` is accepted as-is; `identity < 1` triggers the miniprot-chaining / ORF-rescue path (`run_liftoff.py:171-183`).

#### 6.4.5 Threshold/parameter cross-reference (all defaults)

| Parameter | Default | CLI flag | Used by | Source |
|---|---|---|---|---|
| `mm2_options` | `-a --end-bonus 5 --eqx -N 50 -p 0.5` | `-mm2_options` | minimap2 (Liftoff) | `lifton.py:66` |
| `mp_options` | `""` | `-mp_options` | miniprot | `lifton.py:70` |
| `threads` | `1` | `-t/--threads` | minimap2, Step-7 pool | `lifton.py:107` |
| `overlap` | `0.1` | `-overlap` | `check_ovps_ratio` (miniprot overlap reject) | `lifton.py:136` |
| `min_miniprot` | `0.9` | `-min_miniprot` | `process_miniprot` length-ratio (exclusive lower) | `lifton.py:82` |
| `max_miniprot` | `1.5` | `-max_miniprot` | `process_miniprot` length-ratio (exclusive upper) | `lifton.py:86` |
| `a` (coverage) | `0.5` | `-a` | Liftoff mapping coverage threshold | `lifton.py:73` |
| `s` (sequence id) | `0.5` | `-s` | Liftoff sequence-identity threshold | `lifton.py:77` |
| `sc` | `1.0` | `-sc` | Liftoff copy-number score | `lifton.py:132` |
| `merge_strategy` | `create_unique` | `--merge-strategy` | FeatureDB ingest of Liftoff/miniprot GFFs | `lifton.py:29` |
| recursion-limit floor | `10000` | (internal) | `run_liftoff` Liftoff traversal | `run_liftoff.py:57` |
| stream drain chunk | `65536` (1<<16) | (internal) | `_drain_stream_chunks` | `run_miniprot.py:6,171` |
| Popen pipe bufsize | `1048576` (1<<20) | (internal) | streaming/native miniprot | `run_miniprot.py:168` |
| parasail protein gap open/extend | `11` / `1` | (internal) | protein alignment | `align.py:96-97` |
| parasail DNA gap open/extend | `5` / `2` | (internal) | DNA alignment | `align.py:149-150` |
| parasail DNA match/mismatch | `+1` / `−3` | (internal) | DNA matrix | `align.py:148` |
| mappy `best_n` | `50` | (internal, matches `-N 50`) | `MinimapAligner` | `minimap_facade.py:45` |

**Gotcha (the `=STR` metavar oddity):** `mm2_options` and `mp_options` use single-dash long flags (`-mm2_options`, `-mp_options`) with metavar `=STR`. Users must pass them as `-mm2_options="-a --eqx ..."`; the quoting matters because the value is later naively `.split(" ")`-tokenised.
