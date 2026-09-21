from lifton import (align, coding, coreutils, drop_ledger, logger, lifton_class, lifton_utils,
                    orf_completion)
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
import subprocess, os, sys
from intervaltree import Interval, IntervalTree
from lifton.tool_execution import (
    describe_returncode,
    detect_miniprot_stage,
    execution_record,
    record_execution,
    run_with_bounded_stderr,
)


def _drain_stream_chunks(proc, *, chunk_size: int = 65536):
    """Phase 15c (V3.10) — drain a Popen's stdout/stderr in bounded
    chunks instead of `proc.communicate()`.

    ``communicate()`` allocates ``bytes(stdout)`` + ``bytes(stderr)`` and
    holds both simultaneously. Chunked reads keep each live pipe read bounded
    by ``chunk_size``. ``BytesIO.getvalue()`` then exposes each accumulated
    buffer without the explicit second full-size allocation required by
    ``b"".join(chunks)``.

    Returns ``(stdout_bytes, stderr_bytes, returncode)``.
    """
    import threading

    stdout_buffer = BytesIO()
    stderr_buffer = BytesIO()

    def _drain(stream, sink):
        if stream is None:
            return
        try:
            while True:
                buf = stream.read(chunk_size)
                if not buf:
                    return
                sink.write(buf)
        finally:
            try:
                stream.close()
            except Exception:
                pass

    t_out = threading.Thread(
        target=_drain, args=(proc.stdout, stdout_buffer), daemon=True,
    )
    t_err = threading.Thread(
        target=_drain, args=(proc.stderr, stderr_buffer), daemon=True,
    )
    t_out.start()
    t_err.start()
    rc = proc.wait()
    t_out.join()
    t_err.join()
    return stdout_buffer.getvalue(), stderr_buffer.getvalue(), rc


def _resolve_miniprot_threads(threads):
    """Resolve the miniprot thread count, honouring the LIFTON_MINIPROT_THREADS
    escape hatch (Iteration 17).

    Returns the integer thread count to request via ``-t``, or ``None`` to
    request none (leaving miniprot at its own hard-coded binary default of 4 —
    the pre-Iteration-17 behaviour). Semantics:

      * env unset                 -> use ``threads`` (LiftOn's -t/--threads
                                     budget), but gated to ``> 1`` so the
                                     default ``-t 1`` run emits NO ``-t`` and is
                                     byte-identical to the pre-Iteration-17
                                     command.
      * LIFTON_MINIPROT_THREADS=0 -> None (never add ``-t``; reproduces the old
                                     fixed default — the clean A/B baseline /
                                     opt-out).
      * LIFTON_MINIPROT_THREADS=k -> k for k>=1 (explicit override that bypasses
                                     the ``> 1`` gate so an A/B or user can pin
                                     any count, including 1).

    A non-integer env value is ignored (falls back to ``threads``) so a typo
    cannot crash a run.
    """
    raw = os.environ.get("LIFTON_MINIPROT_THREADS")
    if raw is not None and raw != "":
        try:
            forced = int(raw)
        except (TypeError, ValueError):
            forced = None
        if forced is not None:
            if forced <= 0:
                return None          # =0 -> reproduce the old fixed default
            return forced            # =k -> explicit override (bypasses gate)
    # env unset / unparseable: use the LiftOn budget, gated to >1 so the
    # default -t1 invocation is byte-identical (no -t emitted -> miniprot's
    # own default of 4 threads, exactly as before Iteration 17).
    try:
        n = int(threads)
    except (TypeError, ValueError):
        n = 1
    return n if n > 1 else None


def _build_miniprot_command(miniprot_path, tgt_genome, ref_proteins_file,
                            mp_options, threads, table=None):
    """Build the miniprot subprocess command list (Iteration 17).

    Plumbs LiftOn's -t/--threads into miniprot's own ``-t`` so it scales past
    its hard-coded binary default of 4 threads. Previously miniprot ran at 4
    regardless of LiftOn's -t, becoming the serial tail of the concurrent
    Step 4 on cross-species data (where miniprot does heavy protein alignment).

    The thread flag is appended only when BOTH:
      * a thread count is resolved (see :func:`_resolve_miniprot_threads`), AND
      * the user has not already supplied an explicit ``-t`` in ``mp_options``
        (matched by an EXACT-token / glued-token check ``opt.startswith("-t")``
        — NOT a bare substring, so ``--trans`` does not false-positive and
        silently suppress the flag, and ``-t8`` / ``-t 8`` are both respected).

    The default -t1 run resolves to ``None`` -> no ``-t`` added -> byte-identical
    to the pre-Iteration-17 command.
    """
    opts = [opt for opt in mp_options.split(" ") if opt]
    command = [miniprot_path, "--gff-only", tgt_genome, ref_proteins_file] + opts
    n = _resolve_miniprot_threads(threads)
    user_set_threads = any(opt.startswith("-t") for opt in opts)
    if n is not None and not user_set_threads:
        command += ["-t", str(n)]
    if table is not None and table != coding.DEFAULT_TRANSL_TABLE:
        # miniprot speaks one genetic code per invocation. Emitting -T only for
        # a non-standard group keeps the ordinary command byte-identical, and
        # -P keeps the extra group's MP ids from colliding with the first
        # group's (both would otherwise start at MP000001).
        if not any(opt.startswith("-T") for opt in opts):
            command += ["-T", str(table)]
        if not any(opt.startswith("-P") for opt in opts):
            command += ["-P", f"MPT{table}"]
    return command


def transl_table_groups(ref_proteins_file):
    """Split the reference proteins into one group per declared genetic code.

    Returns ``[(table, fasta_path), ...]`` with the standard-code group first.
    An annotation that declares nothing unusual -- which is almost all of them
    -- yields exactly ``[(1, ref_proteins_file)]``, so the file handed to
    miniprot and the command built from it are unchanged.

    Splitting matters because miniprot scores against one code per run: on the
    13 human mitochondrial proteins, running them under the standard code
    instead of table 2 costs every one of them identity (mean 1.0000 -> 0.9302,
    up to 16 spurious in-frame stops on COX1) and moves 4 of the 13 hits'
    coordinates.
    """
    from lifton import extract_sequence

    overrides = extract_sequence.read_transl_table_sidecar(ref_proteins_file)
    if not overrides:
        return [(coding.DEFAULT_TRANSL_TABLE, ref_proteins_file)]
    grouped = {}
    with open(ref_proteins_file) as handle:
        table, buffered = coding.DEFAULT_TRANSL_TABLE, None
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].split()[0] if len(line) > 1 else ""
                table = overrides.get(identifier, coding.DEFAULT_TRANSL_TABLE)
                buffered = grouped.setdefault(table, [])
            if buffered is not None:
                buffered.append(line)
    if len(grouped) <= 1 and coding.DEFAULT_TRANSL_TABLE in grouped:
        return [(coding.DEFAULT_TRANSL_TABLE, ref_proteins_file)]
    groups = []
    for table in sorted(grouped, key=lambda t: (t != coding.DEFAULT_TRANSL_TABLE, t)):
        if table == coding.DEFAULT_TRANSL_TABLE:
            path = ref_proteins_file + ".table1.faa"
        else:
            path = f"{ref_proteins_file}.table{table}.faa"
        with open(path, "w") as handle:
            handle.writelines(grouped[table])
        groups.append((table, path))
    return groups


@dataclass(frozen=True)
class MiniprotArtifact:
    """Published database and bounded diagnostics from streamed miniprot."""

    database_path: str
    byte_count: int
    feature_count: int
    stderr_tail: str
    returncode: int


class _BoundedStderr:
    def __init__(self, limit: int = 64 << 10):
        self.limit = max(1, int(limit))
        self._tail = bytearray()
        self._error_carry = b""
        self.error_seen = False

    def feed(self, chunk: bytes) -> None:
        upper = self._error_carry + bytes(chunk).upper()
        if b"ERROR" in upper:
            self.error_seen = True
        self._error_carry = upper[-4:]
        self._tail.extend(chunk)
        if len(self._tail) > self.limit:
            del self._tail[:-self.limit]

    @property
    def text(self) -> str:
        return bytes(self._tail).decode("utf-8", errors="replace")


def _with_stderr_tail(message: str, stderr_state: _BoundedStderr) -> str:
    """Attach bounded miniprot diagnostics without retaining full stderr."""

    tail = stderr_state.text.strip()
    if not tail:
        return message
    return f"{message}\nminiprot stderr (bounded tail):\n{tail}"


def _stop_process(proc) -> None:
    """Best-effort process teardown used by all streaming failure paths."""
    try:
        running = proc.poll() is None
    except Exception:
        running = True
    if running:
        try:
            proc.terminate()
        except Exception:
            pass
    try:
        proc.wait(timeout=5)
        return
    except Exception:
        pass
    try:
        proc.kill()
    except Exception:
        pass
    try:
        proc.wait(timeout=5)
    except Exception:
        pass


def _remove_duckdb_staging(path: str) -> None:
    """Remove a staged DuckDB file and any crash-left WAL sidecar."""

    for candidate in (path, path + ".wal"):
        try:
            os.unlink(candidate)
        except FileNotFoundError:
            pass


def _refuse_destination_wal(path: str) -> None:
    """Fail closed instead of publishing a database beside an old WAL."""

    wal_path = path + ".wal"
    if os.path.lexists(wal_path):
        raise RuntimeError(
            "refusing to replace the miniprot database while a destination "
            f"WAL exists: {wal_path}; recover or remove the stale database "
            "and WAL together before retrying"
        )


def run_miniprot_streaming_db(
    outdir, args, tgt_genome, ref_proteins_file, *, chunk_size=1 << 20,
) -> MiniprotArtifact:
    """Stream miniprot stdout directly into a staged gffbase database.

    The published database is replaced atomically only after miniprot exits
    successfully, stderr contains no ``ERROR`` token, and at least one mRNA
    can be queried from a checkpointed database.
    """
    import threading
    from lifton.gffbase import ingest as _ingest

    miniprot_outdir = os.path.join(outdir, "miniprot")
    os.makedirs(miniprot_outdir, exist_ok=True)
    final_path = os.path.join(miniprot_outdir, "miniprot.duckdb")
    partial_path = final_path + ".partial"
    _remove_duckdb_staging(partial_path)
    _refuse_destination_wal(final_path)

    command = _build_miniprot_command(
        "miniprot", tgt_genome, ref_proteins_file, args.mp_options,
        getattr(args, "threads", 1),
    )
    try:
        proc = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            bufsize=1 << 20,
        )
    except (FileNotFoundError, PermissionError, NotADirectoryError,
            OSError, ValueError) as exc:
        record_execution(args, execution_record(
            "miniprot", command=command, stage="process_launch",
            status="launch_error", returncode=None,
            details={"error": str(exc)},
        ))
        raise
    stderr_state = _BoundedStderr()
    execution_was_recorded = False

    def _drain_stderr():
        if proc.stderr is None:
            return
        try:
            while True:
                chunk = proc.stderr.read(1 << 16)
                if not chunk:
                    return
                stderr_state.feed(chunk)
        finally:
            try:
                proc.stderr.close()
            except Exception:
                pass

    stderr_thread = threading.Thread(
        target=_drain_stderr, name="lifton-miniprot-stderr", daemon=True,
    )
    stderr_thread.start()
    decoder = _ingest.GFF3ChunkDecoder(max_chunk_bytes=chunk_size)
    connection = None

    def _records():
        if proc.stdout is None:
            return
        try:
            while True:
                chunk = proc.stdout.read(chunk_size)
                if not chunk:
                    break
                yield from decoder.feed(chunk)
            yield from decoder.finish()
        finally:
            try:
                proc.stdout.close()
            except Exception:
                pass

    try:
        connection, stats = _ingest.from_record_stream(
            _records(), dbfn=partial_path, force=True,
            dialect={"fmt": "gff3"}, directives=decoder.directives,
            build_rtree=False,
        )
        returncode = proc.wait()
        stderr_thread.join()
        if returncode != 0:
            record_execution(args, execution_record(
                "miniprot", command=command,
                stage=detect_miniprot_stage(stderr_state.text),
                status="failed", returncode=returncode,
                stderr_tail=stderr_state.text,
            ))
            execution_was_recorded = True
            raise RuntimeError(_with_stderr_tail(
                f"miniprot {describe_returncode(returncode)}", stderr_state,
            ))
        if stderr_state.error_seen:
            record_execution(args, execution_record(
                "miniprot", command=command,
                stage=detect_miniprot_stage(stderr_state.text),
                status="failed", returncode=returncode,
                stderr_tail=stderr_state.text,
                details={"reason": "ERROR token on stderr"},
            ))
            execution_was_recorded = True
            raise RuntimeError(_with_stderr_tail(
                "miniprot reported ERROR on stderr", stderr_state,
            ))
        if decoder.byte_count == 0 or stats.n_features_raw == 0:
            record_execution(args, execution_record(
                "miniprot", command=command,
                stage=detect_miniprot_stage(stderr_state.text),
                status="failed", returncode=returncode,
                stderr_tail=stderr_state.text,
                details={"reason": "empty GFF3 stream"},
            ))
            execution_was_recorded = True
            raise RuntimeError("miniprot produced an empty GFF3 stream")
        mrna_count = connection.execute(
            "SELECT COUNT(*) FROM features WHERE featuretype = 'mRNA'"
        ).fetchone()[0]
        if not mrna_count:
            record_execution(args, execution_record(
                "miniprot", command=command,
                stage=detect_miniprot_stage(stderr_state.text),
                status="failed", returncode=returncode,
                stderr_tail=stderr_state.text,
                details={"reason": "no mRNA features"},
            ))
            execution_was_recorded = True
            raise RuntimeError("miniprot output contains no mRNA features")
        connection.execute("CHECKPOINT")
        connection.close()
        connection = None
        if os.path.exists(partial_path + ".wal"):
            raise RuntimeError(
                "checkpointed miniprot database retained an unexpected WAL"
            )
        # Repeat the check at publication time in case a concurrent or
        # interrupted opener created a sidecar while miniprot was running.
        _refuse_destination_wal(final_path)
        os.replace(partial_path, final_path)
        from lifton.output_transaction import _fsync_directory
        _fsync_directory(Path(final_path).parent)
        record_execution(args, execution_record(
            "miniprot", command=command,
            stage=detect_miniprot_stage(stderr_state.text),
            status="success", returncode=returncode,
            stderr_tail=stderr_state.text,
        ))
        execution_was_recorded = True
        return MiniprotArtifact(
            database_path=final_path,
            byte_count=decoder.byte_count,
            feature_count=stats.n_features_raw,
            stderr_tail=stderr_state.text,
            returncode=returncode,
        )
    except BaseException:
        if not execution_was_recorded:
            try:
                failed_returncode = proc.poll()
            except Exception:
                failed_returncode = None
            record_execution(args, execution_record(
                "miniprot", command=command,
                stage=detect_miniprot_stage(stderr_state.text),
                status="failed", returncode=failed_returncode,
                stderr_tail=stderr_state.text,
                details={"reason": "streaming ingest or process failure"},
            ))
        _stop_process(proc)
        stderr_thread.join(timeout=5)
        if connection is not None:
            try:
                connection.close()
            except Exception:
                pass
        _remove_duckdb_staging(partial_path)
        raise


def check_miniprot_installed():
    """
        This function checks if miniprot is installed.

        Parameters:
        None

        Returns:
        True if miniprot is installed, False otherwise.
    """
    miniprot_path = "miniprot"
    command = [miniprot_path, "--version"]
    installed = False
    # V1.1b fix: narrow the bare `except:` so a real environmental failure
    # (MemoryError, OSError, KeyboardInterrupt) is not silently misreported
    # as "miniprot is not installed". Only the things that genuinely mean
    # "binary is missing or unrunnable" return False here.
    try:
        # Suppress the child's output and require a successful exit, matching the
        # documented twin run_liftoff.check_minimap2_installed. Without
        # stdout=DEVNULL the version banner is written to the inherited fd 1, which
        # under `-o stdout` injects it straight into the emitted GFF3 stream
        # (redirect_stdout only rebinds the Python-level object, not the fd). A
        # present-but-broken binary also used to report "installed".
        completed = subprocess.run(
            command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        installed = completed.returncode == 0
    except (FileNotFoundError, PermissionError, NotADirectoryError,
            subprocess.SubprocessError):
        pass
    return installed


def probe_miniprot_version():
    """Return miniprot's first version line, or ``None`` when unavailable."""

    try:
        completed = subprocess.run(
            ["miniprot", "--version"], stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, check=False, timeout=5,
        )
    except (FileNotFoundError, PermissionError, NotADirectoryError,
            OSError, subprocess.SubprocessError):
        return None
    output = completed.stdout or ""
    return next(
        (line.strip() for line in output.splitlines() if line.strip()), None,
    )


def run_miniprot(outdir, args, tgt_genome, ref_proteins_file):
    """
    Run miniprot and return the GFF3 result.

    Two return shapes are supported, controlled by ``args.stream``:

    * ``args.stream`` is False (default): write the GFF3 to disk at
      ``<outdir>/miniprot/miniprot.gff3`` and return the **path**.
      This is the legacy Phase 5 contract — preserved byte-for-byte.

    * ``args.stream`` is True: decode miniprot's stdout incrementally into a
      staged DuckDB database and return that database's **path**. No whole
      GFF3 byte blob or intermediate text file is materialized.

    Failure modes return ``None`` so the pipeline can finish a diagnostic
    Liftoff-only partial result. The CLI does not publish that partial GFF3
    unless ``--allow-partial-output`` was explicitly requested.
      1. Non-zero exit code from miniprot.
      2. miniprot prints "ERROR" to stderr (exit 0 + error is a known miniprot quirk).
      3. Output is absent or empty.

    Parameters
    ----------
    outdir : str
        Output directory for intermediate files.
    args : argparse.Namespace
        Must carry ``mp_options``; may carry ``stream`` (default False).
    tgt_genome : str
        Path to target genome FASTA.
    ref_proteins_file : str
        Path to reference proteins FASTA.

    Returns
    -------
    str | None
        Path to the miniprot GFF3 file (legacy mode) or direct DuckDB input
        (streaming mode), or None on failure.
    """
    # os.path.join: the caller passes "<out>/lifton_output" without a trailing
    # separator, and plain concatenation wrote "<out>/lifton_outputminiprot/"
    # (v1.0.10-v1.0.11). The streaming path above already joins correctly.
    miniprot_outdir = os.path.join(outdir, "miniprot", "")
    os.makedirs(miniprot_outdir, exist_ok=True)
    miniprot_output = miniprot_outdir + "miniprot.gff3"
    miniprot_path = "miniprot"
    groups = transl_table_groups(ref_proteins_file)
    commands = [
        _build_miniprot_command(
            miniprot_path, tgt_genome, path,
            args.mp_options, getattr(args, "threads", 1), table=table,
        )
        for table, path in groups
    ]
    command = commands[0]

    stream_mode = bool(getattr(args, "stream", False))
    if stream_mode and len(commands) > 1:
        # --stream ingests ONE process's stdout. Rather than let it and the
        # default path disagree on a multi-code reference -- they are required
        # to produce the same annotation -- this run takes the file path, which
        # runs every group. The only cost is the disk round-trip --stream
        # exists to avoid.
        logger.log_warning(
            f"The reference declares {len(groups)} genetic codes "
            f"({', '.join(str(t) for t, _ in groups)}); --stream runs one "
            f"miniprot invocation, so this run writes miniprot.gff3 instead. "
            f"Output is unaffected."
        )
        stream_mode = False
    if stream_mode:
        try:
            artifact = run_miniprot_streaming_db(
                outdir, args, tgt_genome, ref_proteins_file,
            )
            if artifact.stderr_tail:
                print(artifact.stderr_tail, end="", file=sys.stderr)
            return artifact.database_path
        except Exception as exc:
            print(
                f"\n[LiftOn] miniprot streaming ingest failed: {exc}\n"
                "Miniprot output will be skipped; LiftOn will stage a "
                "Liftoff-only partial result.",
                file=sys.stderr,
            )
            return None
    try:
        # Legacy file-write branch (Phase 5 baseline behaviour). Streaming
        # returned above after publishing its direct DuckDB artifact.
        execution_was_recorded = False
        stderr_text, stderr_error_seen, return_code = "", False, 0
        with open(miniprot_output, "w") as fw:
            # The first group writes straight through, so a single-code
            # reference -- which is nearly every reference -- produces exactly
            # the bytes it always did. A further group goes via its own file so
            # its repeated "##gff-version 3" can be dropped: one GFF3 may carry
            # that directive only once, at the top.
            proc = run_with_bounded_stderr(commands[0], stdout=fw)
            stderr_text += proc.stderr_tail
            stderr_error_seen = bool(
                getattr(proc, "stderr_error_seen", False)
                or proc.stderr_tail.upper().find("ERROR") >= 0
            )
            return_code = proc.returncode
            for table, group_command in zip((t for t, _ in groups[1:]),
                                            commands[1:]):
                if return_code != 0:
                    break
                part = f"{miniprot_output}.table{table}.part"
                with open(part, "w") as handle:
                    proc = run_with_bounded_stderr(group_command, stdout=handle)
                stderr_text += proc.stderr_tail
                stderr_error_seen = stderr_error_seen or bool(
                    getattr(proc, "stderr_error_seen", False)
                    or proc.stderr_tail.upper().find("ERROR") >= 0
                )
                if proc.returncode != 0:
                    # Report against the invocation that actually failed.
                    return_code, command = proc.returncode, group_command
                    break
                with open(part) as handle:
                    for line in handle:
                        if not line.startswith("#"):
                            fw.write(line)
                os.remove(part)
        output_size = (os.path.getsize(miniprot_output)
                       if os.path.exists(miniprot_output) else 0)

        # ── Failure mode 1: non-zero exit code ────────────────────────────
        if return_code != 0:
            record_execution(args, execution_record(
                "miniprot", command=command,
                stage=detect_miniprot_stage(stderr_text), status="failed",
                returncode=return_code, stderr_tail=stderr_text,
            ))
            execution_was_recorded = True
            print(
                f"\n[LiftOn] miniprot {describe_returncode(return_code)}. "
                "Miniprot output will be skipped; LiftOn will stage a "
                "Liftoff-only partial result.",
                file=sys.stderr,
            )
            return None

        # ── Failure mode 2: miniprot printed ERROR on stderr ─────────────
        if stderr_error_seen:
            record_execution(args, execution_record(
                "miniprot", command=command,
                stage=detect_miniprot_stage(stderr_text), status="failed",
                returncode=return_code, stderr_tail=stderr_text,
                details={"reason": "ERROR token on stderr"},
            ))
            execution_was_recorded = True
            print(
                "\n[LiftOn] miniprot reported an ERROR during mapping "
                "(exit code 0 but ERROR seen in output). "
                "Miniprot output will be skipped; LiftOn will stage a "
                "Liftoff-only partial result.",
                file=sys.stderr,
            )
            return None

        # ── Failure mode 3: output is absent or empty ─────────────────────
        if output_size == 0:
            record_execution(args, execution_record(
                "miniprot", command=command,
                stage=detect_miniprot_stage(stderr_text), status="failed",
                returncode=return_code, stderr_tail=stderr_text,
                details={"reason": "empty output"},
            ))
            execution_was_recorded = True
            print(
                "\n[LiftOn] miniprot produced an empty output. "
                "Miniprot output will be skipped; LiftOn will stage a "
                "Liftoff-only partial result.",
                file=sys.stderr,
            )
            return None

        record_execution(args, execution_record(
            "miniprot", command=command,
            stage=detect_miniprot_stage(stderr_text), status="success",
            returncode=return_code, stderr_tail=stderr_text,
        ))
        execution_was_recorded = True
    except Exception as exc:
        if not locals().get("execution_was_recorded", False):
            record_execution(args, execution_record(
                "miniprot", command=command, stage="process_launch",
                status="launch_error", returncode=None,
                details={"error": str(exc)},
            ))
        print(
            f"\n[LiftOn] miniprot failed unexpectedly: {exc}\n"
            "Miniprot output will be skipped; LiftOn will stage a "
            "Liftoff-only partial result.",
            file=sys.stderr,
        )
        return None

    return miniprot_output


def lifton_miniprot_with_ref_protein(
        m_feature, m_feature_db, ref_db, ref_gene_id, ref_trans_id, tgt_fai,
        ref_proteins, ref_trans, tree_dict, ref_features_dict, args,
        state_journal=None):
    """
        This function create a miniprot gene entry with reference protein.

        Parameters:
        - m_feature: miniprot feature
        - m_feature_db: miniprot feature database
        - ref_db: reference database
        - ref_gene_id: reference gene ID
        - ref_trans_id: reference transcript ID
        - tgt_fai: target fasta index
        - ref_proteins: reference protein dictionary
        - ref_trans: reference transcript dictionary
        - tree_dict: tree dictionary
        - ref_features_dict: reference features dictionary
        - args: LiftOn arguments

        Returns:
        lifton_gene: LiftOn gene instance
        lifton_transcript_id: LiftOn transcript ID
        lifton_status: LiftOn status
    """
    mtrans_id = m_feature.attributes["ID"][0]
    # Create LifOn gene instance
    m_gene_feature = coreutils.clone_feature(m_feature)
    m_gene_feature.featuretype = "gene"
    lifton_gene = lifton_class.Lifton_GENE(
        ref_gene_id, m_gene_feature,
        coreutils.clone_attributes(ref_db[ref_gene_id].attributes), tree_dict,
        ref_features_dict, args, state_journal=state_journal,
    )
    lifton_gene.update_gene_info(m_feature.seqid, m_feature.start, m_feature.end)
    # Create LifOn transcript instance
    Lifton_trans = lifton_gene.add_miniprot_transcript(ref_trans_id, coreutils.clone_feature(m_feature), ref_db[ref_trans_id].attributes, ref_features_dict)
    lifton_gene.update_trans_info(Lifton_trans.entry.id, m_feature.seqid, m_feature.start, m_feature.end)
    # Create exon / CDS entries
    cdss = m_feature_db.children(m_feature, featuretype='CDS')  # Replace 'exon' with the desired child feature type
    for cds in list(cdss):
        lifton_gene.add_exon(Lifton_trans.entry.id, cds)
        cds_copy = coreutils.clone_feature(cds)
        lifton_gene.add_cds(Lifton_trans.entry.id, cds_copy)
    # Update LiftOn status
    lifton_status = lifton_class.Lifton_Status()                
    m_entry = m_feature_db[mtrans_id]
    m_lifton_aln = align.lifton_parasail_align(Lifton_trans, m_entry, tgt_fai, ref_proteins, ref_trans_id)
    lifton_status.annotation =  "miniprot"
    lifton_status.lifton_aa = m_lifton_aln.identity
    return lifton_gene, Lifton_trans, Lifton_trans.entry.id, lifton_status


def _miniprot_rescue_band_ok(ratio, args):
    """Iteration 22: is `ratio` inside the (wider) miniprot-only-rescue length
    band? Consulted ONLY when --miniprot-rescue is ON, for candidates that fall
    outside the default -min_miniprot/-max_miniprot band. The rescue quality
    gate is the protein-identity floor, not this band (the band is just a sanity
    bound against catastrophically mis-scaled miniprot hits)."""
    lo, hi = getattr(args, "miniprot_rescue_len", (0.5, 2.0))
    return lo < ratio < hi


def process_miniprot(
        mtrans, ref_db, m_feature_db, tree_dict, tgt_fai, ref_proteins,
        ref_trans, ref_features_dict, fw_score, m_id_2_ref_id_trans_dict,
        ref_features_len_dict, ref_trans_exon_num_dict,
        ref_features_reverse_dict, args, state_journal=None):
    if m_feature_db is None:
        return None
    mtrans_id = mtrans.attributes["ID"][0]
    mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
    is_overlapped = lifton_utils.check_ovps_ratio(mtrans, mtrans_interval, args.overlap, tree_dict)
    lifton_gene = None
    lifton_trans = None
    if not is_overlapped:
        ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_miniprot(ref_features_reverse_dict, mtrans_id, m_id_2_ref_id_trans_dict)
        if lifton_utils.record_unresolved_miniprot_hit(
                ref_gene_id, ref_trans_id, mtrans_id):
            return None
        if ref_trans_id not in ref_proteins or ref_trans_id not in ref_trans:
            drop_ledger.record("reference_protein_sequence", ref_trans_id)
            return None
        # Check if the additional copy is valid:
        # 1. Remove processed pseudogenes: one miniprot CDS but >1 ref exon.
        # 2. Check the transcript ratio against the default miniprot band.
        if (len(list(m_feature_db.children(mtrans, featuretype='CDS'))) == 1
                and ref_trans_exon_num_dict.get(ref_trans_id, 0) > 1):
            return None
        ref_feature_len = lifton_utils.reference_length_or_none(
            ref_features_len_dict, ref_gene_id)
        if not ref_feature_len:
            return None
        miniprot_trans_ratio = (
            (mtrans.end - mtrans.start + 1) / ref_feature_len
        )
        if (miniprot_trans_ratio > args.min_miniprot
                and miniprot_trans_ratio < args.max_miniprot):
            lifton_gene, lifton_trans, transcript_id, lifton_status = \
                lifton_miniprot_with_ref_protein(
                    mtrans, m_feature_db, ref_db.db_connection,
                    ref_gene_id, ref_trans_id, tgt_fai, ref_proteins,
                    ref_trans, tree_dict, ref_features_dict, args,
                    state_journal=state_journal,
                )
            lifton_gene.transcripts[transcript_id].entry.attributes[
                "miniprot_annotation_ratio"
            ] = [f"{miniprot_trans_ratio:.3f}"]
        else:  # Invalid default miniprot transcript; separate rescue may add it.
            # The separate post-Step-8 pass keeps this decision independent of
            # rescue mode and preserves off ⊆ on.
            return None
        lifton_trans_aln, lifton_aa_aln = lifton_gene.orf_search_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)
        # miniprot's CDS ends at the last aligned codon, so the model lacks the
        # stop the reference protein carries. Complete it AFTER the ORF search,
        # which therefore sees exactly the sequence it would have seen.
        orf_completion.complete_and_rescore(
            lifton_trans, mtrans, tgt_fai, ref_proteins, ref_trans_id,
            lifton_status, args)
        lifton_utils.print_lifton_status(transcript_id, mtrans, lifton_status, DEBUG=args.debug)
        lifton_gene.add_lifton_gene_status_attrs("miniprot")
        lifton_gene.add_lifton_trans_status_attrs(transcript_id, lifton_status)
        lifton_utils.write_lifton_status(fw_score, transcript_id, mtrans, lifton_status)
    return lifton_gene
