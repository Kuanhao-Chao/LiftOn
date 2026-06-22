from lifton import align, logger, lifton_class, lifton_utils
import subprocess, os, sys, copy
from intervaltree import Interval, IntervalTree


def _drain_stream_chunks(proc, *, chunk_size: int = 65536):
    """Phase 15c (V3.10) — drain a Popen's stdout/stderr in bounded
    chunks instead of `proc.communicate()`.

    `communicate()` allocates `bytes(stdout)` + `bytes(stderr)` and
    holds both simultaneously; for a multi-GB miniprot run that's an
    8-figure peak-RSS pathology. Chunked reads keep peak buffer at
    `chunk_size` per stream until both are drained, after which we
    `b"".join()` once. The final bytes object is unavoidable because
    gffbase.create_db(from_string=True) needs the full blob — but we
    avoid the doubled-allocation interim.

    Returns ``(stdout_bytes, stderr_bytes, returncode)``.
    """
    import threading

    stdout_chunks: list[bytes] = []
    stderr_chunks: list[bytes] = []

    def _drain(stream, sink):
        if stream is None:
            return
        try:
            while True:
                buf = stream.read(chunk_size)
                if not buf:
                    return
                sink.append(buf)
        finally:
            try:
                stream.close()
            except Exception:
                pass

    t_out = threading.Thread(
        target=_drain, args=(proc.stdout, stdout_chunks), daemon=True,
    )
    t_err = threading.Thread(
        target=_drain, args=(proc.stderr, stderr_chunks), daemon=True,
    )
    t_out.start()
    t_err.start()
    rc = proc.wait()
    t_out.join()
    t_err.join()
    return b"".join(stdout_chunks), b"".join(stderr_chunks), rc


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
                            mp_options, threads):
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
    return command


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
        subprocess.run(command)
        installed = True
    except (FileNotFoundError, PermissionError, NotADirectoryError,
            subprocess.SubprocessError):
        pass
    return installed


def run_miniprot(outdir, args, tgt_genome, ref_proteins_file):
    """
    Run miniprot and return the GFF3 result.

    Two return shapes are supported, controlled by ``args.stream``:

    * ``args.stream`` is False (default): write the GFF3 to disk at
      ``<outdir>/miniprot/miniprot.gff3`` and return the **path**.
      This is the legacy Phase 5 contract — preserved byte-for-byte.

    * ``args.stream`` is True: pipe miniprot's stdout straight into
      RAM via ``Popen.communicate()`` and return the **bytes blob**.
      No disk write. Phase 7 streaming-adapter fast path.

    Failure modes (any of them returns ``None`` so the pipeline
    continues with Liftoff-only results):
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
    str | bytes | None
        Path to the miniprot GFF3 output file, an in-memory bytes
        blob when streaming, or None on failure.
    """
    miniprot_outdir = outdir + "miniprot/"
    os.makedirs(miniprot_outdir, exist_ok=True)
    miniprot_output = miniprot_outdir + "miniprot.gff3"
    miniprot_path = "miniprot"
    command = _build_miniprot_command(
        miniprot_path, tgt_genome, ref_proteins_file,
        args.mp_options, getattr(args, "threads", 1),
    )

    stream_mode = bool(getattr(args, "stream", False))
    native_mode = bool(getattr(args, "native", False))
    try:
        if stream_mode and native_mode:
            # ── Phase 11 native + streaming path: route through the
            # MiniprotIndex facade. Today the facade still uses the
            # subprocess underneath, so the bytes are byte-identical
            # to the streaming branch (verified by
            # tests/test_native_bindings.py::TestMiniprotFacadeStreamingParity).
            # When a real `pyminiprot` PyO3 binding lands the facade
            # will route through it transparently — no caller change.
            from lifton.native_bindings import MiniprotIndex
            try:
                idx = MiniprotIndex(
                    tgt_genome,
                    mp_options=args.mp_options,
                    miniprot_path=miniprot_path,
                    ref_proteins_path=ref_proteins_file,
                    threads=getattr(args, "threads", 1),
                )
                bundle = idx.align_all()
                stdout_bytes = bundle.raw_bytes
                stderr_text = ""
                return_code = 0
                output_size = len(stdout_bytes)
            except RuntimeError as exc:
                # The facade re-raises miniprot's own ERROR / non-zero
                # exit conditions; surface them through the same code
                # path as the legacy branches so logging stays uniform.
                print(f"\n[LiftOn] miniprot (native facade) failed: {exc}",
                      file=sys.stderr)
                return None
        elif stream_mode:
            # ── Streaming branch: stdout -> RAM, bounded-chunk drain ──────
            # Phase 15c: replace proc.communicate() with chunked drain
            # so peak RAM is bounded by chunk_size per stream during
            # collection, not by 2× the full output.
            proc = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                bufsize=1 << 20,
            )
            stdout_bytes, stderr_bytes, return_code = \
                _drain_stream_chunks(proc, chunk_size=1 << 16)
            stderr_text = stderr_bytes.decode("utf-8", errors="replace") \
                if stderr_bytes else ""
            output_size = len(stdout_bytes) if stdout_bytes else 0
        else:
            # ── Legacy file-write branch (Phase 5 baseline behaviour) ─────
            with open(miniprot_output, "w") as fw:
                proc = subprocess.run(
                    command,
                    stdout=fw,
                    stderr=subprocess.PIPE,  # capture stderr so we can scan it
                    text=True,
                )
            stderr_text = proc.stderr or ""
            return_code = proc.returncode
            output_size = (os.path.getsize(miniprot_output)
                           if os.path.exists(miniprot_output) else 0)

        # Print stderr so the user sees miniprot's own log lines
        if stderr_text:
            print(stderr_text, end="", file=sys.stderr)

        # ── Failure mode 1: non-zero exit code ────────────────────────────
        if return_code != 0:
            print(
                f"\n[LiftOn] miniprot exited with code {return_code}. "
                "Miniprot output will be skipped — LiftOn will continue "
                "using Liftoff results only.",
                file=sys.stderr,
            )
            return None

        # ── Failure mode 2: miniprot printed ERROR on stderr ─────────────
        if stderr_text and "ERROR" in stderr_text.upper():
            print(
                "\n[LiftOn] miniprot reported an ERROR during mapping "
                "(exit code 0 but ERROR seen in output). "
                "Miniprot output will be skipped — LiftOn will continue "
                "using Liftoff results only.",
                file=sys.stderr,
            )
            return None

        # ── Failure mode 3: output is absent or empty ─────────────────────
        if output_size == 0:
            print(
                "\n[LiftOn] miniprot produced an empty output. "
                "Miniprot results will be skipped — LiftOn will continue "
                "using Liftoff results only.",
                file=sys.stderr,
            )
            return None

    except Exception as exc:
        print(
            f"\n[LiftOn] miniprot failed unexpectedly: {exc}\n"
            "Miniprot output will be skipped — LiftOn will continue "
            "using Liftoff results only.",
            file=sys.stderr,
        )
        return None

    if stream_mode:
        return stdout_bytes

    return miniprot_output


def lifton_miniprot_with_ref_protein(m_feature, m_feature_db, ref_db, ref_gene_id, ref_trans_id, tgt_fai, 
ref_proteins, ref_trans, tree_dict, ref_features_dict, args):
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
    m_gene_feature = copy.deepcopy(m_feature)
    m_gene_feature.featuretype = "gene"
    lifton_gene = lifton_class.Lifton_GENE(ref_gene_id, m_gene_feature, copy.deepcopy(ref_db[ref_gene_id].attributes), tree_dict, ref_features_dict, args)
    lifton_gene.update_gene_info(m_feature.seqid, m_feature.start, m_feature.end)
    # Create LifOn transcript instance
    Lifton_trans = lifton_gene.add_miniprot_transcript(ref_trans_id, copy.deepcopy(m_feature), ref_db[ref_trans_id].attributes, ref_features_dict)
    lifton_gene.update_trans_info(Lifton_trans.entry.id, m_feature.seqid, m_feature.start, m_feature.end)
    # Create exon / CDS entries
    cdss = m_feature_db.children(m_feature, featuretype='CDS')  # Replace 'exon' with the desired child feature type
    for cds in list(cdss):
        lifton_gene.add_exon(Lifton_trans.entry.id, cds)
        cds_copy = copy.deepcopy(cds)
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


def process_miniprot(mtrans, ref_db, m_feature_db, tree_dict, tgt_fai, ref_proteins, ref_trans, ref_features_dict, fw_score, m_id_2_ref_id_trans_dict, ref_features_len_dict, ref_trans_exon_num_dict, ref_features_reverse_dict, args):
    if m_feature_db is None:
        return None
    mtrans_id = mtrans.attributes["ID"][0]
    mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
    is_overlapped = lifton_utils.check_ovps_ratio(mtrans, mtrans_interval, args.overlap, tree_dict)
    lifton_gene = None
    lifton_trans = None
    if not is_overlapped:
        ref_trans_id = m_id_2_ref_id_trans_dict[mtrans_id]            
        ref_gene_id, ref_trans_id = lifton_utils.get_ref_ids_miniprot(ref_features_reverse_dict, mtrans_id, m_id_2_ref_id_trans_dict)
        if ref_gene_id is None:
            return None
        if ref_trans_id in ref_proteins.keys() and ref_trans_id in ref_trans.keys():
            if ref_trans_id != None:
                # Check if the additional copy is valid
                # 1. Remove processed pseudogenes: 1 CDS in miniprot but >1 CDS in reference
                # 2. Check the trans ratio is in the range of min_minprot and max_minprot
                if len(list(m_feature_db.children(mtrans, featuretype='CDS'))) == 1 and ref_trans_exon_num_dict[ref_trans_id] > 1: # Processed pseudogene
                    return None
                miniprot_trans_ratio = (mtrans.end - mtrans.start + 1) / ref_features_len_dict[ref_gene_id]
                if miniprot_trans_ratio > args.min_miniprot and miniprot_trans_ratio < args.max_miniprot:
                    lifton_gene, lifton_trans, transcript_id, lifton_status = lifton_miniprot_with_ref_protein(mtrans, m_feature_db, ref_db.db_connection, ref_gene_id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, tree_dict, ref_features_dict, args)
                    lifton_gene.transcripts[transcript_id].entry.attributes["miniprot_annotation_ratio"] = [f"{miniprot_trans_ratio:.3f}"]
                else: # Invalid miniprot transcript
                    # Iteration 23: miniprot-only rescue for genuinely-missing
                    # genes (ratio outside the default band) is no longer done
                    # in-loop -- it mutated the shared Step-8 suppression
                    # tree_dict mid-loop (off ⊄ on). It now runs as a separate
                    # post-Step-8 pass (lifton.miniprot_rescue) so this default
                    # decision is independent of --miniprot-rescue (off ⊆ on).
                    return None
            else: # Skip those cannot be found in reference.
                return None
        lifton_trans_aln, lifton_aa_aln = lifton_gene.orf_search_protein(lifton_trans.entry.id, ref_trans_id, tgt_fai, ref_proteins, ref_trans, lifton_status)
        lifton_utils.print_lifton_status(transcript_id, mtrans, lifton_status, DEBUG=args.debug)
        lifton_gene.add_lifton_gene_status_attrs("miniprot")
        lifton_gene.add_lifton_trans_status_attrs(transcript_id, lifton_status)
        lifton_utils.write_lifton_status(fw_score, transcript_id, mtrans, lifton_status)
    return lifton_gene