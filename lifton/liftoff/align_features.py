from multiprocessing import Pool
import gzip
import math
import os
import shlex
import sys
import tempfile
from functools import partial
import numpy as np
from pyfaidx import Fasta, Faidx
import subprocess
import pysam
from lifton.liftoff  import aligned_seg, liftoff_utils
from lifton.exceptions import LiftOnAlignmentError


MMI_MAGIC = b"MMI\x02"
LFS_POINTER_PREFIX = b"version https://git-lfs.github.com/spec/v1"


class _SamTargetMismatch(LiftOnAlignmentError):
    """A minimap2 SAM header does not describe the requested target FASTA."""


class _Minimap2CommandError(LiftOnAlignmentError):
    """A minimap2 subprocess could not be launched or exited unsuccessfully."""


def align_features_to_target(ref_chroms, target_chroms, args, feature_hierarchy, liftover_type, unmapped_features):
    # Phase 11 native path: route Liftoff's minimap2 alignment through the
    # in-process mappy binding instead of per-chromosome subprocesses.
    #
    # Iteration 7 (2026-06-14): this is now OPT-IN via
    # LIFTON_NATIVE_LIFTOFF_ALIGN=1, NOT implied by --native. The A/B
    # (benchmarks/compare/native_liftoff_ab.py) showed the in-process mappy
    # Liftoff alignment is STRICTLY INFERIOR to the subprocess path — slower
    # (a per-query Python loop vs minimap2's batched C threading: up to ~2.2×
    # on mouse_to_rat) and slightly less accurate cross-species (mappy's
    # Aligner can't take --end-bonus/-p). So under plain --native, Liftoff
    # alignment falls back to the proven subprocess path while --native still
    # delivers its real value elsewhere (in-process miniprot facade +
    # Step-7 threading unlock). The mappy Liftoff path stays reachable for
    # environments without a minimap2 binary or future batching work. The
    # polish subcommand always needs a real SAM file, so it never routes here.
    if (
        getattr(args, "native", False)
        and getattr(args, "subcommand", None) != "polish"
        and os.environ.get("LIFTON_NATIVE_LIFTOFF_ALIGN")
    ):
        from lifton.liftoff import native_align as _na
        return _na.align_features_to_target_native(
            ref_chroms, target_chroms, args,
            feature_hierarchy, liftover_type, unmapped_features,
        )
    if args.subcommand == "polish":
        sam_files = [args.directory + "/polish.sam"]
    else:
        target_fasta_dict = split_target_sequence(target_chroms, args.target, args.directory)
        target_lengths = {name: len(sequence) for name, sequence in target_fasta_dict.items()}
        genome_size = get_genome_size(target_fasta_dict)
        threads_per_alignment = max(1, math.floor(int(args.threads) / len(ref_chroms)))
        sam_files = []
        print("aligning features")
        func = partial(align_single_chroms, ref_chroms, target_chroms, threads_per_alignment, args, genome_size,
                       liftover_type, target_lengths=target_lengths)
        # A worker exception (including a checked minimap2 failure) must tear
        # down the remaining workers instead of leaking a Pool while the error
        # propagates to the CLI. Keep explicit close/join calls for compatibility
        # with Liftoff's lightweight Pool test doubles.
        pool = Pool(int(args.threads))
        try:
            for result in pool.imap_unordered(func, np.arange(0, len(target_chroms))):
                sam_files.append(result)
        except BaseException:
            pool.terminate()
            raise
        else:
            pool.close()
        finally:
            pool.join()
    return parse_all_sam_files(feature_hierarchy, unmapped_features, liftover_type, sam_files)


def split_target_sequence(target_chroms, target_fasta_name, inter_files):
    Faidx(target_fasta_name)
    target_fasta_dict = Fasta(target_fasta_name, key_function=lambda x: x.split()[0])
    for chrm in target_chroms:
        if chrm != target_fasta_name:
            with open(inter_files + "/" + chrm + ".fa", 'w') as out:
                out.write(">" + chrm + "\n" + str(target_fasta_dict[chrm]))
    return target_fasta_dict


def get_genome_size(target_fasta_dict):
    genome_size = 0
    for value in target_fasta_dict.values():
        genome_size += len(value)
    return genome_size


def align_single_chroms(ref_chroms, target_chroms, threads, args, genome_size, liftover_type, index,
                        target_lengths=None):
    max_single_index_size = 4000000000
    features_file, features_name = get_features_file(ref_chroms, args, liftover_type, index)
    target_file, output_file = get_target_file_and_output_file(liftover_type, target_chroms, index, features_name, args)
    threads_arg = str(threads)
    minimap2_path = get_minimap_path(args)
    target_prefix = get_target_prefix_name(target_chroms, index, args, liftover_type)
    if target_lengths is None:
        expected_lengths = _read_fasta_lengths(target_file)
    elif target_file == args.target:
        expected_lengths = target_lengths
    else:
        expected_lengths = {target_prefix: target_lengths[target_prefix]}
    if genome_size > max_single_index_size:
        split_prefix = args.directory + "/" + features_name + "_to_" + target_prefix + "_split"
        command = [minimap2_path] + _minimap2_options(args) + [
            "--split-prefix", split_prefix, '-t', threads_arg]
        _run_minimap2_to_sam(
            command, [target_file, features_file], output_file, target_file,
            expected_lengths, index_file=None,
        )
    else:
        index_was_reused = _find_reusable_minimap2_index(target_file, args, target_prefix) is not None
        minimap2_index = build_minimap2_index(
            target_file, args, threads_arg, minimap2_path, target_prefix=target_prefix,
        )
        for attempt in range(2):
            command = [minimap2_path] + _minimap2_options(args) + ['-t', threads_arg]
            try:
                _run_minimap2_to_sam(
                    command, [minimap2_index, features_file], output_file, target_file,
                    expected_lengths, index_file=minimap2_index,
                )
                break
            except (_SamTargetMismatch, _Minimap2CommandError):
                if attempt == 1 or not index_was_reused:
                    raise
                print(
                    "[LiftOn/Liftoff] WARNING: alignment with a cached minimap2 index failed for "
                    f"'{target_file}'; rebuilding the index once in '{args.directory}'.",
                    file=sys.stderr,
                )
                minimap2_index = build_minimap2_index(
                    target_file, args, threads_arg, minimap2_path,
                    target_prefix=target_prefix, force_rebuild=True,
                )
                index_was_reused = False
    return output_file


def get_features_file(ref_chroms, args, liftover_type, index):
    if ref_chroms[index] == args.reference and (liftover_type == "chrm_by_chrm" or liftover_type == "copies"):
        features_name = 'reference_all'
    elif liftover_type == "unmapped":
        features_name = "unmapped_to_expected_chrom"
    elif liftover_type == "unplaced":
        features_name = "unplaced"
    else:
        features_name = ref_chroms[index]
    return args.directory + "/" + features_name + "_genes.fa", features_name


def get_target_file_and_output_file(liftover_type, target_chroms, index, features_name, args):
    if liftover_type != "chrm_by_chrm" or target_chroms[0] == args.target:
        target_file = args.target
        out_file_target = "target_all"
    else:
        target_file = args.directory + "/" + target_chroms[index] + ".fa"
        out_file_target = target_chroms[index]
    output_file = args.directory + "/" + features_name + "_to_" + out_file_target + ".sam"
    return target_file, output_file


def get_minimap_path(args):
    if args.m is None:
        minimap2 = "minimap2"
    else:
        minimap2 = args.m
    return minimap2


def get_target_prefix_name(target_chroms, index, args, liftover_type):
    if liftover_type != "chrm_by_chrm" or target_chroms[0] == args.target:
        prefix = "target_all"
    else:
        prefix = target_chroms[index]
    return prefix


def classify_minimap2_index(index_file):
    """Classify a minimap2 cache without asking minimap2 to interpret it.

    Minimap2 treats an unresolved Git LFS pointer as an empty FASTA and exits
    successfully, so subprocess status is not sufficient for GH #57. The
    binary MMI magic is the reliable boundary between an index and text input.
    """
    if not os.path.exists(index_file):
        return "missing"
    try:
        with open(index_file, "rb") as handle:
            prefix = handle.read(max(len(LFS_POINTER_PREFIX), len(MMI_MAGIC)))
    except OSError:
        return "unreadable"
    if prefix.startswith(LFS_POINTER_PREFIX):
        return "lfs_pointer"
    if prefix.startswith(MMI_MAGIC):
        return "valid"
    return "invalid"


def _minimap2_options(args):
    try:
        return shlex.split(args.mm2_options)
    except ValueError as exc:
        raise LiftOnAlignmentError(f"Invalid -mm2_options value: {exc}") from exc


def _run_checked_minimap2(command, stage, target_file):
    try:
        run_kwargs = {"check": True}
        if stage == "index build":
            # The normal Liftoff options include -a. With `minimap2 -d`, that
            # otherwise emits a SAM header to stdout even though the intended
            # artifact is the binary index; progress remains visible on stderr.
            run_kwargs["stdout"] = subprocess.DEVNULL
        subprocess.run(command, **run_kwargs)
    except subprocess.CalledProcessError as exc:
        raise _Minimap2CommandError(
            f"minimap2 {stage} failed for target '{target_file}' with exit code "
            f"{exc.returncode}. Command: {shlex.join(command)}"
        ) from exc
    except (FileNotFoundError, PermissionError, NotADirectoryError, OSError) as exc:
        raise _Minimap2CommandError(
            f"Unable to run minimap2 during {stage} for target '{target_file}': {exc}. "
            f"Command: {shlex.join(command)}"
        ) from exc


def _warn_unusable_index(index_file, state, replacement):
    if state == "lfs_pointer":
        reason = "is an unresolved Git LFS pointer"
        recovery = " Run 'git lfs pull' to hydrate repository data if desired."
    elif state == "unreadable":
        reason = "cannot be read"
        recovery = ""
    elif state == "stale":
        reason = "is older than the target FASTA"
        recovery = ""
    else:
        reason = "is not a valid minimap2 MMI file"
        recovery = ""
    print(
        f"[LiftOn/Liftoff] WARNING: minimap2 cache '{index_file}' {reason}; "
        f"rebuilding from the target FASTA at '{replacement}'.{recovery}",
        file=sys.stderr,
    )


def _safe_index_prefix(target_prefix):
    safe = "".join(char if char.isalnum() or char in "._-" else "_" for char in str(target_prefix))
    return safe or "target_all"


def _local_minimap2_index(args, target_prefix):
    return os.path.join(args.directory, _safe_index_prefix(target_prefix) + ".mmi")


def _index_state_for_target(index_file, target_file):
    state = classify_minimap2_index(index_file)
    if state != "valid":
        return state
    try:
        if os.stat(index_file).st_mtime_ns < os.stat(target_file).st_mtime_ns:
            return "stale"
    except OSError:
        return "unreadable"
    return "valid"


def _find_reusable_minimap2_index(target_file, args, target_prefix):
    input_index = target_file + ".mmi"
    local_index = _local_minimap2_index(args, target_prefix)
    for candidate in (input_index, local_index):
        if _index_state_for_target(candidate, target_file) == "valid":
            return candidate
    return None


def build_minimap2_index(target_file, args, threads, minimap2_path, target_prefix="target_all",
                         force_rebuild=False):
    """Return a validated index, rebuilding invalid caches in the run directory.

    Valid input-adjacent indexes remain reusable for backward compatibility.
    Missing or unusable caches are never overwritten in the input tree; their
    replacements live with LiftOn's other intermediate artifacts.
    """
    os.makedirs(args.directory, exist_ok=True)
    local_index = _local_minimap2_index(args, target_prefix)
    input_index = target_file + ".mmi"

    if not force_rebuild:
        input_state = _index_state_for_target(input_index, target_file)
        if input_state == "valid":
            return input_index
        if input_state not in ("missing",) and os.path.abspath(input_index) != os.path.abspath(local_index):
            _warn_unusable_index(input_index, input_state, local_index)

        local_state = _index_state_for_target(local_index, target_file)
        if local_state == "valid":
            return local_index
        if local_state != "missing":
            _warn_unusable_index(local_index, local_state, local_index)

    fd, temporary_index = tempfile.mkstemp(
        prefix="." + os.path.basename(local_index) + ".", suffix=".tmp", dir=args.directory,
    )
    os.close(fd)
    try:
        # User options come first so LiftOn's owned output/thread arguments win
        # even if -mm2_options contains -d or -t.
        command = [minimap2_path] + _minimap2_options(args) + [
            '-t', str(threads), '-d', temporary_index, target_file,
        ]
        _run_checked_minimap2(command, "index build", target_file)
        built_state = classify_minimap2_index(temporary_index)
        if built_state != "valid":
            raise LiftOnAlignmentError(
                f"minimap2 reported a successful index build for '{target_file}', but "
                f"'{temporary_index}' is {built_state} instead of an MMI file. "
                "Check that the target FASTA is non-empty and readable."
            )
        try:
            os.replace(temporary_index, local_index)
        except OSError as exc:
            raise LiftOnAlignmentError(
                f"Unable to install rebuilt minimap2 index '{local_index}': {exc}"
            ) from exc
    finally:
        try:
            os.unlink(temporary_index)
        except FileNotFoundError:
            pass
    return local_index


def _read_fasta_lengths(target_file):
    """Read FASTA names/lengths without creating a sidecar in the input tree."""
    try:
        with open(target_file, "rb") as raw_handle:
            compressed = raw_handle.read(2) == b"\x1f\x8b"
        opener = gzip.open if compressed else open
        lengths = {}
        current_name = None
        with opener(target_file, "rt") as handle:
            for line in handle:
                if line.startswith(">"):
                    current_name = line[1:].strip().split()[0]
                    if not current_name:
                        raise ValueError("empty FASTA sequence name")
                    if current_name in lengths:
                        raise ValueError(f"duplicate FASTA sequence name '{current_name}'")
                    lengths[current_name] = 0
                elif current_name is not None:
                    lengths[current_name] += len("".join(line.split()))
    except (OSError, UnicodeError, ValueError) as exc:
        raise LiftOnAlignmentError(f"Unable to inspect target FASTA '{target_file}': {exc}") from exc
    if not lengths or any(length <= 0 for length in lengths.values()):
        raise LiftOnAlignmentError(
            f"Target FASTA '{target_file}' contains no non-empty target sequences."
        )
    return lengths


def validate_sam_target(sam_file, target_file, index_file=None, expected_lengths=None):
    """Require the SAM sequence dictionary to match the requested FASTA."""
    if expected_lengths is None:
        expected_lengths = _read_fasta_lengths(target_file)
    try:
        with pysam.AlignmentFile(sam_file, 'r', check_sq=False, check_header=False) as alignment_file:
            observed_lengths = dict(zip(alignment_file.references, alignment_file.lengths))
    except (OSError, ValueError) as exc:
        raise _SamTargetMismatch(
            f"minimap2 produced an unreadable SAM file '{sam_file}' for target '{target_file}': {exc}"
        ) from exc

    index_description = f" using index '{index_file}'" if index_file else ""
    if not observed_lengths:
        raise _SamTargetMismatch(
            f"minimap2 produced a SAM file with no @SQ target records for '{target_file}'"
            f"{index_description}. The target or cached index is empty, corrupt, incompatible, "
            "or an unresolved Git LFS pointer. Delete the cached .mmi or run 'git lfs pull', then retry."
        )
    if observed_lengths != expected_lengths:
        missing = sorted(set(expected_lengths) - set(observed_lengths))
        extra = sorted(set(observed_lengths) - set(expected_lengths))
        wrong_lengths = sorted(
            name for name in set(expected_lengths) & set(observed_lengths)
            if expected_lengths[name] != observed_lengths[name]
        )
        details = []
        if missing:
            details.append("missing=" + ",".join(missing[:5]))
        if extra:
            details.append("extra=" + ",".join(extra[:5]))
        if wrong_lengths:
            details.append("length_mismatch=" + ",".join(wrong_lengths[:5]))
        raise _SamTargetMismatch(
            f"minimap2 SAM target dictionary does not match FASTA '{target_file}'"
            f"{index_description} ({'; '.join(details)}). The cached index is stale or belongs "
            "to another assembly."
        )
    return None


def _run_minimap2_to_sam(command, positional_inputs, output_file, target_file, expected_lengths, index_file):
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    fd, temporary_sam = tempfile.mkstemp(
        prefix="." + os.path.basename(output_file) + ".", suffix=".tmp", dir=os.path.dirname(output_file) or ".",
    )
    os.close(fd)
    try:
        # Keep LiftOn's owned -o after user-supplied options so the temporary
        # artifact cannot be redirected elsewhere by -mm2_options.
        full_command = command + ['-o', temporary_sam] + positional_inputs
        _run_checked_minimap2(full_command, "alignment", target_file)
        validate_sam_target(
            temporary_sam, target_file, index_file=index_file, expected_lengths=expected_lengths,
        )
        try:
            os.replace(temporary_sam, output_file)
        except OSError as exc:
            raise LiftOnAlignmentError(
                f"Unable to install validated minimap2 SAM file '{output_file}': {exc}"
            ) from exc
    finally:
        try:
            os.unlink(temporary_sam)
        except FileNotFoundError:
            pass


def parse_all_sam_files(feature_hierarchy, unmapped_features, liftover_type, sam_files):
    aligned_segments_dict = {}
    for file in sam_files:
        aligned_segments = parse_alignment(file, feature_hierarchy, unmapped_features, liftover_type)
        aligned_segments_dict.update(aligned_segments)
    return aligned_segments_dict


def parse_alignment(file, feature_hierarchy, unmapped_features, search_type):
    all_aligned_blocks = {}
    aln_id = 0
    name_dict = {}
    align_count_dict = {}
    try:
        with pysam.AlignmentFile(file, 'r', check_sq=False, check_header=False) as sam_file:
            if not sam_file.references:
                raise LiftOnAlignmentError(
                    f"Cannot parse SAM file '{file}': it contains no @SQ target records. "
                    "The minimap2 target/index may be empty or corrupt."
                )
            # Preserve Liftoff's legacy iteration semantics; the context
            # manager adds safe closure without changing which records feed
            # downstream recovery.
            for ref_seq in sam_file.fetch():
                if ref_seq.is_unmapped is False:
                    aln_id = add_alignment(
                        ref_seq, align_count_dict, search_type, name_dict, aln_id,
                        feature_hierarchy, all_aligned_blocks,
                    )
                elif ref_seq.query_name in feature_hierarchy.parents:  # guard (GH #39)
                    unmapped_features.append(feature_hierarchy.parents[ref_seq.query_name])
    except LiftOnAlignmentError:
        raise
    except (OSError, ValueError) as exc:
        raise LiftOnAlignmentError(f"Unable to parse minimap2 SAM file '{file}': {exc}") from exc
    remove_alignments_without_children(all_aligned_blocks, unmapped_features, feature_hierarchy)
    return all_aligned_blocks


def add_alignment(ref_seq, align_count_dict, search_type, name_dict, aln_id, feature_hierarchy,
                  all_aligned_blocks):
    ref_seq.query_name = edit_name(search_type, ref_seq, name_dict)
    aln_id += 1
    if ref_seq.query_name in align_count_dict:
        align_count = align_count_dict[ref_seq.query_name] + 1
    else:
        align_count = 0
    align_count_dict[ref_seq.query_name] = align_count
    aligned_blocks = get_aligned_blocks(ref_seq, aln_id, feature_hierarchy, search_type)
    if ref_seq.query_name in all_aligned_blocks:
        all_aligned_blocks[ref_seq.query_name].extend(aligned_blocks)
    else:
        all_aligned_blocks[ref_seq.query_name] = aligned_blocks
    return aln_id


def edit_name(search_type, ref_seq, name_dict):
    if search_type != "copies":
        return ref_seq.query_name + "_0"
    else:
        if ref_seq.query_name not in name_dict:
            name_dict[ref_seq.query_name] = 0
        name_dict[ref_seq.query_name] += 1
        return ref_seq.query_name + "_" + str(name_dict[ref_seq.query_name])


_WARNED_MISSING_PARENT = set()


def _known_parent_key(feature_hierarchy, query_name):
    """Return ``convert_id_to_original(query_name)`` iff it names a known parent,
    else ``None`` (GitHub issue #39). Some RefSeq gene-family copy/fragment query
    names convert to an id that is not in ``feature_hierarchy.parents``; the old
    code indexed it unconditionally and a single such alignment KeyError-crashed
    the WHOLE chromosome. Returning None lets the caller skip just that alignment
    so every other feature on the chromosome still lifts. Warns once per key."""
    key = liftoff_utils.convert_id_to_original(query_name)
    if key in feature_hierarchy.parents:
        return key
    if key not in _WARNED_MISSING_PARENT:
        _WARNED_MISSING_PARENT.add(key)
        print("[LiftOn/Liftoff]   WARNING: aligned query '{0}' -> '{1}' is not a "
              "known parent feature; skipping it (GitHub issue #39).".format(
                  query_name, key), file=sys.stderr)
    return None


def get_aligned_blocks(alignment, aln_id, feature_hierarchy, search_type):
    cigar_operations = get_cigar_operations()
    cigar = alignment.cigar
    parent_key = _known_parent_key(feature_hierarchy, alignment.query_name)
    if parent_key is None:
        return []
    parent = feature_hierarchy.parents[parent_key]
    query_start, query_end = get_query_start_and_end(alignment, cigar, cigar_operations)
    children = feature_hierarchy.children[parent_key]
    end_to_end = is_end_to_end_alignment(parent, query_start, query_end)
    if search_type == "copies" and end_to_end is False:
        return []
    reference_block_start, reference_block_pos = alignment.reference_start, alignment.reference_start
    query_block_start, query_block_pos = query_start, query_start
    new_blocks, mismatches = [], []
    merged_children_coords = liftoff_utils.merge_children_intervals(children)
    for operation, length in cigar:
        if base_is_aligned(operation, cigar_operations):
            query_block_pos, reference_block_pos = add_aligned_base(operation, query_block_pos, reference_block_pos,
                                                                    length, cigar_operations, mismatches)
            if query_block_pos == query_end:
                add_block(query_block_pos, reference_block_pos, aln_id, alignment, query_block_start,
                          reference_block_start, mismatches, new_blocks, merged_children_coords, parent)
                break
        elif is_alignment_gap(operation, cigar_operations):
            add_block(query_block_pos, reference_block_pos, aln_id, alignment, query_block_start, reference_block_start,
                      mismatches, new_blocks, merged_children_coords, parent)
            mismatches, query_block_start, reference_block_start, query_block_pos, reference_block_pos = \
                end_block_at_gap(
                    operation, query_block_pos, reference_block_pos, length, cigar_operations)
    return new_blocks


def get_cigar_operations():
    return {"insertion": 1, "deletion": 2, "hard_clip": 5, "match": 7, "mismatch": 8}


def get_query_start_and_end(alignment, cigar, cigar_operations):
    query_start = alignment.query_alignment_start
    query_end = alignment.query_alignment_end
    if cigar[0][0] == cigar_operations["hard_clip"]:
        query_start += cigar[0][1]
        query_end += cigar[0][1]
    return query_start, query_end


def is_end_to_end_alignment(parent, query_start, query_end):
    return parent.end - parent.start + 1 == query_end - query_start


def base_is_aligned(operation, cigar_operations):
    return operation == cigar_operations["match"] or operation == cigar_operations["mismatch"]


def add_aligned_base(operation, query_block_pos, reference_block_pos, length, cigar_operations, mismatches):
    if operation == cigar_operations["mismatch"]:
        for i in range(query_block_pos, query_block_pos + length):
            mismatches.append(i)
    query_block_pos, reference_block_pos = adjust_position(operation, query_block_pos, reference_block_pos,
                                                           length, cigar_operations)
    return query_block_pos, reference_block_pos


def adjust_position(operation, query_block_pos, reference_block_pos, length, cigar_operations):
    if operation == cigar_operations["match"] or operation == cigar_operations["mismatch"] or operation == \
            cigar_operations["insertion"]:
        query_block_pos += length
    if operation == cigar_operations["match"] or operation == cigar_operations["mismatch"] or operation == \
            cigar_operations["deletion"]:
        reference_block_pos += length
    return query_block_pos, reference_block_pos


def add_block(query_block_pos, reference_block_pos, aln_id, alignment, query_block_start, reference_block_start,
              mismatches, new_blocks, merged_children_coords, parent):
    query_block_end = query_block_pos - 1
    reference_block_end = reference_block_pos - 1
    new_block = aligned_seg.aligned_seg(aln_id, alignment.query_name, alignment.reference_name, query_block_start,
                                        query_block_end,
                                        reference_block_start, reference_block_end, alignment.is_reverse,
                                        np.array(mismatches).astype(int))
    overlapping_children = find_overlapping_children(new_block, merged_children_coords, parent)
    if overlapping_children != []:

        new_blocks.append(new_block)




def find_overlapping_children(aln, children_coords, parent):
    overlapping_children = []
    for child_interval in children_coords:
        relative_start = liftoff_utils.get_relative_child_coord(parent, child_interval[0], aln.is_reverse)
        relative_end = liftoff_utils.get_relative_child_coord(parent, child_interval[1], aln.is_reverse)
        child_start, child_end = min(relative_start, relative_end), max(relative_start, relative_end)
        overlap = liftoff_utils.count_overlap(child_start, child_end, aln.query_block_start, aln.query_block_end)
        if overlap > 0:
            overlapping_children.append(child_start)
            overlapping_children.append(child_end)
    return overlapping_children


def is_alignment_gap(operation, cigar_operations):
    return operation == cigar_operations["insertion"] or operation == cigar_operations["deletion"]


def end_block_at_gap(operation, query_block_pos, reference_block_pos, length, cigar_operations):
    mismatches = []
    query_block_pos, reference_block_pos = adjust_position(operation, query_block_pos, reference_block_pos,
                                                           length, cigar_operations)
    query_block_start = query_block_pos
    reference_block_start = reference_block_pos
    return mismatches, query_block_start, reference_block_start, query_block_pos, reference_block_pos


def remove_alignments_without_children(all_aligned_blocks, unmapped_features, feature_hierarchy):
    features_to_remove = []
    for seq in all_aligned_blocks:
        if all_aligned_blocks[seq] == []:
            features_to_remove.append(seq)
            key = _known_parent_key(feature_hierarchy, seq)   # guard (GH #39)
            if key is not None:
                unmapped_features.append(feature_hierarchy.parents[key])
    for feature in features_to_remove:
        del all_aligned_blocks[feature]
    return all_aligned_blocks
