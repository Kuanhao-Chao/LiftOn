### 6.5 Vendored Liftoff internals — Part A: lift pipeline & mapping

This section specifies the vendored Liftoff subtree's DNA-level lift-over pipeline at reimplementation depth, covering four files: the top-level driver `liftoff_main.py`, the reference-feature extractor `extract_features.py`, the minimap2-driving aligner `align_features.py`, and the placement-scoring DP `find_best_mapping.py`. Liftoff is invoked **as a library** (not a subprocess) by LiftOn via `lifton/run_liftoff.py`, which calls `liftoff_main.run_all_liftoff_steps(args, ref_db)`. Liftoff itself shells out to `minimap2` (or, on the `--native` path, drives `mappy` in-process — see §6.5 Part B for `native_align.py`).

All line citations are against the four primary files plus the four supporting data-structure/utility modules: `aligned_seg.py`, `new_feature.py`, `feature_hierarchy.py`, `liftoff_utils.py`, and the calling module `liftover_types.py`.

#### 6.5.1 Supporting data structures

These small classes are the universal currency of the lift pipeline; every downstream algorithm consumes or emits them.

##### `aligned_seg` (`aligned_seg.py:1-12`)

A single contiguous aligned block (one ungapped run inside one minimap2 alignment record). Constructor signature:

```python
aligned_seg(aln_id, query_name, reference_name, query_block_start, query_block_end,
            reference_block_start, reference_block_end, is_reverse, mismatches)
```

| Field | Type | Meaning |
|---|---|---|
| `aln_id` | int (or str sentinel `"start"`/`"end"`) | Identifier shared by all blocks belonging to one minimap2 alignment record; distinct records get distinct ids. |
| `query_name` | str | Reference-feature id with copy-tag suffix, e.g. `gene1_0` (or `"start"`/`"end"` for sentinels). |
| `reference_name` | str | Target-genome seqid the block maps onto (a.k.a. the alignment "reference"; confusingly this is the *target* genome in LiftOn terms). |
| `query_block_start` | int | 0-based start of the block in the **gene-relative query** coordinate system (position within the extracted reference-gene sequence). |
| `query_block_end` | int | 0-based inclusive end of the block in query coordinates. |
| `reference_block_start` | int | 0-based start of the block in the target-genome coordinate (minimap2 `reference_start`). |
| `reference_block_end` | int | 0-based inclusive end of the block in target-genome coordinates. |
| `is_reverse` | bool / int | True if the alignment record is reverse-strand. Sentinels use `-1` (start) and `0` (end). |
| `mismatches` | `numpy.ndarray[int]` (or `[]`) | Query-coordinate positions of every mismatched base inside the block. |

Gotcha: the `reference_*` fields name the **target** genome, while `query_*` fields name positions inside the extracted reference gene. This inversion of "reference" vs "query" relative to LiftOn's own vocabulary persists throughout these four files.

##### `new_feature` (`new_feature.py:1-11`)

A lightweight, mutable replacement for a `gffutils.Feature`. Constructor: `new_feature(id, featuretype, seqid, source, strand, start, end, frame, attributes)`.

| Field | Type | Meaning |
|---|---|---|
| `id` | str | Feature ID. |
| `featuretype` | str | e.g. `gene`, `mRNA`, `exon`, `CDS`. |
| `seqid` | str | Sequence/chromosome name. |
| `source` | str | GFF column 2 (e.g. `Liftoff`). |
| `strand` | str | `+`, `-`, or `.`. |
| `start` | int | 1-based start (GFF convention). |
| `end` | int | 1-based inclusive end. |
| `frame` | str/int | CDS phase (GFF column 8). |
| `attributes` | dict | GFF column-9 attribute dict (values are lists). |

##### `feature_hierarchy` (`feature_hierarchy.py:1-5`)

Holds the three-level reference annotation decomposition. `feature_hierarchy(parents, intermediates, children)`:

| Field | Type | Meaning |
|---|---|---|
| `parents` | `{parent_id: new_feature}` | Top-level (highest) parents whose `featuretype` is in the lift list (default `["gene"]`). |
| `intermediates` | `{feature_id: new_feature}` | Features that are both a parent and a child (e.g. `mRNA` under `gene`). |
| `children` | `{parent_id: [new_feature,...]}` | Lowest-level children (e.g. exon/CDS) keyed by their *top-level* parent id. |

#### 6.5.2 `run_all_liftoff_steps` / `run_all_liftoff_steps_inmemory` — pipeline driver

Both entrypoints delegate to a single shared body `_run_liftoff_pipeline(args, ref_db, *, polish_intermediate_write)` (`liftoff_main.py:60-108`). The only difference is whether the final on-disk `write_new_gff` fires.

- `run_all_liftoff_steps(args, ref_db)` (`liftoff_main.py:10-23`): calls `_run_liftoff_pipeline(..., polish_intermediate_write=True)`, then **always** writes the final GFF3 to `args.output` via `write_new_gff.write_new_gff(lifted_feature_list, args, feature_db)`. This is LiftOn's Phase-5 baseline / `-L` disk path.
- `run_all_liftoff_steps_inmemory(args, ref_db)` (`liftoff_main.py:26-57`): calls `_run_liftoff_pipeline(..., polish_intermediate_write=False)` and returns the 4-tuple `(lifted_feature_list, feature_db, ref_parent_order, unmapped_features)` **without** writing to `args.output`. Used by LiftOn's `--inmemory-liftoff` path; the caller serialises `lifted_feature_list` to GFF3 bytes itself.

Return artefacts (both paths produce these four objects):

| Artefact | Type | Meaning |
|---|---|---|
| `lifted_feature_list` | `{copy_id: [list of lifted child new_feature]}` | The mapped output, keyed by per-copy gene id (e.g. `gene1_0`). Each list begins with the top-level parent. |
| `feature_db` | `gffutils.FeatureDB` | Reference DB; carries the `dialect` the GFF emitter uses to pick gff3 vs gtf output. |
| `ref_parent_order` | `numpy.ndarray` of `[id, feature]` rows | Parents sorted by `(seqid, start)` (`liftoff_utils.find_parent_order`, `:49-51`). |
| `unmapped_features` | `list[new_feature]` | Features that failed to map (already written to `args.u` inside the pipeline). |

##### `_run_liftoff_pipeline` — the exact step sequence (`liftoff_main.py:60-108`)

1. **Resolve chromosome lists** (`:66-70`). If `args.chroms is not None`, call `parse_chrm_files(args.chroms)` (`:244-253`) → `(ref_chroms, target_chroms)` parsed from a comma-separated 2-column file (column 0 = ref chrom, column 1 = target chrom if present). Else `ref_chroms = [args.reference]` and `target_chroms = [args.target]` (i.e. the whole-genome FASTA paths act as sentinel chrom names).
2. **Determine parent feature types** (`:71`): `parent_features_to_lift = get_parent_features_to_lift(args.features)` (`:256-262`). Always seeds `["gene"]`; appends one type per non-empty line of the `-f` file if given.
3. **Initialise outputs** (`:72-73`): `lifted_feature_list = {}`, `unmapped_features = []`.
4. **Lift the original annotation** (`:74-79`): call `liftover_types.lift_original_annotation(ref_chroms, target_chroms, lifted_feature_list, args, unmapped_features, parent_features_to_lift, ref_db)`. This is the main pass — it extracts ref features, aligns, lifts (via `find_best_mapping` → `lift_features`), and fixes overlaps. Returns `(feature_db, feature_hierarchy, ref_parent_order)`. Details in §6.5.6.
5. **Recover unmapped features against whole genome** (`:81-84`): `map_unmapped_features(...)` (`:265-274`). Only fires when `len(unmapped_features) > 0 and target_chroms[0] != args.target` (i.e. we were doing per-chromosome restriction and some features fell off). Re-runs the lift with `ref_chroms=[args.reference]`, `target_chroms=[args.target]`, liftover_type `"unmapped"`.
6. **Map from unplaced sequences** (`:85-88`): `map_features_from_unplaced_seq(...)` (`:277-284`). Fires only when `args.unplaced is not None and args.chroms is not None`; liftover_type `"unplaced"`.
7. **Write unmapped-features file** (`:89`): `write_unmapped_features_file(args.u, unmapped_features)` (`:287-291`) — one `gene.id` per line. (Under LiftOn defaults `args.u = "unmapped_features.txt"`.)
8. **Map extra copies** (`:90`): `map_extra_copies(...)` (`:294-300`). Fires only when `args.copies` is set; liftover_type `"copies"`, with `min_cov=0`, `min_seqid=args.sc` (default `1.0`).
9. **CDS status annotation** (`:92-93`): if `args.cds and args.polish is False` → `check_cds(lifted_feature_list, feature_hierarchy, args)` (`:339-345`). Annotates each CDS's status (partial / missing start / missing stop / inframe stop). `args.cds` defaults **True** (`:219`).
10. **Polish branch** (`:94-107`): if `args.polish` (default False, `:218`):
    a. print `"polishing annotations"`; run `check_cds(...)`.
    b. if `polish_intermediate_write` (legacy disk path only): write the pre-polish GFF3 to disk so polish can consume it.
    c. `find_and_polish_broken_cds(...)` (`:303-336`) — re-aligns each lifted feature with `args.d` forced to `100000000` (effectively disables the distance-scaling cutoff so any chain is allowed), re-lifts into `polish_lifted_features`, re-runs `check_cds`, then **replaces** the original lift with the polished one when the polished version has strictly more `valid_ORFs`, or equal `valid_ORFs` with strictly higher `sequence_ID`, or (else) strictly higher `coverage` (`:327-336`).
    d. if `args.output != 'stdout'`: append `"_polished"` to `args.output`.
11. **Return** `(lifted_feature_list, feature_db, ref_parent_order, unmapped_features)`.

Gotcha (LiftOn integration): LiftOn always invokes Liftoff with `args.polish=False` and `args.copies=False` in the default golden path, so steps 5–8 and 10 are typically no-ops; the load-bearing work is steps 1–4 + 9.

##### `parse_args` and minimap2-option normalisation (`liftoff_main.py:112-241`)

Liftoff defines its own argparse parser. The fields that govern this section's algorithms:

| Arg | Default | Meaning / use |
|---|---|---|
| `-mm2_options` | `'-a --end-bonus 5 --eqx -N 50 -p 0.5'` | Space-delimited minimap2 flags (see §6.5.4). |
| `-a` (coverage) | `0.5` | Min alignment coverage to call a feature mapped. |
| `-s` (seq-id) | `0.5` | Min child sequence identity to call a feature mapped. |
| `-d` (distance) | `2.0` | Distance scaling factor; chain edges rejected when actual target distance exceeds `D × expected` (§6.5.7). |
| `-flank` | `0` | Fraction of gene length to add as flanking sequence on extraction. |
| `-p`/`--threads` | `1` | Parallel alignment processes. |
| `-copies` | False | Look for extra gene copies. |
| `-sc` | `1.0` | Min identity to call a copy (must be ≥ `-s`). |
| `-overlap` | `0.1` | Max fraction of overlap allowed between two lifted features. |
| `-mismatch` | `2` | Per-mismatch node-weight penalty in the chaining DP. |
| `-gap_open` | `2` | Gap-open penalty in the chaining DP. |
| `-gap_extend` | `1` | Gap-extend penalty in the chaining DP. |
| `-cds` | True | Annotate CDS status. |
| `-polish` | False | Run the polish pass. |

After parse, the parser **idempotently re-injects defaults into `-mm2_options`** if the user removed them (`:227-236`): appends ` -a` if `-a` absent, ` --eqx` if absent, ` -N 50` if `-N` absent, ` -p 0.5` if `-p` absent, and `--end-bonus 5` if `--end-bonus` absent. Gotcha: the `--end-bonus` append at `:236` concatenates without a leading space (`+= "--end-bonus 5"`), so if it triggers immediately after another append it can fuse onto the prior token — this matters only when a user passes a custom `mm2_options` lacking `--end-bonus`; the default already contains it so the branch is normally skipped. Validation: `-sc` must be ≥ `-s` (`:237-238`); `-unplaced` requires `-chroms` (`:239-240`).

#### 6.5.3 `extract_features.py` — feature hierarchy & reference sequence extraction

Entry: `extract_features_to_lift(ref_chroms, liftover_type, parents_to_lift, ref_db, args)` (`:19-27`):

1. print `"extracting features"`.
2. create `args.directory` (default `intermediate_files`) if absent (`:21-22`).
3. `feature_db = ref_db.db_connection` (`:24`) — reuse the gffutils DB connection already opened by LiftOn (the local `build_database` machinery, `:59-241`, is the standalone-Liftoff fallback; in LiftOn the DB comes from `ref_db`).
4. `feature_hierarchy, parent_order = seperate_parents_and_children(feature_db, parents_to_lift)` (`:25`).
5. `get_gene_sequences(feature_hierarchy.parents, ref_chroms, args, liftover_type)` (`:26`) — writes the per-gene FASTA(s).
6. return `(feature_hierarchy, feature_db, parent_order)`.

##### `seperate_parents_and_children` (`:260-278`) — the hierarchy decomposition

This is a pure SQL-driven graph partition over gffutils' `features` and `relations` tables.

1. Query all parent→child relations where parent≠child (`:262-264`):
   `SELECT * FROM relations join features as a on a.id=relations.parent join features as b on b.id=relations.child` filtered to rows where `feature[0] != feature[1]` (drop self-relations).
2. `all_ids` = every feature id (`SELECT * FROM features`, `:265`).
3. `all_children_ids` = the child column of each relation; `all_parent_ids` = the parent column (`:266-267`).
4. **Classify each id by set difference** (`numpy.setdiff1d`):
   - `lowest_children = setdiff1d(all_ids, all_parent_ids)` — ids that are never a parent (`:268`).
   - `highest_parents = setdiff1d(all_ids, all_children_ids)` — ids that are never a child (`:269`).
   - `intermediates = set(all_children_ids) ∩ set(all_parent_ids)` — ids that are both (`:270`).
5. Populate dicts via `add_parents`, `add_children`, `add_intermediates` (`:272-274`).
6. `parent_order = liftoff_utils.find_parent_order([non-None parents])` (`:275-276`) → array of `[id, feature]` sorted by `(seqid, start)`.
7. return `feature_hierarchy(parent_dict, intermediate_dict, child_dict)` and `parent_order`.

`add_parents` (`:281-291`): for each highest-parent id, build a `new_feature` from the SQL row tuple — column mapping `(id=tup[0], featuretype=tup[3], seqid=tup[1], source=tup[2], strand=tup[7], start=tup[4], end=tup[5], frame=tup[8], attributes=json.loads(tup[9]))`. Only retain it if `parent.featuretype in parent_types_to_lift`; then `parent_dict[id]=parent` and seed `child_dict[id]=[]`.

`add_children` (`:294-317`): join relations to features for the lowest children; for each, the parent is `tup[0]` and the child is built from the joined feature columns (offset by the relation columns: `id=tup[3], featuretype=tup[6], seqid=tup[4], source=tup[5], strand=tup[10], start=tup[7], end=tup[8], frame=tup[11], attributes=json.loads(tup[12])`). Skip `intron` features (`:308`). If the child lacks a `Parent` attribute, call `add_parent_tag` (`:309-310`, `:321-330`) to backfill it from gffutils' `parents()`. Append to `child_dict[parent]`. Finally, any lowest-child that was never added but **is itself a parent** (a single-level feature with no children) becomes its own child: `child_dict[feature] = [parent_dict[feature]]` (`:313-316`).

`add_intermediates` (`:334-345`): build a `new_feature` for each intermediate id; backfill `Parent` via `add_parent_tag` if missing.

##### `get_gene_sequences` & FASTA writing (`:348-398`)

`get_gene_sequences(parent_dict, ref_chroms, args, liftover_type)`:
1. open the reference FASTA with `pyfaidx.Fasta(args.reference)`.
2. if `liftover_type == "unplaced"`, truncate `args.directory + "/unplaced_genes.fa"` first (it is later appended to).
3. for each `chrom` in `ref_chroms`: open the output FASTA via `get_fasta_out` (naming rules below), sort parents by `seqid`, `sys.exit` if there are zero parents (message: "GFF does not contain any gene features…"), then `write_gene_sequences_to_file`.

`get_fasta_out` (`:363-376`) FASTA naming: when chrom == the reference-FASTA sentinel name and liftover_type ∈ {`chrm_by_chrm`, `copies`} → `reference_all`; `unmapped`→`unmapped_to_expected_chrom`; `unplaced`→`unplaced` (opened in **append** mode); else the chrom name. File path is `{directory}/{name}_genes.fa`. These names must exactly match the consumer side in `align_features.get_features_file` (§6.5.4).

`write_gene_sequences_to_file` (`:379-398`) — the per-gene sequence extraction with flanking:
1. Determine `current_chrom`: if `chrom_name == reference_fasta_name` use `parents[0].seqid`, else `chrom_name`; load `chrom_seq = reference_fasta_idx[current_chrom][:].seq`.
2. For each `parent` (in seqid-sorted order):
   - skip if `parent.seqid not in reference_fasta_idx.keys()` (`:387-388`).
   - if the parent's seqid differs from `current_chrom` and we're in whole-genome mode, reload `chrom_seq` for the new chrom (`:390-392`).
   - proceed only if `parent.seqid == chrom_name` or whole-genome mode (`:393`).
   - **Flank expansion** (`:394-396`): `gene_length = parent.end - parent.start + 1`; mutate `parent.start = round(max(1, parent.start - args.flank * gene_length))` and `parent.end = round(min(parent.end + args.flank * gene_length, len(chrom_seq)))`. With default `flank=0` the coordinates are unchanged.
   - extract `parent_seq = chrom_seq[parent.start - 1 : parent.end]` (1-based→0-based slice) and write `>{parent.id}\n{parent_seq}\n`.

Gotcha: `parent.start`/`parent.end` are **mutated in place** by the flank expansion. All downstream gene-relative coordinate math (`get_relative_child_coord`) uses these possibly-expanded parent bounds, so the query coordinate system is "position within the (possibly flank-expanded) extracted gene sequence." With `flank=0` (LiftOn default) this is a no-op and the query system is exactly the gene span.

#### 6.5.4 `align_features.py` — minimap2 invocation, SAM parsing, block model

Entry: `align_features_to_target(ref_chroms, target_chroms, args, feature_hierarchy, liftover_type, unmapped_features)` (`:12-43`):

1. **Native dispatch** (`:19-27`): if `getattr(args,"native",False)` is True **and** `getattr(args,"subcommand",None) != "polish"`, delegate to `native_align.align_features_to_target_native(...)` (the `mappy` in-process path — §6.5 Part B) and return. The legacy subprocess path below is the default and the fallback.
2. **Polish special case** (`:28-29`): if `args.subcommand == "polish"`, the SAM is the single pre-existing `args.directory + "/polish.sam"`.
3. **Normal path** (`:30-42`):
   a. `split_target_sequence(target_chroms, args.target, args.directory)` (`:46-53`) — build a `.fai` via `Faidx`, open `Fasta(... key_function=lambda x: x.split()[0])`, and for each target chrom that isn't the whole-FASTA sentinel write `{directory}/{chrm}.fa`. Returns the `Fasta` dict.
   b. `genome_size = get_genome_size(target_fasta_dict)` (`:56-60`) — sum of all sequence lengths.
   c. `threads_per_alignment = max(1, floor(int(args.threads) / len(ref_chroms)))` (`:33`).
   d. **Parallel alignment** (`:35-42`): `Pool(int(args.threads))`, `pool.imap_unordered(func, arange(0, len(target_chroms)))` where `func = partial(align_single_chroms, ref_chroms, target_chroms, threads_per_alignment, args, genome_size, liftover_type)`. Each result is a SAM path; collect into `sam_files`.
4. `return parse_all_sam_files(feature_hierarchy, unmapped_features, liftover_type, sam_files)`.

##### `align_single_chroms` — the exact minimap2 argv (`:63-80`)

1. `max_single_index_size = 4_000_000_000` (`:64`).
2. `features_file, features_name = get_features_file(...)` (`:83-92`) — the reference-gene FASTA written by §6.5.3; naming mirror of `get_fasta_out`: ref-FASTA+`chrm_by_chrm`/`copies`→`reference_all`; `unmapped`→`unmapped_to_expected_chrom`; `unplaced`→`unplaced`; else chrom name; path `{directory}/{name}_genes.fa`.
3. `target_file, output_file = get_target_file_and_output_file(...)` (`:95-103`): if liftover_type≠`chrm_by_chrm` or whole-genome, `target_file=args.target` and out-target name `target_all`; else per-chrom `{directory}/{chrom}.fa`. `output_file = {directory}/{features_name}_to_{out_file_target}.sam`.
4. `minimap2_path = get_minimap_path(args)` (`:106-111`) — `args.m` if set, else literal `"minimap2"` (resolved from `PATH`).
5. **Branch on genome size** (`:70-79`):
   - **Large genome** (`genome_size > 4e9`): use minimap2's split-index mode. argv:
     `[minimap2_path, '-o', output_file, target_file, features_file] + args.mm2_options.split(" ") + ['--split-prefix', split_prefix, '-t', threads_arg]` where `split_prefix = {directory}/{features_name}_to_{target_prefix}_split`.
   - **Normal genome**: first build/reuse a minimap2 index via `build_minimap2_index` (`:122-127`) — runs `[minimap2_path, '-d', target_file+'.mmi', target_file] + mm2_options + ['-t', threads]` only if the `.mmi` doesn't already exist; returns `target_file + ".mmi"`. Then align:
     `[minimap2_path, '-o', output_file, minimap2_index, features_file] + args.mm2_options.split(" ") + ['-t', threads_arg]`.
6. `subprocess.run(command)`; return `output_file`.

The **effective default minimap2 flag set** (from `-mm2_options` after normalisation) is: `-a --end-bonus 5 --eqx -N 50 -p 0.5 -t {threads}`. Meaning: `-a` = SAM output (required so pysam can read it); `--eqx` = emit `=`/`X` CIGAR ops (required — the parser distinguishes match op 7 from mismatch op 8); `--end-bonus 5` = soft bonus for aligning to ends; `-N 50` = report up to 50 secondary alignments per query (so multiple placements feed `find_best_mapping`); `-p 0.5` = report secondaries scoring ≥ 0.5 of the primary.

Gotcha (byte-identity): `--eqx` is structurally required. Without it, minimap2 emits `M` ops (a single op covering both matches and mismatches), but `get_cigar_operations` only recognises match=7/mismatch=8, so block-splitting and mismatch counting would silently break. The normaliser at `liftoff_main.py:229-230` re-adds `--eqx` defensively.

##### SAM parsing pipeline (`:130-302`)

`parse_all_sam_files` (`:130-135`): for each SAM file, call `parse_alignment` and merge the resulting `{query_name: [aligned_seg,...]}` dicts.

`parse_alignment(file, feature_hierarchy, unmapped_features, search_type)` (`:138-152`):
1. open with `pysam.AlignmentFile(file, 'r', check_sq=False, check_header=False)`; iterate `.fetch()`.
2. For each record `ref_seq`: if `ref_seq.is_unmapped is False`, call `add_alignment(...)`; else append `feature_hierarchy.parents[ref_seq.query_name]` to `unmapped_features`.
3. `remove_alignments_without_children(...)` (`:294-302`) — any query whose block list ended up empty is moved to `unmapped_features` and deleted from the dict.
4. return `all_aligned_blocks` (`{query_name: [aligned_seg,...]}`).

`add_alignment` (`:155-169`):
1. **Rename query** via `edit_name(search_type, ref_seq, name_dict)` (`:172-179`): for non-`copies` searches the query name becomes `{original}_0` (copy-tag `_0`); for `copies` it becomes `{original}_{N}` with a per-name incrementing counter starting at 1. This copy-tag suffix is later stripped by `liftoff_utils.convert_id_to_original`.
2. increment a global `aln_id` (one per SAM record, `:158`).
3. track per-query alignment count in `align_count_dict`.
4. `aligned_blocks = get_aligned_blocks(ref_seq, aln_id, feature_hierarchy, search_type)`.
5. extend/insert `all_aligned_blocks[query_name]`.

##### `get_aligned_blocks` — CIGAR walk into ungapped blocks (`:182-209`)

CIGAR op codes (`get_cigar_operations`, `:212-213`): `insertion=1, deletion=2, hard_clip=5, match=7, mismatch=8`.

1. `parent = feature_hierarchy.parents[convert_id_to_original(query_name)]` (`:185`).
2. `query_start, query_end = get_query_start_and_end(...)` (`:216-222`): start from pysam `query_alignment_start`/`query_alignment_end`; **if the first CIGAR op is a hard clip, add its length to both** (so coordinates are in full-query space).
3. `children = feature_hierarchy.children[convert_id_to_original(query_name)]`.
4. `end_to_end = is_end_to_end_alignment(parent, query_start, query_end)` = `parent.end - parent.start + 1 == query_end - query_start` (`:225-226`). For `search_type == "copies"`, if **not** end-to-end, return `[]` (copies must align the whole gene).
5. `merged_children_coords = liftoff_utils.merge_children_intervals(children)` (`:194`; `liftoff_utils.py:17-29`) — sort child `[start,end]` intervals and merge overlapping/abutting (`current[0] <= previous[1]`) ones.
6. **Walk the CIGAR** building blocks. Maintain `query_block_pos`/`reference_block_pos` (running positions, init at the alignment starts) and `query_block_start`/`reference_block_start` (current block origin) and a `mismatches` list:
   - For an **aligned op** (match or mismatch, `base_is_aligned`, `:229-230`): `add_aligned_base` (`:233-239`) — if mismatch, append every position `query_block_pos .. query_block_pos+length-1` to `mismatches`; then advance positions via `adjust_position` (`:242-249`): query advances on match/mismatch/insertion; reference advances on match/mismatch/deletion. If `query_block_pos == query_end`, emit the final block via `add_block` and `break`.
   - For a **gap op** (insertion or deletion, `is_alignment_gap`, `:281-282`): emit the current block via `add_block`, then `end_block_at_gap` (`:285-291`) resets `mismatches=[]`, advances over the gap, and sets the new block start. Gaps thus split blocks.
7. return `new_blocks`.

`add_block` (`:252-264`): set `query_block_end = query_block_pos - 1`, `reference_block_end = reference_block_pos - 1`; construct an `aligned_seg(aln_id, query_name, reference_name, query_block_start, query_block_end, reference_block_start, reference_block_end, is_reverse, np.array(mismatches).astype(int))`. **Only append the block if it overlaps at least one child** (`find_overlapping_children` returns non-empty, `:260-263`). `find_overlapping_children` (`:268-278`) converts each merged child interval to query-relative coordinates via `get_relative_child_coord` (reverse-aware: `parent.end - coord` if reverse, else `coord - parent.start`; `liftoff_utils.py:9-14`) and tests `count_overlap > 0` against the block's query span.

Gotcha: blocks that fall entirely in intronic/intergenic flank (no child overlap) are silently dropped, which is why `remove_alignments_without_children` can leave a query with zero blocks → unmapped.

#### 6.5.5 `find_best_mapping.py` — best-placement selection (chaining DP)

This module turns the set of `aligned_seg` blocks for one query into a single best placement by building a DAG of blocks, weighting nodes/edges by mismatch and gap penalties, and finding the minimum-weight source→sink path (a shortest-path DP, **not** greedy). It then converts the chosen chain back into lifted child coordinates and computes coverage + sequence identity.

Entry: `find_best_mapping(alignments, query_length, parent, feature_heirarchy, previous_feature_start, previous_feature_ref_start, previous_gene_seq, inter, lifted_features_list, args)` (`:6-25`):

1. `children = feature_heirarchy.children[parent.id]`; `children_coords = merge_children_intervals(children)`.
2. `node_dict, aln_graph = intialize_graph()` (`:28-33`): a `networkx.DiGraph`; node `0` is a sentinel `aligned_seg("start","start","start",-1,-1,-1,-1,-1,[])`.
3. `head_nodes = add_single_alignments(...)` — add one graph node per block, chained intra-alignment (§below).
4. `chain_alignments(head_nodes, ...)` — add inter-alignment edges (§below).
5. `add_target_node(...)` — append the sink node and connect all terminal nodes to it.
6. `shortest_path_nodes = find_shortest_path(node_dict, aln_graph)`. If empty, return `({}, 0, 0)`.
7. `mapped_children, alignment_coverage, seq_id = convert_all_children_coords(shortest_path_nodes, children, parent)`.
8. return `(mapped_children, alignment_coverage, seq_id)`.

##### `add_single_alignments` — node creation + intra-alignment edges (`:36-62`)

1. If `inter is not None`, drop whole alignments that overlap already-placed features: `remove_alignments_with_overlap(...)` (`:81-95`) groups blocks by `aln_id`, computes the alignment's `[min_ref_start, max_ref_end]` span and strand (`get_strand`, `:109-117`, flips when `is_reverse`), and keeps the group only if `liftoff_utils.find_overlaps(...)` returns empty (no disallowed overlap with prior lifts).
2. Sort blocks by `query_block_start`, then re-sort with `sort_alignments` (`:65-78`): primary key prefers same-chromosome placements and minimal `distance_difference`; `distance_difference` (`:102-106`) = `|expected_distance - actual_distance|` where `expected = feature_ref_start - previous_feature_ref_start` and `actual = feature_target_start - query_start - previous_feature_target_start`. This biases placement toward the location consistent with the previously-lifted neighbour. A stable secondary sort orders blocks within each alignment by `query_block_start`.
3. Walk the sorted blocks, building nodes (`:50-61`):
   - When `aln.aln_id` changes, reset `previous_node = 0` (start fresh chain at the source) and update `previous_node_id`.
   - `is_valid_alignment(previous_node, aln, ...)` (`:120-127`): reject if `aln.reference_block_start == -1`; if `previous_node != 0`, reject if `spans_overlap_region(...)` (the gap between the previous block and this one overlaps a disallowed prior lift).
   - If valid: if `previous_node == 0` this block starts a new chain → record `node_num` in `head_nodes`. `add_to_graph(...)` adds the node with its weight and an edge from `previous_node` with the edge cost. Advance `previous_node = node_num`, `node_num += 1`.
4. return `head_nodes` (the per-alignment chain heads).

##### Node and edge weights (the scoring model)

`get_node_weight(aln, children_coords, parent, args)` (`:161-168`): sum over child intervals of `(# mismatches inside the child's query span) × args.mismatch`. Only mismatches that fall inside a child (exon/CDS) region count; intronic mismatches are free. Default `mismatch=2`.

`get_edge_weight(from_node, to_node, children_coords, parent, args)` (`:171-193`): penalises the **unaligned exon bases** in the query gap between two chained blocks.
1. `node_overlap = get_node_overlap(from_node, to_node)` (`:144-148`) = `max(0, from_node.query_block_end - to_node.query_block_start + 1)` (0 if either node is a sentinel) — the query overlap to be trimmed.
2. `unaligned_range = [from_node.query_block_end + 1, to_node.query_block_start + node_overlap - 1]`.
3. For each child interval, convert to query-relative coords (reverse-aware; uses `to_node.is_reverse` when `from_node` is the start sentinel, else `from_node.is_reverse`) and `overlap = count_overlap(child_span, unaligned_range)`.
4. **Deletion special case** (`:185-188`): if `overlap == 1` **and** `unaligned_range[0] == unaligned_range[1] + 1` (a zero-width gap) **and** same chromosome, the penalty is the *reference-side* gap (a genomic deletion) scaled by `gap_extend`: `((to.reference_block_start + node_overlap) - from.reference_block_end - 1) * args.gap_extend`.
5. Else `unaligned_exon_bases += max(0, overlap) * args.gap_extend` (`:190`).
6. If any unaligned exon bases accrued, add the gap-open premium once: `+= (args.gap_open - args.gap_extend)` (`:191-192`). Defaults `gap_open=2`, `gap_extend=1`.

##### `chain_alignments` — inter-alignment edges (`:196-235`)

For every head node, `add_edges` tries to draw an edge from every other node to it, gated by `is_valid_edge(from, to, ...)` (`:214-235`):

| # | Condition (reject the edge if…) | Line |
|---|---|---|
| 1 | `from_node.aln_id == to_node.aln_id` (same alignment — already chained) | `:217-218` |
| 2 | `from_node.query_block_end >= to_node.query_block_end` (no forward progress in query) | `:219-220` |
| 3 | `from_node.is_reverse != to_node.is_reverse` (strand mismatch) | `:221-222` |
| 4 | `from_node.reference_name != to_node.reference_name` (different target chrom) | `:223-224` |
| 5 | `to_node.reference_block_start < from_node.reference_block_end` (target order violated / would imply overlap) | `:228-229` |
| 6 | `actual_distance > args.d * expected_distance` where `expected = to.query_block_end - from.query_block_start`, `actual = to.reference_block_end - from.reference_block_start` (target span too dilated vs query span) | `:226-231` |
| 7 | `spans_overlap_region(from, to, ...)` (the inter-block gap hits a disallowed prior lift) | `:232-234` |

Otherwise the edge is added with `cost = get_edge_weight(...)`. Constant `args.d` defaults `2.0`; in polish mode it is forced to `100000000` so condition 6 never fires.

##### `add_target_node` & shortest path (`:238-272`)

`add_target_node` (`:238-247`): create sink node `max(node_dict)+1` as `aligned_seg("end","end","end", query_length, query_length, query_length, query_length, 0, [])`; for every **terminal** node (no successors, `is_terminal_node`, `:250-254`) other than the sink, add an edge to the sink with `get_edge_weight` cost.

`find_shortest_path` (`:257-265`): `nx.shortest_path(aln_graph, source=0, target=len(node_dict)-1, weight=<lambda>)`. The combined weight function `get_weight(u, v, d, G)` (`:268-272`) is **half of each endpoint node's weight plus the edge cost**: `node_u_wt/2 + node_v_wt/2 + edge_wt` (edge default cost 1 if missing). Because every interior node is entered once and exited once, the half-weights sum to one full node weight per traversed interior node — so the total path cost = Σ(interior node weights) + Σ(edge costs). Then strip the source/sink (keep indices `1 .. len-2`) into `shortest_path_nodes` and `trim_path_boundaries` (`:275-281`): for consecutive nodes, advance the later node's `query_block_start` and `reference_block_start` by `node_overlap` to remove the double-counted overlap region.

Gotcha (order/determinism): `nx.shortest_path` is a Dijkstra-style search; ties between equal-cost paths are broken by node insertion order, which is governed by the deterministic `sort_alignments` ordering. Any change to block ordering, weight arithmetic, or the half-weight convention can alter which equal-cost chain wins and break LiftOn's byte-identity matrix. The sink/source are at fixed ids `0` and `len(node_dict)-1`.

##### `convert_all_children_coords` — lift coords + coverage + identity (`:284-311`)

For each child (exon/CDS), map its reference-gene coordinates onto the target genome through the chosen chain, then aggregate alignment statistics.

1. `shortest_path_nodes.sort(key=query_block_start)`.
2. `total_bases, mismatches, insertions, deletions, matches = 0,0,0,0,0`.
3. For each child:
   - `total_bases += child.end - child.start + 1`.
   - `nearest_start_coord, nearest_end_coord, relative_start, relative_end = find_nearest_aligned_start_and_end(...)` (`:314-320`) — convert child start/end to gene-relative coords (reverse-aware) then snap to the nearest aligned position inside the chain (`find_nearest_aligned_start` `:323-333`, `find_nearest_aligned_end` `:336-349`). Unaligned overhang at the child's ends becomes deletions.
   - If both nearest coords found (`!= -1`):
     - `lifted_start, start_node = convert_coord(nearest_start_coord, ...)`; `lifted_end, end_node = convert_coord(nearest_end_coord, ...)` (`:352-361`): linear lift inside the containing node: `lifted = node.reference_block_start + (relative_coord - node.query_block_start)`.
     - `deletions += find_deletions(start_node, end_node, ...)` (`:364-370`, query gaps between consecutive nodes) `+ (nearest_start_coord - relative_start) + (relative_end - nearest_end_coord)` (unaligned child-end overhang).
     - `mismatches += find_mismatched_bases(child.start, child.end, ...)` (`:382-392`, count node mismatches inside the child's relative span).
     - `insertions += find_insertions(start_node, end_node, ...)` (`:373-379`, reference gaps between consecutive nodes).
     - `strand = get_strand(shortest_path_nodes[0], parent)`.
     - Build a lifted `new_feature(child.id, child.featuretype, reference_name, 'Liftoff', strand, min(lifted_start,lifted_end)+1, max(lifted_start,lifted_end)+1, child.frame, dict(child.attributes))` (the `+1` converts 0-based lifted coords to 1-based GFF). Insert into `mapped_children[new_child.id]`.
   - Else (child entirely unaligned): `deletions += child.end - child.start + 1`.
4. `alignment_length = total_bases + insertions`.
5. **Return** `mapped_children`, **coverage** `= (total_bases - deletions) / total_bases`, **sequence identity** `= (alignment_length - insertions - mismatches - deletions) / alignment_length` (`:309-311`).

These two ratios are exactly the `-a` (coverage) and `-s` (sequence identity) quantities thresholded by the caller in `lift_features` to accept/reject a mapping; `convert_id_to_original`/`get_copy_tag` (`liftoff_utils.py:54-64`) strip/parse the `_0`/`_N` copy-tag and `_frag` suffixes that key the lifted output.

#### 6.5.6 How `lift_original_annotation` drives the above (`liftover_types.py:4-36`)

LiftOn's main lift pass calls `lift_original_annotation(ref_chroms, target_chroms, lifted_features_list, args, unmapped_features, parents_to_lift, ref_db)`:
1. `liftover_type = "chrm_by_chrm"`.
2. **Threshold selection** (`:6-9`): if `target_chroms[0] == args.target and args.exclude_partial == False` → `min_cov, min_seqid = 0.05, 0.05` (permissive — partial mappings retained and later flagged). Else `min_cov, min_seqid = args.a, args.s` (the strict `0.5`/`0.5` defaults).
3. `extract_features_to_lift(...)` (§6.5.3) → hierarchy, db, parent_order.
4. `align_and_lift_features(...)` (`:21-36`): `align_features_to_target` (§6.5.4) → `aligned_segments`; then `lift_features.lift_all_features(aligned_segments, min_cov, feature_db, feature_hierarchy, unmapped_features, lifted_features_list, min_seqid, None, args, ref_parent_order)` (which internally calls `find_best_mapping` per gene and thresholds on `min_cov`/`min_seqid`); then `fix_overlapping_features.fix_incorrectly_overlapping_features(...)` resolves remaining overlaps using `args.overlap` (default `0.1`).
5. return `(feature_db, feature_hierarchy, ref_parent_order)`.

Gotcha (threshold inversion): in the common LiftOn case (whole-genome target, `exclude_partial` unset) the **effective** alignment thresholds during the first pass are `0.05/0.05`, not the nominal `0.5/0.5`. Partial/low-identity mappings are kept and annotated downstream rather than dropped at this stage. The strict `-a`/`-s` values only apply when `exclude_partial` is set or when a per-chromosome `-chroms` restriction is in effect.
