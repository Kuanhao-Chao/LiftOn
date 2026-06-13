### 6.5.B Vendored Liftoff internals — Part B: feature resolution, polish, output & native path

This part continues §6.5. It documents the second half of the vendored Liftoff fork (`lifton/liftoff/`): how a chosen alignment is expanded into lifted child features, how parents are reconstructed, how cross-feature overlaps are resolved, how copy-number / unmapped / unplaced recovery passes are orchestrated, the GFF3/GTF serializer (file and in-memory variants), the optional `--polish` protein-aware refinement, the shared helpers, the core data structures, and the `--native` in-process mappy alignment path. All line citations are real, taken from the files as they exist on the `devel` branch.

Throughout this part, "Feature" means a `new_feature.new_feature` instance (§6.5.B.7) unless qualified. Coordinates inside `aligned_seg` records (§6.5.B.7) are **0-based, query-relative**; coordinates inside `new_feature` objects are **1-based, inclusive** GFF3 coordinates. The bridge between the two conventions happens in `convert_all_children_coords` (`find_best_mapping.py:284`), which adds `+1` when materializing lifted starts/ends.

#### 6.5.B.0 Where this fits in the call graph

The vendored Liftoff entrypoint is `_run_liftoff_pipeline` (`liftoff_main.py:60`). It runs these stages in order, all sharing one `lifted_feature_list` dict, one `feature_hierarchy`, one `feature_db`, and one accumulating `unmapped_features` list:

1. `liftover_types.lift_original_annotation` (`liftover_types.py:4`) — the primary chromosome-by-chromosome pass.
2. `map_unmapped_features` → `liftover_types.map_unmapped_genes_agaisnt_all` (`liftover_types.py:39`) — retry unmapped genes against the whole genome.
3. `map_features_from_unplaced_seq` → `liftover_types.map_unplaced_genes` (`liftover_types.py:62`).
4. `write_unmapped_features_file` (`liftoff_main.py:287`) — dump unmapped IDs to `args.u`.
5. `map_extra_copies` → `liftover_types.map_extra_copies` (`liftover_types.py:85`) — extra-copy detection, only when `args.copies`.
6. Optional `check_cds` / `--polish` refinement.

The final output is written by `write_new_gff.write_new_gff` (`write_new_gff.py:13`) for the legacy/disk path, or produced as a byte string by `inmemory_emitter.lifted_features_to_gff3_bytes` (`inmemory_emitter.py:30`) for the `--inmemory-liftoff` path. Both produce **byte-identical** output (the in-memory emitter reuses the same `finalize_parent_features` / `build_parent_dict` / `write_feature` functions).

---

#### 6.5.B.1 `lift_features.py` — placing children within the chosen mapping

`lift_all_features(alns, threshold, feature_db, feature_hierarchy, unmapped_features, lifted_feature_list, seq_id_threshold, feature_locations, args, ref_parent_order)` (`lift_features.py:4`) is the driver that turns the `{query_name: [aligned_seg, ...]}` dict produced by alignment into entries in `lifted_feature_list`.

##### Driver: `lift_all_features` (`lift_features.py:4`)

Numbered algorithm:

1. `features_to_lift = feature_hierarchy.parents` (a `{ref_parent_id: Feature}` dict).
2. `feature_order = get_feature_order(feature_db)` (see below). This is a `{featuretype: int}` rank used later for output ordering.
3. `alignments = sort_alignments(features_to_lift, alns)` — produces a list of per-parent alignment groups sorted by genomic position of the reference parent (see below).
4. For each `alignment` (one parent's full alignment group, a `list[aligned_seg]`):
   1. `previous_feature_start, previous_feature_seq, previous_feature_ref_start = find_neighbor_location(...)` — locate the already-lifted upstream neighbor of this parent (used to bias mapping selection toward synteny).
   2. `lifted_features, parent_name = lift_single_feature(...)`.
   3. If `lifted_features != []`: `lifted_feature_list[parent_name] = lifted_features`. **Gotcha:** `parent_name` is `aligned_feature[0].query_name`, i.e. the **copy-tagged** name like `geneA_0` (suffix added by `edit_name`, §6.5.B.10), not the bare ref id. The `lifted_feature_list` is therefore keyed by copy id.

##### `get_feature_order(gene_db)` (`lift_features.py:28`)

Builds a deterministic featuretype→rank map. Steps:
1. `feature_types = list(gene_db.featuretypes())`, `index = 0`, `feature_order = {}`.
2. If `'exon'` in `feature_types`: `feature_order['exon'] = 0`, `index = 1`.
3. If `'CDS'` in `feature_types`: `feature_order['CDS'] = index`, `index += 1`.
4. For every remaining `feature_type` not yet present: assign the next `index`.

**Gotcha:** `exon` and `CDS` are *forced* to ranks 0 and 1 (in that order) regardless of their order in the DB. This rank is consumed by `process_final_features_list` to sort children below their parent. It is byte-identity-affecting: children are emitted exon-first, then CDS, then everything else.

##### `sort_alignments(parent_dict, alignments)` (`lift_features.py:45`)

Sorts the alignment groups so parents are lifted in genome order (important: `find_neighbor_location` requires upstream neighbors to already be lifted).
1. Build `parent_list` = `[parent_dict[convert_id_to_original(group[0].query_name)] for each group]`.
2. Sort `parent_list` by `(parent.seqid, parent.start)`.
3. Assign each parent an incrementing `order` into `order_dict[parent.id]`.
4. Sort the alignment-group values by `order_dict[convert_id_to_original(group[0].query_name)]`.

**Gotcha:** `convert_id_to_original` (§6.5.B.6) strips both the `_<copy>` and `_frag…` suffix to recover the bare reference parent id used to key `parent_dict`.

##### `find_neighbor_location(ref_parents, alignment, lifted_feature_list, ref_parent_order)` (`lift_features.py:60`)

Returns `(previous_feature_start, previous_feature_seq, previous_feature_ref_start)` for the **upstream non-overlapping reference neighbor** of the parent currently being lifted. Steps:
1. `ref_feature = ref_parents[convert_id_to_original(alignment[0].query_name)]`.
2. `ref_neighbor_name = liftoff_utils.find_nonoverlapping_upstream_neighbor(ref_parent_order, ref_feature.id)`.
3. If `ref_neighbor_name is not None`:
   1. `ref_neighbor_key = ref_neighbor_name + "_0"` (the copy-0 key in `lifted_feature_list`).
   2. If that key is present: return `(lifted_feature_list[key][0].start, lifted_feature_list[key][0].seqid, ref_parents[ref_neighbor_name].start)`.
4. Otherwise return `(0, "", 0)`.

These three values feed `find_best_mapping`'s `sort_alignments` (in `find_best_mapping.py`, §6.5.B.5) to prefer placements that preserve the reference distance between this gene and its left neighbor.

##### `lift_single_feature(...)` (`lift_features.py:73`)

The per-parent worker. Signature parameters in order: `(threshold, feature_order, features_to_lift, feature_hierarchy, previous_feature_start, previous_feature_ref_start, previous_gene_seq, unmapped_features, aligned_feature, seq_id_threshold, feature_locations, lifted_features_list, args)`.

Steps:
1. `new_parent_name = aligned_feature[0].query_name` (copy-tagged).
2. `original_parent_name = convert_id_to_original(new_parent_name)`.
3. `parent = features_to_lift[original_parent_name]`.
4. If `len(aligned_feature) > 0`:
   1. `lifted_children, alignment_coverage, seq_id = find_best_mapping.find_best_mapping(aligned_feature, parent.end - parent.start + 1, parent, feature_hierarchy, previous_feature_start, previous_feature_ref_start, previous_gene_seq, feature_locations, lifted_features_list, args)`. The second argument `parent.end - parent.start + 1` is the **query length** (the parent's reference span in bp). `alignment_coverage` ∈ [0,1] is fraction of child bases retained; `seq_id` ∈ [0,1] is sequence identity over the aligned length.
   2. `lifted_features = merge_lifted_features.merge_lifted_features(lifted_children, parent, unmapped_features, threshold, new_parent_name, feature_order, feature_hierarchy, alignment_coverage, seq_id, seq_id_threshold)`.
5. Else: `unmapped_features.append(parent)` (note: `lifted_features` is then unbound — but step 4's `len(aligned_feature) > 0` is effectively always true because `lift_single_feature` is only called with non-empty groups).
6. Return `(lifted_features, aligned_feature[0].query_name)`.

The coordinate transform itself — converting query-relative aligned-block coordinates into 1-based target genome coordinates per child — lives in `find_best_mapping.convert_all_children_coords` (§6.5.B.5), not in `lift_features.py`. `lift_features.py` is purely the orchestration layer.

---

#### 6.5.B.5 `find_best_mapping.py` — chaining and the coordinate transform

`find_best_mapping(alignments, query_length, parent, feature_heirarchy, previous_feature_start, previous_feature_ref_start, previous_gene_seq, inter, lifted_features_list, args)` (`find_best_mapping.py:6`) selects the best chain of alignment blocks via a shortest-path search over a DAG and emits lifted child Features.

##### Top-level algorithm

1. `children = feature_heirarchy.children[parent.id]`; `children_coords = liftoff_utils.merge_children_intervals(children)` (merged, sorted `[start,end]` intervals in 1-based ref coords).
2. `node_dict, aln_graph = intialize_graph()` — node 0 is a sentinel "start" `aligned_seg` with all coords `-1` (`find_best_mapping.py:31`).
3. `head_nodes = add_single_alignments(...)` — add one node per aligned block, connecting consecutive blocks of the same `aln_id` into linear sub-chains; record the first node of each sub-chain in `head_nodes`.
4. `chain_alignments(head_nodes, ...)` — add cross-alignment edges between compatible nodes (different `aln_id`, same strand, same contig, monotone query/ref order, distance within `args.d` factor, no overlap with already-lifted homologues).
5. `add_target_node(...)` — append a sentinel "end" node at `query_length` and connect every terminal node to it.
6. `shortest_path_nodes = find_shortest_path(node_dict, aln_graph)`; if empty, return `({}, 0, 0)`.
7. `mapped_children, alignment_coverage, seq_id = convert_all_children_coords(shortest_path_nodes, children, parent)`; return them.

##### Node/edge weights (`find_best_mapping.py:161-193`)

- **Node weight** `get_node_weight` (`:161`): for each child interval, count this block's `mismatches` falling inside the child's relative span and multiply by `args.mismatch` (default **2**). Sum over children.
- **Edge weight** `get_edge_weight` (`:171`): penalizes unaligned exonic bases in the gap between two chained blocks at `args.gap_extend` (default **1**) per base; if any gap exists, add `(args.gap_open - args.gap_extend)` once (default `2 - 1 = 1`) as the gap-open surcharge.
- **Path weight** `get_weight` (`:268`): `node_u_wt/2 + node_v_wt/2 + edge_wt`. `nx.shortest_path` minimizes the sum of these.

##### Edge validity `is_valid_edge` (`find_best_mapping.py:214`)

A directed edge from→to is **invalid** (returns False) if any holds:
- `from_node.aln_id == to_node.aln_id` (same alignment),
- `from_node.query_block_end >= to_node.query_block_end` (not strictly progressing in query),
- `from_node.is_reverse != to_node.is_reverse`,
- `from_node.reference_name != to_node.reference_name`,
- `to_node.reference_block_start < from_node.reference_block_end` (reference overlap / backwards),
- `actual_distance > args.d * expected_distance` where `expected_distance = to.query_block_end - from.query_block_start`, `actual_distance = to.reference_block_end - from.reference_block_start` (`args.d` default **2.0**),
- `spans_overlap_region(...)` is true (the gap between the two blocks would overlap an already-lifted higher-priority homologue).

##### The coordinate transform `convert_all_children_coords` (`find_best_mapping.py:284`)

This is where query-relative block coordinates become 1-based target genome Features. For each `child` in `children`:
1. `total_bases += child.end - child.start + 1`.
2. Find nearest aligned start/end relative coords inside the path (`find_nearest_aligned_start_and_end`, `:314`).
3. If both found:
   1. `lifted_start, start_node = convert_coord(nearest_start_coord, path)`; `lifted_end, end_node = convert_coord(nearest_end_coord, path)`. `convert_coord` (`:352`) maps a query-relative coord to a 0-based reference coord by `block.reference_block_start + (relative_coord - block.query_block_start)`.
   2. Accumulate `deletions`, `mismatches`, `insertions` across blocks spanned (`:294-297`).
   3. `strand = get_strand(path[0], parent)` — combines the alignment's `is_reverse` flag with the parent's strand (§6.5.B.6 `get_strand`).
   4. If `"ID"` not in `child.attributes`: set `child.attributes["ID"] = [child.id]`.
   5. Build `new_child = new_feature.new_feature(child.id, child.featuretype, path[0].reference_name, 'Liftoff', strand, min(lifted_start, lifted_end) + 1, max(lifted_start, lifted_end) + 1, child.frame, dict(child.attributes))`. **The `+1` converts the 0-based ref coord to 1-based GFF3.** Source is hardcoded `'Liftoff'`.
   6. `mapped_children[new_child.id] = new_child`.
4. Else (child not covered by the path): `deletions += child.end - child.start + 1`.
5. After the loop: `alignment_length = total_bases + insertions`. Return:
   - `mapped_children`,
   - `alignment_coverage = (total_bases - deletions) / total_bases`,
   - `seq_id = (alignment_length - insertions - mismatches - deletions) / alignment_length`.

**Gotcha:** coverage and identity are computed over **child (exon/CDS) bases only**, not the whole gene span. A gene whose introns map badly but whose exons map cleanly still scores high.

---

#### 6.5.B.2 `merge_lifted_features.py` — reconstructing the parent and intermediates

`merge_lifted_features(mapped_children, parent, unmapped_features, aln_cov_threshold, copy_id, feature_order, feature_hierarchy, aln_cov, seq_id, seq_id_threshold)` (`merge_lifted_features.py:4`) takes the flat `{child_id: Feature}` dict from `find_best_mapping` and rebuilds the full parent→intermediate→child hierarchy, then applies the coverage/identity acceptance test.

Numbered algorithm:
1. Init `feature_list = {}`, `final_features = []`, `non_parents = []`, `top_target_feature = None`.
2. If `len(mapped_children) == 0`: `unmapped_features.append(parent)`; return `[]`.
3. For each `child_feature` in `mapped_children.values()`:
   1. `feature_list[child_feature.id] = child_feature`.
   2. If `is_top_parent(child_feature, parent)` is False (i.e. `child.id != parent.id`): append `(child_feature, child_feature.attributes["Parent"][0])` to `non_parents`.
   3. Else: `top_target_feature = child_feature` (the top-level gene was itself a "child" in the lift; this happens when the gene has no sub-features).
4. While `non_parents` is non-empty: `non_parents, top_target_feature = create_parents(non_parents, parent, feature_hierarchy, feature_list)` — iteratively manufacture each missing intermediate/parent Feature whose bounds enclose its children.
5. `final_features = process_final_features_list(feature_list, top_target_feature, seq_id, seq_id_threshold, aln_cov, aln_cov_threshold, unmapped_features, parent, feature_order, copy_id)`.
6. Return `final_features`.

##### `create_parents` (`merge_lifted_features.py:31`) and `make_new_parent` (`:45`)

For each `(child, parent_id)` pair whose `parent_id` hasn't been created yet this round:
1. `target_parent_feature = make_new_parent(feature_list, parent_id, feature_hierarchy)`:
   1. `children = [f for f in feature_list.values() if "Parent" in f.attributes and f.attributes["Parent"][0] == parent_id]`.
   2. `ref_parent = get_ref_parent(parent_id, feature_hierarchy)` — looks first in `feature_hierarchy.parents`, then `feature_hierarchy.intermediates` (`:57`).
   3. New Feature: `new_feature(ref_parent.id, ref_parent.featuretype, children[0].seqid, 'Liftoff', children[0].strand, min(child starts), max(child ends), ref_parent.frame, dict(ref_parent.attributes))`. **The manufactured parent's span is the min/max of its already-lifted children**; its strand/seqid come from the children, but its featuretype, frame, and attributes come from the reference.
   4. Insert into `feature_list`.
2. If the manufactured parent is the top parent, set `top_target_feature`; else re-queue it as a non-parent (its own parent may still need creating).

This produces, on the second loop iteration, the grandparent above an intermediate (e.g. `gene` above `mRNA` above `exon`/`CDS`), handling arbitrarily deep hierarchies.

##### `process_final_features_list` (`merge_lifted_features.py:65`) — ordering + acceptance test

1. `final_features = [f for f in feature_list.values() if f != top_target_feature]`.
2. Sort by `(seqid, start)`.
3. Then sort by `feature_order[featuretype]` (stable; exon=0, CDS=1, others after — see §6.5.B.1).
4. Insert `top_target_feature` at index 0.
5. **Acceptance test:** if `aln_cov < aln_cov_threshold` OR `seq_id < seq_id_threshold`:
   - `final_features = []`; `unmapped_features.append(parent)`. The gene is rejected entirely.
6. Else (accepted):
   - `top_target_feature.score = 1 - seq_id`,
   - `top_target_feature.attributes["copy_num_ID"] = [copy_id]` (e.g. `["geneA_0"]`),
   - `top_target_feature.attributes["coverage"] = [str(aln_cov)[0:5]]` (truncated to 5 chars, e.g. `"0.987"`),
   - `top_target_feature.attributes["sequence_ID"] = [str(seq_id)[0:5]]`.
7. Return `final_features` (parent first, then children in exon→CDS→other / position order).

**Gotcha:** the truncation `str(x)[0:5]` keeps 5 characters of the raw `repr` of the float — e.g. `1.0` → `"1.0"`, `0.99999…` → `"0.999"`. These string forms become the `coverage`/`sequence_ID` GFF3 attributes downstream and are re-parsed with `float(...)` in `write_new_gff.add_attributes`. This truncation is byte-identity-affecting.

---

#### 6.5.B.3 `liftover_types.py` — orchestration, copy-number & unmapped tracking

Five public functions, all converging on `align_and_lift_features` with different `(min_cov, min_seqid, liftover_type)` settings.

##### `lift_original_annotation(ref_chroms, target_chroms, lifted_features_list, args, unmapped_features, parents_to_lift, ref_db)` (`liftover_types.py:4`)

The primary pass. Steps:
1. `liftover_type = "chrm_by_chrm"`.
2. Thresholds: if `target_chroms[0] == args.target` AND `args.exclude_partial == False`: `(min_cov, min_seqid) = (0.05, 0.05)` (permissive — accept partial maps, defer rejection to later); else `(args.a, args.s)` (defaults **0.5, 0.5**).
3. `feature_hierarchy, feature_db, ref_parent_order = extract_features.extract_features_to_lift(...)`.
4. `align_and_lift_features(...)` with `max_overlap = args.overlap` (default **0.1**).
5. Return `(feature_db, feature_hierarchy, ref_parent_order)`.

##### `align_and_lift_features(...)` (`liftover_types.py:21`)

The shared core every pass calls:
1. `aligned_segments = align_features.align_features_to_target(ref_chroms, target_chroms, args, feature_hierarchy, liftover_type, unmapped_features)` — dispatches to subprocess minimap2 or the native mappy path (§6.5.B.11).
2. `feature_locations = None`.
3. `lift_features.lift_all_features(aligned_segments, min_cov, feature_db, feature_hierarchy, unmapped_features, lifted_features_list, min_seqid, feature_locations, args, ref_parent_order)`.
4. `fix_overlapping_features.fix_incorrectly_overlapping_features(lifted_features_list, lifted_features_list, aligned_segments, unmapped_features, min_cov, feature_hierarchy, feature_db, ref_parent_order, min_seqid, args, max_overlap)` (§6.5.B.4).

##### `map_unmapped_genes_agaisnt_all(...)` (`liftover_types.py:39`)

Whole-genome retry for genes that failed chromosome-by-chromosome lifting.
1. `liftoff_utils.clear_scores(lifted_features_list, feature_hierarchy.parents)` — resets each already-lifted top feature's `score` to `-1` (§6.5.B.6).
2. `unmapped_dict = get_unmapped_genes(unmapped_features)` — `{feature.id: feature}`.
3. Thresholds: if `args.exclude_partial`: `(args.a, args.s)`; else `(0.05, 0.05)`. **Note the inversion vs the primary pass.**
4. `liftover_type = "unmapped"`.
5. `extract_features.get_gene_sequences(unmapped_dict, ...)` re-extracts just those genes.
6. Reset `unmapped_features = []` (local), call `align_and_lift_features`, return the fresh `unmapped_features`.

##### `map_unplaced_genes(...)` (`liftover_types.py:62`)

For genes on reference scaffolds listed in `args.unplaced`. `liftover_type = "unplaced"`. `unplaced_dict` = every parent in `feature_hierarchy` whose `seqid in ref_chroms` (`get_features_from_unplaced_seq`, `:76`). Same threshold logic as the unmapped pass.

##### `map_extra_copies(...)` (`liftover_types.py:85`)

Detects additional paralogous copies. `liftover_type = "copies"`. Re-extracts **all** parent sequences (`feature_hierarchy.parents`). Thresholds: `min_cov = 0`, `min_seqid = args.sc` (the copy-specific sequence-identity threshold, default **1.0** — i.e. only near-perfect extra copies). With `min_cov = 0` the coverage gate is effectively disabled; only `args.sc` matters. Only invoked when `args.copies` is set (`liftoff_main.py:295`); `liftoff_main.py:237` guards `args.s > args.sc`.

**Copy-number bookkeeping summary:** The `_<n>` suffix (copy index) is assigned by `align_features.edit_name` (§6.5.B.10) at SAM-parse time: copy 0 (`_0`) for the primary search types, incrementing copy ids for the `"copies"` search type. The suffix is what `lifted_feature_list` is keyed by. The final `extra_copy_number` attribute and `copy_num_ID` rewrite happen in `write_new_gff.finalize_parent_features` (§6.5.B.8), which counts how many times each bare parent id appears and assigns `extra_copy_number` = 0, 1, 2, … in `id`-sorted order.

---

#### 6.5.B.4 `fix_overlapping_features.py` — overlap resolution

After lifting, two lifted homologues may incorrectly occupy overlapping target coordinates (e.g. a gene placed where a paralog already sits). This module detects such collisions and re-lifts the loser elsewhere, removing it if no valid alternative exists.

##### Entry: `fix_incorrectly_overlapping_features(...)` (`fix_overlapping_features.py:6`)

1. `features_to_remap, feature_locations = check_homologues(all_lifted_features, features_to_check, feature_hierarchy.parents, ref_parent_order, max_overlap)`.
2. `resolve_overlapping_homologues(...)`.

##### `check_homologues(...)` (`fix_overlapping_features.py:18`)

1. `all_feature_list = liftoff_utils.get_parent_list(all_lifted_features)` (top feature of each lifted entry).
2. `target_parent_order = liftoff_utils.find_parent_order(all_feature_list)` — an Nx2 numpy array of `[id, Feature]` sorted by `(seqid, start)`.
3. `feature_locations = build_interval_list(all_feature_list)` — an `interlap.InterLap` of `[start-1, end-1, [copy_id, Feature]]` (0-based half-open for InterLap).
4. For each `feature` in `features_to_check_list`:
   1. `overlaps = liftoff_utils.find_overlaps(feature.start-1, feature.end-1, feature.seqid, feature.strand, feature.attributes["copy_num_ID"][0], feature_locations, parent_dict, all_lifted_features, max_overlap)` (§6.5.B.6).
   2. For each `overlap` with a *different* `copy_num_ID`: `compare_overlapping_feature(...)` decides which of the two to remap and adds its `copy_num_ID` to the `remap_features` set.
5. Return `(remap_features, feature_locations)`.

##### Which feature loses — `find_feature_to_remap(...)` (`fix_overlapping_features.py:52`)

Decision cascade (first matching rule wins, top to bottom):

| # | Condition | Loser (returned to remap) |
|---|---|---|
| 1 | `feature` is a copy (id not ending `_0`) and overlap is not | `feature` |
| 2 | overlap is a copy and `feature` is not | `overlap` |
| 3 | `feature` already queued in `remap_features` | `feature` |
| 4 | overlap already queued | `overlap` |
| 5 | `feature.score < overlap.score` (`has_greater_seq_id(feature, overlap)`) | `overlap` |
| 6 | `overlap.score < feature.score` | `feature` |
| 7 | out-of-order neighbor heuristic returns a feature (`find_out_of_order_feature`) | that feature |
| 8 | `feature` is shorter (`feature.end-feature.start < overlap…`) | `feature` |
| 9 | overlap is shorter | `overlap` |
| 10 | fallback | `overlap` |

**Gotcha:** `score` is `1 - seq_id` (set in `process_final_features_list`). So "greater seq id" means **lower** `score`; `has_greater_seq_id(f1, f2)` returns True when `f1.score < f2.score` (`:87`). Rule 5/6: the feature with the *worse* identity is remapped. `is_copy` (`:79`) checks `copy_num_ID[0][-2:] != '_0'`.

##### Out-of-order heuristic — `find_out_of_order_feature` (`:93`) / `check_order` (`:111`)

For each non-copy participant, compute the index distance, in `target_parent_order`, between the feature and its reference upstream neighbor (`find_distance_between_features`, `:120`). If the *other* overlapping feature's id lies between them in the target order, subtract 1 from the distance. The feature whose neighbor distance is *larger* (more displaced from its expected position) is the loser. Returns None if neither has a determinable neighbor.

##### `resolve_overlapping_homologues(...)` (`fix_overlapping_features.py:138`)

Iterative remap loop:
1. `iter = 0`, `max_iter = 10 * len(features_to_remap)`.
2. While `features_to_remap` non-empty:
   1. `iter += 1`; if `iter > max_iter`: **break** (safety valve against oscillation).
   2. `aligned_segs_for_remap = remove_features_and_get_alignments(...)` — `del lifted_feature_list[loser]` for each, and collect their original alignment groups.
   3. `lift_features.lift_all_features(aligned_segs_for_remap, ...)` — re-lift them (now with the overlap-occupied loci masked via `feature_locations`).
   4. `features_to_check = get_successfully_remapped_features(...)` — those that re-appeared in `lifted_feature_list`.
   5. `features_to_remap, feature_locations = check_homologues(...)` again on the re-lifted subset.
3. `remove_unresolved_features(...)`: any id still in `features_to_remap` is appended (its bare ref parent) to `unmapped_features` and `del`-eted from `lifted_feature_list`.

`build_interval_list` (`:36`) and `is_shorter` (`:134`) are the helpers; `find_distance_between_features` uses `np.where` on the `target_parent_order[:,0]` id column.

---

#### 6.5.B.8 `write_new_gff.py` — output GFF3/GTF serialization (disk path)

`write_new_gff(lifted_features, args, feature_db)` (`write_new_gff.py:13`) writes the final annotation.

##### Top-level flow

1. Open `args.output` for writing (or `sys.stdout` when `args.output == 'stdout'`).
2. `out_type = feature_db.dialect['fmt']` — `'gff3'` or `'gtf'`.
3. `write_header(f, out_type)` (`:6`): if gff3, emit `##gff-version 3\n`; always emit `# LiftOn v<__version__>\n` (`__version__ = 'v1.0.8'`, so line is `# LiftOn vv1.0.8`) and `# <space-joined sys.argv>\n`.
4. `parents = liftoff_utils.get_parent_list(lifted_features)` (top feature of each entry).
5. `parents.sort(key=lambda x: x.id)` — **first** sorted by bare id (stable basis for copy numbering).
6. `final_parent_list = finalize_parent_features(parents, args)` (assigns copy numbers and attributes; see below).
7. `final_parent_list.sort(key=lambda x: (x.seqid, x.start))` — **then** re-sorted by genome position for output.
8. For each `final_parent`:
   1. `child_features = lifted_features[final_parent.attributes["copy_id"][0]]`.
   2. `parent_child_dict = build_parent_dict(child_features, final_parent)` (`:66`) — `{parent_id: [child Feature, …]}`, and stamps each child with `extra_copy_number` copied from the parent.
   3. `write_feature([final_parent], f, child_features, parent_child_dict, out_type)` — recursive depth-first writer.

##### `finalize_parent_features` / `add_attributes` (`write_new_gff.py:30`, `:48`) — copy numbering

For each parent (in id-sorted order):
1. `add_to_copy_num_dict` (`:41`): first time a given `parent.id` is seen → `copy_num = 0`; each subsequent occurrence increments. This is the `extra_copy_number`.
2. `add_attributes(parent, copy_num, args)`:
   1. `parent.score = "."`.
   2. If `"copy_id"` not yet in attributes: `parent.attributes["copy_id"] = parent.attributes["copy_num_ID"]` (preserves the original copy-tagged key for the `lifted_features[...]` lookup in step 8.1 above).
   3. Delete keys `copy_num_ID`, `partial_mapping`, `low_identity`, `extra_copy_number` (so they re-add in a deterministic position at the end).
   4. `parent.attributes["extra_copy_number"] = [str(copy_num)]`.
   5. `parent.attributes["copy_num_ID"] = [parent.id + "_" + str(copy_num)]`.
   6. If `float(parent.attributes["coverage"][0]) < args.a`: set `partial_mapping = ["True"]`.
   7. If `float(parent.attributes["sequence_ID"][0]) < args.s`: set `low_identity = ["True"]`.

**Gotcha:** the acceptance during lifting used the *permissive* `0.05` thresholds (in the common `chrm_by_chrm` case); the **final** partial/low-identity flags here use the *strict* `args.a`/`args.s` (defaults 0.5). So a feature can be emitted but flagged. The `coverage`/`sequence_ID` values being compared are the 5-char-truncated strings written in §6.5.B.2.

##### `build_parent_dict` (`:66`)

For each child with a `Parent`: copy `final_parent.attributes["extra_copy_number"]` into the child, and bucket the child under `child.attributes["Parent"][0]`.

##### Recursive writer `write_feature` (`:78`)

DFS: for each `child` in the current level, `write_line(child, ...)`, then if `child.id` is itself a parent key in `parent_dict`, recurse into its children. Starting from the single top parent, this emits parent → intermediates → leaves in hierarchy order.

##### `write_line` / `make_gff_line` / `edit_copy_ids` (`:87`, `:100`, `:118`)

1. If `feature.attributes["extra_copy_number"][0] != '0'`: `attr_dict = edit_copy_ids(feature)`; else `attr_dict = feature.attributes`.
   - `edit_copy_ids` (`:118`): for every attribute key ending in `_id`, plus `ID` and `Parent`, append `"_" + copy_num` to the value (e.g. extra copy 1 of `geneA` exon becomes `ID=exonX_1`, `Parent=geneA_1`). Wrapped in `try/except (KeyError, IndexError)` printing an error and `continue`-ing on failure.
2. `make_gff_line(attr_dict, feature)` (`:100`):
   - Start `attributes_str`. If `"ID"` present and non-empty: prepend `ID=<value[0]>;`.
   - For each attr (skipping `copy_id`): build a comma-joined value string; if the attr isn't `ID`, append `attr=values;`. (ID is emitted first only, never again.)
   - Line = `seqid \t source \t featuretype \t start \t end \t "." \t strand \t frame \t attributes_str[:-1]` (the trailing `;` is stripped). **Score column is hardcoded `"."`.**
3. `make_gtf_line` (`:147`) is the GTF analogue: `key "value";` syntax, score `"."`, phase/frame column `"."`.
4. `write_line` writes `line` then `"\n"`.

**Gotcha:** the ID-first ordering in `make_gff_line` (`:104`) is deliberate and byte-identity-affecting — `ID` is always the first attribute, then the remaining attributes in dict-insertion order, with `copy_id` always suppressed.

---

#### 6.5.B.9 `inmemory_emitter.py` — `lifted_features_to_gff3_bytes`

`lifted_features_to_gff3_bytes(lifted_features, args, feature_db, *, sys_argv=None) -> bytes` (`inmemory_emitter.py:30`) is the byte-identical RAM counterpart of `write_new_gff.write_new_gff`, used by `--inmemory-liftoff`.

Algorithm (mirrors §6.5.B.8 exactly, writing to an `io.StringIO` instead of a file):
1. `buf = io.StringIO()`; `out_type = feature_db.dialect["fmt"]`.
2. `_write_header_to(buf, out_type, sys_argv=sys_argv)` (`:20`) — same three header lines as `write_header`; if `sys_argv is None` it defaults to `sys.argv`.
3. `parents = liftoff_utils.get_parent_list(lifted_features)`; `parents.sort(key=lambda x: x.id)`.
4. `final_parent_list = write_new_gff.finalize_parent_features(parents, args)` — **reuses the disk path's function**, so copy numbering and attribute side-effects are identical.
5. `final_parent_list.sort(key=lambda x: (x.seqid, x.start))`.
6. For each parent: `child_features = lifted_features[final_parent.attributes["copy_id"][0]]`; `parent_child_dict = write_new_gff.build_parent_dict(...)`; `write_new_gff.write_feature([final_parent], buf, ...)`.
7. Return `buf.getvalue().encode("utf-8")`.

**Why it is byte-identical:** every transformation that touches feature content (`finalize_parent_features`, `build_parent_dict`, `write_feature`, `make_gff_line`) is the *same function object* imported from `write_new_gff`. Only the sink differs (StringIO vs file handle) and the final `.encode("utf-8")`. The header is reproduced line-for-line by `_write_header_to`. The one degree of freedom is `sys_argv`: the caller passes the real `sys.argv` to make the `# <command line>` provenance line match.

**Gotcha:** `finalize_parent_features` **mutates** the Feature objects in place (deletes/re-adds attributes, sets `score = "."`). If both the disk write and the in-memory emit ran over the same `lifted_features` dict, the second run would see already-mutated attributes (e.g. `copy_num_ID` already rewritten to `id_0`). The pipeline avoids this by choosing exactly one of the two paths.

---

#### 6.5.B.10 `liftoff_utils.py` — shared helpers

| Function (`liftoff_utils.py:`) | Signature | Behaviour |
|---|---|---|
| `count_overlap` (`:4`) | `(start1,end1,start2,end2)` | `min(end1,end2) - max(start1,start2) + 1`. Inclusive overlap length; negative when disjoint. |
| `get_relative_child_coord` (`:9`) | `(parent, coord, is_reverse)` | If `is_reverse`: `parent.end - coord`; else `coord - parent.start`. Converts a genomic coord to a 0-based query-relative offset honoring strand. |
| `merge_children_intervals` (`:17`) | `(children)` | Sort `[start,end]` by start; merge intervals where `current[0] <= previous[1]` (inclusive, 1-based-friendly merge). Returns `[]` for empty. |
| `get_parent_list` (`:32`) | `(feature_list)` | For each entry, take `feature_set[0]` (the top feature) when the list is non-empty. |
| `clear_scores` (`:41`) | `(feature_list, parent_dict)` | Set `feature.score = -1` for every feature whose `id` is a key in `parent_dict`. Resets scores between passes. |
| `find_parent_order` (`:49`) | `(parents)` | Sort by `(seqid,start)`; return `np.array([[id, Feature], …])`. |
| `convert_id_to_original` (`:54`) | `(id)` | Strip `_frag…` suffix, then strip the `_<copytag>` suffix: `frag_split[:-copy_tag_len-1]`. Recovers the bare reference parent id. |
| `get_copy_tag` (`:61`) | `(id)` | Return the `_<n>` suffix (e.g. `_0`). |
| `get_strand` (`:67`) | `(aln, parent)` | If `aln.is_reverse`: flip parent strand; else parent strand. |
| `find_overlaps` (`:78`) | `(start,end,chrm,strand,feature_name,intervals,parent_dict,lifted_features_list,max_overlap_non_copies)` | Query the InterLap for overlapping intervals on the same `(chrm, strand)`; keep an overlap only if it is an **incorrect** one. |
| `overlaps_in_ref_annotation` (`:100`) | `(ref1,ref2)` | True iff same seqid, same strand, different id, and `count_overlap > 0`. Used to *exclude* legitimately-overlapping reference features from the collision set. |
| `find_nonoverlapping_upstream_neighbor` (`:113`) | `(parent_order, feature_name)` | Return the id of the immediately-preceding feature in `parent_order` if on the same seqid; else None. |

##### `find_overlaps` acceptance test in detail (`:78`)

For each candidate overlap (already filtered to same chrm & strand):
1. `shortest_feature_length = min(end-start, overlap[1]-overlap[0]) + 1`.
2. `overlap_amount = count_overlap(start, end, overlap[0], overlap[1])`.
3. `ref_feature = parent_dict[convert_id_to_original(feature_name)]`; `ref_overlap_feature = parent_dict[convert_id_to_original(overlap[2][0])]`.
4. `max_overlap = 0` if either feature is a copy (copy tag != `_0`), else `max_overlap_non_copies`.
5. Keep the overlap (mark incorrect) iff **all**:
   - `overlaps_in_ref_annotation(ref_feature, ref_overlap_feature) is False` (the two are *not* expected to overlap in the reference),
   - `overlap[2][0] in lifted_features_list` (the other feature is currently placed),
   - `overlap_amount / shortest_feature_length > max_overlap`.

**Gotcha:** copies get `max_overlap = 0`, meaning *any* overlap of a copy is incorrect; non-copies tolerate up to `args.overlap` (default 0.1) fractional overlap.

---

#### 6.5.B.7 Core data structures: `new_feature`, `aligned_seg`, `feature_hierarchy`

##### `new_feature.new_feature` (`new_feature.py:1`)

A plain mutable record used everywhere a lifted feature is represented (it intentionally mirrors the read fields of a gffutils Feature so the rest of the code is backend-agnostic).

| Field | Type | Meaning |
|---|---|---|
| `id` | str | Feature ID (bare, no copy tag at construction). |
| `featuretype` | str | `gene` / `mRNA` / `exon` / `CDS` / … |
| `seqid` | str | Target (or reference) sequence name. |
| `source` | str | GFF source column; lifted features use `'Liftoff'`. |
| `strand` | str | `'+'` / `'-'`. |
| `start` | int | 1-based inclusive start (GFF3 convention). |
| `end` | int | 1-based inclusive end. |
| `frame` | str | CDS phase column. |
| `attributes` | dict[str, list[str]] | GFF3 attributes; each value is a list (gffutils convention). |

`score` is **not** a constructor argument; it is set as an attribute later (`process_final_features_list` sets `score = 1 - seq_id`; `add_attributes` sets `score = "."`). Accessing `.score` before it is set raises `AttributeError` — but every code path sets it before reading.

##### `aligned_seg.aligned_seg` (`aligned_seg.py:1`)

One contiguous gap-free alignment block (between indels).

| Field | Type | Meaning |
|---|---|---|
| `aln_id` | int/str | Groups blocks belonging to the same SAM record/alignment. Sentinels use `"start"`/`"end"`. |
| `query_name` | str | Copy-tagged query (reference feature) name. |
| `reference_name` | str | Target contig the block maps to. |
| `query_block_start` | int | 0-based query offset, block start. |
| `query_block_end` | int | 0-based query offset, block end (inclusive). |
| `reference_block_start` | int | 0-based target coord, block start. |
| `reference_block_end` | int | 0-based target coord, block end (inclusive). |
| `is_reverse` | bool | Whether the alignment is reverse-strand. |
| `mismatches` | np.ndarray[int] | Query-relative offsets of mismatched bases inside this block. |

##### `feature_hierarchy.feature_hierarchy` (`feature_hierarchy.py:1`)

A three-field container produced by `extract_features.extract_features_to_lift`:

| Field | Type | Meaning |
|---|---|---|
| `parents` | dict[str, Feature] | Top-level features to lift (genes), keyed by bare id. |
| `intermediates` | dict[str, Feature] | Mid-level features (e.g. mRNA) that are neither top parents nor leaf children, keyed by id. |
| `children` | dict[str, list[Feature]] | Leaf children (exon/CDS) keyed by the **bare parent gene id**. |

`make_new_parent`/`get_ref_parent` consult `parents` then `intermediates` to recover reference attributes when reconstructing the hierarchy.

---

#### 6.5.B.6 `polish.py` — optional `--polish` protein-aware refinement

Activated by `--polish`. The pipeline first runs `check_cds` (`liftoff_main.py:339`), which for every lifted feature calls `polish.find_and_check_cds` to annotate ORF validity, then `find_and_polish_broken_cds` (`liftoff_main.py:303`) re-aligns and re-lifts genes whose ORF is broken, keeping the better of the original vs polished version.

##### `find_and_check_cds(...)` (`polish.py:33`) and `count_good_cds` (`:53`)

For each lifted feature's CDS set:
1. `grouped_cds = group_cds_by_tran(cds_features)` — group CDS by `Parent` transcript id.
2. For each CDS group:
   1. `cds_seq = get_seq(cds_group, target_faidx)` — concatenate `chrom[start-1:end]` for each CDS (sorted by start), reverse-complement if strand `-`, uppercase (`get_seq`, `:126`).
   2. `protein = cds_seq.translate()` (Biopython, warnings suppressed).
   3. `transcript.attributes["matches_ref_protein"]` = `'True'`/`'False'` from `matches_ref(...)` (`:152`) which translates the corresponding **reference** CDS group and compares proteins exactly.
   4. `longest_ORF, longest_ORF_coords = get_longest_ORF(cds_seq)` (`:88`).
   5. Re-translate `longest_ORF`; classify the transcript:

| Condition (checked in order) | Attribute set | `valid_ORF` |
|---|---|---|
| `len(protein) < 3` | `partial_ORF = ['True']` | `'False'` |
| `missing_start(protein)` i.e. `protein[0] != 'M'` | `missing_start_codon = ['True']` | `'False'` |
| `missing_stop(protein)` i.e. `protein[-1] != '*'` | `missing_stop_codon = ['True']` | `'False'` |
| `inframe_stop(protein)` i.e. `'*' in protein[:-1]` | `inframe_stop_codon = ['True']` | `'False'` |
| else | — | `'True'`, `good_cds_count += 1` |

   6. If `longest_ORF != cds_seq`: `adjust_cds_coords(cds_group, longest_ORF_coords, feature_list)` — trims/removes CDS features so they span only the longest ORF (`:103`).
3. `target_sub_features[0].attributes["valid_ORFs"] = [str(num_good_cds)]`.

##### `get_longest_ORF(cds_seq)` (`polish.py:88`)

1. Regex `ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)` finds all in-frame, stop-terminated ORFs as `[start,end]`.
2. Sort by length; take the longest `possible_ORF = ORFs[-1]`.
3. Use it **only if** its length differs from the full `cds_seq` length **and** is `>= 180` bp. Otherwise keep the full CDS.

##### `adjust_cds_coords(...)` (`polish.py:103`)

1. `coords = np.hstack(...)` — flattened array of every genomic position covered by the CDS group.
2. For each CDS, find its `relative_start`/`relative_end` indices within `coords`. Then:
   - If `relative_end < longest_ORF_coords[0]`: remove this CDS entirely (it's upstream of the ORF).
   - If `relative_start > longest_ORF_coords[1]`: remove it (downstream).
   - If the ORF start falls inside this CDS: `cds.start += (longest_ORF_coords[0] - relative_start)`.
   - If the ORF end falls inside: `cds.end -= (relative_end - longest_ORF_coords[1] + 1)`.

##### `polish_annotations(...)` (`polish.py:11`) and the SAM realignment

When at least one transcript has `matches_ref_protein` but not `valid_ORF == True`:
1. Open `args.directory + "/polish.sam"`, write a SAM `@SQ` header for every target sequence (`write_sam_header`, `:44`).
2. `polish_annotation(ref_gene, target_gene, ref_children, target_children, ref_fa, target_fa, output_sam)` (`:164`):
   1. Collect ref/target exons (falling back to CDS if no exons). Merge ref exon and CDS intervals.
   2. `add_splice_sites(ref_exons, ref_gene)` (`:192`): extend each exon 2 bp on each side (within the gene bounds) and record those flanking 2-bp regions as splice sites — these get uppercased in the reference sequence to bias alignment toward exact splice boundaries.
   3. Build a custom parasail scoring matrix (`make_scoring_matrix(3)`, `:223`): non-coding match reward **2**, mismatch penalty **−4**, coding (uppercase) match reward **3**, over alphabet `"ACGTacgt*"`.
   4. For each merged ref interval: determine the corresponding target interval (`get_target_interval`, `:238`), compute a `flank` = `max(0, ref_len - target_len) + 10`, fetch sequences, uppercase CDS+splice-site segments (`cds_and_splice_sites_to_upper`, `:270`), reverse-complement if strands differ, align with `parasail.sg_dx_trace_scan_sat(ref, target, 10, 1, matrix)` (gap-open 10, gap-extend 1), and write a SAM record (`write_sam_file`, `:283`).
   5. `remove_splice_sites(ref_exons, ref_gene)` (`:203`): undo the 2-bp extension.

##### CIGAR construction in polish (`make_cigar`, `:305`; `condense_cigar_string`, `:324`)

From the parasail traceback, build an expanded per-position CIGAR: `=` (match), `X` (mismatch), `I` (target gap), `D` (ref gap, i.e. `reference_traceback[i] == '-'`). The leading run starts at the first non-gap reference position. `condense_cigar_string` run-length-encodes via `itertools.groupby`, drops a leading `D` run, and prepends `<hard_clip_start>H` when there is a hard-clip prefix. `write_sam_file` sets bit flag 16 for reverse alignments and computes the SAM `POS` from the flanked target start plus the alignment offset.

##### Keep-the-better decision (`liftoff_main.py:323-336`)

After re-lifting, for each polished feature: replace the original iff
- the original has no `valid_ORFs`, OR polished `valid_ORFs > original valid_ORFs`; OR
- equal `valid_ORFs` and polished `sequence_ID > original` (string comparison); OR
- equal `valid_ORFs` and polished `coverage > original` (string comparison).

**Gotcha:** the tie-breaks compare `sequence_ID`/`coverage` as **strings**, not floats (`'0.987' > '0.98'` lexicographically). This is the legacy Liftoff behaviour preserved verbatim. Output suffix becomes `args.output += "_polished"`.

---

#### 6.5.B.11 `native_align.py` — the `--native` in-process mappy path

`align_features.align_features_to_target` (`align_features.py:12`) dispatches to the native path when `getattr(args, "native", False)` is truthy **and** `getattr(args, "subcommand", None) != "polish"` (`align_features.py:19-27`). The polish subcommand always needs a real SAM file, so it stays on the subprocess path.

##### `is_mappy_available()` (`minimap_facade.py:22`)

`try: import mappy; return True / except ImportError: return False`. The single switch the native path consults.

##### `align_features_to_target_native(...)` (`native_align.py:155`)

1. If `not is_mappy_available()`: write a warning to `stderr` ("`--native` requested but `mappy` is not installed; falling back to the subprocess minimap2 path…") and return `_af.align_features_to_target(...)` (the legacy path). **Gotcha:** the warning is emitted **once per call** to `align_features_to_target_native` — and that function is called once per liftover pass (primary, unmapped, unplaced, copies), so a missing binding can print the warning several times in one run.
2. Build one `MinimapAligner(target_fasta, mm2_options=getattr(args,"mm2_options",""), threads=int(getattr(args,"threads",1) or 1))` for the whole target genome (`minimap_facade.py:30`).
3. `hits_by_query = {}`. For each `chrom_idx`:
   1. `features_file, _ = _af.get_features_file(ref_chroms, args, liftover_type, chrom_idx)` — reuses the *same* per-chromosome reference-feature FASTA the subprocess path consumes.
   2. For each `(query_name, query_seq)` from `_iter_query_fasta(features_file)` (`:127`, mappy `fastx_read` with a pure-Python FASTA fallback): for each `hit` in `aligner.map(query_name, query_seq)`: append to `hits_by_query[query_name]`.
4. `return parse_alignment_from_hits(hits_by_query, feature_hierarchy, unmapped_features, liftover_type)`.

##### `MinimapAligner.map` (`minimap_facade.py:74`)

Wraps `mappy.Aligner.map(query_seq)` (constructed with `best_n=50` by default) and projects each `mappy.Alignment` into a frozen `MinimapHit` dataclass (`types.py:11`) carrying `query_name, ctg, r_st, r_en, q_st, q_en, strand (+1/-1), mapq, NM, cigar_str, is_primary`. `mappy.Aligner.map` releases the GIL during the C kernel, which is what makes the native path safe under a `ThreadPoolExecutor`.

##### `_PysamShim` (`native_align.py:57`) — adapting mappy hits to the pysam consumer

The legacy SAM parser `align_features.add_alignment` / `get_aligned_blocks` reads pysam-style attributes. `_PysamShim` (a `__slots__` class) exposes exactly those from a `MinimapHit`:

| pysam attr | source | note |
|---|---|---|
| `query_name` | `hit.query_name` | reference feature name |
| `reference_name` | `hit.ctg` | target contig |
| `reference_start` | `int(hit.r_st)` | **both mappy and pysam are 0-based — no shift** |
| `query_alignment_start` | `int(hit.q_st)` | |
| `query_alignment_end` | `int(hit.q_en)` | |
| `cigar` | `cigar_str_to_pysam_tuples(hit.cigar_str)` | list of `(op_int, length)` |
| `is_reverse` | `int(hit.strand) == -1` | mappy strand −1 ⇒ reverse |
| `is_unmapped` | `False` | hits are by definition mapped |

`cigar_str_to_pysam_tuples` (`:46`) tokenizes a CIGAR string with regex `(\d+)([MIDNSHP=X])` and maps op chars to pysam op codes via `_CIGAR_OP_TO_PYSAM` (`:38`): `M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8`. **Gotcha:** `get_aligned_blocks` (`align_features.py:182`) only treats op 7 (`=`) and op 8 (`X`) as aligned bases and op 5 (`H`) for the hard-clip query-start offset. minimap2 must therefore be run with `--eqx` (in `args.mm2_options`) for the `=`/`X` operators to appear — a plain `M` (op 0) CIGAR would yield no aligned blocks. This matches the subprocess path's requirement.

##### `parse_alignment_from_hits(...)` (`native_align.py:82`)

In-memory mirror of `align_features.parse_alignment`:
1. `all_aligned_blocks = {}`, `aln_id = 0`, `name_dict = {}`, `align_count_dict = {}`.
2. For each `query_name` in **`sorted(hits_by_query)`** (deterministic order so `aln_id` values are reproducible): for each `hit`: wrap in `_PysamShim`, call `_af.add_alignment(ref_seq, align_count_dict, search_type, name_dict, aln_id, feature_hierarchy, all_aligned_blocks)` (returns the next `aln_id`).
3. For each `(query_name, parent)` in `feature_hierarchy.parents.items()`: if the query produced no hits, append `parent` to `unmapped_features` (mirrors the SAM iterator's `is_unmapped` branch).
4. `_af.remove_alignments_without_children(all_aligned_blocks, unmapped_features, feature_hierarchy)`.
5. Return `all_aligned_blocks`.

Because steps 2–4 call the **same** `add_alignment` / `get_aligned_blocks` / `remove_alignments_without_children` functions as the subprocess parser, the resulting `aligned_seg` records are byte-equivalent — which is what keeps the `--native` cells of the 24-cell matrix green. The sole ordering precaution is the explicit `sorted(...)` at step 2, replacing the SAM-file iteration order.

**Gotcha — `aln_id` determinism:** in the subprocess path, `aln_id` increments in SAM-record order; in the native path it increments in `sorted(query_name)` order then hit order. For byte-identical block construction the `aln_id` values must coincide with what the SAM parser would have produced for the same input. This holds because (a) `aln_id` is an opaque grouping key (not emitted in output), and (b) the downstream `find_best_mapping` chaining treats blocks of the same `aln_id` identically regardless of the integer value — only equality/inequality of `aln_id` matters, never its magnitude.
