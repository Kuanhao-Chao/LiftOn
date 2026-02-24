"""
annotation.py  —  GFF3 / GTF annotation file loading & gffutils database management.

Builds a gffutils SQLite database from an annotation file with:
  • Pre-flight validation (format detection, duplicate ID scan, orphan-parent check)
  • 3-level DB-build fallback (create_unique → unique-ID transform → merge)
  • Rich, actionable error messages at every failure point
  • Auto GTF→GFF3 conversion via gffread or agat
"""

import gffutils
import sys
import os
import subprocess
from collections import defaultdict
import re

from lifton import extract_sequence, logger
from lifton.annotation_validator import (
    validate_annotation_file,
    print_validation_report,
    print_db_build_error,
    print_db_build_success,
)


# ─────────────────────────────────────────────────────────────────────────────
# Annotation class
# ─────────────────────────────────────────────────────────────────────────────

class Annotation:
    """
    Wraps a gffutils FeatureDB for a GFF3 or GTF annotation file.

    The constructor:
      1. Detects the file format (GTF vs GFF3).
      2. Optionally converts GTF to GFF3.
      3. Runs pre-flight validation (duplicate IDs, orphan parents, …).
      4. Builds a gffutils database with up to 3 fallback strategies.
    """

    def __init__(
        self,
        file_name: str,
        infer_genes: bool,
        infer_transcripts: bool,
        merge_strategy: str = "create_unique",
        id_spec=None,
        force: bool = False,
        verbose: bool = False,
        auto_convert_gtf: bool = True,
    ):
        self.file_name      = file_name
        self.merge_strategy = merge_strategy
        self.id_spec        = id_spec
        self.force          = force
        self.verbose        = verbose
        self.auto_convert_gtf = auto_convert_gtf

        # ── 0. Check if the file is already a gffutils database ───────────────
        if os.path.exists(self.file_name) and os.path.getsize(self.file_name) > 0:
            try:
                # If this succeeds, it's a valid SQLite DB matching gffutils schema
                self._db_connection = gffutils.FeatureDB(self.file_name)
                if self.verbose:
                    logger.log_success(f"Loaded existing gffutils DB directly from {self.file_name}")
                return
            except Exception:
                pass  # It's not a DB file; proceed to parse as GFF3/GTF

        # ── 1. Detect file format ─────────────────────────────────────────────
        file_format = self._detect_file_format()

        if file_format == "GTF format":
            self._handle_gtf_input(file_name, infer_genes, infer_transcripts)
        elif file_format == "Unknown format":
            if self.verbose:
                logger.log_warning(
                    f"Could not determine file format for {file_name!r}. "
                    "Assuming GFF3."
                )
            self.infer_genes        = infer_genes
            self.infer_transcripts  = infer_transcripts
        else:
            # GFF3 (or any other explicitly identified format)
            self.infer_genes        = infer_genes
            self.infer_transcripts  = infer_transcripts

        # ── 2. Pre-flight validation ──────────────────────────────────────────
        self._run_preflight_validation()

        # ── 3. Build gffutils database ────────────────────────────────────────
        self._get_db_connection()


    # ─────────────────────────────────────────────────────────────────────────
    # Public property
    # ─────────────────────────────────────────────────────────────────────────

    @property
    def db_connection(self):
        return self._db_connection


    # ─────────────────────────────────────────────────────────────────────────
    # Step 1: format detection & GTF handling
    # ─────────────────────────────────────────────────────────────────────────

    def _detect_file_format(self) -> str:
        """Return 'GTF format', 'GFF format', or 'Unknown format'."""
        try:
            return extract_sequence.determine_file_format(self.file_name)
        except Exception as exc:
            if self.verbose:
                logger.log_warning(f"Format detection failed for {self.file_name!r}: {exc}")
            return "GFF format"   # safe default

    def _handle_gtf_input(self, original_file: str, infer_genes: bool, infer_transcripts: bool):
        """Configure inference flags and attempt GTF→GFF3 conversion."""
        print(
            f"\n[LiftOn] Detected GTF format: {original_file}",
            file=sys.stderr,
        )
        print(
            "[LiftOn] GTF support is experimental. "
            "For best results, convert to GFF3 first:\n"
            f"    gffread -E {original_file} -o {original_file}.gff3\n"
            f"    or: agat_sp_gtf2gff.pl --gtf {original_file} -o {original_file}.gff3",
            file=sys.stderr,
        )
        # GTF files always need inference
        self.infer_genes       = True
        self.infer_transcripts = True

        if self.auto_convert_gtf:
            converted = self._convert_gtf_to_gff3()
            if converted and os.path.exists(converted) and os.path.getsize(converted) > 0:
                converted_format = self._detect_file_format_for(converted)
                if converted_format == "GFF format":
                    logger.log_success(f"GTF converted to GFF3: {converted}")
                    self.file_name = converted
                else:
                    logger.log_warning(
                        f"Converted file {converted!r} may not be valid GFF3. "
                        "Using original GTF."
                    )
            else:
                if self.verbose:
                    logger.log_warning(
                        "GTF→GFF3 conversion produced no output. Using original GTF."
                    )
        else:
            self.infer_genes       = infer_genes
            self.infer_transcripts = infer_transcripts

    def _detect_file_format_for(self, path: str) -> str:
        try:
            return extract_sequence.determine_file_format(path)
        except Exception:
            return "GFF format"


    # ─────────────────────────────────────────────────────────────────────────
    # Step 2: pre-flight validation
    # ─────────────────────────────────────────────────────────────────────────

    def _run_preflight_validation(self):
        """
        Validate the annotation file before touching gffutils.
        Prints a boxed report when issues are found; exits on fatal errors only.
        """
        result = validate_annotation_file(
            self.file_name,
            max_duplicate_examples=20,
            check_orphan_parents=True,
        )

        # Always print the report if there are warnings or errors
        if result.errors or result.warnings:
            print_validation_report(result)

        # Hard stop on truly fatal conditions (file missing, empty, unreadable)
        fatal_keywords = [
            "not found",
            "not readable",
            "empty",
            "no valid 9-column",
            "not a regular file",
        ]
        truly_fatal = any(
            any(kw in err.lower() for kw in fatal_keywords)
            for err in result.errors
        )
        if truly_fatal:
            logger.log_error(
                f"Cannot build annotation database. "
                f"Please fix the issues in {self.file_name!r} and re-run."
            )
            sys.exit(1)
        # Duplicate IDs and orphan parents are handled by the DB-build fallback.


    # ─────────────────────────────────────────────────────────────────────────
    # Step 3: database building
    # ─────────────────────────────────────────────────────────────────────────

    def _get_db_connection(self):
        """
        Try to open an existing gffutils DB; if that fails, build one.
        Sets self._db_connection.
        """
        db_path = self.file_name + "_db"
        try:
            feature_db = gffutils.FeatureDB(db_path)
            if self.verbose:
                logger.log_success(f"Opened existing gffutils DB: {db_path}")
            self._db_connection = feature_db
        except Exception:
            self._db_connection = self._build_database()

    def _build_database(self) -> gffutils.FeatureDB:
        """
        Build a gffutils database with 3 progressively more permissive strategies.

        Strategy 1: user-specified merge_strategy (default: create_unique)
        Strategy 2: create_unique + unique-ID transform (handles duplicate IDs)
        Strategy 3: merge strategy + unique-ID transform (maximum permissiveness)

        Raises SystemExit on total failure.
        """
        disable_genes       = not self.infer_genes
        disable_transcripts = not self.infer_transcripts

        # ── Strategy 1: user's chosen merge_strategy ─────────────────────────
        print(
            f"\n[LiftOn] Building annotation database: {self.file_name}",
            file=sys.stderr,
        )
        print(
            f"[LiftOn]   strategy='{self.merge_strategy}'"
            f"  infer_genes={self.infer_genes}"
            f"  infer_transcripts={self.infer_transcripts}",
            file=sys.stderr,
        )
        try:
            db = gffutils.create_db(
                self.file_name,
                self.file_name + "_db",
                merge_strategy=self.merge_strategy,
                id_spec=self.id_spec,
                force=self.force,
                verbose=self.verbose,
                disable_infer_transcripts=disable_transcripts,
                disable_infer_genes=disable_genes,
                transform=self._get_transform_func(),
            )
            print_db_build_success(self.file_name, self.merge_strategy)
            return db
        except Exception as exc1:
            exc1_str = str(exc1)
            print_db_build_error(self.file_name, self.merge_strategy, exc1)

            # Only fall through if it's a known recoverable error
            if not _is_duplicate_id_error(exc1_str):
                # Unrecognised error — give up immediately with full detail
                self._fatal_db_error(exc1, "Strategy 1 (merge_strategy={!r})".format(
                    self.merge_strategy))

        # ── Strategy 2: create_unique + unique-ID transform ───────────────────
        print(
            "[LiftOn]   Retrying with unique-ID transformation "
            "(handles duplicate feature IDs in RefSeq/GENCODE files) …",
            file=sys.stderr,
        )
        try:
            db = gffutils.create_db(
                self.file_name,
                self.file_name + "_db",
                merge_strategy="create_unique",
                id_spec=self.id_spec,
                force=True,
                verbose=self.verbose,
                disable_infer_transcripts=disable_transcripts,
                disable_infer_genes=disable_genes,
                transform=self._get_unique_id_transform(),
            )
            print_db_build_success(self.file_name, "create_unique + unique-ID transform")
            logger.log_warning(
                "Database built using unique-ID transformation. "
                "Duplicate IDs were renamed (e.g. cds-XM_1 → cds-XM_1_dup1). "
                "This is safe for most use cases."
            )
            return db
        except Exception as exc2:
            print_db_build_error(self.file_name, "create_unique + unique-ID transform", exc2)

        # ── Strategy 3: merge + unique-ID transform ────────────────────────────
        print(
            "[LiftOn]   Retrying with merge strategy + unique-ID transformation …",
            file=sys.stderr,
        )
        try:
            db = gffutils.create_db(
                self.file_name,
                self.file_name + "_db",
                merge_strategy="merge",
                id_spec=self.id_spec,
                force=True,
                verbose=self.verbose,
                disable_infer_transcripts=disable_transcripts,
                disable_infer_genes=disable_genes,
                transform=self._get_unique_id_transform(),
            )
            print_db_build_success(self.file_name, "merge + unique-ID transform")
            logger.log_warning(
                "Database built using merge strategy. "
                "Overlapping features with the same ID were merged."
            )
            return db
        except Exception as exc3:
            print_db_build_error(self.file_name, "merge + unique-ID transform", exc3)
            self._fatal_db_error(exc3, "all 3 strategies")

    def _fatal_db_error(self, exc: Exception, strategy_desc: str):
        """Print a comprehensive fatal error block and exit."""
        logger.log_section(
            "❌  ANNOTATION DATABASE BUILD FAILED — CANNOT CONTINUE",
            [
                f"File        : {self.file_name}",
                f"Strategy    : {strategy_desc}",
                f"Error       : {exc}",
                "",
                "POSSIBLE CAUSES:",
                "  1. Duplicate feature IDs not handled by any strategy",
                "  2. Corrupted or malformed GFF3/GTF file",
                "  3. GTF file submitted without conversion",
                "",
                "SUGGESTED FIXES:",
                "  • Convert GTF → GFF3:",
                f"      gffread -E {self.file_name} -o {self.file_name}.gff3",
                "  • Run LiftOn on the converted file instead.",
                "  • If the file is GFF3, check for structural issues:",
                f"      grep -n 'ID=' {self.file_name} | sort -t= -k2 | uniq -D -f1 | head -20",
            ],
            kind="error",
        )
        sys.exit(1)


    # ─────────────────────────────────────────────────────────────────────────
    # Transform functions for gffutils
    # ─────────────────────────────────────────────────────────────────────────

    def _get_transform_func(self):
        if not self.infer_genes:
            return None
        return _transform_func_gtf

    def _get_unique_id_transform(self):
        """
        Return a transform function that makes every feature ID unique by
        appending _dup1, _dup2, … to repeated IDs.

        The closure also applies the GTF transcript_id transform when needed.
        """
        seen_ids: dict = {}
        infer_genes = self.infer_genes

        def unique_id_transform(feature):
            # Apply base GTF transform if needed
            if infer_genes:
                feature = _transform_func_gtf(feature)

            # Determine the current ID
            original_id = None
            if hasattr(feature, "id") and feature.id:
                original_id = feature.id
            elif "ID" in feature.attributes and feature.attributes["ID"]:
                original_id = feature.attributes["ID"][0]

            if not original_id:
                return feature

            if original_id in seen_ids:
                seen_ids[original_id] += 1
                new_id = f"{original_id}_dup{seen_ids[original_id]}"
                if "ID" in feature.attributes:
                    feature.attributes["ID"] = [new_id]
                if hasattr(feature, "id"):
                    feature.id = new_id
            else:
                seen_ids[original_id] = 0

            return feature

        return unique_id_transform


    # ─────────────────────────────────────────────────────────────────────────
    # GTF→GFF3 conversion
    # ─────────────────────────────────────────────────────────────────────────

    def _convert_gtf_to_gff3(self) -> "str | None":
        """
        Attempt conversion via gffread (preferred) then agat.
        Returns the output file path or None on failure.
        """
        base = os.path.splitext(self.file_name)[0]
        output_file = base + "_converted.gff3"

        gffread_ok = self._tool_available("gffread")
        agat_cmd   = self._find_agat_command()

        if not gffread_ok and agat_cmd is None:
            if self.verbose:
                logger.log_warning(
                    "Neither gffread nor agat is available. "
                    "GTF will be used directly (may fail).\n"
                    f"  Install gffread:  conda install -c bioconda gffread\n"
                    f"  Install agat:     conda install -c bioconda agat"
                )
            return None

        # Try gffread first
        if gffread_ok:
            try:
                print(
                    f"[LiftOn]   Converting GTF → GFF3 with gffread: {output_file}",
                    file=sys.stderr,
                )
                result = subprocess.run(
                    ["gffread", "-E", self.file_name, "-o", output_file],
                    capture_output=True, text=True, check=True,
                )
                if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                    return output_file
                logger.log_warning("gffread conversion produced an empty file.")
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                logger.log_warning(f"gffread conversion failed: {exc}")

        # Try agat as fallback
        if agat_cmd:
            try:
                print(
                    f"[LiftOn]   Converting GTF → GFF3 with agat ({agat_cmd}): {output_file}",
                    file=sys.stderr,
                )
                if "gtf2gff" in agat_cmd:
                    cmd = [agat_cmd, "--gtf", self.file_name, "-o", output_file]
                else:
                    cmd = [agat_cmd, "--gff", self.file_name, "-o", output_file]
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                    return output_file
                logger.log_warning("agat conversion produced an empty file.")
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                logger.log_warning(f"agat conversion failed: {exc}")

        return None

    def _tool_available(self, name: str) -> bool:
        try:
            subprocess.run(["which", name], capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    def _find_agat_command(self) -> "str | None":
        for cmd in ["agat_sp_gtf2gff.pl", "agat_convert_sp_gff2gtf.pl", "agat"]:
            if self._tool_available(cmd):
                return cmd
        return None


    # ─────────────────────────────────────────────────────────────────────────
    # Public query helpers (unchanged from original)
    # ─────────────────────────────────────────────────────────────────────────

    def get_protein_coding_features(self, feature_types):
        protein_coding_genes = []
        for f_type in feature_types:
            for feature in self._db_connection.features_of_type(featuretype=f_type):
                if self.is_highest_parent(feature):
                    CDS_list = [
                        child for child in self._db_connection.children(feature)
                        if child.featuretype == "CDS"
                    ]
                    if len(CDS_list) > 0:
                        protein_coding_genes.append(feature.id)
        return protein_coding_genes

    def get_noncoding_features(self, feature_types):
        noncoding_genes = []
        for f_type in feature_types:
            for feature in self._db_connection.features_of_type(featuretype=f_type):
                if self.is_highest_parent(feature):
                    CDS_list = [
                        child for child in self._db_connection.children(feature)
                        if child.featuretype == "CDS"
                    ]
                    if len(CDS_list) == 0:
                        noncoding_genes.append(feature.id)
        return noncoding_genes

    def get_novel_protein_coding_features(self, ref_genes, feature_types):
        novel_protein_coding = []
        all_protein_coding = self.get_protein_coding_features(feature_types)
        for feature in all_protein_coding:
            if feature not in ref_genes:
                novel_protein_coding.append(feature)
        return novel_protein_coding

    def get_novel_noncoding_features(self, ref_genes, feature_types):
        novel_noncoding = []
        all_noncoding = self.get_noncoding_features(feature_types)
        for feature in all_noncoding:
            if feature not in ref_genes:
                novel_noncoding.append(feature)
        return novel_noncoding

    def get_all_parent_feature_ids(self, feature_types):
        feature_ids = []
        for f_type in feature_types:
            feature_ids += [
                feature.id
                for feature in self._db_connection.features_of_type(featuretype=f_type)
                if self.is_highest_parent(feature)
            ]
        return feature_ids

    def make_parent_to_child_dict(self, protein_coding, gene):
        child_dict = defaultdict(list)
        child_list = self._db_connection.children(gene)
        child_count = 0
        for child in child_list:
            child_count += 1
            if (protein_coding and child.featuretype == "CDS") or (
                protein_coding is False and self.is_lowest_child(child)
            ):
                parent = list(self._db_connection.parents(child, level=1))[0]
                child_dict[parent.id].append(child)
        if child_count == 0:
            child_dict[gene] = [self._db_connection[gene]]
        return child_dict

    def get_features_of_type(self, feature_types):
        features_of_type = []
        for f_type in feature_types:
            features_of_type += list(self._db_connection.features_of_type(f_type))
        return features_of_type

    def get_feature_dict(self, feature_types):
        id_to_feature = {}
        features = self.get_features_of_type(feature_types)
        for feature in features:
            id_to_feature[feature.id] = feature
        return id_to_feature

    def get_source_name(self, feature_name):
        feature = self._db_connection[feature_name]
        if "extra_copy_number" in feature.attributes:
            if feature.attributes["extra_copy_number"][0] == feature.id.split("_")[-1]:
                return "_".join(feature.id.split("_")[:-1])
        return feature.id

    def is_lowest_child(self, feature_name):
        return len(list(self._db_connection.children(feature_name))) == 0

    def is_highest_parent(self, feature_name):
        return len(list(self._db_connection.parents(feature_name))) == 0

    def get_paralog_name(self, feature_name):
        gene_attributes = self._db_connection[feature_name].attributes
        if "extra_copy_number" in gene_attributes:
            paralog_name = "_".join(gene_attributes["ID"][0].split("_")[:-1])
            return paralog_name
        return ""

    def get_num_levels(self, feature_name):
        level = 1
        new_children = list(self._db_connection.children(feature_name, level=level))
        while len(new_children) != 0:
            level += 1
            new_children = list(self._db_connection.children(feature_name, level=level))
        return level

    # Keep old attribute name for backward compatibility
    # (some callers access .db_connection directly as an attribute, not property)
    def __getattr__(self, name):
        if name == "db_connection":
            # This is reached only when _db_connection is not yet set
            # (i.e. during __init__ before _get_db_connection is called).
            raise AttributeError(name)
        raise AttributeError(f"'Annotation' object has no attribute {name!r}")


# ─────────────────────────────────────────────────────────────────────────────
# Module-level helpers
# ─────────────────────────────────────────────────────────────────────────────

def _is_duplicate_id_error(exc_str: str) -> bool:
    """Return True when the exception indicates a duplicate/unique constraint."""
    lower = exc_str.lower()
    return any(
        kw in lower
        for kw in (
            "unique constraint failed",
            "uniqueconstraintviolation",
            "duplicate",
        )
    )


def _transform_func_gtf(x):
    """Suffix transcript_id with '_transcript' to avoid ID collisions in GTF files."""
    if "transcript_id" in x.attributes:
        x.attributes["transcript_id"][0] += "_transcript"
    return x


# ─────────────────────────────────────────────────────────────────────────────
# Standalone helpers (kept for backward-compatibility with other modules that
# import them directly from annotation.py)
# ─────────────────────────────────────────────────────────────────────────────

def transform_func(x):
    """Alias for _transform_func_gtf (backward compat)."""
    return _transform_func_gtf(x)


def get_perc_id(feature):
    if "sequence_ID" in feature.attributes:
        return float(feature.attributes["sequence_ID"][0])
    return 0.0


def merge_children_intervals(children):
    if len(children) == 0:
        return []
    intervals = [[child.start, child.end] for child in children]
    intervals.sort(key=lambda interval: interval[0])
    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged