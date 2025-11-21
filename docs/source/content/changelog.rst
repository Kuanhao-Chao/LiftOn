
|

Changelog
===========

v1.0.7
-------

**Bug Fixes:**

- **Fixed gffutils UNIQUE constraint errors**: Enhanced duplicate feature ID handling with automatic recovery strategies. The system now automatically handles duplicate IDs in polished liftoff output and RefSeq annotations with non-overlapping CDS features, preventing silent failures during database creation.

- **Fixed GTF file format processing**: Added automatic GTF format detection and proper handling. GTF files are now correctly processed with automatic gene/transcript inference, and optional automatic conversion to GFF3 format using gffread or agat tools.

- **Fixed ID parsing for IDs ending with numbers**: Improved `get_ID_base()` function to safely handle feature IDs that naturally end with underscore and number (e.g., `FMUND_1`). The function now only removes suffixes when confirmed to be copy numbers, preventing silent failures.

- **Fixed CDS ID preservation**: CDS features now preserve their IDs in output files, complying with GFF3 specification. CDS features from the same mRNA can now share the same ID as required by the standard.

**Improvements:**

- **Enhanced biotype attribute support**: Added support for generic `biotype` attribute as fallback when `gene_biotype` (RefSeq) or `gene_type` (GENCODE/ENSEMBL/CHESS) are not present. This ensures protein-coding features are correctly identified regardless of annotation source.

- **Automatic GTF to GFF3 conversion**: Added optional automatic conversion of GTF files to GFF3 format for better compatibility. Conversion uses gffread (preferred) or agat tools if available, with graceful fallback to direct GTF processing.

- **Improved error handling**: Enhanced error messages and recovery strategies for database creation failures, providing better user guidance and automatic problem resolution.

- **Better format detection**: Improved file format detection logic that checks multiple lines and patterns to reliably distinguish between GTF and GFF3 formats, with GFF3 as the safe default.

v1.0.0
-------

- Initial release of LiftOn
- Release via the documentation (http://ccb.jhu.edu/lifton)
- Released via the paper (bioRxiv coming soon!)


|
|
|
|
|



.. image:: ../_images/jhu-logo-dark.png
   :alt: My Logo
   :class: logo, header-image only-light
   :align: center

.. image:: ../_images/jhu-logo-white.png
   :alt: My Logo
   :class: logo, header-image only-dark
   :align: center

