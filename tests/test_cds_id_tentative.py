"""Tier-2 audit fix — a rejected candidate must not burn stable CDS IDs.

`reserve_explicit` / `allocate_synthetic` mutate `CdsIdAllocator._claims` with no
release. All three rescue paths (miniprot_rescue._stage_gene,
cross_locus_rescue._build_replacement_text, miniprot_pipeline's parallel Step 8)
serialize a candidate gene into a THROWAWAY buffer and only afterwards decide whether to
keep it -- so a candidate that was rejected still permanently claimed its identifiers.
The transcript that legitimately owns e.g. ``cds-NP_9.1`` was then pushed to
``cds-NP_9.1-1``, making the output depend on how many candidates happened to be
rejected and defeating the stable-CDS-ID guarantee those paths exist to protect.

`CdsIdAllocator.tentative()` claims provisionally and rolls back unless committed.
"""
from lifton.cds_id_allocator import CdsIdAllocator


def test_rejected_scope_releases_synthetic_claims():
    allocator = CdsIdAllocator()
    with allocator.tentative():
        assert allocator.allocate_synthetic("cds-X", "rna-rejected") == "cds-X"
        # no commit -> rolled back

    # The legitimate owner still gets the un-suffixed name.
    assert allocator.allocate_synthetic("cds-X", "rna-keeper") == "cds-X"


def test_committed_scope_keeps_claims():
    allocator = CdsIdAllocator()
    with allocator.tentative() as scope:
        assert allocator.allocate_synthetic("cds-X", "rna-first") == "cds-X"
        scope.commit()

    # A different parent must now be pushed off the taken name.
    assert allocator.allocate_synthetic("cds-X", "rna-second") == "cds-X-1"


def test_rejected_scope_releases_explicit_claims():
    allocator = CdsIdAllocator()
    with allocator.tentative():
        assert allocator.reserve_explicit("cds-NP_9.1", "rna-rejected") is True

    # The rightful owner can still claim its own stable ID.
    assert allocator.reserve_explicit("cds-NP_9.1", "rna-keeper") is True


def test_rollback_restores_exactly_the_prior_state():
    allocator = CdsIdAllocator()
    allocator.reserve_explicit("cds-kept", "rna-A")
    before = dict(allocator._claims)

    with allocator.tentative():
        allocator.reserve_explicit("cds-tentative", "rna-B")
        allocator.allocate_synthetic("cds-other", "rna-B")

    assert allocator._claims == before


def test_nested_commit_and_rollback_are_independent():
    allocator = CdsIdAllocator()
    with allocator.tentative() as outer:
        allocator.allocate_synthetic("cds-outer", "rna-outer")
        with allocator.tentative():
            allocator.allocate_synthetic("cds-inner", "rna-inner")
            # inner not committed -> only its claim is dropped
        outer.commit()

    assert allocator.allocate_synthetic("cds-inner", "rna-late") == "cds-inner"
    assert allocator.allocate_synthetic("cds-outer", "rna-late2") == "cds-outer-1"


def test_exception_inside_scope_still_rolls_back():
    allocator = CdsIdAllocator()
    try:
        with allocator.tentative():
            allocator.allocate_synthetic("cds-X", "rna-boom")
            raise RuntimeError("staging blew up")
    except RuntimeError:
        pass
    assert allocator.allocate_synthetic("cds-X", "rna-keeper") == "cds-X"
