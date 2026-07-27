"""Explicit repeated fault exercises excluded from normal focused runs."""

from __future__ import annotations

import os

import pytest

from lifton.output_transaction import OutputTransaction, OutputTransactionState


@pytest.mark.fault_stress
@pytest.mark.skipif(
    os.environ.get("LIFTON_RUN_FAULT_STRESS") != "1",
    reason="set LIFTON_RUN_FAULT_STRESS=1 or run make test-fault-stress",
)
def test_repeated_abort_never_replaces_last_good_destination(tmp_path):
    destination = tmp_path / "result.gff3"
    destination.write_text("known-good\n")

    for index in range(100):
        transaction = OutputTransaction(destination)
        transaction.open().write(f"interrupted-{index}\n")
        transaction.abort()
        assert transaction.state is OutputTransactionState.ABORTED
        assert destination.read_text() == "known-good\n"

    assert transaction.partial_path.read_text() == "interrupted-99\n"
