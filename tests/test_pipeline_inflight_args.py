import pytest

from lifton import lifton


BASE = ["target.fa", "reference.fa", "-g", "reference.gff3"]


def test_all_pipeline_inflight_limits_parse():
    args = lifton.parse_args(BASE + [
        "--step7-max-inflight", "3",
        "--step8-max-inflight", "4",
        "--evaluation-max-inflight", "5",
    ])

    assert args.step7_max_inflight == 3
    assert args.step8_max_inflight == 4
    assert args.evaluation_max_inflight == 5


@pytest.mark.parametrize("option", [
    "--step7-max-inflight",
    "--step8-max-inflight",
    "--evaluation-max-inflight",
])
def test_pipeline_inflight_limits_must_be_positive(option):
    with pytest.raises(SystemExit):
        lifton.parse_args(BASE + [option, "0"])
