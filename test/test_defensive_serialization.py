from lifton.lifton_class import LiftOn_FEATURE


class DummyEntry:
    def __init__(self, id):
        self.id = id
        self.seqid = "chr1"
        self.source = "test"
        self.featuretype = "region"
        self.start = 1
        self.end = 10
        self.score = "."
        self.strand = "+"
        self.frame = "."
        self.attributes = {
            "ID": [id],
            "numerical_list": [1, 2.5, 3],
            "mixed_list": ["text", 4],
            "single_int": 42
        }

def test_defensive_serialization(tmp_path):
    # Numeric attributes must serialize safely without mutating the model.
    
    # Fake entry with numerical properties
    gff_entry = DummyEntry("test_id")
    
    # Initialize FEATURE (copy_num "0" usually leaves ID unchanged)
    feature = LiftOn_FEATURE("parent_id", gff_entry, "0")
    
    # Create fake file
    out_file = tmp_path / "out.gff"
    
    with open(out_file, "w") as fw:
        assert feature.write_entry(fw) is True

    text = out_file.read_text()
    assert "numerical_list=1,2.5,3" in text
    assert "mixed_list=text,4" in text
    assert "single_int=42" in text

    # Canonical serialization is non-destructive.
    attrs = feature.entry.attributes
    assert attrs["numerical_list"] == [1, 2.5, 3]
    assert attrs["mixed_list"] == ["text", 4]
    assert attrs["single_int"] == 42
