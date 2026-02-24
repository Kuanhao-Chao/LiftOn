import pytest
from lifton.lifton_class import LiftOn_FEATURE

class DummyFeatureDB:
    pass

class DummyEntry:
    def __init__(self, id):
        self.id = id
        self.attributes = {
            "ID": [id],
            "numerical_list": [1, 2.5, 3],
            "mixed_list": ["text", 4],
            "single_int": 42
        }

def test_defensive_serialization(tmp_path):
    # This test verifies that we explicitly stringify numerical parameters
    # before they are exported to gffutils to avoid concatenation crashes.
    
    # Fake entry with numerical properties
    gff_entry = DummyEntry("test_id")
    
    # Initialize FEATURE (copy_num "0" usually leaves ID unchanged)
    feature = LiftOn_FEATURE("parent_id", gff_entry, "0")
    
    # Create fake file
    out_file = tmp_path / "out.gff"
    
    # Call the write method
    # Since DummyEntry lacks a proper `__str__` for gffutils, we just care
    # that the attributes were successfully mutated into strings before output strings were formulated.
    with open(out_file, "w") as fw:
        # We manually inject a dummy `__str__` so `str(self.entry)` doesn't fail
        gff_entry.__class__.__str__ = lambda self: "dummy_string"
        feature.write_entry(fw)
        
    # Check that attributes are now all strings
    attrs = feature.entry.attributes
    assert attrs["numerical_list"] == ["1", "2.5", "3"]
    assert attrs["mixed_list"] == ["text", "4"]
    assert attrs["single_int"] == "42"
