import subprocess

with open("test.fa", "w") as f:
    f.write(">target\nACGTACGTACGT\n")

with open("test_prot.fa", "w") as f:
    f.write(">prot\nMMMAAA***\n")

proc = subprocess.run(["miniprot", "--gff-only", "test.fa", "test_prot.fa"], capture_output=True, text=True)
print("STDOUT:", proc.stdout)
print("STDERR:", proc.stderr)
