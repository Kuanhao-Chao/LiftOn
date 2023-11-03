import subprocess, os

def check_miniprot_installed():
    miniprot_path = "miniprot"
    command = [miniprot_path, "--version"]
    installed = False
    try:
        res = subprocess.run(command)
        print(res)
        installed = True
    except: 
        pass
    return installed

def run_miniprot(args, tgt_genome, ref_proteins_file):
    # print(">> run_miniprot")

    miniprot_outdir = os.path.dirname(args.output) + "/miniprot/"
    os.makedirs(miniprot_outdir, exist_ok=True)
    miniprot_output = miniprot_outdir + "/miniprot.gff3"
    
    miniprot_path = "miniprot"
    command = [miniprot_path, "--gff-only", tgt_genome, ref_proteins_file]
    try:
        fw = open(miniprot_output, "w")
        res = subprocess.call(command, stdout=miniprot_output)
        fw.close()
        print(command)
    except: 
        pass
    return miniprot_output
