import subprocess, os, sys

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
    miniprot_output = miniprot_outdir + "miniprot.gff3"
    
    miniprot_path = "miniprot"
    command = [miniprot_path, "--gff-only", tgt_genome, ref_proteins_file]
    print("miniprot command: ", command)
    try:
        fw = open(miniprot_output, "w")
        subprocess.run(command, stdout=fw)
        fw.close()
    except: 
        print("failed to run miniprot")
        sys.exit(1)
    return miniprot_output
