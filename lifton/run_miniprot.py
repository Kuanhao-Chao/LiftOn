import subprocess

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

def run_miniprot():
    print(">> run_miniprot")
    pass