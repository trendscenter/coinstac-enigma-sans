import glob
import json
import os
import subprocess
import sys

from utils import listRecursive


def remote_1(args):
    scriptDir = "/computation/enigma_scripts"
    scriptName = "metaanalysis_thickness_SANSTOT_asis_105ROIs.R"
    RScriptDir = "/usr/bin/Rscript"

    regr_args = [
        RScriptDir,
        os.path.join(scriptDir, scriptName), args["state"]["baseDirectory"],
        args["state"]["transferDirectory"], args["state"]["transferDirectory"]
    ]
    subprocess.call(regr_args,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL)
    os.system("cp -rf " + args["state"]["baseDirectory"] + "/* " + args["state"]["transferDirectory"]) 
    site_dict = dict()
    for site in args["input"]:
        logPath = os.path.join(args["state"]["baseDirectory"], site, '*.log')
        file_dict = dict()
        for file in glob.glob(logPath):
            file_name = os.path.split(file)[-1]
            with open(file, 'r') as f:
                fileContent = f.read()
            file_dict[file_name] = fileContent
        site_dict[site] = file_dict

    computation_output = {
        "output": site_dict,  # should be a list of files created -ross
        "success": True
    }
    return json.dumps(computation_output)


if __name__ == '__main__':

    parsed_args = json.loads(sys.stdin.read())
    phase_key = list(listRecursive(parsed_args, 'computation_phase'))

    if 'local_1' in phase_key:
        computation_output = remote_1(parsed_args)
        sys.stdout.write(computation_output)
    else:
        raise ValueError("Error occurred at Remote")
