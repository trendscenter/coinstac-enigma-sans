import glob
import json
import os
import subprocess
import sys

from utils import listRecursive


def remote_1(args):
    scriptDir = "/computation/enigma_scripts"
    scriptName = "metaanalysis_AllComps2020.R"
    RScriptDir = "/usr/bin/Rscript"

    # numSites = len(args["input"])
    
    baseDir = args["state"]["baseDirectory"]
    transferDir = args["state"]["transferDirectory"]
    outputDir = args["state"]["outputDirectory"]
    
    site_list = os.path.join(outputDir, "site_list.txt")

    with open(site_list, 'w') as fh:
        for site in args["input"]:
            fh.write('%s\n' % site)
        
    regr_args = [
        RScriptDir,
        os.path.join(scriptDir, scriptName), baseDir, outputDir,
        site_list
    ]

    subprocess.call(regr_args,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL)

    # Copying local results to transfer directory
    os.system("cp -rf " + baseDir + "/* " + transferDir)
    
    # Copying meta analysis results to transfer directory
    os.system("cp -rf " + outputDir + "/* " + transferDir)

    computation_output = {
        "output": {},
        # should be a list of files created -ross
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
