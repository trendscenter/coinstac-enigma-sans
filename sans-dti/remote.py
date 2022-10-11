import json
import os
import subprocess
import sys

from utils import listRecursive


def remote_1(args):
    scriptDir = "/computation/enigma_scripts"
    scriptName = "metaanalysis_SANS_AllComps2022_v1.R"
    RScriptDir = "/usr/bin/Rscript"

    baseDir = args["state"]["baseDirectory"]
    transferDir = args["state"]["transferDirectory"]
    outputDir = args["state"]["outputDirectory"]

    regr_args = [
        RScriptDir,
        os.path.join(scriptDir, scriptName), baseDir, outputDir, json.dumps(args["input"])
    ]

    result = subprocess.run(regr_args,
                text=True,
                capture_output=True)

    if result.returncode != 0:
        raise Exception("R script failed: " + result.stderr + "\n" + result.stdout)

    # Copying local results to transfer directory
    os.system("cp -rf " + baseDir + "/* " + transferDir)

    # Copying meta analysis results to transfer directory
    os.system("cp -rf " + outputDir + "/* " + transferDir)

    computation_output = {
        "output": {},
        # should be a list of files created -ross
        "success": True
    }
    return computation_output


def start(parsed_args):
    phase_key = list(listRecursive(parsed_args, 'computation_phase'))

    if 'local_1' in phase_key:
        computation_output = remote_1(parsed_args)
        return computation_output
    else:
        raise ValueError("Error occurred at Remote")
