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

    computation_output = {
        "output": "Results files sent to remote", # should be a list of files created -ross
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
