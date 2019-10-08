import json
import os
import subprocess
import sys
from utils import listRecursive


def local_1(args):
    # Running the first script
    scriptDir = "/computation/enigma_scripts"
    scriptName = "SZ_CortReg_20160420.R"
    RScriptDir = "/usr/bin/Rscript"

    regr_args = [RScriptDir, os.path.join(scriptDir, scriptName)]
    subprocess.call(regr_args, cwd=args["state"]["transferDirectory"],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL)

    output_dict = {
        "computation_phase": "local_1"
    }
    
    computation_output = {"output": output_dict}
    
    return json.dumps(computation_output)


if __name__ == '__main__':

    parsed_args = json.loads(sys.stdin.read())
    phase_key = list(listRecursive(parsed_args, 'computation_phase'))

    if not phase_key:
        computation_output = local_1(parsed_args)
        sys.stdout.write(computation_output)
    else:
        raise ValueError("Error occurred at Local")
