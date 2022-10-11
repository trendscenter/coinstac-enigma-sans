import json
import os
import subprocess
import sys
import csv

from utils import listRecursive


def local_1(args):
    scriptDir = "/computation/enigma_scripts"
    scriptName = "SZ_SANSReg_AllComps_04122022.R"
    RScriptDir = "/usr/bin/Rscript"

    fileKeys = [
        'metr_FA', 'metr_AD','metr_RD', 'metr_MD',
        'SANS', 'Covariates', 'CohortInfo'
    ]

    fileMap = {}
    for i in fileKeys:
        for j in args["input"]["data"]:
            if j.find(i) > -1:
                fileMap[i] = j

    file1 = fileMap[fileKeys[0]]
    file2 = fileMap[fileKeys[1]]
    file3 = fileMap[fileKeys[2]]
    file4 = fileMap[fileKeys[3]]
    file5 = fileMap[fileKeys[4]]
    file6 = fileMap[fileKeys[5]]
    file7 = fileMap[fileKeys[6]]
    regr_args = [
        RScriptDir,
        os.path.join(scriptDir, scriptName), args["state"]["baseDirectory"],
        args["state"]["transferDirectory"], file1, file2, file3, file4, file5, file6, file7
    ]

    result = subprocess.run(
        regr_args,
        text=True,
        capture_output=True)

    if result.returncode != 0:
        raise Exception("R script failed: " + result.stderr + "\n" + result.stdout)

    with open(os.path.join(args["state"]["baseDirectory"], file7)) as f:
        reader = csv.reader(f)
        cohortDirectory = "output_sz_sans_factors_" + next(reader)[0].replace(" ", "_")

    output_dict = {
      "cohortDirectory": cohortDirectory,
      "computation_phase": "local_1"
    }

    computation_output = {"output": output_dict}

    return computation_output


def start(parsed_args):
    phase_key = list(listRecursive(parsed_args, 'computation_phase'))

    if not phase_key:
        computation_output = local_1(parsed_args)
        return computation_output
    else:
        raise ValueError("Error occurred at Local")
