import json
import sys
from utils import listRecursive


def remote_1(args):

    computation_output = {
        "output": "Results files sent to remote",
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
