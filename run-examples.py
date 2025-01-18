import os
import subprocess
import sys


def run_example(example_name):
    print("====> Running {}".format(example_name))
    command = ["cargo", "run", "--example", example_name, "--features", "v2_7"]
    child = subprocess.Popen(command)
    child.communicate()
    if child.returncode != 0:
        print("Command `{}` failed with return code `{}`...".format(command, child.returncode))
        return False
    return True


def run_examples():
    ret = 0
    for example in [f for f in os.listdir('examples')]:
        if not example.endswith('.rs'):
            continue
        if not run_example(example[:-3]):
            ret = 1
    return ret


if __name__ == "__main__":
    sys.exit(run_examples())
