#!/usr/bin/python3

'''
Script for making experimentor's life easier.

Example:
{
    "executable" : "/home/fritz/programming/mco/build/cli/mco_cli",
    "output" : "test_output",
    "parameter" : [
        { "short" : "-cst", "type" : "switch" },
        { "short" : "ep-tsaggouris", "type" : "switch" },
        { "short" : "-e", "type" : "constants", "values": [2, 1.5, 1.25, 1.1] },
        { "short" : "-d", "type" : "switch" },
        { "short" : "-F", "type" : "instance_attached", "extension" : "front" }
    ]
    "instances" : [
            { "extension" : "graph", "folder" : "Grid/" },
            { "extension" : "graph", "folder" : "Random/" }
        ]
}
'''

from sys import argv
import json
from pathlib import Path

# keyword definition
EXECUTABLE = "executable"
OUTPUT = "output"

PARAMETERS = "parameters"
SHORT = "short"
TYPE = "type"
SWITCH = "switch"
CONSTANTS = "constants"
VALUES = "values"
INSTANCE_ATTACHED = "instance_attached"

INSTANCES = "instances"
EXTENSION = "extension"
FOLDER = "folder"

def expand_parameters(parameters):
    


def find_instances(instance_defs):
    instances = []
    for instance_def in instance_defs:
        folder = instance_def[FOLDER]
        extension = instance_def[EXTENSION]

        p = Path(folder)
        if not p.exists():
            raise Exception("Instance folder", folder, "does not exist.")

        instances += list(p.glob('*.' + str(extension)))

    return instances


def perform_experiment(experiment_def):
    executable = experiment_def[EXECUTABLE]
    instance_defs = experiment_def[INSTANCES]
    output_file_name = experiment_def[OUTPUT]
    parameter_list = []

    if PARAMETERS in experiment_def.keys():
        parameters = experiment_def[PARAMETERS]
        parameter_list = expand_parameters(parameters)

    instance_list = find_instances(instance_defs)


def process_json(filename):
    json_file = open(filename, "r")
    json_content = json.loads(json_file)

    # Sanity checking
    if EXECUTABLE not in json_content.keys():
        raise Exception("No executable stated.")

    if INSTANCES not in json_content.keys():
        raise Exception("No instances specified.")

    if OUTPUT not in json_content.keys():
        raise Exception("No output filename specified.")

    instances = json_content["instances"]

    for instance in instances:

        if EXTENSION not in instance.keys():
            raise Exception("No instance extension.")

        if FOLDER not in instance.keys():
            raise Exception("No instance folder specified.")

    json_file.close()

    return json_content


if __name__ == '__main__':
    if len(argv) != 2:
        print("Usage:", argv[0], "<json file>")
    else:
        json_content = process_json(argv[1])
        perform_experiment(json_content)