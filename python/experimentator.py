#!/usr/bin/python3

'''
Script for making experimentor's life easier.

Example:
{
    "executable" : "/home/fritz/programming/mco/build/cli/mco_cli",
    "threads" : 15
    "time_limit" : 3600
    "memory_limit" : 16000
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
THREADS = "threads"
TIME_LIMIT = "time_limit"       # seconds
MEMORY_LIMIT = "memory_limit"   # MB

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
    idxs = {}
    for parameter in parameters:
        parameter_string = ""

        if parameter[TYPE] == SWITCH:
            parameter_string += parameter[SHORT] + " "

        elif parameter[TYPE] == CONSTANTS:

            if parameter[SHORT] not in idxs.keys():
                parameter[SHORT] = 0

            idxs[SHORT] += 1

            if idxs[SHORT] == len(parameter[VALUES]):
                idxs[SHORT] = 0

            parameter_string += parameter[SHORT] + " " + str(VALUES[idxs[SHORT]]) + " "


        elif parameter[TYPE] == INSTANCE_ATTACHED:

        yield parameter_string


def expand_instances(instance_defs):
    for instance_def in instance_defs:
        folder = instance_def[FOLDER]
        extension = instance_def[EXTENSION]

        p = Path(folder)
        if not p.exists():
            raise Exception("Instance folder" + str(folder) + "does not exist.")

        for filename in list(p.glob('*.' + str(extension))):
            yield filename


def perform_experiment(experiment_def):
    executable = experiment_def[EXECUTABLE]
    instance_defs = experiment_def[INSTANCES]
    output_file_name = experiment_def[OUTPUT]
    parameter_list = []

    if PARAMETERS in experiment_def.keys():
        parameters = experiment_def[PARAMETERS]

    for instance_filename in expand_instances(instance_defs):
        for parameter_string in expand_parameters(parameters):
            None


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

    if THREADS not in json_content.keys():
        json_content[THREADS] = 1

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