#!/usr/bin/python3

'''
Script for making experimentor's life easier.

Example:
{
    "executable" : "/home/fritz/programming/mco/build/cli/mco_cli",
    "threads" : 15,
    "time_limit" : 3600,
    "memory_limit" : 16000,
    "output" : "test_output",
    "parameter" : [
        { "short" : "-cst", "type" : "switch" },
        { "short" : "ep-tsaggouris", "type" : "switch" },
        { "type" : "instance" },
        { "short" : "-e", "type" : "constants", "values": [2, 1.5, 1.25, 1.1] },
        { "short" : "-d", "type" : "switch" },
        { "short" : "-F", "type" : "instance_attached", "extension" : "front" }
    ],
    "instances" : [
            { "extension" : "graph", "folder" : "Grid/" },
            { "extension" : "graph", "folder" : "Random/" }
        ]
}
'''

from sys import argv
import json
from pathlib import Path
from itertools import product

# keyword definition
EXECUTABLE = "executable"
OUTPUT = "output"
THREADS = "threads"
TIME_LIMIT = "time_limit"       # seconds
MEMORY_LIMIT = "memory_limit"   # MB

PARAMETERS = "parameter"
SHORT = "short"
TYPE = "type"
INSTANCE = "instance"
SWITCH = "switch"
CONSTANTS = "constants"
VALUES = "values"
INSTANCE_ATTACHED = "instance_attached"

INSTANCES = "instances"
EXTENSION = "extension"
FOLDER = "folder"

def expand_parameters(parameters, instance_filename):
    parameter_strings = []

    for parameter in parameters:
        if parameter[TYPE] == SWITCH:
            parameter_strings += [[parameter[SHORT]]]

        elif parameter[TYPE] == CONSTANTS:
            values = []
            for value in parameter[VALUES]:
                values += [parameter[SHORT] + " " + str(value)]

            parameter_strings += [values]

        elif parameter[TYPE] == INSTANCE:
            parameter_strings += [[str(instance_filename)]]

        elif parameter[TYPE] == INSTANCE_ATTACHED:
            parameter_strings += [[parameter[SHORT] + " " + str(instance_filename) + "." + parameter[EXTENSION]]]

    for parameter_set in product(*parameter_strings):
        yield " ".join(parameter_set)


def expand_instances(instance_defs):
    for instance_def in instance_defs:
        folder = instance_def[FOLDER]
        extension = instance_def[EXTENSION]

        p = Path(folder)
        if not p.exists():
            raise Exception("Instance folder " + str(folder) + " does not exist.")

        for filename in list(p.glob('*.' + str(extension))):
            yield filename


def perform_experiment(experiment_def):
    executable = experiment_def[EXECUTABLE]
    instance_defs = experiment_def[INSTANCES]
    output_file_name = experiment_def[OUTPUT]
    parameters = []

    if PARAMETERS in experiment_def.keys():
        parameters = experiment_def[PARAMETERS]

    for instance_filename in expand_instances(instance_defs):
        for parameter_string in expand_parameters(parameters, instance_filename):
            print(executable, parameter_string, ">>", output_file_name)


def process_json(filename):
    json_file = open(filename, "r")
    json_content = json.load(json_file)

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