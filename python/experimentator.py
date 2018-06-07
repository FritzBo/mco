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
        { "short" : "-F", "type" : "instance_attached", "extension" : "front", "attach" : "tsag", depth = 3}
    ],
    "instances" : [
            { "extension" : "graph", "folder" : "Grid/K*/" },
            { "extension" : "graph", "folder" : "Random/K*/" }
        ]
}
'''

from sys import argv
import json
from pathlib import Path
from itertools import product
from multiprocessing import Pool
from subprocess import Popen
import resource
import functools

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
ATTACH = "attach"
DEPTH = "depth"

INSTANCES = "instances"
EXTENSION = "extension"
FOLDER = "folder"


def limit_thread(time_limit, memory_limit):
    if time_limit is not None:
        resource.setrlimit(resource.RLIMIT_CPU, (time_limit, time_limit))
    if memory_limit is not None:
        resource.setrlimit(resource.RLIMIT_RSS, (memory_limit, memory_limit)) # VMEM needs Posix
        # system


def execute_worker(parameter_tuple, time_limit, memory_limit):
    executable = parameter_tuple[0]
    parameters = parameter_tuple[1]
    output_file_name = parameter_tuple[2]

    limit_thread_wrapper = functools.partial(limit_thread, time_limit=time_limit, memory_limit=memory_limit)

    output_file = open(output_file_name, "w")

    print([executable] + parameters)

    # cmd = Popen([executable] + parameters, preexec_fn=limit_thread_wrapper, stdout=output_file)
    #
    # cmd.wait()


def expand_parameters(parameters, filename_iterator):
    # visit parameters to collect index sets

    index_sets = []

    i = 0
    for parameter in parameters:
        if parameter[TYPE] == CONSTANTS:
            index_sets.append(range(len(parameter[VALUES])))
            parameter["index"] = i
            i += 1

    for indices in product(*index_sets, filename_iterator):

        parameter_string = []

        instance_filename = ""
        instance_attachement = ""

        for parameter in parameters:
            if parameter[TYPE] == SWITCH:
                parameter_string += [parameter[SHORT]]

            elif parameter[TYPE] == CONSTANTS:
                parameter_string += [parameter[SHORT] + " " + str(parameter[VALUES][indices[parameter[
                    "index"]]])]

                instance_attachement += parameter[SHORT].strip("-") + str(parameter[VALUES][
                                                                               indices[
                    parameter[
                    "index"]]]) + "."

            elif parameter[TYPE] == INSTANCE:
                instance_file_url = indices[-1]
                parameter_string += [str(instance_filename)]

            elif parameter[TYPE] == INSTANCE_ATTACHED:
                instance_filename = "/".join(instance_file_url.parts[-(parameter[DEPTH] + 1):])
                parameter_string += [parameter[SHORT] + " " + str(instance_filename) + "." +
                                     parameter[ATTACH] + "." +
                                     instance_attachement +
                                     parameter[EXTENSION]]

        yield list(parameter_string)


def expand_instances(instance_defs):
    for instance_def in instance_defs:
        folder = instance_def[FOLDER]
        extension = instance_def[EXTENSION]

        if "*" in folder:
            paths = Path(folder[:folder.find('*') - 1]).glob(folder[folder.rfind('*'):])
        else:
            paths = [Path(folder)]

        # if not paths.first.is_dir():
        #     raise Exception("Instance folder " + str(folder) + " does not exist.")

        print([x for x in paths])
        for path in paths:
            for filename in list(path.glob('/*.' + str(extension))):
                yield filename


def perform_experiment(experiment_def):
    executable = experiment_def[EXECUTABLE]
    instance_defs = experiment_def[INSTANCES]
    output_file_name = experiment_def[OUTPUT]
    parameters = []

    if PARAMETERS in experiment_def.keys():
        parameters = experiment_def[PARAMETERS]

    time_limit = None
    if TIME_LIMIT in experiment_def.keys():
        time_limit = experiment_def[TIME_LIMIT]

    memory_limit = None
    if MEMORY_LIMIT in experiment_def.keys():
        memory_limit = experiment_def[MEMORY_LIMIT] * 1000 * 1000

    execute_worker_wrapper = functools.partial(execute_worker, time_limit=time_limit,
                                        memory_limit=memory_limit)

    thread_pool = Pool(processes=experiment_def[THREADS])

    for x in thread_pool.imap_unordered(execute_worker_wrapper,
                    product([executable],
                            expand_parameters(parameters, expand_instances(instance_defs)),
                            [output_file_name]
                            )
                    ):
        None

    thread_pool.close()

    thread_pool.join()

    # for parameter_tuple in product([executable], expand_parameters(parameters, expand_instances(
    #         instance_defs)), [output_file_name]):
    #     execute_worker_wrapper(parameter_tuple)


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