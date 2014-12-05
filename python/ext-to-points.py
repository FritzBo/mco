#!/usr/bin/python2.7
'''
Created on 3.12.2014

@author: fritz
'''

from sys import argv

def read_ext(file_name, sep):
    input_file = open(file_name, "r")

    if input_file.readline() != "V-representation\n":
        print "Not a valid cdd ext file."
        return

    if input_file.readline() != "begin\n":
        print "Not a valid cdd ext file."
        return

    descriptor = input_file.readline()
    descriptor = descriptor.split(sep)
    rows = int(descriptor[0])
    cols = int(descriptor[1])

    print rows - cols + 1, "points."

    for line in input_file:
        line.strip();
        elements = line.split(sep)

        if elements[0] == "1":
            str = elements[1]
            for i in xrange(2, len(elements)):
                str += sep + elements[i]
            print str,


if __name__ == '__main__':
    if len(argv) != 3:
        print "Syntax:", argv[0], "<input file> <separator>"
        print "<separator> can be one of: t s"
    else:
        input_file_name = argv[1]
        separator = argv[2]

        if separator == "t":
            separator = "\t"
        if separator == "s":
            separator = " "

        read_ext(input_file_name, separator);