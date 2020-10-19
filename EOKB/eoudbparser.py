#!/usr/bin/env python3

# Written by Steven blanachard & Aneeta Uppal 
# used for parsing EOUdb compounds and oils from database scrape

import argparse
import io
import json
import logging
import gzip
from os import path

_LOGGER = logging.getLogger("%(prog)s")

class StoreList(argparse.Action):
    """
    This is a class implementing argparse.Action and is used to store a comma separate list into
    python array

    Attributes:
        parser ():
        args ():
        values ():
        option_string ():
    """
    # Ignore warning 'Too few public methods (1/2)' since we are implementing class 
    #pylint: disable=R0903
    # def __init__(self, option_strings, dest, nargs=None, **kwargs):
    #     super(StoreList, self).__init__(option_strings, dest kwargs)

    def __call__(self, parser, args, values, option_string=None):
        """
        Constructor

        Parameters:
            parser ():
            args ():
            values ():
            option_string ():
        """
        setattr(args, self.dest, values.split(','))

def handle_args():
    """
    Defines and processes command line arguments

    Returns:
        argparse.Namespace: processed arguments
    """
    parser = argparse.ArgumentParser(description="Boilerplate python script")

    # Defaults for command line options
    parser.set_defaults(verbosity=0)

    parser.add_argument("-v", "--verbose", action="count", dest="verbosity", \
                        help="increase output (can be specified multiple times)")
    parser.add_argument("-f", "--file", type=argparse.FileType('r'), dest="filelist", \
                        metavar="FILE", help="list of input files, one per line")
    parser.add_argument("-l", "--list", action=StoreList, dest="list", metavar="ITEM,...", \
                        help="comma separated list of values")
    parser.add_argument('--version', action='version', version='%(prog)s <version>')
    parser.add_argument("files", type=argparse.FileType('rb'), nargs='*', metavar="FILE", \
                        help="input file(s) to process.  Use '-' for stdin")

    args = parser.parse_args()

    # Process filelist, opening every file and adding file handle args.files.  We also expand
    # tilde (~) to the user's home directory in file paths
    if args.filelist is not None:
        if args.files is None:
            args.files = []
        for line in args.filelist:
            handle = open(path.expanduser(line), 'wb')
            args.files.append(handle)
        args.filelist.close()
        args.filelist = args.filelist.name

    # Check to see if any input files are gzip compressed and wrap the input stream to decompress
    # them on the fly
    if args.files is not None:
        for handle in args.files:
            if handle.name.endswith(".gz"):
                handle = gzip.open(handle)

    return args

def main():
    """Main Function"""
    args = handle_args()
    f= open("outputchem.txt","w+")
    d=open("outputoils.txt", "w+")

    if args.verbosity == 1:
        logging.basicConfig(level=logging.INFO)
    elif args.verbosity > 1:
        logging.basicConfig(level=logging.DEBUG)

    for handle in args.files:
        for chemical in json.load(io.TextIOWrapper(handle)):
            if 'name' in chemical:
                #print("%s\t%s" % (chemical['name'], chemical['id']))
                f.write("%s\t%s\t\n" % (chemical['name'], chemical['id']))
                for oil in chemical['oil']:
                    #print("\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (oil['name'], oil['id'], oil['name_botanical'], oil['abstract_author'], oil['abstract_title'], oil['abstract_publication'], oil['pivot']['compound_id'], oil['pivot']['percentage_average']))
                    d.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" % (oil['name'], oil['id'], oil['name_botanical'], oil['abstract_author'], oil['abstract_title'], oil['abstract_publication'], oil['pivot']['compound_id'], oil['pivot']['percentage_average']))
    # Close open file handles.  Realistically, it is better to close these when finished with them
    for handle in args.files:
        handle.close()

if __name__ == "__main__":
    main()
