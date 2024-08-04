#! /usr/bin/env python3

import os

#
# This class represents the content of a configuration
# file used by mfannot; that file stores a list of external
# programs (with their arguments) that are used to add
# further annotations; the output of the programs are
# expected to be a series of XML documents of type
# AnnotPairCollection.
#

class MFannotGeneCommandSet:
    def __init__(self):
        self.debug       = False
        self.genename    = None
        self.percontig   = None
        self.filerank    = None
        self.commandset  = {}

    def set_debug(self, debug_value):
        if not debug_value:
            return self.debug
        commandset = self.commandset
        for genecom in commandset.values():
            genecom.debug = debug_value
        self.debug = debug_value

    def execute(self, cwd, tmp_dir, substitutions):
        # cwd           == working directory
        # tmp_dir       == where we can create temp files
        # substitutions == hash VAR => val to substitute in place of the string %VAR% in bashcommands

        # if cwd is not a valid directory, raise an error
        if not cwd or not os.path.isdir(cwd):
            raise ValueError(f"Error: no valid CWD supplied for Execute()... got '{cwd}'.")
        # if tmp_dir is not a valid directory, raise an error
        if not tmp_dir or not os.path.isdir(tmp_dir):
            raise ValueError(f"Error: no valid TMPDIR supplied for Execute()... got '{tmp_dir}'.")

        debug = self.debug or False

        commandset  = self.commandset
        commandkeys = sorted([key for key in commandset.keys() if key.isdigit()])

        outfilename = substitutions["OUTFILE"] or ""
        foundoutfile = False
        outcapt     = "(none)"
        errcapt     = "(none)"

        for commandkey in ["pre", *commandkeys, "post"]:
            genecom = commandset[commandkey] if commandkey in commandset else None
            if not genecom:
                continue
            if foundoutfile and commandkey.isdigit():
                continue

            commandlabel = ""
            if commandkey == "pre":
                commandlabel = "A-Pre"
            elif commandkey == "post":
                commandlabel = "C-Post"
            else:
                commandlabel = f"B-{commandkey}"

            outcapt, errcapt = genecom.execute(cwd, tmp_dir, substitutions, commandlabel)
            if os.path.isfile(outfilename) and os.path.getsize(outfilename) > 0:
                foundoutfile = True

        return outcapt, errcapt
