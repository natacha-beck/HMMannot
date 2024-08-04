# !/usr/bin/env python3

import sys
import os
import subprocess

#####################################################################
# Utility functions                                                 #
#####################################################################

def check_path(bashscript, searchfor, val):
    each_string = bashscript.split()
    for string in each_string:
        if searchfor not in string:
            continue
        string = string.replace(searchfor, val)
        if not os.path.exists(string) or not os.path.isfile(string):
            raise ValueError(f"Model file '{string}' doesn't exist or isn't readable.\nCheck installation of mfannot_models (see INSTALL.txt).")

#####################################################################
# MFannotGeneCommands class                                         #
#                                                                   #
# This class represents a single group of commands.                 #
#                                                                   #
# See also the encapsulating object, <MfAnnotCommandSet>            #
#####################################################################

class MFannotGeneCommands:
    def __init__(self):
        self.debug        = False
        self.genename     = None
        self.bashcommands = []

    def execute(self, cwd, tmp_dir, substitutions, commandlabel):
        if not cwd:
            raise ValueError("Error: no valid CWD supplied for Execute()... got '%s'." % cwd)
        if not tmp_dir:
            raise ValueError("Error: no valid TMP_DIR supplied for Execute()... got '%s'." % tmp_dir)

        debug        = self.debug
        genename     = substitutions["GENENAME"] or self.genename or "unk"

        bashcommands = self.bashcommands
        bashscript   = "\n".join(bashcommands) + "\n"

        # make copy; we will add some values right now
        moresubs             = substitutions.copy()
        if "GENENAME" not in moresubs:
            moresubs["GENENAME"] = genename
        if "DEBUG" not in moresubs:
            moresubs["DEBUG"] = "#" if debug else ""

        for var, val in moresubs.items():
            val = str(val)

            if not var.isalpha():
                raise ValueError("Not a proper variable name '%s'." % var)

            searchfor = "%" + var + "%"
            if bashscript.find(searchfor) >= 0 and var == "MODPATH":
                check_path(bashscript, searchfor, val)

            bashscript = bashscript.replace(searchfor, val)

        pid     = os.getpid()
        script  = "%s/bashcom.%s.%s.MfGC.%s" % (tmp_dir, genename, commandlabel, pid)
        outcapt = "%s/out.%s.%s.MfGC.%s"     % (tmp_dir, genename, commandlabel, pid)
        errcapt = "%s/err.%s.%s.MfGC.%s"     % (tmp_dir, genename, commandlabel, pid)

        # Write the script
        with open(script, "w") as ofh:
            ofh.write("#!/bin/bash\n")
            ofh.write("\n")
            ofh.write("# This script created automatically for gene '%s' (block '%s').\n" % (genename, commandlabel))
            ofh.write("\n")
            ofh.write(bashscript)
        ofh.close()

        command = "cd '%s';/bin/bash '%s' >'%s' 2>'%s'" % (cwd, script, outcapt, errcapt)

        if debug:
            debug_text  = "# DEBUG: Executing external script block '%s' for '%s'.\n" % (commandlabel, genename)
            debug_text += "# --- COMMAND : %s\n" % command
            debug_text += "# --- SCRIPT %s START ---\n" % script
            debug_text += bashscript
            debug_text += "\n"
            debug_text += "# --- SCRIPT %s END ---\n" % script

            sys.stderr.write(debug_text)

        subprocess.run(command, shell=True, check=True)

        return (outcapt, errcapt)

