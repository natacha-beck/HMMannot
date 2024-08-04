#!/usr/bin/env python3

import re
from mfannot_gene_command_set import MFannotGeneCommandSet
from mfannot_gene_commands    import MFannotGeneCommands

#
# This class represents the content of a configuration
# file used by mfannot; that file stores a list of external
# programs (with their arguments) that are used to add
# further annotations; the output of the programs are
# expected to be a series of XML documents of type
# AnnotPairCollection.
#
# Sample format of text file
#
# # comment
# GeneName=rnl
# PerContig=true        # optional line; not yet implemented
# PreCommands           # optional block
#   blah blah blah line 1
#   blah blah blah line 2
# EndCommands
# Commands                 # one or many such blocks
#   blah blah blah line 1
#   blah blah blah line 2
# EndCommands
# Postcommands          # optional block
#   blah blah blah line 1
#   blah blah blah line 2
# EndCommands
#

class MFannotExternalProgs:
    def __init__(self):
        self.filename  = None
        self.geneprogs = {}

    def import_from_text_file(self, filename):
        self.filename = filename

        lines = []
        with open(filename, 'r') as f:
            line = f.readlines()
            f.close()

        filerank = 0
        genecoms = {}

        lines    = []
        with open(filename, 'r') as file:
            lines = file.readlines()  # Read all lines into a list


        while lines:
            line = lines.pop(0)
            if not line.strip() or line.strip().startswith('#'):
                continue

            # Expects "genename"
            genename_match = re.match(r'^\s*genenames?\s*=\s*(\w+(\s*,\s*\w+)*)\s*$', line, re.IGNORECASE)
            if not genename_match:
                raise ValueError(f"Error: unparsable line in '{filename}' (expected \"genename=\"), got:\n{line}")

            # potentially a list, like "abc,def"
            genename = genename_match.group(1)
            genename = re.sub(r'\s+', '', genename)
            if genename in genecoms:
                raise ValueError(f"Error: genename '{genename}' seen more than once in file '{filename}'.")

            genecom            = MFannotGeneCommandSet()
            genecom.genename   = genename
            genecom.filerank   = filerank
            genecom.commandset = {}

            commandset         = genecom.commandset

            match_ignore = re.match(r'^\s*$|^\s*#', line)
            while lines and match_ignore:
                lines.pop(0)

            # Parse other optional fields (right now, only "percontig")
            optional_match = re.match(r'^\s*(percontig|otherkeyTODO)\s*=\s*(\S+)\s*$', line, re.IGNORECASE)
            while lines and optional_match:
                key = optional_match.group(1)
                val = optional_match.group(2)
                # boolean fields
                bool_match = re.match(r'^(|f|false|0|no|n)$', val, re.IGNORECASE)
                if bool_match:
                    val = 0
                else:
                    val = 1
                key = key.lower()
                genecom.__dict__[key] = val

                match_ignore = re.match(r'^\s*$|^\s*#', file[0])
                while lines and match_ignore:
                    lines.pop(0)

            command_counter = 0
            while lines:
                line = lines.pop(0)
                while lines and re.match(r'^\s*$|^\s*#', line):
                    line = lines.pop(0)
                if not lines:
                    break

                blocktype_match = re.match(r'^\s*(pre|post)?commands\s*$', line, re.IGNORECASE)
                if not blocktype_match:
                    continue
                blocktype = blocktype_match.group(1) or ''
                blocktype = blocktype.lower()
                line = lines.pop(0)

                coms = []
                # read commands until "endcommands"
                while lines and not re.match(r'^\s*endcommands\s*$', line, re.IGNORECASE):
                    line = line.rstrip()
                    # make sure continuation lines don't have trailing blanks
                    if re.match(r'#\\\s+$', line):
                        line = re.sub(r'\s+$', '', line)
                    coms.append(line)
                    line = lines.pop(0)

                # Create command object
                command_obj              = MFannotGeneCommands()
                command_obj.genename     = genename
                command_obj.bashcommands = coms

                command_key = None
                if blocktype:
                    command_key = blocktype
                else:
                    command_key      = str(command_counter)
                    command_counter += 1

                if command_key in commandset:
                    raise ValueError(f"Error in '{filename}': duplicate command block '{command_key}' for gene '{genename}'.")

                if not coms or all(re.match(r'^\s*$|^\s*#', com) for com in coms):
                    raise ValueError(f"Error in '{filename}', gene '{genename}', command block '{command_key}': all commands are blank?!?")

                # Exception if all items in coms are blank
                # XXXXX TODO: need to be tested XXXXX
                if not coms or all(re.match(r'^\s*$|^\s*#', com) for com in coms):
                    raise ValueError(f"Error in '{filename}', gene '{genename}', command block '{command_key}': all commands are blank?!?")

                commandset[command_key] = command_obj

                if not lines:
                    lines.insert(0, "(EOF)\n") # for error messages

                if not re.match(r'^\s*endcommands\s*$', line, re.IGNORECASE):
                    raise ValueError(f"Error: unparsable line in '{filename}' (expected \"endcommands\" keyword), got:\n{line}")

                if lines and not "0" in commandset:
                    raise ValueError(f"Error in '{filename}': unparsable line in '{genename}':\n{line}\n")

                genecoms[genename] = genecom
                break

        self.geneprogs = genecoms
