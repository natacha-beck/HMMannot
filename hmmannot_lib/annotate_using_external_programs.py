# TODO: document

import os
from PirModels.mfannot_external_progs import MFannotExternalProgs


#####################################################################
# Utility functions                                                 #
#####################################################################

# Create a file for each contig
def create_contig_file(tmp_dir, pirmaster):
    # Create a directory to store the contigs
    contigs_dir = f"{tmp_dir}/Contigs"
    if not os.path.isdir(contigs_dir):
        os.mkdir(contigs_dir, 0o700)

    contigs = pirmaster.contigs
    for contig in contigs:
        seq       = contig.sequence.upper()
        seq       = seq.replace("!", "")
        name      = contig.name

        contig_file = f"{contigs_dir}/{name}.fna"
        with open(contig_file, "w") as CF:
            CF.write(f">{name}\n{seq}\n")

        CF.close()

# Sort the external programs by their filerank
def sorted_external_progs(mfannot_external_progs):
    # Retrieve all keys from the geneprogs dictionary
    keys = mfannot_external_progs.geneprogs.keys()

    # Define the sorting key using a lambda function
    sorting_key = lambda x: mfannot_external_progs.geneprogs[x].filerank

    # Sort the keys based on the filerank attribute
    sorted_keys = sorted(keys, key=sorting_key)

    # Return the sorted list of keys
    return sorted_keys

# Execute an external program
def execute_external_program(tmp_dir, genename, model_path, plainmasterfile, gencode, debug, commandobject, outfile):

    substitutions = {
        "OUTFILE":        outfile,                  # What will be created, a series of AnnotPairCollections
        "MODPATH":        model_path,
        "PLAINFASTAFILE": plainmasterfile,
        "TMPDIR":         tmp_dir,                  # Mfannot's tmp dir
        "GENENAME":       genename,                 # Optional; the Execute command use this but you can override it yourself)
        "GENCODE":        gencode,
        "DEBUG":          "#" if debug else "",     # See the config file text
    }

    commandobject.set_debug(debug)
    out_and_err_dir = f"{os.path.dirname(outfile)}/err_out"
    if not os.path.isdir(out_and_err_dir):
        os.mkdir(out_and_err_dir)
    outfile, errfile = commandobject.execute(out_and_err_dir, out_and_err_dir, substitutions)

    # TODO: Read back the AnnotPairCollections
    # Read back AnnotPairCollections;
    # with open(outfile_gene, "r") as infh:
        # annotpaircollections = PirObject.file_handle_to_object(infh)
    annotpaircollections = []

    return annotpaircollections

#####################################################################
# Main function                                                     #
#####################################################################

# Annotate the contigs using external programs
# This function is called from HMMannot.py
#
# The function reads the configuration file and parses the external programs
# to be executed. It then runs each external program command set.
#
def annotate_using_external_programs(conf_file, tmp_dir, pirmaster, masterfile, gencode, model_path, debug, *ext_selected):
    ext_selected = ext_selected or None

    mfannot_external_progs = MFannotExternalProgs()
    mfannot_external_progs.import_from_text_file(conf_file)

    # Here we build a list of external programs to run; by default all
    # of them are executed, but this can be modified by the option --ext_select, which
    # we parse here.

    # Create one file per contig
    work_dir = f"{tmp_dir}/Contigs"
    if not os.path.isdir(work_dir):
        create_contig_file(tmp_dir, pirmaster)

    # Default: all of them
    allprogs     = [prog for prog in mfannot_external_progs.geneprogs.keys()]
    allprogs     = {prog.lower(): 1 for prog in allprogs}
    progs_to_run = {prog.lower(): 1 for prog in allprogs}

    # What do we have in --ext_select? Parse it
    if ext_selected == "all":
        ext_selected = ""

    no_ext = [ext for ext in ext_selected if ext.startswith("no")]
    if no_ext:
        # There are some 'no', which means ( ALL minus the 'no's )
        for noext in no_ext:
            noext = noext[2:]
            del progs_to_run[noext.lower()]
    elif ext_selected:
        # There are no 'nos', which means --ext_select is the full list
        progs_to_run = {ext.lower(): 1 for ext in ext_selected}

    # Print warnings about unknown names in --ext_select
    for prog in progs_to_run:
        if prog == "none":
            continue
        if prog not in allprogs:
            print(f"Warning: unknown external program '{prog}' in '{conf_file}'")

    # Run each external program command set.
    # Now sorted; these are keys with possible MULTIPLE names, e.g. "rns,rnl"
    extcomsets = sorted_external_progs(mfannot_external_progs)

    dir_dict = {filename: os.path.join(work_dir, filename) for filename in os.listdir(work_dir)}
    for comsetnames in extcomsets:
        for genename in comsetnames.split(","):
            if genename.lower() not in progs_to_run:
                continue
            if genename.lower() != "intronI".lower() and genename.lower() != "intronII".lower():
                if not debug:
                    print(f"    '{genename}'...\n")
            for filename in dir_dict:
                if filename.startswith(".") or filename.endswith(".xml") or not os.path.isfile(os.path.join(work_dir, filename)):
                    continue
                commandobject = mfannot_external_progs.geneprogs[comsetnames]
                mf            = os.path.join(work_dir, filename)
                out           = f"{mf}-{genename}.xml"

                APC           = execute_external_program(tmp_dir, genename, model_path, mf, gencode, debug, commandobject, out)
                # TODO: Implement the following functions
                # &ReconstructFullRNA($genename,$APC) if $genename eq "rns" || $genename eq "rnl";
                # &AnnotateFromExternalAPC($genename,$APC);