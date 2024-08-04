#!/usr/bin/env python3

from pprint import pprint
import pdb


import functools
import re

from masterfile_contig import MasterfileContig
from annot_pair        import AnnotPair


#####################################################################
# Utility functions                                                 #
#####################################################################

# This is a Python implementation of the Perl cmp function
def cmp(a, b):
    return (a > b) - (a < b) # return 1 if a > b, -1 if a < b, 0 otherwise

def sort_annots(a, b):
    # Modified by T. HOELLINGER

    # Simple case, nucleotide position is different:
    if a[0] != b[0]:
        return cmp(a[0], b[0])

    # Simple case, line numbers are available and > 0, so compare them
    if a[1] and b[1]:
        return cmp(a[1], b[1])

    # The new sorting function
    fulllinea   = " "; # G-genname ==> start ;; comment
    fulllineb   = " "; #    ||     ||
    namea       = " "; # the genename
    nameb       = " "; #   ||    ||
    arrowa      = " "; # ==> or <==
    arrowb      = " "; #  ||    ||
    startorenda = " "; # as it said start or end
    startorendb = " "; #         ||      ||

    if a[2] == "S":
        fulllinea = a[3].startline
    if a[2] == "E":
        fulllinea = a[3].endline
    if b[2] == "S":
        fulllineb = b[3].startline
    if b[2] == "E":
        fulllineb = b[3].endline

    #  Cut the line were there are ;;
    cutlinea = fulllinea.split(";;")
    cutlineb = fulllineb.split(";;")

    #  Cut the line and identify pattern
    # Only take the part before ;;, Delete comments
    linea = cutlinea[0] or ""
    #     ||        ||
    lineb = cutlineb[0] or ""

    linea = linea.strip()  # Delete the spaces before, at the beginning of the description
    lineb = lineb.strip()  #     ||        ||
    linea = linea.rstrip() # Delete the spaces at the end of the descriptio
    lineb = lineb.rstrip() #     ||        ||

    linea       = re.match(r'([\w|\_\-\d\(\)\?]+)\s*(==>|<==)\s*(start|end)', linea)
    namea       = linea.group(1) or "" # Get the genename
    arrowa      = linea.group(2) or "" # Ge                                            IsDef = 1 if not char else 0 "" # Get the start or the end

    lineb = re.match(r'([\w|\_\-\d\(\)\?]+)\s*(==>|<==)\s*(start|end)', lineb)
    nameb = lineb.group(1)       or ""
    arrowb = lineb.group(2)      or ""
    startorendb = lineb.group(3) or ""

    if arrowa == "==>" and arrowb == "==>"   and startorenda == "start" and startorendb == "start":
        return cmp(namea, nameb)
    elif arrowa == "==>" and arrowb == "==>" and startorenda == "end"   and startorendb == "end":
        return cmp(nameb, namea)
    elif arrowa == "==>" and arrowb == "==>" and startorenda == "start" and startorendb == "end":
        return cmp(a[2], b[2])  # end before the start
    elif arrowa == "==>" and arrowb == "==>" and startorenda == "end"   and startorendb == "start":
        return cmp(a[2], b[2])  # end before the start
    elif arrowa == "<==" and arrowb == "<==" and startorenda == "start" and startorendb == "start":
        return cmp(nameb, namea)
    elif arrowa == "<==" and arrowb == "<==" and startorenda == "end"   and startorendb == "end":
        return cmp(namea, nameb)
    elif arrowa == "<==" and arrowb == "<==" and startorenda == "start" and startorendb == "end":
        return cmp(b[2], a[2])  # start before the end
    elif arrowa == "<==" and arrowb == "<==" and startorenda == "end"   and startorendb == "start":
        return cmp(b[2], a[2])  # start before the end
    else:
        return cmp(fulllinea, fulllineb)


# Print the sequence in fasta format X characters per line
def fasta_block_to_fh(seq, pos, fh):
    while seq:
        if len(seq) > 60:
            subseq = seq[:60]
            fh.write(f"{pos:6d}  {subseq}\n")
            pos += len(re.findall(r'[acgtACGTnN]', subseq))
            seq = seq[60:]
            continue
        fh.write(f"{pos:6d}  {seq}\n")
        break

# Remove AP (annot_pair) from masterfile annot_pair list
def remove_AP(AP_to_rm=[],contig=MasterfileContig()):

    all_annots = contig.annotations
    for i in range(len(all_annots) - 1, -1, -1):
        contig_AP = all_annots[i]
        id_contig_AP = id(contig_AP)

        # go over the table containing annotation
        for id_rm_AP in AP_to_rm:
            # It mean start is already defined
            if id_rm_AP == id_contig_AP:
                all_annots.remove(contig_AP)


#####################################################################
# Masterfile class                                                 #
#####################################################################

class Masterfile:
    def __init__(self):
        self.masterfile = None
        self.filename   = None
        self.header     = []
        self.comment    = None
        self.contigs    = None


    def clean_pirmaster(self, tmpdir=None):
        AP_to_rm = []

        # Check each annot push all annot to remove on $AP_to_rm, changed the other one
        isUnique = {}
        count    = 0
        for contig in self.contigs:
            annotations = contig.annotations
            contigname  = contig.name
            comments    = contig.namecomments or ""
            header      = contigname + comments

            if header in isUnique:
                isUnique[header] += 1
            else:
                isUnique[header]  = 1

            if isUnique[header] > 1:
                raise ValueError(f"Two contigs have same header '{header}'\n")

            count += 1

            # Clean name and comments
            contig.name         = "contig" + str(count)
            contig.namecomments = ""

            seq                 = contig.sequence
            seq                 = seq.replace("!", "")
            seq                 = seq.upper()
            contig.sequence     = seq
            for annot_pair in annotations:
                type          = annot_pair.type
                startline     = annot_pair.startline or ""
                endline       = annot_pair.endline   or ""

                new_startline = None
                new_endline   = None

                id_AP = id(annot_pair)

                if type == "C":
                    if startline and not startline.startswith(";; mfannot:"):
                        continue
                    if endline and not endline.startswith(";; mfannot:"):
                        continue

                annot_pair.type = "C"

                if endline != None and startline != None:
                    if not startline.startswith(";; mfannot:"):
                        startline_match = re.match(r'(.+)\s*(;;.+)$', startline)
                        if startline_match:
                            if not re.search(r'mfannot:', startline_match.group(2)):
                                new_startline        = startline_match.group(2)
                                annot_pair.startline = new_startline
                                if new_startline == "":
                                    annot_pair.startpos = None

                    if not endline.startswith(";; mfannot:"):
                        endline_match = re.match(r'(.+)\s*(;;.+)$', endline)
                        if endline_match:
                            if not re.search(r'mfannot:', endline_match.group(2)):
                                new_endline        = endline_match.group(2)
                                annot_pair.endline = new_endline
                                if new_endline == "":
                                    annot_pair.endpos = None

                    if new_startline == None and new_endline == None:
                        AP_to_rm.append(id_AP)
                        continue

                elif endline == None:
                    if startline.startswith(";; mfannot:"):
                        AP_to_rm.append(id_AP)
                        continue
                    else:
                        new_startline = ";$startline ;; mfannot: no end found"

                    annot_pair.startline = new_startline
                    continue
                else:
                    if endline.startswith(";; mfannot:"):
                        AP_to_rm.append(id_AP)
                        continue
                    else:
                        new_endline = ";$endline ;; mfannot: no start found"

                    annot_pair.endline = new_endline
                    continue

        remove_AP(AP_to_rm, contig)
        self.object_to_masterfile(f"{tmpdir}/Masterfile_copy")


    def object_from_masterfile(self, filename, RemoveIupac=0):
        self.filename = filename

        with open(filename, 'r') as fh:
            lines = fh.readlines()
        lines = [line.strip() for line in lines]

        row            = 0
        header         = []
        comment_header = []
        after_header   = 0

        while row < len(lines) and not lines[row].startswith('>'):
            line = lines[row]
            if line.startswith(';; end mfannot'):
                row         += 1
                after_header = 1
                continue
            if after_header == 0:
                header.append(line)
                row += 1
            else:
                comment_header.append(line)
                row += 1

        self.header   = header
        MFheaderisDef = len(comment_header)
        if MFheaderisDef == 0:
            self.comment = header
        else:
            self.comment = comment_header

        contigs = []
        # Parse each contigs
        while True:
            if row >= len(lines):
                break

            contig = MasterfileContig()

            line = lines[row]

            row += 1
            match = re.match(r'^>\s*(\S+)(.*)', line)
            if not match:
                raise ValueError(f"Can't parse header line of masterfile '{filename}'. Line is:\n{line}\n")
            else:
                faname       = match.group(1)
                contig.name  = faname
                namecomments = match.group(2)
                if namecomments != "" and re.match(r'/trans\s*=\s*(\d+)/i', namecomments):
                    contig.geneticcode = namecomments.group(1)
                if namecomments != "":
                    contig.namecomments = namecomments

            seq        = ""
            annotexp   = {}
            allannots  = []
            seqpos     = 0 # in computer coordinates, not biological coordinates!
            linenumber = 0

            while row < len(lines):
                linenumber += 1
                line        = lines[row]
                row        += 1

                if not line.strip():
                    continue
                if line.startswith('>'):
                    row -= 1
                    break

                match = re.match(r'^\s*\d*\s*([^;].*)', line)
                if match:
                    dna = match.group(1)
                    if not dna.strip():
                        continue  # ignore blank lines
                    dna = dna.replace(" ", "").replace("\r", "").replace("\n", "").replace("\t", "")
                    if RemoveIupac:
                        if not re.match(r'^[acgtnN!]+$', dna):
                            raise ValueError(f"Bad characters in sequence line?!??\nLine: {line}\n")
                        seq    += dna  # includes special characters such as '!'
                        seqpos += len(re.findall(r'[acgtACGTnN]', dna))  # count without the '!'s.
                    else:
                        if not re.match(r'^[uUyrwskmbdhvxUYRWSKMBDHVXACGTacgtnN!]+$', dna):
                            raise ValueError(f"Bad characters in sequence line?!??\nLine: {line}\n")
                        seq    += dna  # includes special characters such as '!'
                        seqpos += len(re.findall(r'[uUyrwskmbdhvxUYRWSKMBDHVXacgtACGTnN]', dna))  # count without the '!'s.
                        seq     = seq.translate(str.maketrans('uUyrwskmbdhvxUYRWSKMBDHVX', 'TTNNNNNNNNNNNNNNNNNNNNNNN'))
                    continue

                if re.match(r'^;\s*([G])-(\S+)\s+(<==\*?|\*?==>)\s+(start|end|point)(.*)', line):
                    type, name, arrow, startend, comment = re.match(r'^;\s*([G])-(\S+)\s+(<==\*?|\*?==>)\s+(start|end|point)(.*)', line).groups()
                    intron_type = re.search(r'/group=(\S+)', line).group(1) if re.search(r'/group=(\S+)', line) else ""

                    multicomment = []
                    while row < len(lines):
                        if not line.startswith(';') or not lines[row-1].endswith('\\'):
                            break
                        multicomment.append(lines[row])
                        row += 1

                    pos = None
                    if arrow == '>':
                        if startend == 'start' or startend == 'point':
                            pos = seqpos + 1
                        else:
                            pos = seqpos
                    else:
                        if startend == 'end':
                            pos = seqpos + 1
                        else:
                            pos = seqpos

                    annotkey = f"{type}-{name}".lower()

                    ### structure of date : a hash called annotexp, contains a
                    ### table. This one contains an other hash for storing annotation
                    ### The last hash array has 3 keys :
                    ### 1 -> start : if it's defined, start pos is defined in annotation
                    ### 2 -> end : if it's defined, end pos is defined in annotation
                    ### 3 -> annot : contains the annotation itself

                    if startend == "point":
                        # Create a new AnnotPair
                        annot = AnnotPair(
                            type       = type,
                            introntype = intron_type,
                            genename   = name,
                            direction  = arrow
                        )
                        annot.startpos          = pos
                        annot.endpos            = pos
                        annot.startline         = line
                        annot.startmulticomment = multicomment
                        annot.startlinenumber   = linenumber
                        allannots.append(annot)
                    else:
                        if annotkey in annotexp:
                            table = annotexp[annotkey]
                            if startend == "start":
                                count      = 0
                                hasdefined = 0
                                # go over the table containing annotation
                                while count < len(table):
                                    infos = table[count]
                                    # here we get if it's start or stop
                                    state = infos['startend']
                                    # It mean start is already defined
                                    if state == 'start' or state == 'finish':
                                        count += 1
                                        continue
                                    # Means that end
                                    else:
                                        annot = infos['annotation']

                                        annot.startpos          = pos
                                        annot.startline         = line
                                        annot.startmulticomment = multicomment
                                        annot.startlinenumber   = linenumber

                                        allannots.append(annot)
                                        infos['startend'] = 'finish'
                                        hasdefined = 1
                                        break
                                # It means that annotation array has been running without finding an empty
                                # case, a corresponding annotation in existing annotation
                                if hasdefined == 0:
                                    newinfos = {}
                                    # It's just to define the thing
                                    newinfos['startend'] = 'start'
                                    annot = AnnotPair(
                                        type       = type,
                                        introntype = intron_type,
                                        genename   = name,
                                        direction  = arrow
                                    )
                                    annot.startpos          = pos
                                    annot.startline         = line
                                    annot.startmulticomment = multicomment
                                    annot.startlinenumber   = linenumber

                                    newinfos['annotation'] = annot
                                    table.append(newinfos)
                            else:
                                count      = 0
                                 # a variable that allows to
                                hasdefined = 0
                                # go over the table containing annotation
                                while count < len(table):
                                    infos = table[count]
                                    # here we get if it's start or stop
                                    state = infos['startend']
                                    # It mean start is already defined
                                    if state == 'end' or state == 'finish':
                                        count += 1
                                        continue
                                    # Means that end
                                    else:
                                        annot = infos['annotation']
                                        annot.endpos          = pos
                                        annot.endline         = line
                                        annot.endmulticomment = multicomment
                                        annot.endlinenumber   = linenumber

                                        allannots.append(annot)
                                        infos['startend'] = 'finish'
                                        hasdefined = 1
                                        break
                                # It means that annotation array has been running without finding an empty
                                # case, a corresponding annotation in existing annotation
                                if hasdefined == 0:
                                    newinfos = {}
                                    # It's just to define the thing
                                    newinfos['startend'] = 'end'
                                    annot = AnnotPair(
                                        type       = type,
                                        introntype = intron_type,
                                        genename   = name,
                                        direction  = arrow
                                    )
                                    annot.endpos          = pos
                                    annot.endline         = line
                                    annot.endmulticomment = multicomment
                                    annot.endlinenumber   = linenumber

                                    newinfos['annotation'] = annot
                                    table.append(newinfos)
                        else:
                            annot =  AnnotPair(
                                type       = type,
                                introntype = intron_type,
                                genename   = name,
                                direction  = arrow,
                            )
                            if startend == "start":
                                annot.startpos          = pos
                                annot.startline         = line
                                annot.startmulticomment = multicomment
                                annot.startlinenumber   = linenumber
                            if startend == "end":
                                annot.endpos          = pos
                                annot.endline         = line
                                annot.endmulticomment = multicomment
                                annot.endlinenumber   = linenumber
                            # hash array containing, start stop, and the annot
                            infos = {}
                            # It's just to define the thing
                            infos['startend']   = startend
                            infos['annotation'] = annot
                            table = []
                            # we put the hashing table in the table
                            table.append(infos)
                            annotexp[annotkey] = table
                    continue

                if line.startswith(';'):
                    multicomment = []
                    while row < len(lines):
                        if not lines[row].startswith(';') or not lines[row-1].endswith('\\'):
                            row += 1
                            break
                        multicomment.append(lines[row])

                    annot = AnnotPair(
                        type              = "C",  # a comment
                        startpos          = seqpos+1,
                        startline         = line,
                        startmulticomment = multicomment,
                        startlinenumber   = linenumber,
                    )

                    allannots.append(annot)
                    continue

                raise ValueError(f"Unexpected line:\n{line}\n")

            #  checking annotations
            for annot in allannots:
                if annot.type == 'G':
                    # For the genename, deleting the _X : the number of each copy
                    newgename = annot.genename
                    newgename = re.sub(r'\_\d+', '', newgename)

                    if re.match(r'-I\d+-(\S+)$', newgename):
                        # special case for G-cox1_2-I3-orf232
                        newgename = re.match(r'-I\d+-(\S+)$', newgename).group(1)

                    annot.genename = newgename

                    # If it's a gene
                    if re.match(r'Sig-(.+)$', annot.genename):
                        annot.genename = re.match(r'Sig-(.+)$', annot.genename).group(1)
                        annot.type = "S"
                    # If it's an exon
                    elif re.match(r'-E\d+$', annot.genename) or re.match(r'-E\d+-\S+$', annot.genename):
                        genenamecut = annot.genename.split("-")
                        annot.genename = genenamecut[0] if genenamecut[0] != "" else annot.genename
                        annot.type = "E"
                    # If it's an intron
                    elif re.match(r'-I\d+$', annot.genename) or re.match(r'-I\d+-\S+$', annot.genename):
                        genenamecut = annot.genename.split("-")
                        annot.genename = genenamecut[0] if genenamecut[0] != "" else annot.genename
                        annot.type = "I"
                    # If it's a tRNA
                    elif re.match(r'trn[\w|?]*\([\w|?]*\)', annot.genename):
                        annot.genename = re.match(r'trn([\w|?]*)\([\w|?]*\)', annot.genename).group(1)

            annotexp = None
            contig.sequence = seq

            contig.annotations    = allannots
            seqlen                = len(re.findall(r'[acgtACGTnN]', seq))
            contig.sequencelength = seqlen

            contigs.append(contig)

        self.contigs = contigs
        return self

    def object_to_masterfile(self, filename):
        # Print masterfile header
        fh     = open(filename, 'w')
        header = self.header or []

        # Remove trailing blank lines
        while len(header) and re.match(r'^\s*$', header[-1]):
            header.pop()

        if len(header):
            fh.write("\n".join(header) + "\n")

        contigs = self.contigs

        for contig in contigs:
            name         = contig.name
            namecomments = contig.namecomments or ""

            fh.write("\n\n")
            fh.write(f">{name}{namecomments}\n")

            posannots = []
            annots    = contig.annotations
            fullseq   = contig.sequence
            start_ac  = {}
            stop_ac   = {}

            for annot in annots:
                type  = annot.type
                if type != 'AC':
                    continue

                start = annot.startpos
                start_ac[start] = 1

                stop = annot.endpos
                stop_ac[stop] = 1

            # Num will be bio coords
            num = 0
            for i in range(len(fullseq)):
                c = fullseq[i]
                if c == '!':
                    continue

                # Num is bio coords
                num += 1
                if start_ac.get(num) or stop_ac.get(num):
                    if fullseq[i-1] != '!' and start_ac.get(num):
                        fullseq = fullseq[:i] + '!' + fullseq[i:]
                        i += 1
                    if fullseq[i] != '!' and stop_ac.get(num):
                        fullseq = fullseq[:i+1] + '!' + fullseq[i+1:]
                        i += 1

            for annot in annots:
                type  = annot.type
                start = annot.startpos
                if type == 'AC':
                    continue

                end   = annot.endpos
                # <== or ==> or <==* or *==> or undef
                arrow = annot.direction or ">"

                if arrow.startswith('>'):
                    if start:
                        posannots.append([start-1, annot.startlinenumber or 0, "S", annot])
                    if end:
                        posannots.append([end, annot.endlinenumber or 0, "E", annot])
                else:
                    if start:
                        posannots.append([start, annot.startlinenumber or 0, "S", annot])
                    if end:
                        posannots.append([end-1, annot.endlinenumber or 0, "E", annot])

            key_sort_annots = functools.cmp_to_key(sort_annots)
            posannots = sorted(posannots, key=key_sort_annots)

            seqpos  = 0
            charpos = 0
            for annotinfo in posannots:
                atpos, se, annot = annotinfo[0], annotinfo[2], annotinfo[3]
                name = annot.genename

                block    = ""
                blockpos = seqpos
                while seqpos < atpos:
                    char     = fullseq[charpos]
                    block   += char
                    charpos += 1

                    if char != '!':
                        seqpos += 1

                if fullseq[charpos] == "!" and len(block) >= 4 and "!" in block[-4:]:
                    block += "!"
                    charpos += 1
                if block:
                    fasta_block_to_fh(block, blockpos+1, fh)

                multicomment = [];
                if se == "S":
                    fh.write(f"{annot.startline}\n")
                    multicomment = annot.startmulticomment
                else:
                    fh.write(f"{annot.endline}\n")
                    multicomment = annot.endmulticomment
                if multicomment:
                    fh.write("\n".join(multicomment) + "\n")

            block = ""
            if charpos < len(fullseq):
                block = fullseq[charpos:]
            else:
                block = ""

            if block:
                fasta_block_to_fh(block, seqpos+1, fh)

        fh.close()

    # Utility routine that precomputes an internal hash the first
    # time it's called.
    def get_contig_by_name(self, contigname):
        if not self['!contigbyname!']:
            cache   = self['!contigbyname!'] = {}
            contigs = self.contigs or []
            for contig in contigs:
                name = contig.name
                if name:
                    cache[name] = contig

        self['!contigbyname!'][contigname]  # returns undef if not found.


    def sort_masterfile(self, direction):
        contigs = self.contigs
        contigs = sorted(contigs, key=lambda x: len(x.sequence), reverse=(direction == 1))
