#!/usr/bin/env python3

import pdb
from pprint import pprint

import argparse
import os
import shutil
import subprocess

from masterfile import Masterfile

##############################################################################################
# Utilities Methods                                                                          #
##############################################################################################

def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

##############################################################################################
# Set up the command line arguments                                                          #
##############################################################################################

parser = argparse.ArgumentParser(
    prog='HMMannot',
    description='HMMannot: A pipeline for the annotation of genes in a genome.',
    epilog='That\'s all folks!',
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument('masterfile', help='The masterfile to be annotated.')
parser.add_argument('-g', '--genetic',
                    type=int,
                    default=1,
                    choices=[1, 2, 3, 4, 4, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23],
                    help="""Genetic code used.
---------------------------------
1 =>  Standard (default)
2 =>  Vertebrate Mitochondrial                  AGA=Ter(*),AGG=Ter(*),AUA=Met(M),UGA=Trp(W)
3 =>  Yeast Mitochondrial                       ATA=Met(M),CTN=Thr(T),TGA=Trp(W)
4 =>  Mold Mitochondrial                        TGA=Trp(W)
5 =>  Invertebrate Mitochondrial                AGA=Ser(S),AGG=Ser(S),ATA=Met(M),TGA=Trp(W)
6 =>  Ciliate Dasycladacean Hexamita Nuclear    TAA=Gln(Q),TAG=Gln(Q)
9 =>  Echinoderm Flatworm Mitochondrial         AAA=Asn(N),AGA=Ser(S),AGG=Ser(S),TGA=Trp(W)
10 => Euplotid Nuclear                          TGA=Cys(C)
11 => Bacterial and Plant Plastid
12 => Alternative Yeast Nuclear                 CTG=Ser(S)
13 => Ascidian Mitochondrial                    AGA=Gly(G),AGG=Gly(G),ATA=Met(M),TGA=Trp(W)
14 => Alternative Flatworm Mitochondrial        AAA=Asn(N),AGA=Ser(S),AGG=Ser(S),TAA=Tyr(Y),TGA=Trp(W)
15 => Blepharisma Macronuclear                  TAG=Gln(Q)
16 => Chlorophycean Mitochondrial               TAG=Leu(L)
21 => Trematode Mitochondrial                   TGA=Trp(W),ATA=Met(M),AGA=Ser(S),AGG=Ser(S)
22 => Scenedesmus Obliquus Mitochondrial        TCA=Stop(*),TAG=Leu(L)
23 => Thraustochytrium Mitochondrial            TTA=Stop(*)""")

parser.add_argument('-T', '--tmpdir',type=dir_path, help='Temporary directory to be used.')

parser.add_argument('-d', '--debug', action='store_true', help='Print debugging information.')

parser.add_argument('--minorflen',
                    type=int,
                    default=40,
                    help="""This option allows the user to choose the size of the minumum ORFs
(in amino acids) that are produced using Esl-translate. This value must be
an integer.""")

parser.add_argument('--eslstrand',
                    choices=['watson', 'crick'],
                    help="""The strand of the sequence.""")

args = parser.parse_args()

##############################################################################################
# Should run before any step                                                                 #
##############################################################################################

#--------------------------------#
# Create the temporary directory #
#--------------------------------#

# If args.tmpdir is not None, then set TMPDIR to args.tmpdir
TMPDIR = None
if args.tmpdir:
    TMPDIR = args.tmpdir
else:
    TMPDIR = f"/tmp/HMMannot.{os.getpid()}"
if args.debug:
    print(f"TMPDIR: {TMPDIR}\n")
os.mkdir(TMPDIR, 0o700)

#---------------------#
# Read the masterfile #
#---------------------#

pirmaster = Masterfile()
pirmaster.object_from_masterfile(args.masterfile)
pirmaster.clean_pirmaster(TMPDIR)

#-----------------------------#
# Create the raw contigs file #
#-----------------------------#

raw_contigs_file = f"{TMPDIR}/contigs.fna"
# Write the contigs to a file
contigs   = pirmaster.contigs
for contig in contigs:
    contigname = contig.name
    seq        = contig.sequence
    with open(raw_contigs_file, "w") as ofh:
        ofh.write(f">{contigname}\n{seq}\n\n")

ofh.close()


##############################################################################################
# Main steps                                                                                 #
##############################################################################################


# TODO


##############################################################################################
# Should run after all step                                                                  #
##############################################################################################

# Remove the temporary directory
# if not in debug mode
if not args.debug:
    print(f"Removing temporary directory: {TMPDIR}\n")
    os.rmtree(TMPDIR)
