#!/usr/bin/env python3

from annot_pair import AnnotPair

class MasterfileContig:
    def __init__(self):
        self.name         = None
        self.uniq_name    = None
        self.genetic_code = None
        self.namecomments = None
        self.annotations  = []

    def add_annot_pair(self, annot_pair):
        self.annot_pairs.append(annot_pair)

