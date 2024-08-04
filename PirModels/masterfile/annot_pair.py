#!/usr/bin/env python3

class AnnotPair:
    # __init__ with all the attributes optional
    def __init__(self, type=None, genename=None, startpos=None, endpos=None, direction=None, startline=None, endline=None, startmulticomment=None, endmulticomment=None, startlinenumber=None, endlinenumber=None, introntype=None, globaltype=None):
        self.type              = type
        self.genename          = genename
        self.startpos          = startpos
        self.endpos            = endpos
        self.direction         = direction
        self.startline         = startline
        self.endline           = endline
        self.startmulticomment = startmulticomment
        self.endmulticomment   = endmulticomment
        self.startlinenumber   = startlinenumber
        self.endlinenumber     = endlinenumber
        self.introntype        = introntype
        self.globaltype        = globaltype
