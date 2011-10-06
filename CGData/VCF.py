#!/usr/bin/python

import re

import CGData

#TODO: probably shouldn't be using assert statements

class VCFFormatError(Exception):

    def __init__(self, text):
        Exception.__init__(self, text)

class Filter:
    def __init__(self, id, description):
        self.id = id
        self.description = description

class Info:
    def __init__(self, id, number, type, description):
        self.id = id
        self.number = number
        self.type = type
        self.description = description

def parse_info(info):
    infoDict = dict()
    for i in info.split(";"):
        if "=" in i:
            id,value = i.split("=")
            infoDict[id]=value
        else:
            infoDict[i]=True
    return infoDict

class Data:
    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, genotype_data = None, headers = None):
        self.chrom = chrom
        self.pos = int(pos)
        self.id = id
        self.ref = ref
        self.alt = alt.split(",")
        self.qual = qual
        self.filter = filter.split(";")
        self.info = parse_info(info)
        
        self.genotype = dict()
        assert(len(genotype_data) == len(headers))
        genotypeInfo = genotype_data[0]
        assert(headers[0] == "FORMAT")
        for sample,genotypeValues in zip(headers[1:],genotype_data[1:]):
            self.genotype[sample] = dict()
            for i,g in zip(genotypeInfo.split(":"), genotypeValues.split(":")):
                #if i == "GT":
                #   if g == ".":
                #       self.genotype[sample][i] = ()
                #   elif "/" in g:
                #       self.genotype[sample][i] = tuple(map(int, g.split("/")))
                #   elif "|" in g:
                #       self.genotype[sample][i] = tuple(map(int, g.split("|")))
                #   else:
                #       self.genotype[sample][i] = g
                #else:
                #   self.genotype[sample][i] = g
                self.genotype[sample][i] = g


class VCF(CGData.CGDataSetObject):
    def __init__(self):
        self.meta = []
        self.infos = []
        self.filters = []
        self.formats = []
        self.headers = []
        self.data = []
        

    def add_info(self, value):
        #make re match spec
        match = re.match(r"""<ID=(.*),Number=(.*),Type=(.*),Description=['"](.*)['"]""", value)
        if match:
            self.infos.append(Info(*match.groups()))
        else:
            raise VCFFormatError("Improperly formatted INFO metadata:" + value)
    
    def add_filter(self, value):
        #make re match spec
        match = re.match(r"""<ID=(.*),Description=['"](.*)['"]>""", value)
        if match:
            self.filters.append(Filter(*match.groups()))
        else:
            raise VCFFormatError("Improperly formatted FILTER metadata: " + value)
        #self.filters.append(Filter(values[0], values[1]))

    def set_headers(self, header_list):
        self.headers = header_list

    def make_data(self, data_list):
        baseData = data_list[:8]
        genotypeData = data_list[8:]
        data = baseData+[genotypeData]+[self.headers[8:]]
        return Data(*data)
        
    def add_data(self, data_list):
        self.data.append(self.make_data(data_list))
        
    def read(self, handle):
        line = ""
        
        #read in header information
        #start with metadata
        for line in handle:
            if not line.startswith("##"):
                break
            key,value = line[2:].split("=",1)
            if key == "INFO":
                self.add_info(value.strip())
            elif key == "FILTER":
                self.add_filter(value.strip())
            else:
                self.meta.append(line[2:])
                #print >> sys.stderr, "WARN: unsupported meta-information line, skipping\n" + line
            #TODO: read in FORMAT information
            #TODO: read in other information?
    
        #if necessary (maybe always necessary?) read FORMAT lines
        headers = []
        #read in column labels
        if line.startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"):
            headers = line[1:].strip().split("\t")
            self.set_headers(headers)
        else:
            raise VCFFormatError("Improperly formatted header: " + headers)
    
        for line in handle:
            data_as_list = line.strip().split("\t")
            if len(data_as_list) == len(headers):
                self.data.append(self.make_data(data_as_list))
    
    def __iter__(self):
        for d in self.data:
            yield d
