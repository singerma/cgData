#!/usr/bin/python

import cgData, operator

gafHeaders = [
    "Entry Number", "FeatureID", "FeatureType", "FeatureDBSource",
    "FeatureDBVersion", "FeatuerDBDate", "FeatureSeqFileName", "Composite",
    "CompositeType", "CompositeDBSource", "CompositeDBVersion",
    "CompositeDBDate", "AlignmentType", "FeaturedCoordinates",
    "CompositeCoordinates", "Gene", "GeneLocus", "FeatureAliases",
    "FeatureInfo"]

gafCols = [
    "entryNumber", "featureID", "featureType", "featureDBSource",
    "featureDBVersion", "featureDBDate", "featureSeqFileName", "composite",
    "compositeType", "compositeDBSource", "compositeDBVersion",
    "compositeDBDate", "alignmentType", "featureCoordinates",
    "compositeCoordinates", "gene", "geneLocus", "featureAliases",
    "featureInfo"]


def parseFeatureCoordinates(featureCoordinates):
    featureCoordinatesList = featureCoordinates.split(";")
    accumulator = []
    for i in featureCoordinatesList:
        i = i.split(",")
        i = map(lambda x: tuple(x.split("-")), i)
        accumulator.append(i)
    return accumulator

def parseChromCoordinates(chromCoordinates):
    chromCoordinatesList = chromCoordinates.split(";")
    chromCoordinatesList = map(lambda x: x.split(":"), chromCoordinatesList)
    unzippedList = zip(*chromCoordinatesList)
    unzippedList[1] = map(parseFeatureCoordinates, unzippedList[1])
    assert(reduce(operator.mul, map(len, unzippedList[1])) == 1) #make sure all the lists have len 1
    unzippedList[1] = unzippedList[1][0]
    return zip(*unzippedList)


class gafLine:
    def __init__(
        self, entryNumber, featureID, featureType, featureDBSource,
        featureDBVersion, featureDBDate, featureSeqFileName, composite,
        compositeType, compositeDBSource, compositeDBVersion, compositeDBDate,
        alignmentType, featureCoordinates, compositeCoordinates, gene,
        geneLocus, featureAliases, featureInfo):
        
        self.entryNumber = entryNumber
        self.featureID = featureID
        self.featureType = featureType
        self.featureDBSource = featureDBSource
        self.featureDBVersion = featureDBVersion
        self.featureDBDate = featureDBDate
        self.featureSeqFileName = featureSeqFileName
        self.composite = composite
        self.compositeType = compositeType
        self.compositeDBSource = compositeDBSource
        self.compositeDBVersion = compositeDBVersion
        self.compositeDBDate = compositeDBDate
        self.alignmentType = alignmentType
        self.featureCoordinates = parseFeatureCoordinates(featureCoordinates)
        self.compositeCoordinates = parseChromCoordinates(compositeCoordinates)
        self.gene = gene
        self.geneLocus = geneLocus
        self.featureAliases = featureAliases
        self.featureInfo = featureInfo


class gaf(cgData.baseObject):
    def __init__(self):
        cgData.baseObject.__init__(self)
        self.gafData = []
    
    def read(self, handle):
        assert(handle.readline()[:-1].split("\t") == gafHeaders)
        for line in handle:
            line = line.rstrip("\n")
            splitLine = line.split("\t")
            assert(len(splitLine) == len(gafCols))
            gafData.append(gafLine(**dict(zip(gafCols, splitLine))))

    def __iter__(self):
        for i in self.gafData:
            yield i
