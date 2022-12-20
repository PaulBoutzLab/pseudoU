#!/usr/bin/env python
# encoding: utf-8
#PLB 10/26/22


"""

Input format:
pkl from Annotator_v6, modified bedgraph

"""
"""
Modules

"""

import sys
import pickle
import copy
import pysam
from bedmodules import *
from bioTools import *
import subprocess
import numpy as np

#hg19 = pysam.Fastafile ("/home/paul/Desktop/SPLICE_GTF/hg19.fa")
#mm10 = pysam.Fastafile ("/home/paul/Desktop/SPLICE_GTF/mm10.fa")
#mm9 = pysam.Fastafile ("/home/paul/Desktop/SPLICE_GTF/mm9.fa")
hg38 = pysam.Fastafile ("/home/paul/Downloads/alpha_0.2.2_source/Genome/GRCh38.primary_assembly.genome.fa")


"""
Classes
"""

class Gene:
    def __init__(self,chrom,name,strand,geneintronlist):
        self.chrom = chrom
        self.name = name
        self.strand = strand
        self.geneintronlist = geneintronlist
        self.genetranscriptlist = []

    def __str__(self):
        return("GENE"+"\t"+self.name+" "+self.locus)


class ncRNA:
    def __init__(self,chrom, start, end, uid,counts,strand,genename,geneuid,geneclass,genestart,geneend,controlcounts):
        self.chrom = chrom
        self.strand = strand
        if self.strand == "+":
            self.start = int(start) + 1
            self.end = int(end) + 1
        else:
            self.start = int(start)
            self.end = int(end)
        self.uid = uid
        self.counts = counts
        self.genename = genename
        self.geneuid = geneuid
        self.geneclass = geneclass
        self.genestart = int(genestart)
        self.geneend = int(geneend)
        self.controlcounts = controlcounts
        self.genelength = self.geneend - self.genestart        
        self.locus = self.chrom+":"+str(self.start)+"-"+str(self.end)
        self.genelocus = self.chrom+":"+str(self.genestart)+"-"+str(self.geneend)
        #filter out end reads
        if self.strand == "+":
            if self.end - self.genestart > 5:
                self.endread = False
            else:
                self.endread = True
        elif self.strand == "-":
            if self.geneend - self.end > 5:
                self.endread = False
            else:
                self.endread = True
        #get sequences  
        if self.strand == "+":
                self.seqstart = self.start - 5
                self.seqend = self.start + 2
                self.cminusone = hg38.fetch(self.chrom,self.start - 3,self.start - 2)
                self.c = hg38.fetch(self.chrom,self.start - 2,self.start -1)
                self.cplusone = hg38.fetch(self.chrom,self.start - 1, self.start)
        else:
                self.seqstart = self.start - 2
                self.seqend = self.start + 5
                self.cminusone = hg38.fetch(self.chrom,self.start,self.start + 1)
                self.c = hg38.fetch(self.chrom,self.start + 1,self.start + 2)
                self.cplusone = hg38.fetch(self.chrom,self.start + 2,self.start + 3)
        self.sequence = hg38.fetch(self.chrom,self.seqstart,self.seqend)
        if self.strand == "-":
            self.sequence = reverseComplement(self.sequence, "DNA")
            self.cminusone = reverseComplement(self.cminusone, "DNA")
            self.c = reverseComplement(self.c, "DNA")
            self.cplusone = reverseComplement(self.cplusone, "DNA")

    def __str__(self):
        return ("ncRNA"+"\t"+self.genename+"\t"+self.geneclass+"\t"+self.uid+"\t"+self.locus+"\t"+self.sequence+"\t"+self.cminusone+"\t"+self.c+"\t"+self.cplusone+"\t"+self.counts+"\t"+self.controlcounts)


    def ncRNAbed(self):
        uslength = self.start - int(peakoffset) - self.genestart
        dslength = self.geneend - (self.end + int(peakoffset))
        if self.genelength <= windowsize:
            print "Short ncRNA:", str(uslength), str(dslength)
            bedchrom = self.chrom
            usbedstart = self.genestart
            usbedend = self.start - int(peakoffset)
            dsbedstart = self.end + int(peakoffset)
            dsbedend = self.geneend
            bedname = self.uid
            bedscore = self.counts
            bedstrand = self.strand
            minibed = open("TEMP"+outstops+".bed", "w")
            writeBed(bedchrom,str(usbedstart),str(usbedend),bedname,bedscore,bedstrand,minibed)
            writeBed(bedchrom,str(dsbedstart),str(dsbedend),bedname,bedscore,bedstrand,minibed)
            minibed.close()
            try:
                self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP"+outstops+".bed",bedgraph,bgcutoff)
                self.controlbgcoveragelist, self.controlbgmedian, self.controlbgpercentile = coverageBed("TEMP"+outstops+".bed",controlbedgraph,bgcutoff)
            except:
                self.bgcoveragelist = "NA"
                self.bgmedian = "NA"
                self.bgpercentile = "NA"
                self.controlbgcoveragelist = "NA"
                self.controlbgmedian = "NA"
                self.controlbgpercentile = "NA"
        elif uslength >= windowsize/2 and dslength >= windowsize/2:
            print "Us good ds good ncRNA:", str(uslength), str(dslength)
            bedchrom = self.chrom
            usbedstart = self.start - int(peakoffset) - windowsize/2
            usbedend = self.start - int(peakoffset)
            dsbedstart = self.end + int(peakoffset)
            dsbedend = self.end + int(peakoffset) + windowsize/2
            bedname = self.uid
            bedscore = self.counts
            bedstrand = self.strand
            minibed = open("TEMP"+outstops+".bed", "w")
            writeBed(bedchrom,str(usbedstart),str(usbedend),bedname,bedscore,bedstrand,minibed)
            writeBed(bedchrom,str(dsbedstart),str(dsbedend),bedname,bedscore,bedstrand,minibed)
            minibed.close()
            try:
                self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP"+outstops+".bed",bedgraph,bgcutoff)
                self.controlbgcoveragelist, self.controlbgmedian, self.controlbgpercentile = coverageBed("TEMP"+outstops+".bed",controlbedgraph,bgcutoff)
            except:
                self.bgcoveragelist = "NA"
                self.bgmedian = "NA"
                self.bgpercentile = "NA"
                self.controlbgcoveragelist = "NA"
                self.controlbgmedian = "NA"
                self.controlbgpercentile = "NA"
        elif uslength >= windowsize/2 and dslength < windowsize/2:
            print "Us good ds short ncRNA:", str(uslength), str(dslength)
            dsdeficit = windowsize/2 - dslength
            bedchrom = self.chrom
            try:
                usbedstart = self.start - int(peakoffset) - dsdeficit - windowsize/2
                usbedend = self.start - int(peakoffset)
            except:
                usbedstart = self.genestart
                usbedend = self.start - int(peakoffset)
            dsbedstart = self.end + int(peakoffset)
            dsbedend = self.geneend
            bedname = self.uid
            bedscore = self.counts
            bedstrand = self.strand
            minibed = open("TEMP"+outstops+".bed", "w")
            writeBed(bedchrom,str(usbedstart),str(usbedend),bedname,bedscore,bedstrand,minibed)
            writeBed(bedchrom,str(dsbedstart),str(dsbedend),bedname,bedscore,bedstrand,minibed)
            minibed.close()
            try:
                self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP"+outstops+".bed",bedgraph,bgcutoff)
                self.controlbgcoveragelist, self.controlbgmedian, self.controlbgpercentile = coverageBed("TEMP"+outstops+".bed",controlbedgraph,bgcutoff)
            except:
                self.bgcoveragelist = "NA"
                self.bgmedian = "NA"
                self.bgpercentile = "NA"
                self.controlbgcoveragelist = "NA"
                self.controlbgmedian = "NA"
                self.controlbgpercentile = "NA"
        elif uslength < windowsize/2 and dslength >= windowsize/2:
            print "Us short ds good ncRNA:", str(uslength), str(dslength)
            usdeficit = windowsize/2 - uslength
            bedchrom = self.chrom
            usbedstart = self.genestart
            usbedend = self.start - int(peakoffset)
            try:
                dsbedstart = self.end + int(peakoffset)
                dsbedend = windowsize/2 + self.end + int(peakoffset) + usdeficit
            except:
                dsbedstart = self.end + int(peakoffset)
                dsbedend = self.geneend
            bedname = self.uid
            bedscore = self.counts
            bedstrand = self.strand
            minibed = open("TEMP"+outstops+".bed", "w")
            writeBed(bedchrom,str(usbedstart),str(usbedend),bedname,bedscore,bedstrand,minibed)
            writeBed(bedchrom,str(dsbedstart),str(dsbedend),bedname,bedscore,bedstrand,minibed)
            minibed.close()
            try:
                self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP"+outstops+".bed",bedgraph,bgcutoff)
                self.controlbgcoveragelist, self.controlbgmedian, self.controlbgpercentile = coverageBed("TEMP"+outstops+".bed",controlbedgraph,bgcutoff)
            except:
                self.bgcoveragelist = "NA"
                self.bgmedian = "NA"
                self.bgpercentile = "NA"
                self.controlbgcoveragelist = "NA"
                self.controlbgmedian = "NA"
                self.controlbgpercentile = "NA"
        else:
            print "ncBed Error:", self.uid


class RTStop:
    def __init__(self,chrom, start, end, uid,counts,strand,gene,geneclass):
        self.chrom = chrom
        self.coordinate = start
        self.start = int(self.coordinate) - 2
        self.end = int(self.coordinate) + 5
        self.uid = uid
        self.gene = gene
        self.genename = gene.name
        self.geneclass = geneclass
        self.genestart = gene.start
        self.geneend = gene.end
        self.counts = counts
        self.strand = strand
        #filter out end reads
        if self.strand == "+":
            if self.end - self.genestart > 5:
                self.endread = False
            else:
                self.endread = True
        elif self.strand == "-":
            if self.geneend - self.end > 5:
                self.endread = False
            else:
                self.endread = True

        print self.chrom, str(self.start), str(self.end), str(self.uid), str(self.genename), str(self.geneclass), str(self.counts), str(self.strand)
        self.bgdist = []
        self.locus = self.chrom+":"+str(self.start)+"-"+str(self.start)+"/t"+str(self.strand)
        self.sequence = hg38.fetch(self.chrom,self.start,self.end)
        self.cminusone = hg38.fetch(self.chrom,int(self.coordinate),int(self.coordinate) + 1)
        self.c = hg38.fetch(self.chrom,int(self.coordinate) + 1,int(self.coordinate) + 2)
        self.cplusone = hg38.fetch(self.chrom,int(self.coordinate) + 2,int(self.coordinate) + 3)
        if self.strand == "-":
            self.sequence = reverseComplement(self.sequence, "DNA")
            self.cminusone = reverseComplement(self.cminusone, "DNA")
            self.c = reverseComplement(self.c, "DNA")
            self.cplusone = reverseComplement(self.cplusone, "DNA")

    def __str__(self):
        return ("RTSTOP"+"\t"+self.genename+"\t"+self.geneclass+"\t"+self.uid+"\t"+self.locus+"\t"+self.sequence+"\t"+self.cminusone+"\t"+self.c+"\t"+self.cplusone+"\t"+self.counts)

    def transcriptMaker(self):
        exonorder = sorted(self.gene.geneexonlist, key= lambda e: e.start, reverse = False)
        for e in exonorder:
            if e.consensus == True:
                if int(self.coordinate) >= int(e.start) and int(self.coordinate) <= int(e.end):
                    print "Exonic site", e , self.coordinate
                    i=exonorder.index(e)
                    print "Index", i
                    uslength = int(self.coordinate) - int(e.start)
                    dslength = int(e.end) - int(self.coordinate)
                    if uslength >= windowsize/2 and dslength >= windowsize/2:
                        bedchrom = e.chrom
                        bedstart = int(self.coordinate) - windowsize/2 -int(peakoffset)
                        bedend = int(self.coordinate) - int(peakoffset)
                        bedname = e.name
                        bedscore = "TypeI"
                        bedstrand = e.strand
                        minibed = open("TEMP.bed", "w")
                        bedfileout = open("Coding_RTStops_global.bed","a")
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,minibed)
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                        bedstart = int(self.coordinate) + int(peakoffset)
                        bedend = int(self.coordinate) + windowsize/2 + int(peakoffset)
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,minibed)
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                        minibed.close()
                        bedfileout.close()
                        self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP.bed",bedgraph,bgcutoff)
                        break
                    elif uslength >= windowsize/2 and dslength < windowsize/2:
                        bedchrom = e.chrom
                        bedstart = int(self.coordinate) - windowsize/2 - int(peakoffset)
                        bedend = int(self.coordinate) - int(peakoffset)
                        bedname = e.name
                        bedscore = "TypeII"
                        bedstrand = e.strand
                        minibed = open("TEMP.bed", "w")
                        bedfileout = open("Coding_RTStops_global.bed","a")
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,minibed)
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                        bedstart = int(self.coordinate) + int(peakoffset)
                        bedend = e.end
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,minibed)
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                        dsdeficit = windowsize/2 - dslength
                        while dsdeficit < windowsize/2:
                            i += 1
                            try:
                                nextexon = exonorder[i]
                            except:
                                next
                            if nextexon.consensus == False:
                                next
                            else:
                                exonlength = int(nextexon.end) - int(nextexon.start)
                                dslength += exonlength
                                if dslength >= windowsize/2:
                                    bedchrom = e.chrom
                                    bedstart = nextexon.start
                                    bedend = int(nextexon.start) + dsdeficit
                                    bedname = e.name
                                    bedscore = "TypeII"
                                    bedstrand = e.strand
                                    writeBed(bedchrom,bedstart,str(bedend),bedname,bedscore,bedstrand,minibed)
                                    writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                                    minibed.close()
                                    bedfileout.close()
                                    try:
                                        self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP.bed",bedgraph,bgcutoff)
                                        break
                                    except:
                                        print self
                                else:
                                    bedchrom = e.chrom
                                    bedstart = nextexon.start
                                    bedend = nextexon.end
                                    bedname = e.name
                                    bedscore = "TypeII"
                                    bedstrand = e.strand
                                    writeBed(bedchrom,bedstart,bedend,bedname,bedscore,bedstrand,minibed)
                                    writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                                    dsdeficit -= exonlength
                    elif uslength < windowsize/2 and dslength >= windowsize/2:
                        bedchrom = e.chrom
                        bedstart = e.start
                        bedend = int(self.coordinate) - int(peakoffset)
                        bedname = e.name
                        bedscore = "TypeIII"
                        bedstrand = e.strand
                        minibed = open("TEMP.bed", "w")                        
                        bedfileout = open("Coding_RTStops_global.bed","a")
                        writeBed(bedchrom,bedstart,str(bedend),bedname,bedscore,bedstrand,minibed)
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                        bedstart = int(self.coordinate) + int(peakoffset)
                        bedend = int(self.coordinate) + int(peakoffset) + windowsize/2
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,minibed)
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                        usdeficit = windowsize/2 - uslength
                        while usdeficit < windowsize/2:
                            i += -1
                            try:
                                nextexon = exonorder[i]
                            except:
                                next
                            if nextexon.consensus == False:
                                next
                            else:
                                exonlength = int(nextexon.end) - int(nextexon.start)
                                uslength += exonlength
                                if uslength >= usdeficit:
                                    bedchrom = e.chrom
                                    bedstart = int(nextexon.end) - usdeficit
                                    bedend = nextexon.end
                                    bedname = e.name
                                    bedscore = "TypeIII"
                                    bedstrand = e.strand
                                    writeBed(bedchrom,str(bedstart),bedend,bedname,bedscore,bedstrand,minibed)
                                    writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                                    minibed.close()
                                    bedfileout.close()
                                    try:
                                        self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP.bed",bedgraph,bgcutoff)
                                    except:
                                        print self
                                    break
                                else:
                                    bedchrom = e.chrom
                                    bedstart = nextexon.start
                                    bedend = nextexon.end
                                    bedname = e.name
                                    bedscore = "TypeIII"
                                    bedstrand = e.strand
                                    writeBed(bedchrom,bedstart,bedend,bedname,bedscore,bedstrand,minibed)
                                    writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                                    usdeficit -= exonlength
                    elif uslength < windowsize/2 and dslength < windowsize/2:
                        bedchrom = e.chrom
                        bedstart = e.start
                        bedend = int(self.coordinate) - int(peakoffset)
                        bedname = e.name
                        bedscore = "TypeIV"
                        bedstrand = e.strand
                        minibed = open("TEMP.bed", "w")
                        bedfileout = open("Coding_RTStops_global.bed","a")
                        writeBed(bedchrom,bedstart,bedend,bedname,bedscore,bedstrand,minibed)
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                        bedstart = int(self.coordinate) + int(peakoffset)
                        bedend = e.end
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,minibed)
                        writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                        usdeficit = windowsize/2 - uslength
                        dsdeficit = windowsize/2 - dslength
                        while usdeficit < windowsize/2:
                            i += -1
                            try:
                                nextexon = exonorder[i]
                            except:
                                next
                            if nextexon.consensus == False:
                                next
                            else:
                                exonlength = int(nextexon.end) - int(nextexon.start)
                                uslength += exonlength
                                if uslength >= windowsize/2:
                                    bedchrom = e.chrom
                                    bedstart = int(nextexon.end) - dsdeficit
                                    bedend = nextexon.end
                                    bedname = e.name
                                    bedscore = "TypeIV"
                                    bedstrand = e.strand
                                    writeBed(bedchrom,str(bedstart),bedend,bedname,bedscore,bedstrand,minibed)
                                    writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                                    break
                                else:
                                    bedchrom = e.chrom
                                    bedstart = nextexon.start
                                    bedend = nextexon.end
                                    bedname = e.name
                                    bedscore = "TypeIV"
                                    bedstrand = e.strand
                                    writeBed(bedchrom,bedstart,bedend,bedname,bedscore,bedstrand,minibed)
                                    writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                                    usdeficit -= exonlength
                        while dsdeficit < windowsize/2:
                            i += 1
                            try:
                                nextexon = exonorder[i]
                            except:
                                next
                            if nextexon.consensus == False:
                                next
                            else:
                                exonlength = int(nextexon.end) - int(nextexon.start)
                                dslength += exonlength
                                if dslength >= windowsize/2:
                                    bedchrom = e.chrom
                                    bedstart = nextexon.start
                                    bedend = int(nextexon.start) + dsdeficit
                                    bedname = e.name
                                    bedscore = "TypeIV"
                                    bedstrand = e.strand
                                    writeBed(bedchrom,bedstart,str(bedend),bedname,bedscore,bedstrand,minibed)
                                    writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                                    minibed.close()
                                    bedfileout.close()
                                    try:
                                        self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP.bed",bedgraph,bgcutoff)
                                        break
                                    except:
                                        print self
                                else:
                                    bedchrom = e.chrom
                                    bedstart = nextexon.start
                                    bedend = nextexon.end
                                    bedname = e.name
                                    bedscore = "TypeIV"
                                    bedstrand = e.strand
                                    writeBed(bedchrom,bedstart,bedend,bedname,bedscore,bedstrand,minibed)
                                    writeBed(bedchrom,str(bedstart),str(bedend),bedname,bedscore,bedstrand,bedfileout)
                                    dsdeficit -= exonlength            

class Transcript:
    def __init__(self,gene,exonlist,num,start,end):
        self.gene = gene
        self.start = int(start)
        self.end = int(end)
        self.transcriptid = gene+"_ts_"+num
        self.transcriptexonlist = exonlist

    def __str__(self):
        return("TRANSCRIPT"+"\t"+self.transcriptid)

class Exon:
    def __init__(self,chrom,start,end,genename,strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = genename
        self.strand = strand
        if self.strand == "+":
            self.threepss = self.start
            self.fivepss = self.end
        elif self.strand == "-":
            self.threepss = self.end
            self.fivepss = self.start
        self.locus = self.chrom+":"+str(self.start)+"-"+str(self.end)
        self.subsumed = None

    def __str__(self):
        return ("EXON"+"\t"+self.name+"\t"+self.locus+"\n")

class Intron:
    def __init__(self,chrom,start,end,intronid,gene,number,introntype,spliceunit,DI,pAsites,strand,juncs):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.intronid = intronid
        self.gene = gene
        self.number = int(number)
        self.introntype = introntype
        self.spliceunit = spliceunit
        self.DI = DI
        self.pAsites = pAsites
        self.strand = strand
        self.juncs = int(juncs)
        self.locus = chrom+":"+str(start)+"-"+str(end)
        if self.strand == "+":
            self.fivepss = self.start
            self.threepss = self.end
        elif self.strand == "-":
            self.fivepss = self.end
            self.threepss = self.start

    def __str__(self):
        return ("INTRON"+"\t"+self.intronid+"\t"+self.gene+"\t"+str(self.number)+"\t"+self.introntype+"\t"+self.spliceunit+"\t"+str(self.juncs)+"\t"+self.locus+"\n")



"""
Functions
"""

def reverseComplement(instring,nucleic):
    rev = instring[::-1]
    revcomp = ''
    for letter in rev:
        if letter == 'G':
            revcomp += 'C'
        if letter == 'C':
            revcomp += 'G'
        if letter == 'T':
            revcomp += 'A'
        if letter == 'A':
            if nucleic == 'DNA':
                revcomp += 'T'
            elif nucleic == 'RNA':
                revcomp += 'U'
    return revcomp
 
def pickler(rootname,output):
    print 'Pickling genes...'
    picklefilename = rootname+'.pkl'
    filehandler = open(picklefilename, 'wb')
    pickle.dump(output,filehandler, -1)
    filehandler.close()

def unpickle(rootname):
    print 'Unpickling genes...'
    pkl_file = rootname+'.pkl'
    picklefile = open(pkl_file,'rb')
    output = pickle.load(picklefile)
    picklefile.close()
    return output

def writeBed(chrom,start,end,name,score,strand,bedfile):
    newentry = "%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom,start,end,name,score,strand)
    bedfile.write(newentry)
    return None

def getSequence(chrom,coordinate,strand,mode,peakoffset,windowsize):
    if strand == "+":
        if mode == "-s":
            #single position is -1 for conversion 1 to zero based, plus -1 for the upstream nucleotide from 3' end
            start = int(coordinate) - 2
            end = int(coordinate) - 1
        if mode == "-6":
            start = int(coordinate) - 4
            end = int(coordinate) + 3
        if mode == "-centered":
            start = int(coordinate) -  (windowsize / 2) - 2
            end = int(coordinate) + (windowsize / 2) - 1
        if mode == "-us":
            start = int(coordinate) - windowsize - 2
            end = int(coordinate) - 1
        if mode == "-ds":
            start = int(coordinate) - 1
            end = int(coordinate) + windowsize - 1
##- strand not fixed yet
    elif strand == "-":
        if mode == "-s":
            start = int(coordinate)
            end = int(coordinate) + 1
        if mode == "-6":
            start = int(coordinate) - 4
            end = int(coordinate) + 3
        if mode == "-centered":
            start = int(coordinate) -  (windowsize / 2) - 1
            end = int(coordinate) + (windowsize / 2)
        if mode == "-us":
            start = int(coordinate) - windowsize - 1
            end = int(coordinate) 
        if mode == "-ds":
            start = int(coordinate)
            end = int(coordinate) + windowsize
    
    siteout = hg38.fetch (chrom,start,end)
    return siteout

def readBedin(bedin):
    for line in bedin:
        modline = line.strip().split('\t')
        chrom = modline[0]
        start = modline[1]
        end = modline[2]
        uid = modline[3]
        counts = modline[4]
        strand = modline[11]
        genename = modline[9]     
        geneclass = modline[10]
        if int(counts) >= threshold:
            for g in metagenelist:
                if g.name == genename:
                    gene = g
                    break
            s = RTStop(chrom, start, end, uid,counts,strand,gene,geneclass)
            if not s.endread:
                s.transcriptMaker()
                print s
                stoplist.append(s)

def readncBedin(bedin):
    for line in bedin:
        modline = line.strip().split('\t')
        chrom = modline[0]
        start = modline[1]
        end = modline[2]
        uid = modline[3]
        counts = modline[4]
        strand = modline[11]
        genename = modline[9]
        geneuid = modline[9]+"_"+modline[7]
        geneclass = modline[10]
        genestart = modline[7]
        geneend = modline[8]
        controlcounts = modline[12]
        if int(counts) >= threshold:
            r = ncRNA(chrom, start, end, uid,counts,strand,genename,geneuid,geneclass,genestart,geneend,controlcounts)
            if not r.endread:
                print r
                r.ncRNAbed()
                stoplist.append(r)   

def coverageBed(minibed,bedin,bgcutoff):
    proc = subprocess.Popen(["intersectBed", "-a", minibed, "-b", bedin, "-wao", "-S", "-sorted"],stdout = subprocess.PIPE)
    uidbglist = []
    for line in proc.stdout:
        line = str(line)
        line = line.strip().split("\t")
        depth = int(line[10])
        breadth = int(line[19])
        #print "depth=", depth, "breadth=", breadth
        for n in range(1,breadth+1):
            uidbglist.append(depth)
    bgmedian = np.median(uidbglist)
    bgpercentile = np.percentile(uidbglist,bgcutoff)
    print uidbglist, len(uidbglist), np.median(uidbglist), np.percentile(uidbglist,90)
    return uidbglist, bgmedian, bgpercentile

##def getRegion(bam,region):
##    bamname = bam.split(".")[0]
##    minibam = "TEMP_minibam_"+bamname+"_.bam" 
##    subprocess.call(["samtools", "view", "-bh","-o",minibam,bam,region])
##    subprocess.call(["samtools", "index",minibam])
##    return minibam




"""
Program

usage: python2.7 pseudoSeq_background.py input.pkl input.bed mode  [-s, -6, -centered, -us, -ds] windowsize peakoffset

usage: python2.7 pseudoSeq_.py input.pkl input.bed control.bedgraph mode  [-c, -nc, -a] windowsize peakoffset threshold backgroundpercentilecutoff

-nc: non-coding RNA set
-c:  coding RNA set
-a:  analysis mode, unpickle previous run.  input.pkl is the RTStops.pkl file from -nc or -c.  Opens a list of objects which are individual putative stop sites.

"""

if __name__ == "__main__":

    infile = sys.argv[1] #This is the input.pkl from Annotator
    bedin = sys.argv[2]  #This is the processed bedgraph from which peak nucleotides are extracted (gene_exon_only.bed)
    rootname = infile.split(".")[0] #For unpickling genes
    outstops = bedin.split(".")[0]  #Name for the output pkl file
    bedfilein = open(bedin, "r")
    bedgraph = sys.argv[3] #This is a non-read filtered processed bedgraph to extract background counts from CMC+.  If a threshold-filtered input bed is used, use the unfiltered version for this (wt_C_120_+-3_gene.bed)
    controlbedgraph = sys.argv[4] #As above, this is the unfiltered processed bedgraph for the CMC- counts (wt_x_120_+-3_gene.bed)
    mode = sys.argv[5]  #Temporarily, use -nc for everything--once splicing sensitive bed derivation is configured, switch to -c for spliced genes
    windowsize = int(sys.argv[6]) #Size of the window for extracting background counts, it will be symetrical around the peak site if possible, otherwise extended as far to the border of an exon in either direction to obtain a windowsize
    peakoffset = int(sys.argv[7]) #Number of nucleotides excluded on either side of the peak site in the background bed
    threshold = int(sys.argv[8]) #Number of counts minimum to check a potential RTstop site (if the input bed was already filtered, this will be redundant)
    bgcutoff = int(sys.argv[9])  #This determines which percentile is displayed for the background distribution
    stoplist = []
    rootname = infile.split('.')[0]
    bedfilename = rootname+"_background.bed"
    if mode == "-nc":
        readncBedin(bedfilein)
    elif mode == "-c":
        metagenelist = unpickle(rootname)
        readBedin(bedfilein)            
        for s in stoplist:
            print s
    elif mode == "-a":
        stoplist = unpickle(rootname)
        for s in stoplist:
            print s
        #INSERT ANALYSIS FUNCTIONS HERE
    bedfilein.close()
    pickler(outstops+"_RTstops",stoplist)            

