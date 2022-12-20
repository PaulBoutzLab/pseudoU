#!/usr/bin/env python
# encoding: utf-8
#PLB 11/12/22


"""

Input format:
pkl from pseudoSeq_background.py

"""
"""
Modules

"""

import sys
import pickle
# import copy
# import pysam
# from bedmodules import *
# from bioTools import *
# import subprocess
import numpy as np
from scipy.stats import poisson
from statsmodels.stats.multitest import fdrcorrection

#hg19 = pysam.Fastafile ("/home/paul/Desktop/SPLICE_GTF/hg19.fa")
#mm10 = pysam.Fastafile ("/home/paul/Desktop/SPLICE_GTF/mm10.fa")
#mm9 = pysam.Fastafile ("/home/paul/Desktop/SPLICE_GTF/mm9.fa")
#hg38 = pysam.Fastafile ("/home/paul/Downloads/alpha_0.2.2_source/Genome/GRCh38.primary_assembly.genome.fa")


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
        
        ##statistics
        
        
        
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
        return ("ncRNA"+"\t"+self.genename+"\t"+self.geneclass+"\t"+self.uid+"\t"+self.locus+"\t"+self.sequence+"\t"+self.cminusone+"\t"+self.c+"\t"+self.cplusone+"\t"+self.counts+"\t"+self.controlcounts+"\t"+str(self.bgmedian)+"\t"+str(self.controlbgmedian)+"\t"+str(self.p_value)+"\t"+str(self.p_value_control)+"\t"+str(self.best_fit_mu_bg_count)+"\t"+str(self.best_fit_mu_bg_count_control)+"\t"+self.BH_FDR_5+"\t"+self.BH_FDR_10+"\t"+self.BH_FDR_1+"\t"+self.BH_FDR_5_control+"\t"+self.BH_FDR_10_control+"\t"+self.BH_FDR_1_control+"\t"+self.FDR_5_psedu+"\t"+self.FDR_10_psedu+"\t"+self.FDR_1_psedu)
    
    
    
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
            self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP.bed",bedgraph,bgcutoff)
            self.controlbgcoveragelist, self.controlbgmedian, self.controlbgpercentile = coverageBed("TEMP.bed",controlbedgraph,bgcutoff)
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
            self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP.bed",bedgraph,bgcutoff)
            self.controlbgcoveragelist, self.controlbgmedian, self.controlbgpercentile = coverageBed("TEMP.bed",controlbedgraph,bgcutoff)
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
            self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP.bed",bedgraph,bgcutoff)
            self.controlbgcoveragelist, self.controlbgmedian, self.controlbgpercentile = coverageBed("TEMP.bed",controlbedgraph,bgcutoff)
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
            self.bgcoveragelist, self.bgmedian, self.bgpercentile = coverageBed("TEMP.bed",bedgraph,bgcutoff)
            self.controlbgcoveragelist, self.controlbgmedian, self.controlbgpercentile = coverageBed("TEMP.bed",controlbedgraph,bgcutoff)
        else:
            print "ncBed Error:", self.uid

    def possion_test(self):
        ###initialize
        self.BH_FDR_5 = "N/A"
        self.BH_FDR_10 = "N/A"
        self.BH_FDR_1 = "N/A"
        self.BH_FDR_5_control = "N/A"
        self.BH_FDR_10_control = "N/A"
        self.BH_FDR_1_control = "N/A"
            
        bg_list = [item for item in self.bgcoveragelist if item > 0]
        bg_list_control = [item for item in self.controlbgcoveragelist if item > 0]
        
        if(len(bg_list)>30) and (len(bg_list_control)>30):
            self.best_fit_mu_bg_count = str(int(mle_possion(remove_outlier(5,bg_list))))
            self.best_fit_mu_bg_count_control = str(int(mle_possion(remove_outlier(5,bg_list_control))))

            self.p_value = str(hypothesis_test(int(self.best_fit_mu_bg_count),int(self.counts)))
            self.p_value_control = str(hypothesis_test(int(self.best_fit_mu_bg_count_control),int(self.controlcounts)))
        else:
            self.best_fit_mu_bg_count = "N/A"
            self.best_fit_mu_bg_count_control = "N/A"
            self.p_value = "N/A"
            self.p_value_control = "N/A"
        return

    def FDR_set(self,FDR,rejected,control):
        
        if control == False:
            if FDR == 0.05:
                if rejected == True:
                    self.BH_FDR_5 = "Significant"
                else:
                    self.BH_FDR_5 = "Not Significant"
            if FDR == 0.1:
                if rejected == True:
                    self.BH_FDR_10 = "Significant"
                else:
                    self.BH_FDR_10 = "Not Significant"
            if FDR == 0.01:
                if rejected == True:
                    self.BH_FDR_1 = "Significant"
                else:
                    self.BH_FDR_1 = "Not Significant"
            
        else:
            if FDR == 0.05:
                if rejected == True:
                    self.BH_FDR_5_control = "Significant"
                else:
                    self.BH_FDR_5_control = "Not Significant"
            if FDR == 0.1:
                if rejected == True:
                    self.BH_FDR_10_control = "Significant"
                else:
                    self.BH_FDR_10_control = "Not Significant"
            if FDR == 0.01:
                if rejected == True:
                    self.BH_FDR_1_control = "Significant"
                else:
                    self.BH_FDR_1_control = "Not Significant"
        return
    
    def Pseudo_U_calling(self):
        
        if (self.BH_FDR_5 == "Significant") and (self.BH_FDR_5_control == "Not Significant"):
            self.FDR_5_psedu = "True"
        elif (self.BH_FDR_5 == "Significant") and (self.BH_FDR_5_control == "Significant"):
            self.FDR_5_psedu = "High in control"
        else:
            self.FDR_5_psedu = "Not significant"

        if (self.BH_FDR_10 == "Significant") and (self.BH_FDR_10_control == "Not Significant"):
            self.FDR_10_psedu = "True"
        elif (self.BH_FDR_10 == "Significant") and (self.BH_FDR_10_control == "Significant"):
            self.FDR_10_psedu = "High in control"
        else:
            self.FDR_10_psedu = "Not significant"

        if (self.BH_FDR_1 == "Significant") and (self.BH_FDR_1_control == "Not Significant"):
            self.FDR_1_psedu = "True"
        elif (self.BH_FDR_1 == "Significant") and (self.BH_FDR_1_control == "Significant"):
            self.FDR_1_psedu = "High in control"
        else:
            self.FDR_1_psedu = "Not significant"

        return
      


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

 
def pickler(rootname,output):
#     print 'Pickling genes...'
    picklefilename = rootname+'.pkl'
    filehandler = open(picklefilename, 'wb')
    pickle.dump(output,filehandler, -1)
    filehandler.close()

def unpickle(rootname):
#     print 'Unpickling genes...'
    pkl_file = rootname+'.pkl'
    picklefile = open(pkl_file,'rb')
    output = pickle.load(picklefile)
    picklefile.close()
    return output

def remove_outlier(perct,data):
    data_sorted = np.sort(data)
    n = len(data)
    outliers = int(round(n*perct/100))
    trimmed_data = data_sorted[: n-outliers]
    return trimmed_data

def sum_log_prob(data, mu):
    fit = lambda val: poisson.logpmf(val, mu)
    return fit(data).sum()

def mle_possion(trimmed_data):
    guess = np.arange(round(trimmed_data.mean())-round(1.5*trimmed_data.std()), round(trimmed_data.mean())+round(1.5*trimmed_data.std())+1)
    likelihoods={}
    for value in guess:
        likelihoods[value]=sum_log_prob(trimmed_data, value)
    best_fit = max(likelihoods, key=lambda x: likelihoods[x])
    return best_fit

def hypothesis_test(best_fit_mu, target):
    return poisson.sf(target, best_fit_mu) 

def FDR_BH(stoplist):
        all_tested = [item for item in stoplist if (item.p_value != "N/A")]
        all_tested_control = [item for item in stoplist if (item.p_value_control != "N/A") ]
        
        all_test_p_value = [float(item.p_value) for item in all_tested]
        all_test_p_value_control = [float(item.p_value_control) for item in all_tested_control]
        
        for FDR in [0.01,0.05,0.1]:
            rejected,qval=fdrcorrection(all_test_p_value, alpha=FDR, method='indep', is_sorted=False)
            rejected_control,qval_control=fdrcorrection(all_test_p_value_control, alpha=FDR, method='indep', is_sorted=False)
            for i in range(len(all_tested)):
                all_tested[i].FDR_set(FDR,rejected[i],control = False)
            for i in range(len(all_tested_control)):
                all_tested_control[i].FDR_set(FDR,rejected_control[i],control = True)
        return

"""
Program

usage: python2.7 pseudoSeq_stats.py input.pkl 

"""

if __name__ == "__main__":

    infile = sys.argv[1] #This is the pkl results file from running the pseudoSeq_background.py script
    rootname = infile.split(".")[0] #For unpickling genes
    stoplist = unpickle(rootname)
#INSERT ANALYSIS FUNCTIONS HERE
    for s in stoplist:
        s.possion_test()
    FDR_BH(stoplist)
    for s in stoplist:
        s.Pseudo_U_calling()
        print s

    #These are the attributes of each instance that you can call for the functions. For the purpose of this analyis,everything is input as an ncRNA instance with the exon taking place of the full ncRNA gene:
##    s.chrom           chromosome
##    s.strand          gene strand, opposite of read strand
##    s.start           start of peak
##    s.end             end of peak (actual position in 1-based coordinates)
##    s.uid             unique identifier for peak
##    s.genename        gene name
##    s.geneuid         unique gene identifier for ncRNA loci
##    s.geneclass       gene class
##    s.genestart       gene start (in case of spliced RNA here, this is the exon start)
##    s.geneend         gene end (in case of spliced RNA here, this is the exon end)
##    s.genelength      length of gene (exon) 
##    s.locus           chr:start-end for genome browser
##    s.genelocus       gene (exon) chr:start-end
##
##  SEQUENCE DATA
##    s.sequence        7 nucleotide sequence surrounding the center nucleotide, which is the putative pseudoU 
##    s.cminusone       1 nucleotide upstream (genomically, not transcriptionally oriented) of the putative pseudoU
##    s.c               This will be a 'T' if it is a pseudoU and not a different type of strong RT stop
##    s.cplusone        1 nucleotide downstream (genomically, not transcriptionally oriented) of the putative pseudoU


##  THESE ARE THE DATA FOR THE STATS:
##    s.c                       This will be a 'T' if it is a pseudoU and not a different type of strong RT stop
##    s.counts                  3' end counts in the CMC+ (C) sample at peak
##    s.controlcounts           3' end counts in the CMC- (x) sample at peak        
##    s.bgcoveragelist          A list of the 3' end counts at surrounding nucleotides in the CMC+ sample, excluding a distance (peakoffset, in these sets 5 nucleotides) from the peak in each direction
##    s.bgmedian                The median value in bgcoveragelist
##    s.bgpercentile            The value at percentile (in this case, 90%) in the bgcoveragelist
##    s.controlbgcoveragelist   A list of the 3' end counts at surrounding nucleotides in the CMC- sample, excluding a distance (peakoffset, in these sets 5 nucleotides) from the peak in each direction
##    s.controlbgmedian         The median value in controlbgcoveragelist
##    s.controlbgpercentile     The value at percentile (in this case, 90%) in the bgcoveragelist

    #If you want to repickle
#     pickler(rootname,stoplist)            

