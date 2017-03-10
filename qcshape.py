#!/usr/bin/env python

#Input: a directory with the output of the icSHAPE pipeline

import sys
import subprocess
import os
import glob
import re
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
#####

class GtfRec(object):
    def __init__(self,reclist):
         self.seqname=reclist[0]
         self.source=reclist[1]
         self.feature=reclist[2]
         self.start=int(reclist[3])
         self.end=int(reclist[4])
         self.score=reclist[5]
         self.strand=reclist[6]
         self.frame=reclist[7]
         self.attdict={}
         for self.item in reclist[8].strip(';').split('; '):
             self.splitline=self.item.replace('\"','').split(' ')
             self.attdict[self.splitline[0]]=self.splitline[1]

def file_len(fname):
    i=0
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def percReadsMapped(trimmed_fastq,mapct):
    total_reads=file_len(trimmed_fastq)/4
    perc_mapped=mapct/float(total_reads)
    return perc_mapped

#def countBiotypes(gtfdict,samlist):
#    btdict={}
#    for entry in samlist:
#        maptx=entry[2]
#        for rec in gtfdict[maptx]:
#            if rec.feature == 'transcript':
#                biotype=rec.attdict['transcript_biotype']
#                if biotype in btdict:
#                    btdict[biotype] += 1
#                else:
#                    btdict[biotype] = 1
#            break
#    return btdict


def countBiotypes(gtfdict,sam_file):
    btdict={}
    mapped_reads=set()
    with open(sam_file) as infile:
        for line in infile:
            #no headers, no unmapped reads
            if not line.startswith('@'):
                tmpline=line.strip().split('\t')
                if tmpline[1] != '4':
                    maptx=tmpline[2]
                    mapped_reads.add(tmpline[0])
                    for rec in gtfdict[maptx]:
                        if rec.feature == 'transcript':
                            biotype=rec.attdict['transcript_biotype']
                            if biotype in btdict:
                                btdict[biotype]+=1
                            else:
                                btdict[biotype]=1
    return btdict,len(mapped_reads)

def stopPositions(gtfdict,fasta_data,rtfile,min_stops):
    basecounts={}
    with open(rtfile) as infile:
        for count,tx in enumerate(rtfile, start=1):
            #only keep even lines; those are the RT stops (odds are base density)
            if count % 2 == 0:
                txtmp=tx.strip().split('\t')
                txname=txtmp[0]
                #acquire sequence
                exonlist=[[x.seqname,x.start,x.end,x.strand,int(x.attdict['exon_number'])] for x in gtfdict[txname] if x.feature == 'exon']
                exonlist.sort(key=lambda y : y[4])
                fullseq=''
                for exon in exonlist:
                    if exon[3] == '+':
                        range_seq=fasta_data[exon[0]].seq[exon[1]-1:exon[2]]
                    elif exon[3] == '-':
                        range_seq=fasta_data[exon[0]].seq[exon[1]-1:exon[2]].reverse_complement()
                    fullseq=fullseq + range_seq
                #get the base at each stop
                for idx,rts in enumerate(txtmp[3:]):
                    if rts >= min_stops:
                        rtbase=fullseq[idx]
                        #assuming that the base exactly at the stop is what you want.  may have to subtract 1 to get the base 1 5' of this; not currently sure
                        if rtbase in basecounts:
                            basecounts[rtbase]+=1
                        else:
                            basecounts[rtbase]=1
    return basecounts


def genToTx(txname,gtfdict):
    exonlist=[[x.seqname,x.start,x.end,x.strand,int(x.attdict['exon_number'])] for x in gtfdict[txname] if x.feature == 'exon']
    exonlist.sort(key=lambda y : y[4])
    rangelist=[]
    #make exon genome/tx ranges
    tx_start=1
    for exon in exonlist:
        tx_stop=tx_start+(exon[2]-exon[1])
        tmprange=[exon[1],exon[2],tx_start,tx_stop]
        rangelist.append(tmprange)
        tx_start=tx_stop+1
    return rangelist


def shapeByRegion(gtfdict,shapefile):
    #we're going to be agnostic to that the feature type is (unless it's a start codon or a transcript on an exon)
    #just grab shape reactivity for all positions in its range
    #for start codons, we're additionally going to grab 25bp up & down the transcript
    #transcript and exon features will be ignored
    shapeout=[]
    with open(shapefile) as infile:
         for line in infile:
             tmpshape=line.strip().split('\t')
             txname=tmpshape[0]
             #translate the genome positions into tx positions using the magic of exons
             ranges=genToTx(txname,gtfdict)
             feature_list=[x for x in gtfdict[txname] if x.feature not in ['exon','transcript']]
             for feature in feature_list:
                 shapesub=[]
                 #for position in feature
                 for genpos in range(feature.start,feature.end+1):
                     #get corresponding tx pos
                     tmprng=[x for x in ranges if genpos in range(x[0],x[1]+1)][0]
                     txpos=(genpos-tmprng[0])+tmprng[2]
                     #then get the shape at 1-that position in tmpshape[3:]
                     shapeval=tmpshape[3:][txpos-1]
                     shapesub.append(shapeval)
                 shapeout.append([txname,feature.feature,feature.start,feature.end,feature.strand,shapesub])
                 if feature.feature == 'start_codon':
                     shapesub=[]
                     for biggenpos in range(feature.start-25,feature.end+26):
                         #get corresponding tx pos
                         bigtmprng=[x for x in ranges if biggenpos in range(x[0],x[1]+1)]
                         if len(bigtmprng) == 0:
                             shapeval='NA'
                         else:
                             txpos=(biggenpos-bigtmprng[0][0])+bigtmprng[0][2]
                             #then get the shape at 1-that position in tmpshape[3:]
                             shapeval=tmpshape[3:][txpos-1]
                             shapesub.append(shapeval)
                     shapeout.append([txname,feature.feature+'_25',feature.start,feature.end,feature.strand,shapesub])

    return shapeout

def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final

#####

#make this a dict
gtfdict={}
#with open(sys.argv[2] as infile:
with open('../mouse_txome/Mus_musculus.GRCm38.87.gtf') as infile:
    for line in infile:
        if not line.startswith('#'):
            rectmp=GtfRec(line.strip().split('\t'))
            if 'transcript_id' in rectmp.attdict:
                transcript=rectmp.attdict['transcript_id']
                try:
                    gtfdict[transcript].append(rectmp)
                except KeyError:
                    gtfdict[transcript]=[rectmp]



#fasta_dict=SeqIO.index(sys.argv[3], "fasta", alphabet=IUPAC.unambiguous_dna)
fasta_dict=SeqIO.index('../mouse_txome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', "fasta", alphabet=IUPAC.unambiguous_dna)

#this will be sys.argv[1]
#sam_list=glob.glob(workdir+'/*.sam')
id_list=['SRR1534952']
shapefile='icshape.out'
prm=[]
btcount=[]
basestops=[]
for id in id_list:
    sam_file=id+'.sam'
    trimmed_fastq=id+'.trimmed.fastq'
    rt_file=id+'.rt'
    #02 - count reads mapped to each transcript_biotype
    #returns dict of {biotype:count}, plus a count of mapped reads
    bts,mapct=countBiotypes(gtfdict,sam_file)
    btcount.append(bts)
    #01 - percentage reads mapped, period
    #returns float
    prm.append(percReadsMapped(trimmed_fastq,mapct))
    #03 & 05 - count raw stops at each base (also gets you total stops)
    #returns dict of {base:count}
    min_stops=1
    basestops.append(stopPositions(gtfdict,fasta_dict,rt_file,min_stops))
    #temp output
    print prm
    print btcount
    print basestops

#now for the non-file-specific stuff
#06 & 07 shape across RNA regions (5'UTR,3'UTR, orf, start & stop codons), plus a range around start codons
#this will be shape reactivity - i.e. final output, already background subtracted
#returns a list of lists: [tx,genomic_start,genomic_stop,strand,[list of shape values]]
shape_by_reg=shapeByRegion(shapefile,gtfdict)
f=open('shapebyregtest','w')
f.write(flattenList(shape_by_reg))
f.close()

#3.5. Correlation of expression level b/w replicates
#4. Correlation b/w stops in replicates (sample should have more in common than control)
