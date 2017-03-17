#!/usr/bin/env python

#Input: a directory with the output of the icSHAPE pipeline

import sys
import subprocess
import os
import glob
import re
#biopython
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
#scipy
from scipy.stats.stats import pearsonr

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


#def percReadsMapped(trimmed_fastq,mapct):
#    total_reads=file_len(trimmed_fastq)/4
#    perc_mapped=mapct/float(total_reads)
#    return perc_mapped



def countBiotypes(gtfdict,sam_file):
    btdict={}
    mapped_reads=set()
    unmapped_ct=0
    with open(sam_file) as infile:
        for line in infile:
            #no headers, no unmapped reads
            if not line.startswith('@'):
                tmpline=line.strip().split('\t')
                if tmpline[1] == '4':
                    unmapped_ct+=1
                else:
                    maptx=tmpline[2]
                    mapped_reads.add(tmpline[0])
                    for rec in gtfdict[maptx]:
                        if rec.feature == 'transcript':
                            biotype=rec.attdict['transcript_biotype']
                            if biotype in btdict:
                                btdict[biotype]+=1
                            else:
                                btdict[biotype]=1
    return btdict,len(mapped_reads),unmapped_ct



def stopPositions(gtfdict,fasta_data,rtfile,min_stops):
    basecounts={}
    with open(rtfile) as infile:
        for count,tx in enumerate(infile, start=0):
            #only keep even lines; those are the RT stops (odds are base density)
            if count % 2 == 0 and not tx.startswith('#'):
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
                #get the base at each stop.
                #the length of the RT-stop list is one longer than the actual tx
                #there are never stops at the last position (due to 1bp shift)
                #the first would hit the base before the tx - so we don't have the base
                # for now, skipping first (because it has no base) and last (because it will never have any stops)
                for idx,rts in enumerate(txtmp[4:-1]):
                    if rts >= min_stops:
                        #icSHAPE pipeline already shifts the base 1 the left (5')
                        rtbase=fullseq[idx].upper()
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
                     shapeout.append([txname,feature.feature+'_25',feature.start+25,feature.end-25,feature.strand,shapesub])

    return shapeout



def corrRPKM(id_list):
    rpkm_corr_out=[]
    rpkm_list=[]
    for entry in id_list:
        tmpdict={}
        rpkm_file=entry+'.rpkm'
        with open(rpkm_file) as infile:
            for line in infile:
                if not line.startswith('#'):
                    tmpline=line.strip().split()
                    tmpdict[tmpline[0]]=float(tmpline[4])
        rpkm_list.append([entry,tmpdict])
    for i in range(0,len(rpkm_list)):
        for j in range(i,len(rpkm_list)):
            vals1=[]
            vals2=[]
            keyset=set([key for tx in [rpkm_list[i]]+[rpkm_list[j]] for key in tx[1].keys()])
            for key in keyset:
                try:
                    val1=rpkm_list[i][1][key]
                    val2=rpkm_list[j][1][key]
                except KeyError:
                    continue
                vals1.append(val1)
                vals2.append(val2)
            rpkm_r=pearsonr(vals1,vals2)
            rpkm_corr_out.append([rpkm_list[i][0],rpkm_list[j][0],rpkm_r[0],rpkm_r[1]])
    return rpkm_corr_out

def corrRT(id_list,min_stops):
    rt_corr_out=[]
    rt_list=[]
    for entry in id_list:
        tmpdict={}
        rt_file=entry+'.rt'
        with open(rt_file) as infile:
            for count,tx in enumerate(infile, start=0):
                #only keep even lines; those are the RT stops (odds are base density)
                if count % 2 == 0 and not tx.startswith('#'):
                    txtmp=tx.strip().split('\t')
                    #currently calculating this as binary stop presence/absence
                    #maybe later normalize to base density and compare frequency?
                    #as with finding bases, ignoring first and last position of each transcript
                    poi=txtmp[4:-1]
                    bin_list=[]
                    for pos in poi:
                        if float(pos) > min_stops:
                            bin_list.append(1)
                        else:
                            bin_list.append(0)
                    tmpdict[txtmp[0]]=bin_list
        rt_list.append([entry,tmpdict])
    for i in range(0,len(rt_list)):
        for j in range(i,len(rt_list)):
            vals1=[]
            vals2=[]
            keyset=set([key for tx in [rt_list[i]]+[rt_list[j]] for key in tx[1].keys()])
            for key in keyset:
                try:
                    bins1=rt_list[i][1][key]
                    bins2=rt_list[j][1][key]
                except KeyError:
                    continue
                vals1+=bins1
                vals2+=bins2
            rt_r=pearsonr(vals1,vals2)
            rt_corr_out.append([rt_list[i][0],rt_list[j][0],rt_r[0],rt_r[1]])
    return rt_corr_out

def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final

def flattenDict(dictin):
    listin=[[key,val] for key,val in dictin.iteritems()]
    final=flattenList(listin)
    return final

#####

#make this a dict
gtfdict={}
#with open(sys.argv[2] as infile:
#with open('../mouse_txome/Mus_musculus.GRCm38.87.gtf') as infile:
with open('../Mus_musculus.GRCm38.87.gtf') as infile:
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
#fasta_dict=SeqIO.index('../mouse_txome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', "fasta", alphabet=IUPAC.unambiguous_dna)
fasta_dict=SeqIO.index('../Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', "fasta", alphabet=IUPAC.unambiguous_dna)

#this will be sys.argv[1]
#sam_list=glob.glob(workdir+'/*.sam')
id_list=['SRR1534952','SRR1534953','SRR1534954','SRR1534955']
shapefile='icshape.out'
prm=[]
btcount=[]
basestops=[]
map_unmap=[]
for id in id_list:
    sam_file=id+'.sam'
    rt_file=id+'.rt'
    #02 - count reads mapped to each transcript_biotype
    #returns dict of {biotype:count}, plus a count of mapped reads
    bts,mapct,unmapct=countBiotypes(gtfdict,sam_file)
    btcount.append(bts)
    map_unmap.append((mapct,unmapct))
    #01 - percentage reads mapped, period
    #returns float
    #prm.append(percReadsMapped(trimmed_fastq,mapct))
    prm.append(mapct/float(mapct+unmapct))
    #03 & 05 - count raw stops at each base (also gets you total stops)
    #returns dict of {base:count}
    min_stops=1
    basestops.append(stopPositions(gtfdict,fasta_dict,rt_file,min_stops))


e=open('basestoptest','w')
for i in range(0,len(id_list)):
    e.write(id_list[i]+'\n')
    e.write("Mapped, unmapped, percmapped: %s, %s, %s\n" % (map_unmap[i][0],map_unmap[i][1],prm[i]))
    e.write(flattenDict(btcount[0])+'\n')
    e.write(flattenDict(basestops[0]))
e.close()


#now for the non-file-specific stuff
#06 & 07 shape across RNA regions (5'UTR,3'UTR, orf, start & stop codons), plus a range around start codons
#this will be shape reactivity - i.e. final output, already background subtracted
#returns a list of lists: [tx,genomic_start,genomic_stop,strand,[list of shape values]]
shape_by_reg=shapeByRegion(gtfdict,shapefile)
f=open('shapebyregtest','w')
f.write(flattenList(shape_by_reg))
f.close()

#3.5. Correlation of expression level b/w replicates
#4. Correlation b/w stops in replicates (sample should have more in common than control)

#first, expression (RPKM files)

rpkm_corr=corrRPKM(id_list)
g=open('rpkmcorrtest','w')
g.write(flattenList(rpkm_corr))
g.close()

#now RT
rt_corr=corrRT(id_list,min_stops)
h=open('rtcorrtest','w')
h.write(flattenList(rt_corr))
h.close()
