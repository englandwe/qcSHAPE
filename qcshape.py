#!/usr/bin/env python

#Input: a directory with the output of the icSHAPE pipeline

import sys
import os
import re
from multiprocessing import Pool
from collections import Counter
#biopython
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
#scipy
from scipy.stats.stats import pearsonr
#just for timing
from datetime import datetime
#####
#Classes
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

#Multiprocessing
def chunkify(fname,size):
    fileEnd = os.path.getsize(fname)
    with open(fname,'r') as f:
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size,1)
            f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break

def child_initialize(_inData):
     global ggtfdict
     ggtfdict = _inData

def child_initialize_2(_inData,_inData2):
     global ggtfdict, gfasta_data
     ggtfdict = _inData
     gfasta_data = _inData2

#Make pretty output
def flattenList(listin):
    list2=[]
    for item in listin:
        list2.append('\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final

def flattenListAddCol(listin,col1):
    list2=[]
    for item in listin:
        list2.append(str(col1) + '\t' + '\t'.join([str(x) for x in item]))
    final='\n'.join(list2)
    return final

def flattenDict(dictin):
    listin=[[key,val] for key,val in dictin.iteritems()]
    final=flattenList(listin)
    return final

def flattenDictAddCol(dictin,col1):
    listin=[[key,val] for key,val in dictin.iteritems()]
    final=flattenListAddCol(listin,col1)
    return final

#Science goes here
def countBiotypes(sam_file,chunkStart,chunkSize):
    btdict={}
    mapped_ct=0
    unmapped_ct=0
    with open(sam_file) as infile:
        infile.seek(chunkStart)
        lines = infile.read(chunkSize).splitlines()
        for line in lines:
            #no headers, no unmapped reads
            if not line.startswith('@'):
                tmpline=line.strip().split('\t')
                #if tmpline[1] == '4':
                #bitwise is bestwise
                if int(tmpline[1]) & 4:
                    unmapped_ct+=1
                    #with unmapped_ct.get_lock():
                        #unmapped_ct.value+=1
                #else:
                elif not int(tmpline[1]) & 2304:
                    maptx=tmpline[2]
                    mapped_ct+=1
                    #with mapped_ct.get_lock():
                        #mapped_ct.value+=1
                    #mapped_reads.add(tmpline[0])
                    for rec in ggtfdict[maptx]:
                        if rec.feature == 'transcript':
                            biotype=rec.attdict['transcript_biotype']
                            try:
                                btdict[biotype]+=1
                            except KeyError:
                                btdict[biotype]=1
    sys.stderr.write('another chunk bites the dust @ %s' % (str(datetime.now())) + '\n')
    return btdict,mapped_ct,unmapped_ct


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
                    if float(rts) >= min_stops:
                        #icSHAPE pipeline already shifts the base 1 the left (5')
                        rtbase=fullseq[idx].upper()
                        try:
                            basecounts[rtbase]+=1
                        except KeyError:
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

def shapeByRegion(shapefile,chunkStart,chunkSize):
    #we're going to be agnostic to that the feature type is (unless it's a start codon or a transcript on an exon)
    #just grab shape reactivity for all positions in its range
    #for start/stop codons, we're additionally going to grab 25bp up & down the transcript
    #transcript and exon features will be ignored
    shapeout=[]
    with open(shapefile) as infile:
        infile.seek(chunkStart)
        lines = infile.read(chunkSize).splitlines()
        for line in lines:
             tmpshape=line.strip().split('\t')
             txname=tmpshape[0]
             #translate the genome positions into tx positions using the magic of exons
             ranges=genToTx(txname,ggtfdict)
             feature_list=[x for x in ggtfdict[txname] if x.feature not in ['exon','transcript']]
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
                 for idx, item in shapesub:
                     shapeout.append([txname,feature.feature,feature.start,feature.end,feature.strand,idx+1,item])
                 if feature.feature in ['start_codon','stop_codon']:
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
                     for idx, item in shapesub:
                         shapeout.append([txname,feature.feature+'_25',feature.start+25,feature.end-25,feature.strand,idx+1,item])
    sys.stderr.write('another chunk bites the dust @ %s' % (str(datetime.now())) + '\n')
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


#####
#main
sys.stderr.write('importing gtf @ %s' % (str(datetime.now())) + '\n')

gtfdict={}
with open(sys.argv[2]) as infile:
#with open('../../mouse_txome/Mus_musculus.GRCm38.87.gtf') as infile:
#with open('../Mus_musculus.GRCm38.87.gtf') as infile:
    for line in infile:
        if not line.startswith('#'):
            rectmp=GtfRec(line.strip().split('\t'))
            if 'transcript_id' in rectmp.attdict:
                transcript=rectmp.attdict['transcript_id']
                try:
                    gtfdict[transcript].append(rectmp)
                except KeyError:
                    gtfdict[transcript]=[rectmp]


sys.stderr.write('importing fasta @ %s' % (str(datetime.now())) + '\n')

#fasta_dict=SeqIO.index(sys.argv[3], "fasta", alphabet=IUPAC.unambiguous_dna)
#fasta_dict=SeqIO.index('../../mouse_txome/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', "fasta", alphabet=IUPAC.unambiguous_dna)
#fasta_dict=SeqIO.index('../Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', "fasta", alphabet=IUPAC.unambiguous_dna)

#higher memory, hopefully faster
fasta_dict=SeqIO.to_dict(SeqIO.parse(sys.argv[3], "fasta", alphabet=IUPAC.unambiguous_dna))

cores=int(sys.argv[1])

#this will be sys.argv[1]
#sam_list=glob.glob(workdir+'/*.sam')
id_list=['SRR1534952','SRR1534953','SRR1534954','SRR1534955']
shapefile='icshape.out'
prm=[]
btcount=[]
basestops=[]
min_stops=1
for id in id_list:
    sam_file=id+'.sam'
    rt_file=id+'.rt'
    #BIOTYPES - this works
    sys.stderr.write('counting biotypes: %s @ %s' % (id,str(datetime.now())) + '\n')
    btpool = Pool(cores, initializer = child_initialize, initargs = (gtfdict,))
    btjobs = []
    btjobs=[btpool.apply_async(countBiotypes, (sam_file,chunkStart,chunkSize)) for chunkStart,chunkSize in chunkify(sam_file,1024*1024*1024)]
    btpool.close()
    btpool.join()
    btoutput=[btjob.get() for btjob in btjobs]
    biotype_dict=dict(sum([Counter(x[0]) for x in btoutput], Counter()))
    btcount.append(biotype_dict)
    mapct=sum([x[1] for x in btoutput])
    unmapct=sum([x[2] for x in btoutput])
    mapperc=mapct/float(mapct+unmapct)
    prm.append([mapct,unmapct,mapperc])
    sys.stderr.write('done counting biotypes: %s @ %s' % (id,str(datetime.now())) + '\n')
    #STOPS
    sys.stderr.write('counting RT stops: %s @ %s' % (id,str(datetime.now())) + '\n')
    basestops.append(stopPositions(gtfdict,fasta_dict,rt_file,min_stops))
    sys.stderr.write('done counting RT stops: %s @ %s' % (id,str(datetime.now())) + '\n')

#Output reads mapped, biotypes, base stops
mapout=open('reads_mapped.txt','w')
mapout.write("Sample\tMapped\tUnmapped\tPercmapped\n")
for i in range(0,len(id_list)):
    mapout.write("%s\t%s\t%s\t%s\n" % (id_list[i],prm[i][0],prm[i][1],prm[i][2]))
mapout.close()

btout=open('biotypes.txt','w')
btout.write("Sample\tBiotype\tReadsMapped\n")
for i in range(0,len(id_list)):
    btout.write(flattenDictAddCol(btcount[i],id_list[i])+'\n')
btout.close()

rtout=open('stopbases.txt','w')
rtout.write("Sample\tBase\tCount\n")
for i in range(0,len(id_list)):
    rtout.write(flattenDictAddCol(basestops[i],id_list[i])+'\n')
rtout.close()


#now for the non-file-specific stuff
#06 & 07 shape across RNA regions (5'UTR,3'UTR, orf, start & stop codons), plus a range around start codons
#this will be shape reactivity - i.e. final output, already background subtracted
#returns a list of lists: [tx,genomic_start,genomic_stop,strand,[list of shape values]]

sys.stderr.write('shape by region @ %s' % (str(datetime.now())) + '\n')
shapepool = Pool(cores, initializer = child_initialize, initargs = (gtfdict,))
shapejobs = []
shapejobs=[shapepool.apply_async(shapeByRegion, (shapefile,chunkStart,chunkSize)) for chunkStart,chunkSize in chunkify(shapefile,10*1024*1024)]
shapepool.close()
shapepool.join()
shapeoutput=[shapejob.get() for shapejob in shapejobs]
shape_by_reg=[item for sublist in shapeoutput for item in sublist]
regout=open('shaperegions.txt','w')
regout.write(flattenList(shape_by_reg))
regout.close()
sys.stderr.write('Done with shape by region @ %s' % (str(datetime.now())) + '\n')

#3.5. Correlation of expression level b/w replicates
#4. Correlation b/w stops in replicates (sample should have more in common than control)

#first, expression (RPKM files)
sys.stderr.write('RPKM corr @ %s' % (str(datetime.now())) + '\n')
rpkm_corr=corrRPKM(id_list)
rpkmout=open('rpkmcorrtest','w')
rpkmout.write(flattenList(rpkm_corr))
rpkmout.close()

#now RT
sys.stderr.write('RT corr @ %s' % (str(datetime.now())) + '\n')
rt_corr=corrRT(id_list,min_stops)
rtcorrout=open('rtcorrtest','w')
rtcorrout.write(flattenList(rt_corr))
rtcorrout.close()
sys.stderr.write('Done @ %s' % (str(datetime.now())) + '\n')

