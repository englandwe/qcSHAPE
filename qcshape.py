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


def stopPositions(gtfdict,fasta_data,rtfile,min_stops,min_rpkm):
    basecounts={}
    basepercs={}
    with open(rtfile) as infile:
        for count,tx in enumerate(infile, start=0):
            #only keep even lines; those are the RT stops (odds are base density)
            if count % 2 == 0 and not tx.startswith('#'):
                txtmp=tx.strip().split('\t')
                txname=txtmp[0]
                if txtmp[2] >= min_rpkm:
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
                    #skipping first (because it has no base)
                    for idx,rts in enumerate(txtmp[4:]):
                        stopct=float(rts)
                        if stopct >= min_stops:
                            #icSHAPE pipeline already shifts the base 1 the left (5')
                            rtbase=fullseq[idx].upper()
                            try:
                                basecounts[rtbase]+=stopct
                            except KeyError:
                                basecounts[rtbase]=stopct
    for key,val in basecounts.iteritems():
        basepercs[key]=val/float(sum(basecounts.values()))
    return basecounts,basepercs

def enrichBases(gtfdict,fasta_data,shapefile,min_score,min_rpkm):
    basecounts={}
    basepercs={}
    baselist=[]
    with open(shapefile) as infile:
        for tx in infile:
            txtmp=tx.strip().split('\t')
            txname=txtmp[0]
            if float(txtmp[2]) >= min_rpkm:
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
                for idx,rts in enumerate(txtmp[4:]):
                    if rts != 'NULL' and float(rts) >= min_score:
                        #file format is diff from rt, shift by 1
                        rtbase=fullseq[idx+1].upper()
                        baselist.append([rtbase,rts])
                        try:
                            basecounts[rtbase]+=1
                        except KeyError:
                            basecounts[rtbase]=1
    for key,val in basecounts.iteritems():
        basepercs[key]=val/float(sum(basecounts.values()))
    baselist.sort(key=lambda y : y[1])
    baselist2=[[idx,x[0],x[1]] for idx,x in enumerate(baselist)]
    return basecounts,basepercs,baselist2


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
                 for idx, item in enumerate(shapesub):
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
                     for idx, item in enumerate(shapesub):
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
                    #skipping first position
                    poi=txtmp[4:]
                    bin_list=[]
                    for pos in poi:
                        if float(pos) > min_stops:
                            bin_list.append(1)
                            #bin_list.append(float(pos))
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

def focalTx(gtfdict,fasta_data,txid):
    #get the sequence
    exonlist=[[x.seqname,x.start,x.end,x.strand,int(x.attdict['exon_number'])] for x in gtfdict[txid] if x.feature == 'exon']
    exonlist.sort(key=lambda y : y[4])
    fullseq=''
    for exon in exonlist:
        if exon[3] == '+':
            range_seq=fasta_data[exon[0]].seq[exon[1]-1:exon[2]]
        elif exon[3] == '-':
            range_seq=fasta_data[exon[0]].seq[exon[1]-1:exon[2]].reverse_complement()
        fullseq=fullseq + range_seq
    #get the score and position
    with open(shapefile) as infile:
        for line in infile:
            txtmp=line.strip().split('\t')
            txname=txtmp[0]
            if txname == txid:
                scores=[[x,y] for x,y in enumerate(txtmp[3:],start=1)]
    #join them
    outlist=[]
    for i in range(0,len(scores)):
        outlist.append([scores[i][0],fullseq[i],scores[i][1]])

    return outlist
            

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

#sam_list=glob.glob(workdir+'/*.sam')
id_list=['DMSO1','NAI1','DMSO2','NAI2','DMSO3','NAI3']
shapefile='icshape.out'
prm=[]
btcount=[]
basestops=[]
basepercs=[]
stops,rpkm,score=sys.argv[4].split(',')
min_stops=int(stops)
min_rpkm=float(rpkm)
min_score=float(score)

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
    stopct,stopperc=stopPositions(gtfdict,fasta_dict,rt_file,min_stops,min_rpkm)
    basestops.append(stopct)
    basepercs.append(stopperc)
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
rtout2=open('stopbases_perc.txt','w')
rtout.write("Sample\tBase\tCount\n")
rtout2.write("Sample\tBase\tPerc\n")
for i in range(0,len(id_list)):
    rtout.write(flattenDictAddCol(basestops[i],id_list[i])+'\n')
    rtout2.write(flattenDictAddCol(basepercs[i],id_list[i])+'\n')
rtout.close()
rtout2.close()


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
rpkmout=open('rpkmcorr.txt','w')
rpkmout.write(flattenList(rpkm_corr))
rpkmout.close()

#now RT
sys.stderr.write('RT corr @ %s' % (str(datetime.now())) + '\n')
rt_corr=corrRT(id_list,min_stops)
rtcorrout=open('rtcorr.txt','w')
rtcorrout.write(flattenList(rt_corr))
rtcorrout.close()

#shape scores by base
#make a dict (just like raw counts) and also make a list of base, score
sys.stderr.write('Enrich bases @ %s' % (str(datetime.now())) + '\n')
enrichcounts,enrichperc,enrichlist=enrichBases(gtfdict,fasta_dict,shapefile,min_score,min_rpkm)
enout=open('enrichbases.txt','w')
enout.write("Base\tCount\n")
enout.write(flattenDict(enrichcounts))
enout.close()

enout2=open('enrichbases_perc.txt','w')
enout2.write("Base\tPerc\n")
enout2.write(flattenDict(enrichperc))
enout2.close()

enout3=open('enrichbases_list.txt','w')
enout3.write("Base\tScore\n")
enout3.write(flattenList(enrichlist))
enout3.close()
sys.stderr.write('Enrich bases done @ %s' % (str(datetime.now())) + '\n')

# let's yoink 18S
sys.stderr.write('Getting 18S @ %s' % (str(datetime.now())) + '\n')
out18s=focalTx(gtfdict,fasta_dict,'18S')
focout=open('18S_scores.txt','w')
focout.write("Position\tBase\tScore\n")
focout.write(flattenList(out18s))
focout.close()
sys.stderr.write('18S done @ %s' % (str(datetime.now())) + '\n')

sys.stderr.write('Getting 28S @ %s' % (str(datetime.now())) + '\n')
out28s=focalTx(gtfdict,fasta_dict,'28S')
focout2=open('28S_scores.txt','w')
focout2.write("Position\tBase\tScore\n")
focout2.write(flattenList(out28s))
focout2.close()
sys.stderr.write('28S done @ %s' % (str(datetime.now())) + '\n')



sys.stderr.write('Done @ %s' % (str(datetime.now())) + '\n')




