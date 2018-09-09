#!/usr/bin/env python

import sys,glob
#import argparse
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict
from textwrap import fill


def strata(keptFrac,stratFile):
    markDiv=[]
    for line in open(stratFile):
        #print len(line.rstrip('\n').split('\t'))
        if line.startswith('Mark'): continue
        elts=line.rstrip('\n').split('\t')
        gid,fct=elts[0],elts[8]
        if not fct=='ND':
            markDiv.append((gid,float(fct)))
    srtMarkDiv=sorted(markDiv, key=lambda x:x[1],reverse=True)
    kept=int(keptFrac*len(srtMarkDiv))
    #print [elt[0] for elt in srtMarkDiv[0:kept]][0:10]
    return set([elt[0] for elt in srtMarkDiv[0:kept]])


#Usage: concatenate-ext.py <taxlist> <fasta files...>
if not len(sys.argv) >= 4:
    sys.exit('Usage: concatenate-ext.py <taxlist> <suffix> <fasta files...>')


taxfile=sys.argv[1]
suffix=sys.argv[2]
markList=sys.argv[3:]

taxset=[]
for line in open(sys.argv[1]):
    if not line.strip(): continue
    if line.startswith('#'): continue
    cols=line.rstrip().split()
    taxset.append(cols[0])
print len(taxset),'taxa...'

print taxset

#Some parameters
cutoff=0.5
SuperMatrix=defaultdict(str)
partitions=[]
currpos=1
restype='PROT'
#markList=sorted(sys.argv[2:],key=lambda x: int(x.split('.')[0]))
covTax=defaultdict(int)
covClu=defaultdict(int)
nbmar=0
binary=defaultdict(list)
markNames=[]
print "Reading {0} alignments...".format(len(markList))
for i,clusf in enumerate(markList):
    #print clusf
    pid=clusf.split('.')[0]
    alignment = AlignIO.read(open(clusf),"fasta")
    aliLen=alignment.get_alignment_length()
    #covClu[len(alignment)]+=1
    cseq=''
    taxseq={}
    for record in alignment:
        #if '-' in record.id:
        #    taxseq[record.id.split('-')[0]]=record.seq
        # elif '.' in record.id:
        #     taxseq[record.id.split('.')[0]]=record.seq
        #else:
        taxseq[record.id]=record.seq
        cseq+=str(record.seq)
    chset=set(cseq)
    #print len(chset),chset
    if len(chset)==20:
        restype='PROT'
    elif len(chset)==4:
        restype='DNA'

    #trtax=[tax for tax in taxseq.keys() if len(tax)>3]
    trtax=[tax for tax in taxseq.keys() if tax in taxset]
    #if not len(taxseq)>=cutoff*len(taxset):
    #print i,clusf,len(trtax)
    if not len(trtax)>=cutoff*len(taxset):
         continue

    covClu[len(trtax)]+=1

    #taxseq=dict((record.id,record.seq) for record in alignment)
    #print taxseq.keys()
    miss=[]

    markNames.append(clusf.split('.')[0])

    nbmar+=1
    for tax in taxset:
        if tax in taxseq:
            binary[tax].append(1)
            SuperMatrix[tax]+=str(taxseq[tax])
            covTax[tax]+=1
        else:
            binary[tax].append(0)
            SuperMatrix[tax]+="X" * aliLen
            miss.append(tax)
    partitions.append((pid,currpos,currpos+aliLen-1))
    #print i,'Marker:',pid,'Start:',currpos,'End:',currpos+aliLen-1,'Size:',aliLen,'Missing:',len(miss),','.join(miss)
    currpos+=aliLen
print 'Detected format',restype
print "{0} alignments retained!".format(nbmar)
print "{0} total positions".format(currpos)

print "#Taxa\t#Markers"
for cov in covClu:
    print cov,covClu[cov]

print 'Done!'
totlen,totund=0,0
concat=open('Concat-{0}.fa'.format(suffix),'w')
mistats=[]
for tax in SuperMatrix:
    taxseq=SuperMatrix[tax]
    if len(set(taxseq))<4:
        print "Problem with {0}, insufficient alphabet, all missing data!".format(tax)
    totlen+=len(taxseq)
    nbxxx=taxseq.count('X')
    nbgap=taxseq.count('-')
    perdat=(nbxxx+nbgap)
    totund+=perdat
    pexxx=round(100*nbxxx/float(len(taxseq)),2)
    pegap=round(100*nbgap/float(len(taxseq)),2)
    petot=round(100*(nbgap+nbxxx)/float(len(taxseq)),2)
    mistats.append((tax,covTax[tax],pexxx,pegap,petot))
    concat.write(">%s\n%s\n"%(tax,fill(taxseq,width=80,break_on_hyphens=False)))

mistats.sort(key=lambda x: x[4],reverse=True)
outstats=open('Stats-{0}.txt'.format(suffix),'w')
for tax in mistats:
    outstats.write('\t'.join(map(str,tax))+'\n')

print "Alignment length: {0}".format(totlen)
print "%2f total occupancy"%((100-100*totund/float(totlen)))

# print '\t'.join(map(str,[tax,covTax[tax],pexxx,pegap,petot]))

# #print SuperMatrix
# SuperList=[SeqRecord(seq=SuperMatrix[tax],id=tax) for tax in SuperMatrix]
# print 'Superlist Done!'
#
# SuperMatrixAlign=MultipleSeqAlignment(SuperList)
# print 'supermatrix,Done!'
# #print SuperMatrixAlign
#
# concat=open('Concatenate.fa','w')
# AlignIO.write(SuperMatrixAlign, concat, "fasta")
# print 'concat,Done!'
model='LG'
#if restype=='DNA':
 #   model='DNA'

#compl=open('Completness-{0}.txt'.format(suffix),'w')
#compl.write("Species\t{0}\n".format('\t'.join(map(str,markNames))))
#for tax in binary:
#    compl.write("{0}\t{1}\n".format(tax,'\t'.join(map(str,binary[tax]))))




outpart=open("Part-{1}-{0}.txt".format(suffix, model),'w')
for part in partitions:
    idx,start,end=part
    if restype=='DNA':
        #outpart.write("{0}, gene_{1}_p1 = {2}-{4}\\3,{3}-{4}\\3\n".format(model,idx,start,start+1,end))
        #outpart.write("{0}, gene_{1}_p2 = {2}-{3}\\3\n".format(model,idx,start+2,end))
        outpart.write("{0}, gene_{1} = {2}-{3}\n".format(model,idx,start,end))

    else:
        outpart.write("{0}, gene_{1} = {2}-{3}\n".format(model,idx,start,end))
