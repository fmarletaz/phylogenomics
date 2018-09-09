#! /usr/bin/env python

import sys
import subprocess
import argparse
import os
from glob import glob
from datetime import datetime
from itertools import product
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(usage='%(prog)s [options] <CTRL>')
    parser.add_argument('ctrl', help='Control file including names of assembly and paired reads for each library')
    parser.add_argument('-p', type=int, default=8 ,dest='nproc', help='Number of threads (default: 8)')
    parser.add_argument('-f', type=int, default=2, dest='fold',help='Fold-enrichment to discard contig (default: 2)')
    parser.add_argument('-m', type=int, default=2, dest='mincov',help='Minimal coverage of contig by corresp. reads (default: 2)')
    args = parser.parse_args()
    return args

def readfa(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:], [], None
        for l in fp: # read the sequence
            if l[0] in '>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs) # yield a fasta record
            if not last: break
                
            

def kal_quant(rundir,read_1,read_2,asmb_a,asmb_b,nproc):
    if os.path.exists("{0}/{1}-{2}/abundance.tsv".format(rundir,asmb_a,asmb_b)): 
        return "Found!"
    kal_cmd="kallisto quant -t {0} -i {1}/{2}.kix -o {1}/{3}-{2} {4} {5}".format(nproc,rundir,name_b,name_a,read_1,read_2)   
    #print kal_cmd 
    child = subprocess.Popen(kal_cmd,stdout=subprocess.PIPE,shell=True)
    (out,err)=child.communicate()
    return "Done!"
 

def kal_index(rundir,asmb,name):
    kal_cmd="kallisto index -i {0}/{1}.kix {2}".format(rundir,name,asmb) 
    #print kal_cmd    
    child = subprocess.Popen(kal_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    (out,err)=child.communicate()
    #return 'Done!'

def parse_ctrl_0(ctrl):
    samples=[]
    for line in open(args.ctrl):
        if line.startswith('#'): continue
        sp_id,sp_name,sp_ct=line.rstrip().split()
        name=sp_name.split(':')[1].split('-')[0]
        asmb="{}_trinity.fa".format(name)
        read1="{}_1.fastq.gz".format(sp_id)
        read2="{}_2.fastq.gz".format(sp_id)
        samples.append((name,asmb,read1,read2))
    return samples


def trim_asmb(asmb,fold,mincov,asmb_list,rundir):
    Counts=defaultdict(list)
    for asmbb in asmb_list:
        for line in open("{0}/{1}-{2}/abundance.tsv".format(rundir,asmbb,asmb)):
            if line.startswith('target_id'): continue
            ctg,lgth,elgth,e_ct,tpm=line.rstrip().split('\t')
            Counts[ctg].append((asmbb,float(e_ct)))
    kept,lowcov,lowfold=[],[],[]
    for ctg in Counts:
        ctgCt=dict(Counts[ctg])
        refCt=ctgCt.get(asmb,0)
        ctgFilt=sorted([(asb,ct) for asb,ct in ctgCt.items() if not asb==asmb],key=lambda x:x[1],reverse=True)
        altCt=ctgFilt[0][1]
        if refCt < mincov:
            lowcov.append(ctg)
        else:
            ratio=round(refCt/float(altCt),2) if altCt>0 else round(refCt,2)
            if ratio>=fold:
                kept.append(ctg)
            else:
                #print ctg,refCt,altCt,ratio,ctgFilt
                lowfold.append(ctg)
    return kept,lowcov,lowfold
        
    
    


#######################
if __name__ == "__main__":
    
    args=parse_args()
    samples=parse_ctrl_0(args.ctrl)
    #samples=[line.rstrip().split('\t') for line in open(args.ctrl)]
    #set rundir
    #now=datetime.now()
    #rundir="Cross-conta_{0}_{1}:{2}".format(now.date(),now.hour,now.minute)
    rundir="Cross-conta_test"
    if not os.path.exists(rundir):
        os.makedirs(rundir)
    
    #kallisto index
    print "Indexing {} samples:".format(len(samples))
    for i,(name,asmb,read1,read2) in enumerate(samples):
        print "{0}/{1}.kix".format(rundir,name)
        if os.path.exists("{0}/{1}.kix".format(rundir,name)): continue
        print "[{0}/{1}] Indexing {2}...\n".format(i,len(samples),name)
        kal_index(rundir,asmb,name)
    
    total=len(samples)*len(samples)
    for i,(lib_a,lib_b) in enumerate(product(samples,repeat=2)):
        name_a,asmb_a,read_a_1,read_a_2=lib_a
        name_b,asmb_b,read_b_1,read_b_2=lib_b
        print "[{}/{}] {} versus {} ...".format(i,total,name_a,name_b),
        print kal_quant(rundir,read_a_1,read_a_2,name_a,name_b,args.nproc)
    
    stat=open('cross_conta_sum.txt','w')
    stat.write("Sample\t#contigs\t#kept\t\%low_ratio\t\t#low_ratio\t#low_coverage\n")
    assemblies=[elt[0] for elt in samples]
    print "Filtering ..."
    for name,asmb,read1,read2 in samples:
        outfa=open("{}.filt.fa".format(name),'w')
        kept,lowcov,lowfold=trim_asmb(name,args.fold,args.mincov,assemblies,rundir)
        fasta=dict((head.split()[0],seq)for head,seq in readfa(open(asmb)))
        contaFrac=round(len(lowfold)/float(len(fasta)),4)*100
        sumstat="{0}\t{1}\t{2}\t{3}%\t{4}\t{5}\n".format(name,len(fasta),len(kept),contaFrac,len(lowfold),len(lowcov))
        print sumstat
        stat.write(sumstat)
        for ctg in kept:
            outfa.write(">{0}\n{1}\n".format(ctg,fasta[ctg]))
            
    
     #def trim_asmb(asmb,fold,mincov,asmb_list,rundir):
   
        

        
        
        