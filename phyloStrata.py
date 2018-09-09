#! /usr/bin/env python

import sys
from glob import glob
from collections import defaultdict
from collections import Counter
import operator
from ete3 import Tree
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from itertools import combinations

def reject_outliers(datdict, m = 20.):
    data=np.array(datdict.items(),dtype=[('tax', 'S50'),('dist', 'f4')])
    d = data['dist'] - np.median(data['dist'])
    mdev = np.median(np.abs(d))
    s = d/mdev if mdev else 0.
    return data[s>m],data[s<m]

def phylum_mono(t,taxPhyl):
    # Checks monophyly of each phylum
    Phyla=list(set(taxPhyl.values()))
    phylumSt={}
    for phyl in Phyla:
        plf=t.search_nodes(phylum=phyl)
        if len(plf)>1:
            status=t.check_monophyly(values=[phyl], target_attr="phylum")
            mst=status[0]
            #print phyl,[l.name for l in plf],status[0]
        else:
            mst=None
            #print phyl, [l.name for l in plf],'not enough taxa!'
        phylumSt[phyl]=mst
    return phylumSt

def saturation(fafile,trfile):
    #compute pairwise %Id
    aln = AlignIO.read(open(fafile), 'fasta')
    calculator = DistanceCalculator('blosum62')
    dm = calculator.get_distance(aln)
    pwdists,lfpairs=[],[]
    for i,j in combinations(range(len(dm.names)),2):
        lfpairs.append((dm.names[i],dm.names[j]))
        pwdists.append(dm[i][j])

    #Compute patristic Distance from ML Tree
    t=Tree(open(trfile).readline())
    padists=[t.get_distance(lf1,lf2) for lf1,lf2 in lfpairs]
    slope,intersect=np.polyfit(padists,pwdists,1)

    return slope

#Usage: concatenate-ext.py <taxlist> <fasta files...>
if not len(sys.argv) >= 4:
    sys.exit('Usage: phyloStrata.py <taxlist> <suffix> <tree files...>')


#Define phyla and taxlist
TaxPhyl={}
for line in open(sys.argv[1]):
    if line.startswith('#') or line.startswith('$'): continue
    if not line.strip(): continue
    species,phylum=line.rstrip().split()
    TaxPhyl[species]=phylum

out=open('Phylostats-{0}.txt'.format(sys.argv[2]),'w')

out.write("Marker\tntaxa\talpha\ttmed\ttlen\tfast_len\tmonophyletic\tsaturation\tbrL_rej\trej_ntax\n")
nbTaxReject=defaultdict(int)
for trfile in sys.argv[3:]:
    path,fnm=trfile.rsplit('/',1)
    gid=fnm.rsplit('.',2)[0]
    #OG93.al.tr.
    print gid,'...'
    t_str=open(trfile).readline()
    t=Tree(t_str)
    #Add phylum label to branches
    for leaf in t:
        leaf.add_features(phylum=TaxPhyl.get(leaf.name, "none"))
    #print t.get_ascii(attributes=["name", "phylum"], show_internal=False)
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)
    try:
        ancestor = t.get_common_ancestor(t.search_nodes(phylum='Porifera'))
        t.set_outgroup(ancestor)
    except:
        #print "root phylum not monophyletic!"
        #print t.get_ascii(attributes=["name", "phylum"], show_internal=False)
        for pnd in t.search_nodes(phylum='Porifera'):
            try:
                t.set_outgroup(t&pnd)
            except:
                #print 'trying another root...'
                continue
    #Check phylum mononphyly
    Monophyletic=phylum_mono(t,TaxPhyl)
    mct='/'.join(map(str,[Counter(Monophyletic.values())[cat] for cat in [True,False,None]]))
    #Calculate p-dists and reject taxa deviation from median by an order of magnitude
    pdists=dict((leaf.name,round(t.get_distance(leaf),5)) for leaf in t)
    rejected,kept=reject_outliers(pdists)
    #print len(rejected),rejected,len(pdists)
    for tax in rejected['tax']:
        nbTaxReject[tax]+=1

    #generate new alignments with rejected...
    filtered=[]
    fafile='{0}/{1}.al.hc.tr.fa'.format(path,gid)
    for rec in SeqIO.parse(fafile,'fasta'):
        if rec.id in set(kept['tax']):
            filtered.append(rec)
    SeqIO.write(filtered,'{0}/{1}.al.hc.tr.ft.fa'.format(path,gid),'fasta')

    #calculate saturation
    satSlp=saturation(fafile,trfile)
    tlen=sum([node.dist for node in t])
    #grab alpha
    mod=[l for l in open('{0}/{1}.raxml.bestModel'.format(path,gid))]
    alpha=mod[0].split('{')[1].split('}')[0]
    alen=mod[0].split('-')[-1].strip()

    #mean of fasta evolving taxa
    fastPhyl=set(['Rotifera','Gnathostomulida'])
    fastDist=[td['dist'] for td in kept if TaxPhyl[td['tax']] in fastPhyl]
    # print fastDist
    tmed="{0:.4f}".format(np.median(pdists.values()))
    fmean="{0:.4f}".format(np.mean(fastDist))
    outList=[gid,len(pdists),alpha,alen,tlen,tmed,fmean,mct,satSlp,len(rejected),','.join(rejected['tax'])]
    #print outList
    out.write('\t'.join(map(str,outList))+'\n')

for sp,nb in sorted(nbTaxReject.items(), key=operator.itemgetter(1)):
    print sp,nb
