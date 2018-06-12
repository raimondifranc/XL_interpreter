#!/usr/bin/python2.7

import os,sys,operator, math, argparse
import gzip, sqlite3, glob
import numpy as np
import glob, itertools
from Bio import PDB
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import *
from difflib import SequenceMatcher
import timeit

#############################################################################################################################
#### XLinterpreter: a tool for fast structural annotation of restraints from structural proteomics (e.g. primarily MS/XL)####
########################################### Raimondi F. and Russell RB ######################################################

AA=['K','R','H','E','D','N','Q','Y','S','T','C','A','G','P','F','I','L','V','M','W']


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def GetPdbRes(struct, ch):
  resmaps={}
  for m in struct:
    for c in m:
      if c.id == ch:
        ii=1
        for r in c:
          resmaps[ii]=r.id[1]
          ii+=1
  
  return resmaps
  
def GetID(aa):
  labels=""
  resid=""
  for c in aa:
    if c.isdigit() == False:
      labels+=c
    elif c.isdigit() == True:
      resid+=c

  if int(resid)>0:
    return int(resid)



def FindPos(aa, start, seq):
  idx=0
  #print aa, start, seq
  for s in seq:
    if int(start) == int(aa):
      #print "Hello"
      return idx
    #print s, start,idx,aa
    if s != "-" and s != ".":
      start+=1
    idx+=1


def GetPos(idx,start,seq):
  pos=0
  #print idx, start, seq
  for s in seq:
    if int(pos) == int(idx):
      return start,s
    if s != "-" and s != ".":
      start+=1
    pos+=1


def MapVariants(tstart,tseq,rstart,rseq,var,simcutoff):
  offset=5
  #simcutoff=0
  varidx=FindPos(var,tstart,tseq)
  fampos,res=GetPos(varidx,rstart,rseq)
  #Checking on a sequence window equal to offset if the query and reference sequences are identical or not
  #Some tolerance is applied through the similar function to allow small differences in the sequence
  ###Remember to remove the similarity cutoff when performing remote homology searches with either PSI BLAST or HHLIBITS!!!
  #return fampos, res  
  if similar(tseq[varidx:varidx+offset], rseq[varidx:varidx+offset]) >= simcutoff and tseq[varidx:varidx+offset] != "" and rseq[varidx:varidx+offset] != "":# or tseq[varidx:varidx-offset] == rseq[varidx:varidx-offset]:
    #print tseq[varidx:varidx+offset], rseq[varidx:varidx+offset] 
    return fampos, res
  elif  similar(tseq[varidx:varidx-offset], rseq[varidx:varidx-offset]) >= simcutoff  and tseq[varidx:varidx-offset] != "" and rseq[varidx:varidx-offset] != "":
    #print tseq[varidx:varidx-offset], rseq[varidx:varidx-offset] 
    return fampos, res
  else:
    #print tseq[varidx:varidx+offset], rseq[varidx:varidx+offset]
    return 0,""



def calc_dist(resA, resB, atomtype):
  ###Check the distance between specified atom types from two residues
  for atm1 in resA:
    if (atm1.fullname).strip() == atomtype:
      #print "OK", resA, atm1.fullname, atm1.coord
      atm_sele1=atm1
      
  for atm2 in resB:
    if (atm2.fullname).strip() == atomtype:
      #print "OK",resB, atm2.fullname, atm2.coord
      atm_sele2=atm2
  diff  = atm_sele1.coord - atm_sele2.coord
  return np.sqrt(np.sum(diff * diff))


def ExtractFasta(db, idslist, fastaout):
  alias1={}
  alias2={}
  alias3={}
  fastaoutfile=open(fastaout,"w")

  conn = sqlite3.connect(db)

  c = conn.cursor()
  ###We need to retrieve all the pdb2uniprot matches for the input XL
  sql_query=""
  #sql_query = 'select * from swissprot where ' + ' or '.join(("id_ac = '" + str(n)+"'"+" or id_id = '" + str(n)+"'"+" or id_gn = '" + str(n)+"'" for n in idslist)) 
  for n in idslist:
    ###This is a workaround I implemented to process as much data as possible from XlinkAnalyzer
    ###There are a lot of identifiers screwed up, i.e. they are lower case with no correspondence in the Uniprot
    ###while instead Uniprot cotains the corresponding uppercase identifier.
    sql_query = 'select * from swissprot where id_ac = '+"'"+ str(n)+"'"+" or id_id = '" + str(n)+"'"+" or id_gn = '" + str(n)+"'"
    c.execute(sql_query)
    data=c.fetchall()
    if len(data) == 0:
      N=n.upper()
      sql_query = 'select * from swissprot where id_ac = '+"'"+ str(N)+"'"+" or id_id = '" + str(N)+"'"+" or id_gn = '" + str(N)+"'"
      #c.execute(sql_query):
      #data=c.fetchall()
    print sql_query
    
    for line in c.execute(sql_query):
      uniac=line[0]
      unid=line[1]
      gn=line[2]
      fasta=line[3]
      fastaoutfile.write(fasta)
      alias1[uniac]=unid
      alias2[uniac]=gn
      alias3[n]=uniac

  conn.close()
  
  return alias1, alias2, alias3


def GetBlastOut(infile, evalcutoff,nriter,overlap_flag):
  blastout=open(infile,"r")
  pdb_matches={}
  flag1=0
  flag2=0
  flag3=0
  start=0
  nrround=0
  protein_matches=0
  buffer_list={}
  for line in blastout:
    if line.find("Results from round")!= -1:

      start=1
      nrround=int((line.split()[3]).strip("\n"))
      continue
      #if int((line.split()[3]).strip("\n")) == nriter:
      #  #print int((line.split()[3]).strip("\n")), nriter
      #  start=1
      #  continue
      #else:
      #  start=0
      
    if start == 1 and len(line.split()) > 0 and line.split()[0] == "Query=":
      #print "Accession OK"
      flag1=1
      uniac=line.split()[1]
      unid=uniac.split("|")[2]
      uniac=uniac.split("|")[1]
      matches=0

    if start == 1 and line[0] == ">":
      #print "Seq header OK"
      #pdbid=line.split()[0]
      #pdbid=(pdbid.split(":")[0])[1:]
      #pdbid=pdbid[1:]
      chain=line.split("_")[1]
      chain=chain.split()[0]
      pdbid=line.split("_")[0]
      pdbid=(pdbid.lstrip(">")).lstrip()
      flag2=1
      flag3=0
      cont=0
      targ_start=0
      targ_seq=""
      ref_start=0
      ref_seq=""
      targ_stop=0
      ref_stop=0
      f1=0
      f2=0
      matches+=1
      continue

    if start == 1 and flag1 == 1 and flag2 == 1:
      if len(line.split()) > 0 and line.split()[0] == "Score" and line.split()[1] == "=":
        #print "Score OK"
        Evalue=(line.split()[7]).strip(",")
        if Evalue[0] == "e":
          Evalue="1"+Evalue
        Evalue=float(Evalue)          
      if len(line.split()) > 0 and line.split()[0] == "Identities" and line.split()[1] == "=":
        #print "Identity OK"
        identity=(line.split()[3]).strip(",")
        identity=identity.strip("(")
        identity=identity.strip(")")
        identity=identity.strip("%")
        if Evalue < evalcutoff:
          flag3=1
          flag2=0
          continue
        continue
      continue

    if start == 1 and  len(line.split()) > 0 and line.split()[0] == "Query" and flag3 == 1:
      #print "Sequences OK", targ_start, targ_seq
      targ_start=line.split()[1]
      targ_seq=line.split()[2]
      targ_stop=line.split()[3]
      f1=1

    if start == 1 and len(line.split()) > 0 and line.split()[0] == "Sbjct" and flag3 == 1:
      ref_start=line.split()[1]
      ref_seq=line.split()[2]
      ref_stop=line.split()[3]
      f2=1

    if start == 1 and flag3 == 1 and f1 == 1 and f2 == 1:
      ###Using targ_range and ref_range to dynamically monitor whether a given sequence alignment is overlapping with previously aligned regions
      ###If so, discard them to avoid confusing sequence-structure assignments which more often happen with repeat motives
      f1=0
      f2=0
      if buffer_list.has_key(pdbid) == False:
        buffer_list[pdbid]={}
        buffer_list[pdbid][uniac]=[]
        buffer_list[pdbid][uniac].append([Evalue,int(targ_start), int(targ_stop),targ_seq,int(ref_start), int(ref_stop), ref_seq, chain, identity])
        targ_range=[]
        for ii in range(int(targ_start), int(targ_stop)):
          targ_range.append(ii)
        ref_range=[]
        for ii in range(int(ref_start), int(ref_stop)):
          ref_range.append(ii)
      elif buffer_list.has_key(pdbid) and buffer_list[pdbid].has_key(uniac) == False:
        buffer_list[pdbid][uniac]=[]
        buffer_list[pdbid][uniac].append([Evalue,int(targ_start), int(targ_stop),targ_seq,int(ref_start), int(ref_stop), ref_seq, chain, identity])
        targ_range=[]
        for ii in range(int(targ_start), int(targ_stop)):
          targ_range.append(ii)
        ref_range=[]
        for ii in range(int(ref_start), int(ref_stop)):
          ref_range.append(ii)
      else:
        if overlap_flag == 0:
          if int(targ_start) not in targ_range and int(targ_stop) not in targ_range and int(ref_start) not in ref_range and int(ref_stop) not in ref_range:
            buffer_list[pdbid][uniac].append([Evalue,int(targ_start), int(targ_stop),targ_seq,int(ref_start), int(ref_stop), ref_seq, chain, identity])
            for ii in range(int(targ_start), int(targ_stop)):
              targ_range.append(ii)
            for ii in range(int(ref_start), int(ref_stop)):
              ref_range.append(ii)
        elif overlap_flag == 1:
          buffer_list[pdbid][uniac].append([Evalue,int(targ_start), int(targ_stop),targ_seq,int(ref_start), int(ref_stop), ref_seq, chain, identity])
          
      ###Some debugging here - remove
      #if buffer_list.has_key("4Q9U"):
      #  print "Buffer 4Q9U", buffer_list["4Q9U"].keys() 
      #if buffer_list.has_key("4Q9U"):
      #  if buffer_list["4Q9U"].has_key("P20339"):
      #    print "Final 4Q9U", "P20339", buffer_list["4Q9U"]["P20339"], uniac, unid, pdbid, nrround  
      
    if (start == 1 and line.find("Search has CONVERGED!") != -1) or (start == 1 and line.find("Effective search space used:") != -1 and nrround == nriter):
      #print uniac, unid, nrround
      #print buffer_list
      #if buffer_list.has_key("4Q9U"):
      #  if buffer_list["4Q9U"].has_key("P20339"):
      #    print "Final 4Q9U", "P20339", buffer_list["4Q9U"]["P20339"]    
      for pid in buffer_list.keys():
        if pdb_matches.has_key(pid) == False:
          pdb_matches[pid]={}
        for uac in buffer_list[pid].keys():
          #print pdbid, uniac
          pdb_matches[pid][uac]=buffer_list[pid][uac]
      start=0
      buffer_list={}
      protein_matches+=1
      #if pdb_matches.has_key("4Q9U"):
      #  print "Final 4Q9U", pdb_matches["4Q9U"].keys()      
  
  #print protein_matches      
  return pdb_matches

def aaLabel(resname):
  if resname in standard_aa_names:
    aaout=three_to_one(resname)
  else:
    aaout=resname

  return aaout

def GetResIdx(chain, residue):
  cc=0
  for rr in chain:
    if rr.id[1] == residue.id[1]:
      return cc  
    cc+=1

def GetHeader(headline):
  ###I also need to determine here the Field Separator type 
  count=0
  F1=0
  F2=0
  F3=0
  F4=0
  FS=""
  headline=headline.strip("\n")
  headline=headline.strip("\r")
  if headline.find("\t") != -1:
    FS="\t"
  elif headline.find(",") != -1:
    FS=","
  elif headline.find(";") != -1:
    FS=";"
  for field in headline.split(FS):
    if field.find("Protein1") != -1:
      F1=count
    if field.find("Protein2") != -1:
      F2=count
    if field.find("AbsPos1") != -1:
      F3=count
    if field.find("AbsPos2") != -1:
      F4=count
    count+=1
  return F1, F2, F3, F4, FS

def ReadRestr(restraint_list, LOGFILE):
  lc=0
  rest_list=[]
  prot_list=[]
  for l in open(restraint_list, "r"):
    if lc==0:
      f1,f2,f3,f4,fs=GetHeader(l)
      #LOGFILE.write("%s %s %s %s %s\n" % (f1,f2,f3,f4,fs))
    if len(l.split(fs)) >= 3 and lc>0 and l.split(fs)[f1] != "-" and l.split(fs)[f2] != "-": #let's not consider at this stage monolinks
      if l.split(fs)[f1].find("|") != -1:
        prot1=(l.split(fs)[f1]).split("|")[1]
      else:
        prot1=(l.split(fs)[f1])
      if l.split(fs)[f2].find("|") != -1:      
        prot2=(l.split(fs)[f2]).split("|")[1]
      else:
        prot2=(l.split(fs)[f2])
      pos1=l.split(fs)[f3]
      pos2=l.split(fs)[f4]
      if prot1 not in prot_list:
        prot_list.append(prot1)
      if prot2 not in prot_list:
        prot_list.append(prot2)
      rest_list.append([prot1, prot2, pos1, pos2])
    if l[0] != "#":
      lc+=1
  return rest_list, prot_list


def XLScore(distances, threshold, total, ac2gn, ac2xl):
  interfaces={}
  max_sc=len(distances.keys())
  tot_max_sc=total
  d_list=[]
  ###Getting the maximum distance for normalization of the penalty term
  for d in distances.keys():
    for val in distances[d]:
      d_list.append(val)
  max_d=max(d_list)
  sc=0
  count=0
  tot_xl={}
  for d in distances.keys():
    dis=distances[d]
    p1=d.split("|")[0]
    p1=p1.split("/")[0]
    p2=d.split("|")[1]
    p2=p2.split("/")[0]
    gn1=""
    gn2=""
    if ac2gn.has_key(p1):
      gn1=ac2gn[p1]
    if ac2gn.has_key(p2):
      gn2=ac2gn[p2]
    intf=gn1+"|"+gn2
    invintf=gn2+"|"+gn1
    if min(dis) <= threshold:
      part_sc=1
      count+=1
      if interfaces.has_key(intf) == False and interfaces.has_key(invintf) == False:
        interfaces[intf]=1
      elif interfaces.has_key(intf):
        interfaces[intf]+=1
      elif interfaces.has_key(invintf):
        interfaces[invintf]+=1
      
    else:
      dev=min(dis)-threshold
      #print dev, max_d, threshold
      part_sc=1-(float(dev)/float((max_d-threshold)))
    sc+=part_sc
    
    if ac2xl.has_key(p1):
      if tot_xl.has_key(gn1) == False:
        tot_xl[gn1]=ac2xl[p1]
    if ac2xl.has_key(p2):
      if tot_xl.has_key(gn2) == False:
        tot_xl[gn2]=ac2xl[p2]
    
    strout=""
    for intf in interfaces.keys():
      gn1=(intf.split("|")[0])
      gn2=(intf.split("|")[1])
      totxl=[]
      totxl=len(set(tot_xl[gn1]+tot_xl[gn2]))
      fraction=float(interfaces[intf])/float(totxl)
      strout+=intf+":("+str(interfaces[intf])+":"+str(totxl)+"=%3.2f)," % (fraction)

  return count, float(sc)/float(max_sc), float(sc)/float(tot_max_sc), interfaces, strout


def GenPym(xlinfo,outname,inpdb):
  for pdb in xlinfo.keys():
    print outname, pdb
    pymfile=open("%s_%s.pml" % (outname,pdb),"w")
    if inpdb != None:
      #os.popen("cp /net/home.isilon/ds-russell/pdb/%s/pdb%s.ent.gz ." % (pdb[1:3], pdb))
      pymfile.write("load %s\n" % (pdb))    
    else:
      os.popen("cp /net/home.isilon/ds-russell/pdb/%s/pdb%s.ent.gz ." % (pdb[1:3], pdb))
      pymfile.write("load pdb%s.ent.gz\n" % (pdb))

    outstring="""hide everything, all
    #set cartoon_smooth_loops=0
    #set cartoon_fancy_helices=1
    #set cartoon_cylindrical_helices = 1
    set transparency_mode=1
    set ray_shadows=0
    set ray_transparency_shadows=0
    #set ribbon_trace,1
    #set cartoon_transparency=0.5"""
    
    pymfile.write(outstring)
    #pymfile.write("\nshow ribbon\n")
    pymfile.write("\nshow cartoon\n")
    
    cc=0
    for xl in xlinfo[pdb]:
      c1=xl[0]
      c2=xl[1]
      r1=xl[2]
      r2=xl[3]
      idx1=xl[4]
      idx2=xl[5]
      satisfaction=xl[6]
      if inpdb != None:
        sele1="/%s//%s/%s/CA" % (pdb[:-4],c1,idx1)
        sele2="/%s//%s/%s/CA" % (pdb[:-4],c2,idx2)        
      else:
        sele1="/pdb%s//%s/%s/CA" % (pdb,c1,idx1)
        sele2="/pdb%s//%s/%s/CA" % (pdb,c2,idx2)
      
      if satisfaction == "satisfied":
        pymfile.write('cmd.distance("dist%s", "(%s)", "(%s)")\n' % (cc, sele1, sele2))
        pymfile.write('color grey80,dist%s\n' % (cc))
        pymfile.write('set label_outline_color, grey80, dist*\n')
        pymfile.write('et label_color, grey80, dist*\n')
        
      elif satisfaction == "unsatisfied":
        pymfile.write('cmd.distance("undist%s", "(%s)", "(%s)")\n' % (cc, sele1, sele2))
        pymfile.write('color red,undist%s\n' % (cc))
        pymfile.write('set label_outline_color, red, undist*\n')
        pymfile.write('et label_color, red, undist*\n')
      pymfile.write('show spheres, %s or %s\n' % (sele1,sele2))
      pymfile.write('label %s, "%s"\n' % (sele1, r1))
      pymfile.write('label %s, "%s"\n' % (sele2, r2))
      cc+=1

    pymfile.write('util.cbc pdb%s\n' % (pdb))
  
    outstring="""set dash_gap, 0
    #hide label, all
    set dash_width, 5
    set float_labels, on
    set ortho=on
    set depth_cue=off
    set antialias=1

    cmd.bg_color('white')
    set ray_shadows=0
    set ray_transparency_shadows=1

    """
    pymfile.write(outstring)    

  return

def GetNonRedPDB(pdb_site_list):
  ###1)Check the pdb with the highest xl coverage
  ###2)If the coverage between 2 pdb is the same , then take the one with best Evalue from PSI-BLAST
  
  residue_list={}
  eval_list={}
  for pdb in pdb_site_list.keys():
    for match in pdb_site_list[pdb]:
      p1=match[2]
      p2=match[3]
      r1=match[4]
      r2=match[5]
      E1=match[10]
      E2=match[12]
      pair=str(p1)+"|"+str(p2)
      invpair=str(p2)+"|"+str(p1)
      
      setinfo=[[p1,r1,E1],[p2,r2,E2]]
      for si in setinfo:
        p=si[0]
        r=si[1]
        E=si[2]
        #print p
        if residue_list.has_key(pdb) == False:
          residue_list[pdb]={}
          residue_list[pdb][p]=[r]
          eval_list[pdb]={}
          eval_list[pdb][p]=[E]
        elif residue_list.has_key(pdb) and residue_list[pdb].has_key(p) == False:
          residue_list[pdb][p]=[r]
          eval_list[pdb][p]=[E]
        else:
          if r not in residue_list[pdb][p]:
            residue_list[pdb][p].append(r)
          if E not in eval_list[pdb][p]:
            eval_list[pdb][p].append(E)
  
  pdb_matches={}
  for pdb in residue_list.keys():
    totres=0
    for p in residue_list[pdb].keys():
      totres+=len(residue_list[pdb][p])   
    #sorted_pdb=sorted(residue_list[p].items(), key=lambda a: len(a[1]), reverse=True)
    print "INITIAL LIST",pdb, len(residue_list[pdb].keys()), totres, residue_list[pdb]
    pdb_matches[pdb]=totres
  sorted_pdb_matches=sorted(pdb_matches.items(), key=operator.itemgetter(1), reverse=True)

  ###Sorting here according to Evalue, so to consider first matches with better E-value
  ###Remember that in some cases, there might be template matches with same nr of XL satisfied but different E-value. In that case, the one with best E value has to be considered. 
  
  final_pdb=[]
  final_uniac={}  
  for pdb, cc in sorted_pdb_matches:
    for uniac in residue_list[pdb].keys():
      if final_uniac.has_key(uniac) == False:
        final_uniac[uniac]=[]
        for res in residue_list[pdb][uniac]:
          final_uniac[uniac].append(res)
        if pdb not in final_pdb:
          final_pdb.append(pdb)
      elif final_uniac.has_key(uniac):
        for res in residue_list[pdb][uniac]:
          if res not in final_uniac[uniac]:
            final_uniac[uniac].append(res)
            if pdb not in final_pdb:
              final_pdb.append(pdb)

  for pdb in final_pdb:
    print "FINAL LIST", pdb, residue_list[pdb]
  
  return final_pdb

##########################################################

def main():
  parser = argparse.ArgumentParser(description='Tool to map experimental restraints (usually MS/XL) on available 3D structures from the PDB | Basic usage: XL_interprets.py -restraint Input_XL_list.txt')
  parser.add_argument( '-input', action="store", dest='input', required=False, help='List of restraints', default="/net/home.isilon/ag-russell/bq_fraimondi/PROJECTs/XL_interpreter/XlinkAnalyzer/TFIIIC/xlinks/TFIIIC_allxlinks_without_tDNA.csv")
  parser.add_argument( '-inputpdb', action="store", dest='inputpdb', required=False, help='Input PDB: you must also provide a blast DB file with the fasta sequences from the PDB file', default=None)
  parser.add_argument( '-cutoff', action="store", dest='cutoff', required=False, help='Distance cutoff', type=float, default=35 )
  parser.add_argument( '-db', action="store", dest='db', required=False, help='SQL database with PDB2Uniprot blast alignments (default=/net/home.isilon/ag-russell/bq_fraimondi/PROJECTs/XL_interpreter/test/pdbseq)', default="/net/home.isilon/ag-russell/bq_fraimondi/PROJECTs/XL_interpreter/test/pdbseq")
  parser.add_argument( '-output', action="store", dest='output', required=False, help='Output filename', default="test")
  parser.add_argument( '-atom2monitor', action="store", dest='atom2monitor', required=False, help='Atom type whose distance will be monitored in order to model the restraints', default="CA")
  parser.add_argument( '-pymol', action="store_true", dest='pymol', required=False, help='Generates 3D representations of the mapped XL' )
  parser.add_argument( '-mode', action="store", dest='mode', required=False, help='Execution mode fast|thorough \n   -fast=only a non-redundant list of best pdb-query sequence matches is considered\n   -thorough=all pdb-query sequence matches are considered', default="fast" )
  parser.add_argument( '-blast_eval', action="store", dest='Eval', required=False, help='E-value cutoff for Blast results', type=float, default=10 )
  parser.add_argument( '-blast_v', action="store", dest='blast_v', required=False, help='Number of database sequences to show one-line descriptions', type=float, default=500)
  parser.add_argument( '-blast_b', action="store", dest='blast_b', required=False, help='Number of database sequence to show alignments for', type=float, default=250)
  parser.add_argument( '-psiblast_iterations', action="store", dest='psiblast_nriter', required=False, help='Number of PsiBlast iterations', type=int, default=2 )
  parser.add_argument( '-ali_overlap', action="store_true", dest='ali_overlap', required=False, help='Allows overlap between aligned sequence regions' )
  
  
  
  arguments = parser.parse_args()
  restr = arguments.input
  cutoff = float(arguments.cutoff)
  Eval = float(arguments.Eval)
  output = arguments.output
  outfile=open(output+"_XL_matches.txt", "w+t")
  atom2monitor=arguments.atom2monitor
  pymol=arguments.pymol
  blast_v=arguments.blast_v
  blast_b=arguments.blast_b
  psiblast_nriter=int(arguments.psiblast_nriter)
  mode=arguments.mode
  ali_overlap_flag=0
  if arguments.ali_overlap:
    ali_overlap_flag=1
  
  logfile=open(output+"_log.txt", "w+t")
  sysargv_str=""
  for sa in sys.argv:
    sysargv_str += sa + " "
  logfile.write("Command line: "+str(sysargv_str)+"\n")
  
  db=arguments.db
  inputpdb=arguments.inputpdb
  
  start_time=timeit.default_timer()
  
#### 1) Reading the experimental restraint input
  restr_list, id_list=ReadRestr(restr,logfile)
  
  
#### 2) Retrieving corresponding fasta sequences 
  fastaname=(output+"_inputseq.fasta")
  pac2pid, pac2gn,origid2pac = ExtractFasta("/net/home.isilon/ag-russell/bq_fraimondi/DBs/Uniprot/sprot_varsplic_seq_ids.db", id_list, fastaname)

#### 3) PSIBLAST sequences on the PDB
  blastout_file="%s_blast.out" % (output)
#### Remember to uncomment this once the pevolution issue will be fixed and so the pdbseq
  #os.system("/net/home.isilon/ag-russell/install/CentOS-7.2.1511-x86_64/bin/blastall -p blastp -e %f -F F -m 0 -d /net/home.isilon/ds-russell/blastdb/pdbseq -i %s > %s" % (Eval, fastaname, blastout_file))
  os.system("/usr/bin/psiblast -evalue %f -db %s -query %s -num_iterations %d > %s" % (Eval, db, fastaname, psiblast_nriter, blastout_file))

#### 4) Getting the Blast search output for sequences in the query against PDB reference sequences
  pdb2uniprot=GetBlastOut(blastout_file,Eval,psiblast_nriter,ali_overlap_flag)
  elapsed1 = timeit.default_timer() - start_time
  #print pdb2uniprot
#### 5) Determining how many XL sites are found within each PDB
  pdb_w_sites={}
  matched_xl=0
  #logfile.write(str(restr_list)+"\n")
  uniac2xl={}
  for xl in restr_list:
    if origid2pac.has_key(xl[0]):
      p1=origid2pac[xl[0]]
    else:
      p1=xl[0]
    if origid2pac.has_key(xl[1]):
      p2=origid2pac[xl[1]]
    else:
      p2=xl[1]
    r1=int(xl[2])
    r2=int(xl[3])
    xl=p1+"/"+xl[2]+"|"+p2+"/"+xl[3]
    invxl=p2+"/"+xl[3]+"|"+p1+"/"+xl[2]
    flag_match=0
    if uniac2xl.has_key(p1) == False:
      uniac2xl[p1]=[xl]
    elif xl not in uniac2xl[p1] and invxl not in uniac2xl[p1]:
      uniac2xl[p1].append(xl)
    if uniac2xl.has_key(p2) == False:
      uniac2xl[p2]=[xl]
    elif xl not in uniac2xl[p2] and invxl not in uniac2xl[p2]:
      uniac2xl[p2].append(xl)
    
    #print p1, p2, r1, r2
    for pdb in pdb2uniprot.keys():
      if pdb2uniprot[pdb].has_key(p1) and pdb2uniprot[pdb].has_key(p2):
        #print pdb, p1, p2
        ###Checking if the XL sites are within the blast match ranges
        flag1=0
        flag2=0
        m1=[]
        m2=[]
        ###First XL site
        for m in pdb2uniprot[pdb][p1]:
          if int(r1) in range(m[1], m[2]):
            #print p1, r1
            flag1=1
            m1.append(m)
            
        ###Second XL site
        for m in pdb2uniprot[pdb][p2]:
          if int(r2) in range(m[1], m[2]):
            #print p2, r2
            flag2=1
            m2.append(m)
            
        ### We don't know about the stechiometry of the XL, so we need to check all the possible chain combinations here
        if flag1 == 1 and flag2 == 1:
          for match1 in m1:
            for match2 in m2:
              E1=match1[0]
              I1=match1[8]
              E2=match2[0]
              I2=match2[8]
              idx1, aa1=MapVariants(match1[1],match1[3],match1[4],match1[6],int(r1),0.0)
              idx2, aa2=MapVariants(match2[1],match2[3],match2[4],match2[6],int(r2),0.0)
              flag_match=1
              #print "MATCHED_SITE:", pdb, match1[7], match2[7], p1, p2, r1, r2, idx1, idx2, aa1, aa2  
              if pdb_w_sites.has_key(pdb) == False:
                pdb_w_sites[pdb]=[[match1[7], match2[7], p1, p2, r1, r2, idx1, idx2, aa1, aa2, E1, I1, E2, I2]]
              else:
                pdb_w_sites[pdb].append([match1[7], match2[7], p1, p2, r1, r2, idx1, idx2, aa1, aa2, E1, I1, E2, I2])

  print "4Q9U", pdb2uniprot["4Q9U"].keys()
  #loprint pdb_w_sites.keys()  
###Make here pdb_w_sites non redundant if the mode flag is set to fast:
###With fast, we consider the pdb with highest number of xl mapped.
###Pdbs with xl already mapped in other pdb, are just skipped
  if mode == "fast":
    pdb_w_sites_final=GetNonRedPDB(pdb_w_sites)
  elif mode == "thorough":
    pdb_w_sites_final=pdb_w_sites.keys()

#### 6) Checking whether the XL sites found are compatible with XL theoretical distances (35A for DSS, 31A for DSG for example)
#### 7) Sorting pdb template according to restraint-satisfaction score

  outfile.write("MATCHED_SITE:\tPDB ID\tChain Site1\tChain Site2\tIdentity S1\tIdentity S2\tEval S1\tEval S2\tGene Name 1\t Gene Name 2\tUniprot Accession S1\tUniprot Accession S2\tS1 AA position (sequence)\tS2 AA position (sequence)\tS1 AA resseq PDB\tS2 AA resseq PDB\tS1 AA seqpos (structure sequence)\tS2 AA seqpos (structure sequence)\tS1 AA (structure)\tS2 AA (structure)\t%s-%s distance (Angstrom)\n" % (atom2monitor, atom2monitor))
  pdb_stat={}
  pdb4pym={}
  #print pdb_w_sites_final
  XL_satisfied={}
  for pdb in pdb_w_sites_final:
###Calculating the distances on biological assembly only to speed up calculations
### Maybe I can speed up the process by first looping over PDB residues and internally checking
    if inputpdb != None:
      pdbpath=inputpdb
      Pdb=pdb
      #pdb=pdb.lower()
    else:
      Pdb=pdb
      pdb=pdb.lower()
      pdbpath="/net/home.isilon/ds-russell/pdb-biounit/%s/%s.pdb*.gz" % (pdb[1:3],pdb)
    #print pdb, Pdb, pdbpath
    for pdbdir in glob.glob(pdbpath):
      #print pdbdir
      #print pdb2uniprot[Pdb]
      if pdb_stat.has_key(pdb) == False:
        pdb_stat[pdb]={}
      if os.path.isfile(pdbdir):
        #print pdb_w_sites[Pdb]
        for match in pdb_w_sites[Pdb]:
          c1=match[0]
          c2=match[1]
          idx1=match[6]
          idx2=match[7]
          prot1=match[2]
          prot2=match[3]
          aa1=match[4]
          aa2=match[5]
          E1=match[10]
          I1=match[11]
          E2=match[12]
          I2=match[13]
          gn1=""
          gn2=""
          if pac2gn.has_key(prot1):
            gn1=pac2gn[prot1]
          if pac2gn.has_key(prot2):
            gn2=pac2gn[prot2]
          xl=prot1+"/"+str(aa1)+"|"+prot2+"/"+str(aa2)
          invxl=prot2+"/"+str(aa2)+"|"+prot1+"/"+str(aa1)
          if pdbdir[-3:] == ".gz":
            pdbpath=gzip.open(pdbdir, "r")
          else:
            pdbpath=open(pdbdir, "r")
          parser = PDB.PDBParser(QUIET=True)
          struct = parser.get_structure('XXX',pdbpath)
          chain1_check=0
          chain2_check=0
          #for model in struct:
          model = struct[0]
          for chain in model:
            if chain.get_id() == c1:
              chain1_check = 1
            if chain.get_id() == c2:
              chain2_check = 1
          if chain1_check == 1 and chain2_check == 1:
            for ii, res1 in enumerate(model[c1]): 
              if ii+1 == idx1:  ###Checking that pdb residue sequential number (i.e. from PDB Calpha sequence) matches
                for jj, res2 in enumerate(model[c2]):
                  if jj+1 == idx2:
                    distance = calc_dist(res1, res2, atom2monitor)
                    outfile.write("MATCHED_SITE:\t"+pdb+"\t"+c1+"\t"+c2+"\t"+str(I1)+"\t"+str(I2)+"\t"+str(E1)+"\t"+str(E2)+"\t"+gn1+"\t"+gn2+"\t"+str(match[2])+"\t"+str(match[3])+"\t"+str(match[4])+"\t"+str(match[5])+"\t"+str(res1.get_full_id()[3][1])+"\t"+str(res2.get_full_id()[3][1])+"\t"+str(idx1)+"\t"+str(idx2)+"\t"+str(match[8])+"\t"+str(match[9])+"\t"+str(distance)+"\n")
                    if pymol and distance < cutoff:
                      if pdb4pym.has_key(pdb) == False:
                        pdb4pym[pdb]=[[c1,c2,aa1,aa2,res1.get_full_id()[3][1],res2.get_full_id()[3][1], "satisfied"]]
                      else:
                        pdb4pym[pdb].append([c1,c2,aa1,aa2,res1.get_full_id()[3][1],res2.get_full_id()[3][1], "satisfied"])
                    else:
                      if pdb4pym.has_key(pdb) == False:
                        pdb4pym[pdb]=[[c1,c2,aa1,aa2,res1.get_full_id()[3][1],res2.get_full_id()[3][1], "unsatisfied"]]
                      else:
                        pdb4pym[pdb].append([c1,c2,aa1,aa2,res1.get_full_id()[3][1],res2.get_full_id()[3][1], "unsatisfied"])
                      
                    #pdb_stat[pdb].append(distance)
                    if XL_satisfied.has_key(xl) == False and XL_satisfied.has_key(invxl) == False:
                      XL_satisfied[xl]=[distance]
                    elif XL_satisfied.has_key(xl):
                      XL_satisfied[xl].append(distance)
                    elif xl != invxl and XL_satisfied.has_key(invxl):
                      XL_satisfied[invxl].append(distance)
                    if pdb_stat[pdb].has_key(xl) == False and pdb_stat[pdb].has_key(invxl) == False:
                      pdb_stat[pdb][xl]=[distance]
                    elif pdb_stat[pdb].has_key(xl):
                      pdb_stat[pdb][xl].append(distance)
                    elif xl != invxl and pdb_stat[pdb].has_key(invxl):
                      pdb_stat[pdb][invxl].append(distance)

          else:
            continue
          continue
  
  pdb_stat_sorted=sorted(pdb_stat.items(), key=lambda a: len(a[1]), reverse=True)
  outfile.write("SUMMARY:\tPDB ID\tNr. of total XL \tNr. of total XL mapped to the PDB\tNr. of XL mapped to the specific PDB\tNr. of XL distances satisfied\tPartial Score\tTotal Score\tRepresented Interactions(:Nr of Satisfied XL)\tXL sites mapped to PDB chains\n")
#  for pdb, dist in pdb_stat_sorted:
#    if len(dist) > 0:
#      satisfied, score=XLScore(dist, cutoff)
#      outfile.write("SUMMARY:\t"+pdb+"\t"+str(len(restr_list))+"\t"+str(len(XL_satisfied.keys()))+"\t"+str(len(dist))+"\t"+str(satisfied)+"\t"+str(score)+"\n")
  for pdb in pdb_stat.keys():
    dist=pdb_stat[pdb]
    if len(dist) > 0:
      satisfied, partial_score, total_score, intfs, outstr=XLScore(dist, cutoff, len(restr_list), pac2gn, uniac2xl)
      outfile.write("SUMMARY:\t"+pdb+"\t"+str(len(restr_list))+"\t"+str(len(XL_satisfied.keys()))+"\t"+str(len(dist.keys()))+"\t"+str(satisfied)+"\t"+str(partial_score)+"\t"+str(total_score)+"\t"+outstr+"\t"+str(dist.keys())+"\n")

###Finding here the best combination of template structures: i.e. covering the highest fraction of sequence space and restraints
###In other words, I calculate the Jaccard score of XL sites in pairs of chains and sort the combinations with lowest
  
  elapsed2 = timeit.default_timer() - start_time

  logfile.write("Total XL number: %d\n" % (len(restr_list))) 
  logfile.write("Nr. XL sites matched to the PDB: %d\n" % (len(XL_satisfied.keys()))) 
  xl_sat=0
  for k in XL_satisfied.keys():
    if min(XL_satisfied[k]) < cutoff:
      xl_sat+=1
  logfile.write("Nr. XL sites matched to the PDB and within distance threshold: %d\n" % ((xl_sat))) 
  logfile.write("Time to get fasta sequences and blast them on the PDB = %s\n" % (elapsed1))
  logfile.write("Time for the whole analysis = %s\n" % (elapsed2))

  #print pdb4pym
  if pymol:
    GenPym(pdb4pym, output,inputpdb)

  sys.exit(0)
 
if __name__ == '__main__' : main()


