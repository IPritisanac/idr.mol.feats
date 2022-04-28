import io
import os,sys
import math
import random
import re
import time
import numpy as np # IP 14/9/2020
import scipy.stats
from data import ParseInput
from sequence_features import SequenceFeatures

class Alignment():

	def __init__(self):
		self.names = [] #list of names of sequences in the alignment
		self.seq = {} #the actual sequence indexed by their names
		self.info = "" #any other information

	def read_mfa(self,file): #read from a multiple fasta file
		f = open(file, 'r') #open the file
		for line in f:
			line=line.rstrip("\n\r")
			if ('>' in line):
				self.names.append(line) # add this header line to the list
			else:
				seqname=self.names[-1] # this is the most recently found header line
				if (not seqname in self.seq):
					self.seq[seqname] = line
				else:
					self.seq[seqname] += line

		f.close()

	def print_mfa(self): #print fasta format to the screen
		for s in self.names:
			if ('>' in s):
				 print(s)
			else:
				print ('>' + s)
			print(self.seq[s])

	def print_maf(self): #print maf format to the screen
						#assumes the score line has been stored as the info for the alignment
		print(self.info),;
		for s in self.names:
			print ("s "),;print(s),;print(self.seq[s])
		print ('\n')


class MolFeats():
    def __init__(self,inputf):
        self.dirpath=os.path.dirname(os.path.abspath(inputf)) # extract directory name out of the input filename
        if not os.path.exists(self.dirpath+os.sep+"RES"): # if the output results directory does not exist on this path
            os.makedirs(self.dirpath+os.sep+"RES/")   # create the directory
        self.outdir=self.dirpath+os.sep+"RES/"   # set the output directory
        # create instance of the parser class to parse the input file
        self.aastring="ACDEFGHIKLMNPQRSTVWY"
        self.parse_param_files(inputf)
        self.get_features()
        # this should get inherited
        #self.get_functions()
        # //TO DO// Inherit from DATA and from SEQ FEAT - make this a super class
        #- ensure input file data structures are all generated
        # self.repeats, self.motifs, self.motifs_EXP, results_file, sequences_file, run_indels, featnamesorted, maxsize,indelprob

    ## run parser and return only file paths of various parameter files
    def parse_param_files(self,input_file):
        P=ParseInput()
        try:
            P.parse(input_file)
        except Exception as e:
            print("ERROR: Could not parse input file:\n%s"%e)
            sys.exit(1)
        try:
            P.check_variables()
        except Exception as e:
            print("ERROR in input variables!\n")
            print(e)
            sys.exit(1)

        P.read_in_files()
        self.repeats=P.repeats
        #print(self.repeats)
        self.motifs=P.motifs
        self.aln_dir=P.aln_dir
        #self.aa_freq=P.aa_freq

    def get_features(self,):
        SF=SequenceFeatures() # //TO DO// Make a super constructor that inherits INIT of seq features
        self.aafeats=SF.aafeats
        self.aafeats_names=SF.aafeats_names
        self.seqfeats=SF.seqfeats
        self.seqfeats_names=SF.seqfeats_names

    def seq_choosing_heuristic(self,aln): ##try to get dis_tot, but don't use "redundant" sequences
        seqlist = []
        seqnames = []
        for s in aln.names:
            ugseq=aln.seq[s].replace("-","")
            if ('X' in ugseq):
                continue #no bad data
            if len(ugseq)==0:
                continue
            seqlist.append(ugseq) # store the actual sequence
            seqnames.append(s) # store the name associated with sequence

        return seqlist,seqnames

    def write_out_features(self,featdict,outpath,split_each=False):
        if split_each:
            keys=[key for key in featdict.keys()]
            first_entry=keys[0]
            output=outpath.replace(".txt",'').replace('.out','')
            covered=[]
            for seqid,feats in featdict.items():
                uniprot=seqid.split("_")[0]
                if uniprot not in covered:
                    covered.append(uniprot)
                    fout=open(output+uniprot[1:]+".out.txt","a")
                    for entry in featdict[first_entry]:
                        fout.write("\t%s"%(entry))

            for seqid,feats in featdict.items():
                uniprot=seqid.split("_")[0]
                fout=open(output+uniprot[1:]+".out.txt","a")
                fout.write("\n%s"%(seqid))
                for featid,featval in feats.items():
                    fout.write("\t%s"%(featval))

        else:
            fout=open(outpath,"w") # open a file on path
            #if not os.path.isfile(outfile):
            keys=[key for key in featdict.keys()]
            first_entry=keys[0]

            for entry in featdict[first_entry]:
                fout.write("\t%s"%(entry))
            for seqid,feats in featdict.items():
                fout.write("\n%s"%(seqid))
                for featid,featval in feats.items():
                    fout.write("\t%s"%(featval))

        fout.close()

    def compute_feats_dir(self):
        featnamesorted = (self.aafeats_names + self.seqfeats_names +  list(self.repeats.keys()) + list(self.motifs.keys()))
        for mfa in os.listdir(self.aln_dir):
            #start new data structures for each fasta file
            mol_feats = {}
            obs={}
            pred={}
            Z = {}
            if (not mfa.endswith(".fasta")) and (not mfa.endswith(".fa")):
                continue

            output_file=self.outdir+os.sep+mfa[:-4]+"_FEAT.out.txt"
            ALN = Alignment()
            ALN.read_mfa(self.aln_dir+os.sep+mfa)
            print(str(len(ALN.names))+' sequences read in '+mfa);
            
            for f in featnamesorted:
                obs[f]=[]

            sequences_to_use,names_to_use =self.seq_choosing_heuristic(ALN) ### this should do all the quality control
            print(str(len(sequences_to_use))+' sequences passed filtering');

            for x in range(len(sequences_to_use)):
                ugseq=sequences_to_use[x]
                nseq=names_to_use[x]
                mol_feats.setdefault(nseq,{})

                for r in list(self.repeats.keys()):
                    pat = re.compile(self.repeats[r])
                    all_m=[] #IP added 1/11/2019
                    repeatlength=0
                    for m in pat.finditer(ugseq):
                        all_m.append(m)
                        repeatlength += (len(m.group(0)) - 1) #first residue in the repeat isn't counted
                    obs[r].append(float(repeatlength))
                    mol_feats[nseq][r]=float(repeatlength)

                for m in list(self.motifs.keys()):
                    pat = re.compile("".join(self.motifs[m]))
                    mol_feats[nseq][m]=float(len(pat.findall(ugseq)))
                    obs[m].append(float(len(pat.findall(ugseq))))

                obsaa = {aa : ugseq.count(aa) for aa in self.aastring}

                # iterates over a set of functions stored in a list aafeats
                #aafeats[i](obsaa) - a call to a function under the index i
                #obsaa is the argument passed to each function
                for i in range(0,len(self.aafeats)):
                    mol_feats[nseq][self.aafeats_names[i]]=self.aafeats[i](obsaa)
                    if self.aafeats[i](obsaa)==None:
                        print("SEQ: %s"%ugseq)
                    obs[self.aafeats_names[i]].append(self.aafeats[i](obsaa))
                for i in range(0,len(self.seqfeats)):
                    mol_feats[nseq][self.seqfeats_names[i]]=self.seqfeats[i](ugseq)
                ## IF SPLITTING EACH OUTPUT BY UNIPROT ID, add SPLIT=True

            self.write_out_features(mol_feats,output_file,split_each=False)
            print("Computed molecular features for %s"%(mfa))
            print("Written molecular features out to >> %s"%(output_file))
            print("Exiting ...")
        sys.exit(0)

if __name__=="__main__":
    input_file=sys.argv[1]
    MF=MolFeats(input_file)
    MF.compute_feats_dir()
