"""
__author__: Iva Pritisanac

Input files loading and manipulation; relies on class Parser
"""
import os,sys
from myparser import Parser

## wrapper class that inherits all the methods of a class Parser
#the class also inherits the constructor of the Parser
## Parser allows the extraction of keywords from an input text file
#input text file should contain a keyword per line, followed by one or more parameters (separated from the keyword by a tab)
class ParseInput(Parser):

    def __init__(self):
        super(ParseInput,self).__init__()
        self.add('motifs_file','motifs_file','str',"")
        #self.add('exp_motifs_n_file','exp_motifs_file','str',"")
        self.add('repeats_file','repeats_file','str',"")
        self.add('aa_freq_file','aa_freq_file','str',"")
        self.add('align_dir','aln_dir','str',"") # path to directory with alignment files

        #self.add('REF_NUM','ref_n','int',0) # #the position of the reference species in the file
        #self.add('REF_NAME','ref_name','str',"")

    def check_variables(self):

        if os.path.isfile(self.motifs_file) == False:

            raise Exception("Motifs file does not exist or the path to the file is incorrect!")
            sys.exit(0)

        #if os.path.isfile(self.exp_motifs_file) == False:

            #raise Exception("Expected_n_motifs file does not exist or the path to the file is incorrect!")
            #sys.exit(0)

        if os.path.isfile(self.repeats_file) == False:
            raise Exception("Repeats file does not exist or the path to the file is incorrect!")
            sys.exit(0)

        #if os.path.isfile(self.aa_freq_file) == False:
        #    raise Exception("AA frequency file does not exist or the path to the file is incorrect!")
        #    sys.exit(0)

        if not os.path.isdir(self.aln_dir):
            raise Exception("Alignment directory does not exist or the path to the directory is incorrect!")
            sys.exit(0)

    def read_in_files(self):
        #self.read_aa_frequency(self.aa_freq_file) # sets aa_freq object - a dictionary that stores precomputed aa frequencies of long human IDRs
        self.read_repeats_from_file(self.repeats_file) # sets R object - a dictionary that stores repeats and their regex
        self.read_motifs_from_file(self.motifs_file) # sets M object - a dictionary that stores motives and their regex
        #self.get_expected_n_matches(self.exp_motifs_file) # sets self.motifs_expect - a dictionary that stores motives and their expected frequency

    def read_repeats_from_file(self,infile):
        self.repeats = {}

        f = open(infile,'r') #open the file
        for line in f:
            line=line.rstrip("\n\r")
            if ('id=' in line):
                mname=line.split("\t")[0].split("id=")[1]
                mregex=line.split("\t")[1].replace("m=","")
                self.repeats.setdefault(mname,mregex.replace("'","").replace('"',""))

        f.close()
        #mname=""
        #for line in f:
        #    line=line.rstrip("\n\r")
        #    if ('id=' in line):
        #        mname=line.replace("\t"," ")
        #        self.repeats[mname] = [] # add this header line to the list
        #    else:
        #        vals=line.split("\t")
        #        self.repeats[mname].append(vals[1]) ## these are fixed width motifs

    def read_aa_frequency(self,infile="PARAM_FILES"+os.sep+"AA_COMPOSITION.txt"):
        self.aa_freq={}
        f = open(infile,'r') #open the file
        for line in f:
            line=line.rstrip("\n\r")
            splitted=line.split("\t")
            self.aa_freq.setdefault(splitted[0],float(splitted[1]))

        #return aa_freq

    #IP added variable length motifs reading
    #used to be called read_f_motifs_from_file
    def read_motifs_from_file(self,infile):
        self.motifs = {}
        f = open(infile,'r') #open the file
        for line in f:
            line=line.rstrip("\n\r")
            if ('id=' in line):
                mname=line.split("\t")[0].split("id=")[1]
                mregex=line.split("\t")[1].replace("m=","")
                self.motifs.setdefault(mname,mregex.replace("'",""))
        f.close()
        #return M


    def get_expected_n_matches(self,expfile):
        self.motifs_expect = {}
        f = open(expfile, 'r') #open the file
        for line in f:
            line=line.rstrip("\n\r")
            if len(line)==0:
                continue
            mname=line.split("\t")[0]
            mexp=float(line.split("\t")[1])
            self.motifs_expect.setdefault(mname,mexp)
        #return motifs_expect
