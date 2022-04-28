"""
__author__: Matteo Degiacomi
input file parser. Extracted from BioBOx, (c) Matteo Degiacomi
"""

import sys
import numpy as np

##Parser allowing the extraction of keywords from an input file.
#input text file should contain a keyword per line, followed by one or more parameters.
class Parser(object):

    ##initialize parameters and assign them their default values
    def __init__(self):

        ##dictionary storing parameter definitions (key=parameter_name, value=[variable_name,variable_type,default_value]).
        self.parameters={}
        self._add_standard()

    ##insert a new keyword entry in parameters dictionary.
    #@param key name of keyword in input file
    #@param variable name of variable within the code
    #@param vartype variable type. Accepts int, str, float, array int, array str, array float, dictionary
    #@param default default value
    def add(self,key,variable,vartype,default):
        self.parameters[key]=[variable,vartype,default]

    ##parse input file and replace default values.
    #param name of file to parse
    def parse(self,infile):

        self._set_default_values()

        f = open(infile,'r+')
        line = f.readline()
        while line:
            #print(line)
            w = line.split()

            if len(w) > 0 and str(w[0][0])!='#':

                #val=[variable_name,variable_type,default_value]
                try:
                    val=self.parameters[w[0]] # get the default parameter of the value
                except KeyError:
                    print("unrecognised keyword %s"%w[0])
                    sys.exit(1)

                #if type is string
                if val[1].split()[0]=='str': # val[1] = 'str' or 'int'
                    exec('self.%s=%s("%s")'%(val[0],val[1],w[1]))

                #if type is int or float
                elif val[1].split()[0]=='int' or val[1].split()[0]=='float':
                    exec('self.%s=%s(%s)'%(val[0],val[1],w[1]))

                #if type is an array of int, float, or str
                elif val[1].split()[0]=='array' and (val[1].split()[1]=='int' or val[1].split()[1]=='float' or val[1].split()[1]=='str'):
                    exec('self.%s=np.array(%s).astype(%s)'%(val[0],w[1:len(w)],val[1].split()[1]))

                #if type is dictionary
                elif val[1].split()[0]=="dictionary":
                    exec('self.%s[\"%s\"]=%s'%(val[0],w[1],w[2:]))

                else:
                    print("unrecognised type for keyword %s: %s"%(w[0],val[1]))
                    sys.exit(1)

            line = f.readline()

        f.close()


    def _add_standard(self):

        #keyword_name, variable_name, type, default value
        self.add('monomer','monomer','str',"")
        self.add('polyhedron','polyhedron','str',"")

    #set default values for all defined keywords Xx Parsing the self paramater and creating a self of all the values xX
    def _set_default_values(self):
        #print(self.parameters)

        for k,v in self.parameters.items():
            exec('self.%s=v[2]'%v[0]) # -> this is where the self.style comes fro
