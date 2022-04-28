"""
__authors__: Alan Moses, Khaled Elemam & Iva Pritisanac
"""
import math
import numpy as np


class SequenceFeatures:

	def __init__(self):# -> None:
		#pass
		self.get_aafeats()
		self.get_seqfeats()

	def get_aafeats(self):
		#self.aafeats = {"net_charge": self.net_charge, "WF_complexity": self.WF_complexity, "KL_hydropathy": self.KL_hydropathy, "isoelectric_point":self.isoelectric_point, "FCR": self.frac_charged_residues, "ED_ratio": self.ED_ratio, "RK_ratio":self.RK_ratio }
		self.aafeats_names = ["net_charge", "WF_complexity", "KL_hydropathy", "isoelectric_point", "FCR", "ED_ratio", "RK_ratio"]
		self.aafeats = [self.net_charge,self.WF_complexity,self.KL_hydropathy,self.isoelectric_point,self.frac_charged_residues, self.ED_ratio,self.RK_ratio]

		#return aafeats

	def get_seqfeats(self):
		#self.seqfeats = {"SCD": self.sequence_charge_decoration, "my_kappa": self.my_kappa, "my_omega": self.my_omega}

		self.seqfeats_names = ["SCD", "kappa*", "omega*","aromatic_spacing*","omega_aromatic*"]
		self.seqfeats = [self.sequence_charge_decoration, self.my_kappa,self.my_omega,self.my_aromspace,self.my_omega_aromatic]
		#return seqfeats


	def net_charge(self, freq): ##computes the net charge from the frequences of residues
		nc=float(0)			## freq should be a dict with amino acid keys and counts/numbers as values
		for aa in list("RK"):
			nc+=freq[aa]
		for aa in list("DE"):
			nc-=freq[aa]
		return nc

	def frac_charged_residues(self,freq):
		nc=float(0)
		for aa in list("DERK"):
			nc+=freq[aa]
		tot=float(0)
		for aa in freq.keys():
			tot += freq[aa]
		try:
			return nc/tot
		except ZeroDivisionError:
			#print(freq)
			return np.nan

	def WF_complexity(self, freq):
		tot=float(0)
		wfc=float(0)
		for aa in freq.keys():
			tot += freq[aa]
			wfc += math.lgamma(freq[aa]+float(1))/math.log(float(20)) #convert to log20 because lgamma is natural
		try:
			return (math.lgamma(tot+float(1))/math.log(float(20)) - wfc)/tot
		except ZeroDivisionError:
			"""
			// CHECK THE BELOW //
			"""
			#eturn math.log(float(20))
			return np.nan

	def KL_hydropathy(self,freq):
		KLscale={
		'I' : 4.5, 'V' : 4.2, 'L': 3.8, 'F' : 2.8, 'C' : 2.5, 'M' : 1.9, 'A' : 1.8, 'G' : -0.4, 'T' : -0.7,
		'S' : -0.8, 'W' : -0.9, 'Y' : -1.3, 'P' : -1.6, 'H' : -3.2, 'E' : -3.5, 'Q' : -3.5, 'D' : -3.5, 'N' : -3.5,
		'K' : -3.9, 'R' : -4.5
		}
		hyd=float(0)
		tot=float(0)
		for aa in freq.keys():
			hyd += freq[aa]*KLscale[aa]
			tot+=freq[aa]
		try:
			return hyd/tot
		except ZeroDivisionError:
			return np.nan

	def isoelectric_point(self, freq)	:  ###based on An algorithm for isoelectric point estimation David L. Tabb
									###and the bisection method for finding 0s, similar to python SeqUtils
		ph=float(0)
		#find the limits
		lowerph=float(0)
		while (self.charge_at_PH(freq,ph)>float(0)):
			lowerph=ph
			ph+=1.0
		#print (str(lowerph) + "is the lower limit on iso point")
		upperph = float(1) + lowerph

		#refine the value using bisection
		while (upperph - lowerph > 0.001):
			#print ("iso point is between "+str(upperph)+" and "+str(lowerph))
			ph = (upperph+ lowerph)/float(2)
			ch=self.charge_at_PH(freq,float(ph))
			if (ch>float(0)):
				lowerph=ph
			else:
				upperph=ph

		return ph


	def charge_at_PH(self, freq,ph):
		ch=float(0)
		pKaapos={ 'K': 10.0, 'R': 12.0, 'H': 5.98, }
		pKaaneg={ 'D': 4.05, 'E': 4.45, 'C': 9.0, 'Y': 10.0  }
		pKN =float(7.5)
		pKC= float(3.55)
		cr = 10**(pKN-ph)
		ch +=  cr/(cr+float(1))
		for aa in pKaapos.keys():
			cr = 10**(pKaapos[aa]-ph)
			ch += freq[aa]*cr/(cr+float(1))
		cr = 10**(ph-pKC)
		ch -=  cr/(cr+float(1))
		for aa in pKaaneg.keys():
			cr = 10**(ph-pKaaneg[aa])
			ch -= freq[aa]*cr/(cr+float(1))
		return ch

	def ED_ratio(self, freq):
		return math.log1p(float(freq['E'])) - math.log1p(float(freq['D']))

	def RK_ratio(self, freq):
		return math.log1p(float(freq['R'])) - math.log1p(float(freq['K']))


	def count_aa(self, seq):
		return { aa:seq.count(aa) for aa in "ACDEFGHIKLMNPQRSTVWY" }


	def sequence_charge_decoration(self, seq):
		scd=float(0)
		charge = {'D':-1,'E':-1, 'R':1,'K':1}
		positions=[]
		for i in range(1,len(seq)):
			if (seq[i] in charge):
					positions.append(i)
		for i in range(1,len(positions)):
			for j in range(0,i):
				scd+= float(charge[seq[positions[i]]]*charge[seq[positions[j]]])*math.sqrt(float(positions[i]-positions[j]))

		try:
			return scd/float(len(seq))
		except ZeroDivisionError:
			return np.nan

	def my_omega(self, seq):	#das and pappu suggest overlapping blobs of 5 or 6
		blob=5 ### my test is based on comparison of spacing between same residue type.
		### if blob gets too small, then we are just capturing the repeats.
		procharge = 'DERKP'
		positions=[]
		obsdis = []

		for i in range (0,len(seq)):
			if (seq[i] in procharge):
					positions.append(i)
		if (len(positions)>1): #need at least two in the sequence
			if ((positions[0]+len(seq)-positions[-1])<=blob): #circularize the sequence
				obsdis.append(positions[0]+len(seq)-positions[-1])
			for j in range(1,len(positions)):
				dis = positions[j]-positions[j-1]
				if (dis<=blob):
					obsdis.append(dis)
			geompar=float(len(positions))/float(len(seq))
			geomprob=float(0) ## compute the probability of observerations in the blob size
			for i in range(0,blob):
				geomprob += geompar*(float(1) - geompar)**i #python's range is not including blob, so don't need the -1
			#print (len(obsdis)),;print("observed within " + str(blob))
			#print(obsdis)
			#print geomprob*float(len(positions)),;print(" of "+str(len(positions))+" expected")
			##use normal approx. to binomial
			if ((float(1)-geomprob)*geomprob > float(0)):
				return (float(len(obsdis)) - geomprob*float(len(positions)))/math.sqrt((float(1)-geomprob)*geomprob*float(len(positions)))
			else:
				#print("problem with " + str(geomprob) + " as a binomial distribution parameter")
				return(np.nan)
		else:
			return(np.nan)


	def my_kappa(self, seq): #das and pappu suggest overlapping blobs of 5 or 6
		blob=5 ### my test is based on comparison of spacing between same charges.
		### if blob gets too small, then we are just capturing the repeats.
		charge = {'D':-1,'E':-1, 'R':1,'K':1}
		positions=[]
		tot= {1:0 , -1:0}
		dsame=[]
		for i in range (0,len(seq)):
			if (seq[i] in charge):
				if (len(positions)>0):
					dis = i - positions[-1]
					if (dis<=blob):
						if (charge[seq[i]] is charge[seq[positions[-1]]]):
							dsame.append(dis)
				positions.append(i)
				tot[charge[seq[i]]] += 1

		if (len(positions)>1):
			geompr = {  }
			sumpr=float(0)
			sumsq=float(0)
			for c in tot.keys():
				geompr[c] = float(tot[c])/float(len(seq))
				sumpr += geompr[c]
			for c in tot.keys():
				sumsq += geompr[c]*geompr[c]/sumpr
			prsame=float(0) ## use the bivariate geometric distribution to compute the probability of
						##finding two of the same within the a blob
			for d in range(0,blob):			#python's range is not including blob, so don't need the -1
				prsame += sumsq*(float(1)-sumpr)**d
			#print (len(dsame)),;print("observed within " + str(blob))
			#print(dsame)
			#print(obsdis)
			#print prsame*float(len(positions)),;print(" of "+str(len(positions))+" expected")
			##use normal approx. to binomial
			if ((float(1)-prsame)*prsame >0):
				return (float(len(dsame)) - prsame*float(len(positions)))/math.sqrt((float(1)-prsame)*prsame*float(len(positions)))
			else:
				#print("problem with " + str(prsame) + " as a binomial distribution parameter")
				return(np.nan)
		else:
			return(np.nan)

	def my_omega_aromatic(self,seq):	#das and pappu suggest overlapping blobs of 5 or 6
		blob=5 ### my test is based on comparison of spacing between same residue type.
		 ### if blob gets too small, then we are just capturing the repeats.
		proarom = 'FYW'
		positions=[]
		obsdis = []

		for i in range (0,len(seq)):
			if (seq[i] in proarom):
					positions.append(i)
		if (len(positions)>1): #need at least two in the sequence
			if ((positions[0]+len(seq)-positions[-1])<=blob): #circularize the sequence
				obsdis.append(positions[0]+len(seq)-positions[-1])
			for j in range(1,len(positions)):
				dis = positions[j]-positions[j-1]
				if (dis<=blob):
					obsdis.append(dis)
			geompar=float(len(positions))/float(len(seq))
			geomprob=float(0) ## compute the probability of observerations in the blob size
			for i in range(0,blob):
				geomprob += geompar*(float(1) - geompar)**i #python's range is not including blob, so don't need the -1
			#print (len(obsdis)),;print("observed within " + str(blob))
			#print(obsdis)
			#print geomprob*float(len(positions)),;print(" of "+str(len(positions))+" expected")
			##use normal approx. to binomial
			if ((float(1)-geomprob)*geomprob > float(0)):
				return (float(len(obsdis)) - geomprob*float(len(positions)))/math.sqrt((float(1)-geomprob)*geomprob*float(len(positions)))
			else:
				#print("problem with " + str(geomprob) + " as a binomial distribution parameter")
				return(np.nan)
		else:
			return(np.nan)

	# version of my kappa for spacing of aromatic residues Y, F, W
	def my_aromspace(self,seq): #das and pappu suggest overlapping blobs of 5 or 6
		blob=5 ### my test is based on comparison of spacing between same charges.
		 ### if blob gets too small, then we are just capturing the repeats.
		aromatic = {'Y':1,'F':1, 'W':1}
		positions=[]
		tot= {1:0}
		dsame=[]
		for i in range (0,len(seq)):
			if (seq[i] in aromatic):
				if (len(positions)>0):
					dis = i - positions[-1]
					if (dis<=blob):
						if (aromatic[seq[i]] == aromatic[seq[positions[-1]]]):
							dsame.append(dis)
				positions.append(i)
				tot[aromatic[seq[i]]] += 1

		if (len(positions)>1):
			geompr = {  }
			sumpr=float(0)
			sumsq=float(0)
			for c in list(tot.keys()):
				geompr[c] = float(tot[c])/float(len(seq))
				sumpr += geompr[c]
			for c in list(tot.keys()):
				sumsq += geompr[c]*geompr[c]/sumpr
			prsame=float(0) ## use the bivariate geometric distribution to compute the probability of
						##finding two of the same within the a blob
			for d in range(0,blob):			#python's range is not including blob, so don't need the -1
				prsame += sumsq*(float(1)-sumpr)**d
			#print (len(dsame)),;print("observed within " + str(blob))
			#print(dsame)
			#print(obsdis)
			#print prsame*float(len(positions)),;print(" of "+str(len(positions))+" expected")
			##use normal approx. to binomial
			if ((float(1)-prsame)*prsame >0):
				return (float(len(dsame)) - prsame*float(len(positions)))/math.sqrt((float(1)-prsame)*prsame*float(len(positions)))
			else:
				#print("problem with " + str(prsame) + " as a binomial distribution parameter")
				return(np.nan)
		else:
			return(np.nan)
