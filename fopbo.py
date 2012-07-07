#!/usr/bin/env python
#
#       fopBO.py
#       Coding by Raul Mera-Adasme.
#		Chemistry by Raul Mera-Adasme, Fernando Mendizabal
#		Sebastian Miranda-Rojas, Claudio Olea-Azar and Patricio Fuentealba.
#
#       Copyright 2009 Raul Mera-Adasme <rmeraa[at]gmail[dot]com>
#       ESR/NMR Free Radicals and Computational Chemistry Group,
#       Theoretical Inorganic Chemistry Lab, Universidad de Chile        
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#
##
##Dedicated to the long life of the Ven. Khempo Phunstok Tenzin Rinpoche
##
#usage: fopbo.py filename atindex1 atindex2 uhf/rhf [verbosiness] 
#
#
#This uses a default (i.e. no need for any keyword in the FILE.47) NBO output 
#as an input, so you need the NBO program by Weinhold et al. (http://www.chem.wisc.edu/~nbo5/)
#


import sys, copy
from numpy.oldnumeric import array
from numpy.oldnumeric.mlab import std, mean


deloc_dict=dict(id=0,E2=0.0,F=0.0,deltaE=0.0,coef=0.0,donorb=0,aceptorb=0,donor=[],aceptor=[])


#gets a preeliminary bond order consisting in the bonding-antibonding NBOs 
def nbo_bond_order(fin,pair,spin):
	bond_order=0
	ids=[]
	while True:
		i=fin.readline()
		if " NHO DIRECTIONALITY AND BOND BENDING (deviations from line of nuclear centers)" in i: #this is where the bond listing section ends.
			break 
		if i[0:4].replace(" ","").replace("-","").replace(".","").isdigit() and len(i)>=28:
			if "BD" in i and int(i[25:28]) in pair and int(i[31:34]) in pair:
				ids.append(int(i[:4]))
				bond_ponderation=1 #change to include only covalent contributions
				if "BD*" in i:
					continue
				if spin=="rhf":
					bond_order=bond_order+1.0
				else:
					bond_order=bond_order+0.5
	return bond_order,ids

	
#get a correction of the normal NBO single-structure bond order based on the
#first order perturbation theory calculated donor-aceptor interactions. 
#main function of the program.
def deloc_bond_order(fin,pair,ids,spinstate,verbosity):
	if spinstate=="uhf":
		ponderation=0.5
	else:
		ponderation=1
	delocs=[] #every delocalization acepted in principle, goes here.
	special_delocs=[]
	#########
	#In this big loop  we just set aside the orbital interaction which will be taken into 
	#Account int eh next block. We avoid interaction involving previously used orbitals 
	#(in ids) consider only delocalizataions without participation of already analyzed NBOs and
	#with participation of each atom. discard delocalizations in which te same atom is in the 
	#aceptor and donor grups.
	deloc_id=0 #an id for each saved delocalization
	while True:
		antibond=1 #-1 for antibonding orbital, 1 for bonding orbital
		aindex=[[],[]] #atoms in the bond data.
		i=fin.readline()
		if "RAL BOND ORBITALS (Summary)" in i:
			break
		if len(i)<62:
			continue
		if (not i[55:62].replace(" ","").replace("-","").replace(".","").isdigit()): #skip if the delocalization energy is below thres_deloc
			continue
		#this is an ugly, Q&D fix to read atoms > 100 and < 110. I need a nicer fix.
		while  "A" in i: 
			place=i.index("A")
			i=i[:place]+"10"+i[place+1]+i[place+3:]
		if i[6:8]=="BD":
			aindex[0].append(int(i[15:17]))  #the ids of the donating atoms.
			aindex[0].append(int(i[20:22]))
		else:
			aindex[0].append(int(i[15:17]))  #if is not BD or BD*, is a lone pair of core, so it only have one atom id
		if i[33:35]=="BD": #the receiving bond
			aindex[1].append(int(i[42:44])) #the ids of the receiving atoms
			aindex[1].append(int(i[47:49]))
			if i[35]=="*" and aindex[1][0] in pair and aindex[1][1] in pair: #value of antibond allows to substract contributions donated to antibonding orbitals, so is important that the atoms in the antibond belongs to the pair.
				antibond=-1	
		else:
			aindex[1].append(int(i[42:44]))
		if not aindex[0] or not aindex[1]:
			print "empty index"
			continue	
		if (aindex[0]==aindex[1] or aindex[1].reverse()==aindex[0]): #Probably not needed since a donatio from an orbital to the same would have no effect, but is nicer this way.
			continue
		delocs.append(copy.deepcopy(deloc_dict))
		delocs[-1]["E2"]=float(i[54:63])
		delocs[-1]["deltaE"]=float(i[63:70])
		delocs[-1]["F"]=float(i[71:79])*antibond
		delocs[-1]["donorb"]=int(i[0:4])
		delocs[-1]["aceptorb"]=int(i[27:31])
		delocs[-1]["donor"]=aindex[0]
		delocs[-1]["aceptor"]=aindex[1]
		delocs[-1]["id"]=deloc_id
		if len(delocs[-1].keys())!=9:
			raise KeyError # Something would be wrong with the code.
		deloc_id=deloc_id+1
	####
	####  This is the fun part, where we use the data collected upthere to get the delocalization
	####  contriburions to bond order.
	####

	#the 2 previous should go with one more identation level andunder the next for bucle, thus saving the additional cleaning of these values later (lins 207 and 208).
	bond_contribs=[] # contributions to bonding. the sum of this is the return value
	for d in range(len(pair)):  #first 0 donor 1 aceptor, then reversed.
		comp_delocs=[] #competing delocalizataions with those of current pair will be used twice, once for each member of the pair acting as donor.
		bond_delocs=[] #delocalizations between the current pair in the current donor-aceptor order
		a=d-1 #d for donor a for aceptor. if d=0 a=-1, if d=1 a=0
		if verbosity in (1,3):
			print "donor:", pair[d], "aceptor:", pair[a]
		for i in delocs:
			if pair[a] in i["aceptor"] and pair[d] in i["donor"]:
				if  not (i["donorb"] in ids and not i["aceptorb"] in ids): #condition added 20/10/10 to prevent donations from bonds to be counted as contributions
					bond_delocs.append(copy.deepcopy(i))
			elif pair[d] in i["donor"]:
				comp_delocs.append(copy.deepcopy(i))
		if len(bond_delocs)<1: ####no donor aceptor interaction from this pair in this order
			continue
		for j in bond_delocs:
			comp_delocs_orb=[] ##competing delocs for current_orbital
			comp_coefs_sq_sum=0 #sum of squares of all elements in previous list
			coef=(1.0*j["F"])/j["deltaE"]   ##float assured
			for k in comp_delocs:  #competing interaction with other atoms
				if k["donorb"]==j["donorb"] and k["id"]!=j["id"]:  #looking for delocalizations with the same donor, i.e. competing.
					c=(1.0*k["F"])/k["deltaE"] 
					comp_coefs_sq_sum=comp_coefs_sq_sum+(c**2)
			for k in bond_delocs: #competing interactions with different orbitals from the same atom
				if k["donorb"]==j["donorb"] and k["id"]!=j["id"]:
					c=(1.0*k["F"])/k["deltaE"]   ##float assured
					comp_coefs_sq_sum=comp_coefs_sq_sum+(c**2)		
			den=(1+(coef**2)+comp_coefs_sq_sum)**(0.5)
			coef=coef/den
			if coef<0:  #means donation to antibonding. These are handled later.
				coef=0
			coef=(coef**2)
			coef=coef*ponderation #ponderation is 1 for restricted, 0.5 for UHF
			bond_contribs.append(coef)
			if verbosity in (1,3):
				print "coef:",coef,"den:",den,"deltaE",j["deltaE"],"E2",j["E2"],"F",j["F"],"donorb",j["donorb"],"aceptorb",j["aceptorb"],"donor",j["donor"],"aceptor",j["aceptor"], "unperturbed:", (1.0/den)**2 
	additional_interactions=0.0
	for j in bond_contribs:
		additional_interactions=additional_interactions+j    #just add all the coefs
	if verbosity in (1,3):
		print "ids corrections"  #donations received on antibonding orbtitals and donations from bonding orbitals from the interaction in study. i.e. everything that decreases the bond order
	destabilization_count=0
	for j in ids: #in ids we have bonds and antibonds. Donation to antibondings were skipped previously, are handled now.
		bond_delocs=[] #here these 2 lists are declared in a good place so there is not need for cleaning their values
		comp_delocs=[]
		for k in delocs:
			comp_coefs_sq_sum=0
			sign=1
			if k["donorb"]==j or k["aceptorb"]==j:
				for l in delocs:
					if l["donorb"]==k["donorb"] and l["id"]!=k["id"]:
						c=(1.0*l["F"])/l["deltaE"]
						comp_coefs_sq_sum=comp_coefs_sq_sum+(c**2)
				if k["donorb"]==j:
					sign=-1
				else:
					if k["F"]<0:
						sign=-1
						if (verbosity in (1,3)):  #donations to antobonding were tagged previously adding a minus sign to the fock matrix value.
							print "donation to an antibonding"
					else:
						continue  # Donations to orbitals already full at unperturbed level are not accepted.
				coef=(1.0*k["F"])/k["deltaE"]
				den=(1+(coef**2)+comp_coefs_sq_sum)**(0.5)
				coef=coef/den
				coef=sign*(coef**2) 
				coef=coef *ponderation 
				if verbosity in (1,3):
					print "coef:",coef,"den:",den,"deltaE",k["deltaE"],"E2",k["E2"],"F",k["F"],"donorb",k["donorb"],"aceptorb",k["aceptorb"],"donor",k["donor"],"aceptor",k["aceptor"], "contribution of unperturbed:", (1.0/den)**2 ####last part changed 20/10/10##################33	
				#bond_contribs.append(coef)
				destabilization_count=destabilization_count+coef
	total_border=0.0 #total bond order, sum of bond_contribs. The return value
#	for i in bond_contribs:
#		total_border=total_border+i
	total_border=additional_interactions+destabilization_count
	if verbosity in (2,3):
		print "additional contributions to bonding:", additional_interactions
		print "destabilizations:                   ", destabilization_count
	return total_border



				
#
#just an interface function to read both spins.
#
def pair_bond(filename,pair,spinstate,verbosity):
	data=[] #data for one atom. Each element is a list of 2 elements. A float for bond order and a list of IDs of NBO IDs which  
	#are bonding or antibonding the 2 atoms
	fin=open(filename,"r")
	spin="t"
	while True: 
		i=fin.readline()
		if len(i)==0:
			break
		if "****         Alpha spin orbitals         ****" in i:
			spin="a"
		if "****         Beta  spin orbitals         ****" in i:
			spin="b"
		if "(Occupancy)   Bond orbital/ Coefficients/ Hybrids" in i:
			bond_order,ids=nbo_bond_order(fin,pair,spinstate)   #list of at most 3 lists
			data.append([bond_order,ids])
			if verbosity in (2,3):
				print "Unperturbed Bond Order", bond_order
		if "SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS" in i:
			data[-1][0]=data[-1][0]+deloc_bond_order(fin,copy.copy(pair),copy.copy(data[-1][1]),spinstate,verbosity)
	print "Bond Orders"	
	for j in data:
		print j[0], #bond order in each spin
	if spinstate=="uhf":
		print data[0][0]+data[1][0] #total bond order







#usage: program filename atom1 atom2 spinmode verbosity
#
filename=sys.argv[1]
atoms=[]

if len(sys.argv)==5:
	sys.argv.append(0) 

pair_bond(filename,[int(sys.argv[2]),int(sys.argv[3])],sys.argv[4],int(sys.argv[5]))
	
	
