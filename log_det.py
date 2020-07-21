import numpy as np
import math
import fileinput
import sys

DNA_DIST_MAX = 5.0


def SpeciesFileReader(folder):
	species = []

	species_file = open(folder+'/species_names.txt', 'r') 
	for spec in species_file: 
		species.append(spec.strip())
	return species


def GTRFileReader(folder):
	GTRfile = open(folder+'/DistMatrix.txt','r')
	count  = 0
	lines = []
	for line in GTRfile:	#fileinput.FileInput("DistMatrix.txt",inplace=1,backup='.bak'):
		if line.rstrip():
			count += 1
			if count %5 ==1 :
				continue
			else:
				lines.append(line)
	return lines

def Lines_To_Matrix(lines):
	b = []
	a = []
	count = 0
	for line in lines:
		a.append([float(x) for x in line.split()])
		count += 1
		
		if count %4 == 0:
			b.append(a)
			a = []
	return b
	

def stat_reader(folder):
	statdict = {}
	lines =[]
	for line in fileinput.FileInput(folder+"/stats.txt"):
		lines.append(line)
		#print(line)
	j = 0
	for i in range(int(len(lines)/2)):
		statdict[lines[j].split('\n')[0]] = lines[j+1].split('\n')[0]
		j = j+2
	return statdict

def Pi_array_builder(statdict):
	PiDict = {}
	for taxa in statdict:
		counts_taxa = statdict[taxa].split()
		counts_taxa = np.array([int(i) for i in counts_taxa])
		summation = sum(counts_taxa[0:4]) 
		PiDict[taxa] = counts_taxa[0:4]/summation
	return PiDict

#/************************* FASTME DET METHOD ********************************/
def nextPerm(p, index, size, length):
	temp = 0
	if 0 ==  (index % np.math.factorial(size)) :
		temp = p[length-1];
		p[length-1] = p[length-1-size];
		p[length-1-size] = temp;
		return p
	else:
		return nextPerm(p, index, size-1, length)

def permDiagProduct(P, p, d):
	prod = 1.0

	for i in range(d):
		prod = prod * P[i][p[i]]

	return prod

def det(P, d):
	p = []
	signum = 1;
	det = 0;

	p = [i for i in range(d)]
	numPerms = np.math.factorial(d)

	for i in range(numPerms):
		det += signum * permDiagProduct (P, p, d)
		p = nextPerm (p, i+1, d-1, d)
		signum = -1 * signum
	return det
	
#/*********************************************************/

def Dist_Matrix_Builder(nospecies, species, b, statdict, PiDict):
	DistMatrix = []
	for i in range(nospecies):
		DistMatrix.append([0 for j in range(nospecies)])
	
	i = 0
	j = 1
	matrixno = 0
	count = 0
	
	for matrix in b:
		if j ==nospecies:
			i = i+1
			j = i+1	
			
		######    x_ii clc    #############
		taxa1 = species[i]
		taxa2 = species[j]
		counts_of_taxa1 = statdict[taxa1].split()
		counts_of_taxa2 = statdict[taxa2].split()
		no_base = []
		
		for c in range(4):
			no_base.append(int(counts_of_taxa1[c]) + int(counts_of_taxa2[c]))  # total no of A, C, G, T count
		seqlen = max(int(counts_of_taxa1[4]), int(counts_of_taxa2[4]))
		
		matrix[0][0] = ((no_base[0] - (matrix[0][1] + matrix[0][2] + matrix[0][3])*seqlen)/(seqlen*2)) # ((#A's-(x_AC+x_AG+x_AT)*seqlen)/(2))/seqlen
		matrix[1][1] = ((no_base[1] - (matrix[1][0] + matrix[1][2] + matrix[1][3])*seqlen)/(seqlen*2)) # ((#C's-(x_CA+x_CG+x_CT)*seqlen)/(2))/seqlen
		matrix[2][2] = ((no_base[2] - (matrix[2][0] + matrix[2][1] + matrix[2][3])*seqlen)/(seqlen*2)) # ((#G's-(x_GA+x_GC+x_GT)*seqlen)/(2))/seqlen
		matrix[3][3] = ((no_base[3] - (matrix[3][0] + matrix[3][1] + matrix[3][2])*seqlen)/(seqlen*2)) # ((#T's-(x_TA+x_TC+x_TG)*seqlen)/(2))/seqlen
		###
		
		
		Pi1 = PiDict[taxa1]
		Pi2 = PiDict[taxa2]
		
		print(taxa1,taxa2)
		
		######### Normalization of F matrix #####################
		summation = 0.0
		for k in range(4):
			for l in range(4):
				summation += matrix[k][l]

		for k in range(4):
			for l in range(4):
				if summation != 0.0:
			 		matrix[k][l] = matrix[k][l]/summation
			 		
		detP = det(matrix,4)
		
		print("Determinant : ",detP)
		
		if detP < 0.0:
			DistMatrix[i][j] = DNA_DIST_MAX
		else:
			DistMatrix[i][j] = -0.25 * np.log(detP)					
			print("Log_det : ",DistMatrix[i][j])
			for base in range(4):
				if Pi1[base] <= 0 or Pi2[base] <= 0:
					exit()
				print(np.log(Pi1[base]), np.log(Pi2[base]))
				DistMatrix[i][j] += (np.log(Pi1[base])+np.log(Pi2[base]))/16.0
		print("Final Log_det : ",DistMatrix[i][j])
		DistMatrix[j][i] = DistMatrix[i][j]
		count +=1
		j = j+1
		matrixno +=1 
	return DistMatrix
	
	
def FileWriter(DistMatrix, nospecies, species,folder):

	writefile = open(folder+"/Dist_Matrix.txt","w+")

	writefile.write("%d\n" %(nospecies))
	for i in range(nospecies):
		writefile.write("%s " %species[i])
		for j in range(nospecies):
			writefile.write("%g " %(DistMatrix[i][j]))
		writefile.write("\n")
	writefile.close()	

def main():
	folder = sys.argv[1]
	#print(folder)
	species = SpeciesFileReader(folder) #species_names.txt file, where the species names are written in the order of Skmer ref_dir folder.
	nospecies = len(species)
	
	lines = GTRFileReader(folder)		#F matrix file, generated from Our GTR Approach. n * (n-1) matrices.
	b = Lines_To_Matrix(lines)			
	statdict = stat_reader(folder)		#Statistics file, where A,T,C,G counts are saved.
	
	PiDict = Pi_array_builder(statdict) #Pi Arrays, similiar to FastME.
	DistMatrix = Dist_Matrix_Builder(nospecies, species, b, statdict, PiDict)	#Final Distance Matrix
	FileWriter(DistMatrix, nospecies, species,folder)
    
if __name__== "__main__" :
	main()












