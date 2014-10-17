from rosetta import *
rosetta.init()
import numpy
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from numpy import linalg


def find_CN_lens(CN_len,pose):
	'''Takes in a pose and list of CN_lengths for other poses. Finds CN lengths in pose
	and updates list with new values. Returns updated list. '''
	#cycles through every residue in pose
	for i in range(1,pose.total_residue()):
		#ensures that residues on different chains don't skew results
		if pose.pdb_info().chain(i)==pose.pdb_info().chain(i+1):
			#find atom coordinates
			C=pose.residue(i).xyz('C')
			N=pose.residue(i+1).xyz('N')
		
			C_coord=numpy.array((C[0],C[1],C[2]))
			N_coord=numpy.array((N[0],N[1],N[2]))
			#find bond length
			bond_len=numpy.linalg.norm(C_coord-N_coord)
			CN_len.append(bond_len)

	return CN_len
	
	
def make_plots(CN_len):
	'''Takes in list of CN bond lengths, creates a histogram of the data and then overlays
	the CHARMm parameter equation onto it'''

	#finds endpoints of data
	mx=max(CN_len)
	mn=min(CN_len)
	
	#finds number of bins needed
	r=(mx-mn)/.01
	#creates bins of appropriate size
	bins=[]
	for i in range(int(r)+1):
		bins.append(mn+i*.01)
	bins.append(mx)
	
	#labels plot
	f=plt.figure(1)	
	#plots histogram
	#n is probabilities (y-axis values)
	n=plt.hist(CN_len,bins,normed=True)
	plt.xlabel('C-N Bond Length (Angstoms)')
	plt.ylabel('Probability (%) of Occurance')
	plt.title('Occurance of CN Bond Lengths')
	f.show()

		
		
	p=[]
	e=[]
	kT=.6
	
	#calculate energy from probability
	for i in range(len(n[0])):
		e.append(-kT*numpy.log(n[0][i]/100.0))
	#removes max value from bins, for purpose of making bar graph
	bins.remove(mx)

	g=plt.figure(2)	
	Kb=471
	b0=1.33
	
	K=600
	#range of plotted function
	x=numpy.linspace(mn-.01,mx+.01,100)
	#plots energy function
	plt.plot(x,Kb*(x-b0)**2,color='red')
	plt.plot(x,K*(x-b0)**2+e[int(.5*len(e))]-1,color='green')
	plt.xlabel('C-N Bond Length (Angstroms)')
	plt.ylabel('Energy')
	plt.title('Bond Length v. Bond Energy')
	plt.bar(bins,e,.01)
	g.show()
	
	#keeps figures visible
	raw_input()
	
def main():
	CN_len=[]
	from toolbox import cleanATOM
	#Note that I could not access the WHATIF set so I used a set of structures from the PDB
	pdb_ID=['EGFR','centuximab','1BKR','5P21','1E6K','1F21','1R9H','2HDA','2O72','2IT6']
	
	for i in range(len(pdb_ID)):
		cleanATOM(pdb_ID[i]+'.pdb')
		#create pose from PDB
		pose=pose_from_pdb(pdb_ID[i]+'.clean.pdb')
		
		#find CN lengths for all bonds in PDB
		#updates list that contains bond lengths of all poses
		CN_len=find_CN_lens(CN_len,pose)
		
		
	#write CN lengths to a file
	write_file=open('CA_N_Bond_Lens.txt','w')
	for i in range(len(CN_len)):
		write_file.write(str(CN_len[i])+'\n')
	
	#generate plots
	make_plots(CN_len)
	
main()