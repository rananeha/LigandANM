#! /usr/bin/python
import sys, os, os.path, shutil, getopt, math, glob, time
sys.path.append(os.getcwd())
sys.path.append('/home/nrana/Desktop/Software/ENM_concoord_3.0')
import scipy
import numpy
from numpy import *
#import scipy.weave
#import scipy.weave.converters
#from scipy.sparse import linalg
from scipy import weave
from scipy.weave import converters
from prody import *
import numpy as np
import scipy.linalg as linalg
from inspect import currentframe, getframeinfo
import operator
import pp
import pulchra
import computeGOAP 
import computeRMSD3
import computeQMEAN
import computeJerniganPotential
#import pd2_ca2
import scwrl_run
import datetime
#from computeEnergy1 import computeEnergy_allmodels
#from computeRMSD import align2EachStructure 
DSSP_exe = "/home/shared/Dock_Neha/ENM_concoord_3.0/dssp-2.0.4-linux-i386"
#DSSP_exe = "/usr/local/dssp-2-linux-i386"
Null = 0.0001
###################################################################################################
###################################################################################################
def	prepareProtein4ENM(filepdb_in, filepdb_out):
	global map_resid_to_PDB
	global map_chain_to_PDB
	global num_lig_atoms, num_prot_atoms
	global lig_coords

	residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HIE", "HID", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
	fin = open(filepdb_in, "r")
	fout = open(filepdb_out, "w")
	
	map_resid_to_PDB = []
	map_chain_to_PDB = []

	ci = 1 
	cr = 0
	old_resnam = "999"
	old_resid = -999
	for line in fin:
		label = line[0:6]
		if label.find("ATOM") == 0:
			resnam = line[17:20].strip()
			chainid = line[21:22]
			resid = int(line[22:26].strip())
			#atmnam = line[13:16].strip()
			alternative = line[16:17]
			if alternative == 'A' or alternative == ' ':
				for res_list in residues:
					if resnam == res_list:
						if resid != old_resid or resnam != old_resnam:
							cr += 1
							map_resid_to_PDB.append(resid)
							map_chain_to_PDB.append(chainid)
							old_resid = resid
							old_resnam = resnam
						line_out = "ATOM % 6d" % ci
						line_out += line[11:16]
						line_out += ' '
						line_out += line[17:21]
						line_out += "% 5d " % cr
						line_out += line[27:]
						fout.write(line_out)
						ci += 1

	num_prot_atoms = ci - 1
	#print num_prot_atoms
	fin.close()
	fout.close()

	# If ligand is present, add ligand.
	lig_coords = []
	num_lig_atoms = 0 
	all_lig_atoms = 1
	for z in range(len(lig_names)):
		L = Ligand(lig_names[z], int(lig_nums[z]))
		num_lig_atoms = L.numAtoms(filepdb_in, num_prot_atoms)
		all_lig_atoms = L.addLigand(filepdb_in, filepdb_out, num_prot_atoms, num_lig_atoms, all_lig_atoms)
		lig_coors = L.findLigCoords(filepdb_in, num_lig_atoms)
		for elements in lig_coors:
			lig_coords.append(elements)
###################################################################################################
###################################################################################################
class Ligand:
	def __init__(self, ligname, lignum):
		self.ligname = ligname
		self.lignum = int(lignum)
#		self.lig_coords = []

	def numAtoms(self, filepdb_in, num_prot_atoms):
		fin = open(filepdb_in, "r")
		num_lig_atoms = 0
		for line in fin:
			label = line[0:6]
			if label.find("HETATM") == 0:
				resnam = line[17:20].strip()
				chainid = line[21:22]
				resid = int(line[22:26].strip())
				alternative = line[16:17]
				if alternative == 'A' or alternative == ' ':
					if resnam == self.ligname and resid == int(self.lignum):
						num_lig_atoms += 1
		fin.close()
		return num_lig_atoms
						
	def addLigand(self, filepdb_in, filepdb_out, num_prot_atoms, num_lig_atoms, all_lig_atoms):
		fin = open(filepdb_in, "r")
		#fout = open(filepdb_out, "a+")
		ligfile = open("lig.pdb","w")
		print num_prot_atoms, num_lig_atoms, all_lig_atoms
		ci = num_prot_atoms + all_lig_atoms
		for line in fin:
			label = line[0:6]
			if label.find("HETATM") == 0:
				resnam = line[17:20].strip()
				chainid = line[21:22]
				resid = int(line[22:26].strip())
				alternative = line[16:17]
				if alternative == 'A' or alternative == ' ':
					if resnam == self.ligname and resid == int(self.lignum):
#						print resnam, resid
						line_out = "HETATM% 5d" % ci
#						line_out = "ATOM % 6d" % ci
						line_out += line[11:16]
						line_out += ' '
						line_out += line[17:21]
						line_out += "% 5d " % int(self.lignum)
						line_out += line[27:]
#						fout.write(line_out)
						ligfile.write(line_out)
						ci += 1
		all_lig_atoms = ci - num_prot_atoms
		print all_lig_atoms
		fin.close()
		#fout.close()
		ligfile.close()
		os.popen("cat %s lig.pdb > %s.tmp"%(filepdb_out, filepdb_out)).read()
		os.popen("cp %s.tmp %s"%(filepdb_out,filepdb_out))
		return all_lig_atoms
	#*************************************************************************************************
	def findLigCoords(self, filepdb, num_lig_atoms):
		fin = open(filepdb, "r")
		tmp_lig_coords = scipy.zeros((num_lig_atoms,3))
		z = 0
		for lines in fin:
			tmp_lines = lines.split()
			if tmp_lines[0][0:6] == "HETATM":
				res_name = lines[16:20].strip()
				chainid = lines[21:22].strip()
				alternative = lines[16:17]
#				print res_name, chain_name, lig_name
				if self.ligname == res_name and self.lignum == int(lines[22:26]):
					if alternative == 'A' or alternative == ' ':
#						print lines
	 					tmp_lig_coords[z][0] = float(lines[30:38].strip())
						tmp_lig_coords[z][1] = float(lines[38:46].strip())
						tmp_lig_coords[z][2] = float(lines[46:54].strip())
						z += 1

		fin.close()
#		print lig_coords
		return tmp_lig_coords
	#*************************************************************************************************
	def findFuncGroups(self):
		func_grps = []
		#os.popen("babel -ipdb lig.pdb -osdf lig.sdf -h &>/dev/null")
		os.popen("babel -ipdb lig.pdb -osdf lig.sdf -h ")
		os.popen("/home/nrana/Downloads/checkmol-0.5-linux-x86_64 -p lig.sdf &> abc.txt")
		for lines in open("abc.txt"):
			# lines.strip()[1:4] = functype, lines.strip().split(':')[1] = # of same functype, atoms involved in functype
			func_grps.append((lines.strip()[1:4],lines.strip().split(':')[1],lines.strip().split(':')[-1]))
		print func_grps
		return func_grps
        #*************************************************************************************************
	def userdefined(self, lig_coords, txtfile):
		centroid_coords = []
		labels = []
		tfile = open("%s"%txtfile)
		for lines in tfile:
			lines=int(lines.strip().split()[0])
			print lines
			centroid_coords.append((lig_coords[lines-1][0],lig_coords[lines-1][1],lig_coords[lines-1][2]))
		tfile.close()
		#tmp_labels = np.zeros((len(centroid_coords),),dtype=np.int)
		#labels = [x for x in range(len(tmp_labels))]
		labels = [5,5,0,3,1,3,2,3,3,3,3,5,5,5,4,4,4,4,4,4,5]
		print centroid_coords, labels
		return centroid_coords, labels		
	#*************************************************************************************************
	def KMeansClusterCentroids_Ligand(self, file_pdb, ligname, lignum, lig_coords):
		from sklearn.cluster import KMeans
		L = Ligand(ligname, lignum)
		num_lig_atoms = L.numAtoms(file_pdb, num_prot_atoms)
		lig_coords_subset = L.findLigCoords(file_pdb, num_lig_atoms)
		k = int(math.ceil(float(num_lig_atoms)/8))
		print "Num of ligand centroids ", k
		
		X = KMeans(n_clusters=k).fit(lig_coords_subset)
		centroid_coords = X.cluster_centers_
		labels = X.labels_
		return centroid_coords, labels

	#*************************************************************************************************
	def KMedoids_Ligand(self, file_pdb, ligname, lignum, lig_coords):
		import pyclust
		L = Ligand(ligname, lignum)
		num_lig_atoms = L.numAtoms(file_pdb, num_prot_atoms)
		lig_coords_subset = L.findLigCoords(file_pdb, num_lig_atoms)
		k = int(math.ceil(float(num_lig_atoms)/8))
		print "Num of ligand centroids ", k

		X = pyclust.KMedoids(n_clusters=k, n_trials=50)
		X.fit(lig_coords_subset)
		print X
		centroid_coords = X.centers_
		labels = X.labels_
		return centroid_coords, labels

	#*************************************************************************************************
	def AffinityPropagation_Ligand(self, lig_coords):
		from sklearn.cluster import AffinityPropagation
		af = AffinityPropagation().fit(lig_coords)
		centroid_coords = af.cluster_centers_
		labels = af.labels_
		return centroid_coords, labels

	#*************************************************************************************************
	def HierarchialClusterCentroids_Ligand(self, lig_coords, num_lig_atoms):
		centroid_coords = []
		from sklearn.cluster import AgglomerativeClustering
		from sklearn.neighbors import kneighbors_graph
		k = int(num_lig_atoms/8)
		connectivity = kneighbors_graph(lig_coords, n_neighbors=10, include_self=True)
		print connectivity
		X = AgglomerativeClustering(n_clusters=k, connectivity=connectivity, linkage='ward').fit(lig_coords)
		labels = X.labels_
		n_clusters_ = len(set(labels))-(1 if -1 in labels else 0)
		for m in range(n_clusters_):
			centroid_jx=0.0
			centroid_jy=0.0
			centroid_jz=0.0
			for n in range(len(labels)):
				if m == n:
					centroid_jx += lig_coords[n][0]
					centroid_jy += lig_coords[n][1]
					centroid_jz += lig_coords[n][2]
			centroid_coords.append((centroid_jx,centroid_jy,centroid_jz))
		return centroid_coords, labels
	#*************************************************************************************************
	def DBSCANClusterCentroids_Ligand(self, lig_coords):
		centroid_coords = []
		from sklearn.cluster import DBSCAN
		db = DBSCAN(eps=1.5,min_samples=1).fit(lig_coords)
		labels = db.labels_
		n_clusters_ = len(set(labels))-(1 if -1 in labels else 0)
		for m in range(n_clusters_):
			centroid_jx=0.0
			centroid_jy=0.0
			centroid_jz=0.0
			for n in range(len(labels)):
				if m == n:
					centroid_jx += lig_coords[n][0]
					centroid_jy += lig_coords[n][1]
					centroid_jz += lig_coords[n][2]
			centroid_coords.append((centroid_jx,centroid_jy,centroid_jz))
		return centroid_coords, labels

###################################################################################################
###################################################################################################
def	MapCentroidsToLigGroups(hbond_file, filepdb, labels, num_residues, ligname, lignum):

	global centroids_mass
	atom_masses = {'C':12.01,'N':14.0,'O':15.999,'F':18.998,'P':30.973,'S':32.066,'CL':35.452,'Si':28.085,'Br':79.904,'I':126.90447,'H':1.00794}
	new_labels = []
	for i in range(len(labels)):
#		new_labels.append(labels[i] + num_prot_atoms)
		new_labels.append(labels[i] + num_residues + 1)
	cent_set = set(new_labels)
	centroids = list(cent_set)
	centroids.sort()
	print labels, new_labels, centroids
	centroids_mass = scipy.zeros(len(centroids))

	atom_mass_type = []
	fin = open("%s"%filepdb)
	fout = open("lig.pdb","w")
	ci = num_residues
	z = 0
	map_type_to_index = []
	for line in fin:
		label = line[0:6]
		if label.find("HETATM") == 0:
			resnam = line[17:20].strip()
			chainid = line[21:22]
			resid = int(line[22:26].strip())
			alternative = line[16:17]
			if alternative == 'A' or alternative == ' ':
				if resnam == ligname and resid == int(lignum):
					atom_mass_type.append(line[70:80].strip())
					line_out = "HETATM% 5d" % ci
#					line_out = "ATOM % 6d" % ci
					line_out += line[11:16]
					line_out += ' '
					line_out += line[17:21]
					line_out += "% 5d " % int(lignum)
					line_out += line[27:]
					fout.write(line_out)
#					ligfile.write(line_out)
					map_type_to_index.append((z,line[11:16].strip()))
					ci += 1
					z += 1
	fin.close()
	fout.close()
	for l in range(len(labels)):
		mass_type = atom_mass_type[l]
		centroids_mass[int(labels[l])] += atom_masses[mass_type]

	print centroids_mass[0]
	#sys.exit(2)
	#Rewrites hb2 to replace ligand's name by GLY and atom's number by centroid
	z = 0
	flag = 0
	skipNextLine = 0
	hb_file = open("%s"%hbond_file)
	hbond_outfile = open("%s_new.hb2"%hbond_file[:-4],"w")
	for line in hb_file:
		"""
		if line.find("[ HYDROPHOBIC CLUSTERS ]") >= 0:
			skipNextLine = 2
			hbond_outfile.write(line)
		if skipNextLine == 2:
			hbond_outfile.write(line)
		if skipNextLine == 1:
			lig_atom_1 = int(line[0:6].strip())
			lig_atom_2 = int(line[26:32].strip())
			res_name_1 = line[17:22].strip()
			res_name_2 = line[43:48].strip()
			if line.find(lig_name) >= 0:
				if res_name_1 == lig_name and res_name_2 != lig_name:
					index = lig_atom_1 - num_prot_atoms - 1
					hbond_outfile.write("%s% 6d  GLY%s"%(line[:11],new_labels[index],line[22:]))
				elif res_name_1 != lig_name and res_name_2 == lig_name:
					index = lig_atom_2 - num_prot_atoms - 1
					print index, line
					hbond_outfile.write("%s% 5d  GLY%s"%(line[:38],new_labels[index],line[48:]))
				else: 
					index_1 = lig_atom_1 - num_prot_atoms - 1
					index_2 = lig_atom_2 - num_prot_atoms - 1
					hbond_outfile.write("%s% 6d  GLY%s 5d  GLY%s"%(line[:11],new_labels[index_1],line[22:38],new_labels[index_2],line[48:]))
			else:
				hbond_outfile.write(line)

		if flag == 1:
			skipNextLine = 1
			flag = 0
		if line.find("[ HBONDS ]") >= 0:
			flag = 1
			hbond_outfile.write(line)
		"""
		if z >= 8:
			lig_atom_1 = int(line[1:5].strip())
			lig_atom_2 = int(line[15:19].strip())
			res_name_1 = line[6:9].strip()
			res_name_2 = line[20:23].strip()
			atom_type_1 = line[10:14].strip()
			atom_type_2 = line[24:28].strip()
			#print lig_atom_1, lig_atom_2, atom_type_1, atom_type_2, res_name_1, res_name_2
			if line.find(ligname) >= 0 and (lig_atom_1 == int(lignum) or lig_atom_2 == int(lignum)):
				if lig_atom_1 == int(lignum) and lig_atom_2 != int(lignum):
					#index = lig_atom_1 - num_prot_atoms - 1
					for items in range(len(map_type_to_index)):
						if atom_type_1 == map_type_to_index[items][1]:
							index = map_type_to_index[items][0]
					hbond_outfile.write("%s%04d-GLY%s"%(line[0:1],new_labels[index],line[9:]))
				elif lig_atom_1 != int(lignum) and lig_atom_2 == int(lignum):
					#index = lig_atom_2 - num_prot_atoms - 1
					for items in range(len(map_type_to_index)):
						if atom_type_2 == map_type_to_index[items][1]:
							index = map_type_to_index[items][0]
					hbond_outfile.write("%s%04d-GLY%s"%(line[:15],new_labels[index],line[23:]))
				else:
			#		print line
					#index_1 = lig_atom_1 - num_prot_atoms - 1
					#index_2 = lig_atom_2 - num_prot_atoms - 1
					for items in range(len(map_type_to_index)):
						if atom_type_1 == map_type_to_index[items][1]:
							index_1 = map_type_to_index[items][0]
						if atom_type_2 == map_type_to_index[items][1]:
							index_2 = map_type_to_index[items][0]
					hbond_outfile.write("%s%04d-GLY%s%04d-GLY%s"%(line[0:1],new_labels[index_1],line[9:15],new_labels[index_2],line[23:]))
			else:
				hbond_outfile.write(line)
			
		z += 1
	hb_file.close()
	hbond_outfile.close()	
	os.popen("cp %s_new.hb2 %s"%(hbond_file[:-4], hbond_file))
	return centroids
###################################################################################################
def	WritePDB_original_ids(filepdb_in, filepdb_out):
	fin = open(filepdb_in, "r")
	fout = open(filepdb_out, "w")

	ci = 1 
	cr = 0
	for line in fin:
		label = line[0:6]
		if label.find("ATOM") == 0:
			resid = int(line[22:26].strip())
			#print resid
			cr = resid-1
			resid = map_resid_to_PDB[cr]
			chainid = map_chain_to_PDB[cr]
			#print cr, resid, chainid
			line_out = line[:21]
			line_out += "%1s" % chainid
			line_out += "% 4d " % resid
			line_out += line[27:]
			fout.write(line_out)

	fin.close()
	fout.close()
	#print map_chain_to_PDB, map_resid_to_PDB
	#sys.exit(0)
###################################################################################################
###################################################################################################
def MinimizeStructure(filein, dirtmp, filemdp):
	# generate pdb file for GROMACS
	fcmd = open("run_pdb2gmx", "w")
#	fcmd.write("pdb2gmx -f tmp.pdb -o gmx.pdb -ignh -renum -ff oplsaa &>/dev/null << EOF\n")
	fcmd.write("pdb2gmx -f %s -o gmx.pdb -ignh -ff oplsaa &>/dev/null << EOF\n"%filein)
	fcmd.write("6\n")
	fcmd.write("EOF")
	fcmd.close()
	os.popen("chmod a+x run_pdb2gmx; ./run_pdb2gmx").read()
	find_flag = 0
	filename = "gmx.pdb"
	for root, dirs, names in os.walk(dirtmp):
		if filename in names:
			find_flag = 1
	if find_flag == 0:
		sys.exit()

	# minimize
	comm1 = "grompp -f %s -c gmx.pdb &>/dev/null" % filemdp
	print comm1
	os.popen(comm1).read()
#	os.popen("grompp -f %s -c gmx.pdb" % filemdp).read()
#	flog.write("grompp -f %s -c gmx.pdb &>/dev/null\n" % filemdp)
	os.popen("mdrun -c em.pdb -nt 1 &>/dev/null").read()
#	os.popen("cp gmx.pdb em.pdb").read()
	os.popen("rm \#*").read()
	os.popen("cp em.pdb %s"%filein).read()
	
###################################################################################################
###################################################################################################
def	assignSecondaryStructure(filepdb, secondary_structure, int_resn, sheet_resn):
	tmp_dir = os.getcwd()
#	flog.write("%s\n" % tmp_dir)
#	flog.write("%s -i %s -o output.dssp &>/dev/null\n" % (DSSP_exe, filepdb))
#	os.popen("%s -i %s -o output.dssp &>/dev/null" % (DSSP_exe, filepdb)).read()
	os.popen("%s -i %s -o output.dssp " % (DSSP_exe, filepdb)).read()
	# read dssp output
	fin = open("output.dssp", "r")
	flag = 0
	while 1:
		line = fin.readline()
		if len(line) < 5:
			break
		if flag == 1:
			if line[13:14] != "!":
				secondary_structure.append(line[16:17])
				int_resn.append(int(line[41:45]))
		if line.find("#  RESIDUE AA") >= 0:
			flag = 1
	fin.close()
	i = 0
	print "length: len resn", len(int_resn)
	for ss,resn in zip(secondary_structure, int_resn):
		i += 1
		if ss == "E":
			resn2 = i + resn
			if secondary_structure[resn2-1] == "E":
				sheet_resn.append("E")
			else:
				sheet_resn.append("A")
		else:
			sheet_resn.append(ss)
	print len(sheet_resn)
	print len(secondary_structure)

###################################################################################################
###################################################################################################
def	tCONCOORD_Analysis(filetdist, filepdb):
	"""
	fcmd = open("run_dist", "w")
	fcmd.write("dist -s %s -inp %s << EOF\n" % (file_pdb, filetdist))
	fcmd.write("1\n")
	fcmd.write("EOF")
	fcmd.close()
	os.popen("chmod a+x run_dist; ./run_dist").read()
	"""
	print "tdist -s %s -inp %s " % (filepdb, filetdist)
#	os.popen("tdist -s %s -inp %s &>/dev/null" % (filepdb, filetdist)).read()
	os.popen("tdist -s %s -inp %s " % (filepdb, filetdist)).read()

###################################################################################################
###################################################################################################
def	hbplus_Analysis(filepdb):
#	os.popen("/usr/local/hbplus/hbplus %s &>/dev/null" % filepdb).read()
	os.popen("/usr/local/hbplus/hbplus %s " % filepdb).read()
###################################################################################################
###################################################################################################
def	IdentifySaltbridgesBwProtLig(func_grps, centroids, centroid_coords, lig_coords, ligname, lignum, num_residues):

#	salt_bridge_atoms = ['N-oxide', 'quaternary ammonium salt', 'carboxylic acid salt', 'diazonium salt', 'sulfuric acid deriv.', 'phosphoric acid deriv.', 'carboxylic acid azide', 'isonitrile', 'nitrate', 'nitro compound', 'azide', 'carbonic acid deriv.']
#	salt_bridge_atoms = ['060','059','077','142','153','181','086','143','152','150','140','117']
	metals = ['Fe','FE','Mg','MG','Ca','CA','Co','CO','Mn','MN','Ni','NI','Mo','MO','Zn','ZN','Cu','CU']
	salt_bridge_atoms_pos = ['048', '059','086','140','142','143']
	salt_bridge_atoms_neg = ['060', '076', '077', '086','117','140','143','150','152','153','181']
	func_grps = []
	salt_bridge_centroids_pos = []
	salt_bridge_centroids_neg = []
	if ligname in metals:
		salt_bridge_centroids_pos.append(centroids[0]-num_residues-1)
		#print centroids[0], num_residues, salt_bridge_centroids_pos, salt_bridge_centroids_neg
		return salt_bridge_centroids_pos, salt_bridge_centroids_neg
	"""
	fin = open("%s"%filepdb)
	fout = open("lig.pdb","w")
	ci = num_residues
	for line in fin:
		label = line[0:6]
		if label.find("HETATM") == 0:
			resnam = line[17:20].strip()
			chainid = line[21:22]
			resid = int(line[22:26].strip())
			alternative = line[16:17]
			if alternative == 'A' or alternative == ' ':
				if resnam == ligname and resid == int(lignum):
					line_out = "HETATM% 5d" % ci
#					line_out = "ATOM % 6d" % ci
					line_out += line[11:16]
					line_out += ' '
					line_out += line[17:21]
					line_out += "% 5d " % int(lignum)
					line_out += line[27:]
					fout.write(line_out)
					#ligfile.write(line_out)
					ci += 1
	fin.close()
	fout.close()
	"""
	#print func_grps
	#check abc.txt file for understanding this section
	for f in range(len(func_grps)):
		closest_cent_dist = 9999.99
		dist = 0.0
		cent = 0
		if func_grps[f][0] in salt_bridge_atoms_pos:
			#print func_grps[f][0], "cation"
			for num_func in range(int(func_grps[f][1])):
				func_atom_num = int(func_grps[f][2].split(',')[num_func].split('-')) - 1
				for cent in range(len(centroid_coords)):
					#dist = math.sqrt(math.pow( (centroid_coords[cent][0] - lig_coords[int(func_grps[f][1])-1][0] ),2) + math.pow( (centroid_coords[cent][1] - lig_coords[int(func_grps[f][1])-1][1] ),2) + math.pow( (centroid_coords[cent][2] - lig_coords[int(func_grps[f][1])-1][2] ),2) )
					dist = math.sqrt(math.pow( (centroid_coords[cent][0] - lig_coords[func_atom_num][0] ),2) + math.pow( (centroid_coords[cent][1] - lig_coords[func_atom_num][1] ),2) + math.pow( (centroid_coords[cent][2] - lig_coords[func_atom_num][2] ),2) )
					if dist < closest_cent_dist:
						closest_cent_dist = dist
						closest_cent = cent
					else:
						continue			
				salt_bridge_centroids_pos.append(closest_cent)

		if func_grps[f][0] in salt_bridge_atoms_neg:
			#print func_grps[f][0], "anion"
			for num_func in range(int(func_grps[f][1])):
				func_atom_num = int(func_grps[f][2].split(',')[num_func].split('-')) - 1
				for cent in range(len(centroid_coords)):
					dist = math.sqrt(math.pow( (centroid_coords[cent][0] - lig_coords[int(func_grps[f][1])-1][0] ),2) + math.pow( (centroid_coords[cent][1] - lig_coords[int(func_grps[f][1])-1][1] ),2) + math.pow( (centroid_coords[cent][2] - lig_coords[int(func_grps[f][1])-1][2] ),2) )
					if dist < closest_cent_dist:
						closest_cent_dist = dist
						closest_cent = cent
					else:
						continue
				salt_bridge_centroids_neg.append(closest_cent)
	
	return salt_bridge_centroids_pos, salt_bridge_centroids_neg				
###################################################################################################	
###################################################################################################
def 	IdentifyAromaticBwProtLig(func_grps, centroid_coords, lig_coords):

	for f in range(len(func_grps)):
		closest_cent_dist = 9999.99
		dist = 0.0
		cent = 0
		func_atom_num = []
		center_of_func = [0.0,0.0,0.0]
		lig_aromatic_cent = []
		if func_grps[f][0] == '201' :
			for num_func in range(int(func_grps[f][1])):
				func_atom_num = func_grps[f][2].split(',')[num_func].split('-')[:-1]
				for x in func_atom_num:
					center_of_func[0] += lig_coords[int(x)-1][0]
					center_of_func[1] += lig_coords[int(x)-1][1]
					center_of_func[2] += lig_coords[int(x)-1][2]
				center_of_func[0] /= len(func_atom_num)
				center_of_func[1] /= len(func_atom_num)
				center_of_func[2] /= len(func_atom_num)
				for cent in range(len(centroid_coords)):
					dist = math.sqrt(math.pow( (centroid_coords[cent][0] - center_of_func[0] ),2) + math.pow( (centroid_coords[cent][1] - center_of_func[1] ),2) + math.pow( (centroid_coords[cent][2] - center_of_func[2] ),2) )
					#print cent, dist, closest_cent_dist
					if dist < closest_cent_dist:
						closest_cent_dist = dist
						closest_cent = cent
					else:
						continue
				lig_aromatic_cent.append(closest_cent)
				func_atom_num = []
				center_of_func = [0.0,0.0,0.0]
				closest_cent_dist = 9999.99
				dist = 0.0
				cent = 0
	print "lig_aromatic_cent", lig_aromatic_cent
	return lig_aromatic_cent
###################################################################################################     
###################################################################################################
def	identifyNumberOfResidues(filepdb):
	residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HIE", "HID", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
	fin = open(filepdb, "r")
      
	# find directly bonded residues in sequence
	global num_residues
	global ENM_res_list

	# find number of residues
	num_residues = 0 
	ENM_res_list = []
	       
	for line in fin:
		label = line[0:6].strip()
		resnam = line[17:20]
		atmnam = line[13:16].strip()
		for res_list in residues:
			if label == 'ATOM' and resnam == res_list and atmnam == 'CA':
				ENM_res_list.append(resnam)
				num_residues += 1
	fin.close()
	return num_residues
###################################################################################################
###################################################################################################
def	identifyBondedRestraints(filepdb, bonded_A, bonded_B, disulfide_A, disulfide_B):
#	global num_residues
	residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HIE", "HID", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
	fin = open(filepdb, "r")
      
	# find directly bonded residues in sequence
	tmp_store_x = []
	tmp_store_y = []
	tmp_store_z = []
	tmp_store_x2 = []
	tmp_store_y2 = []
	tmp_store_z2 = []
	tmp_store_id = []
	tmp_store_id2 = []
	
	flag = 0 
	for line in fin:
		resnam = line[17:20]
		atmnam = line[13:16].strip()
		for res_list in residues:
			if resnam == res_list and atmnam == 'N':
				#print resnam, line[22:27]
				n_x = float(line[30:38])  
				n_y = float(line[38:46])
				n_z = float(line[46:54])
				n_id = line[22:27]  
				flag += 1

			if resnam == res_list and atmnam == 'C':
				o_x = float(line[30:38])  
				o_y = float(line[38:46])
				o_z = float(line[46:54])
				o_id = line[22:27]  
				flag += 1
                
			if flag == 2:
				tmp_store_x.append(n_x)
				tmp_store_y.append(n_y)
				tmp_store_z.append(n_z)
				tmp_store_x2.append(o_x)
				tmp_store_y2.append(o_y)
				tmp_store_z2.append(o_z)
				tmp_store_id.append(n_id)
				tmp_store_id2.append(o_id)
				flag = 0

	cmax = len(tmp_store_id)
	for ci in range(cmax-1):
		o_id = tmp_store_id2[ci]
		o_x = tmp_store_x2[ci]
		o_y = tmp_store_y2[ci]
		o_z = tmp_store_z2[ci]
		n_id = tmp_store_id[ci+1]
		n_x = tmp_store_x[ci+1]
		n_y = tmp_store_y[ci+1]
		n_z = tmp_store_z[ci+1]
		#print math.sqrt(math.pow(o_x - n_x, 2) + math.pow(o_y - n_y, 2) + math.pow(o_z - n_z, 2)), o_id, o_x, o_y, o_z, n_x, n_y, n_z
		if(math.sqrt(math.pow(o_x - n_x, 2) + math.pow(o_y - n_y, 2) + math.pow(o_z - n_z, 2)) < 2.5):		# check distance between C and N
			bonded_A.append(o_id.strip())
			bonded_B.append(n_id.strip())
	#sys.exit(2)
	# find disulfide bridges
	tmp_store_x = []
	tmp_store_y = []
	tmp_store_z = []
	tmp_store_id = []
	tmp_store_name = []

	fin.seek(0)
	for line in fin:
		resnam = line[17:20]
		atmnam = line[13:15]
		chain_id = line[21:22]
		flag = 0 
		if (resnam == 'CYS' or resnam == 'CYX') and atmnam == 'SG':
			o_x = float(line[30:38])  
			o_y = float(line[38:46])
			o_z = float(line[46:54])
			o_id = line[22:27]  
                
			tmp_store_x.append(o_x)
			tmp_store_y.append(o_y)
			tmp_store_z.append(o_z)
			tmp_store_id.append(o_id)
			tmp_store_name.append(resnam)

	ci = 0
	for o_id in tmp_store_id:
		o_x = tmp_store_x[ci]
		o_y = tmp_store_y[ci]
		o_z = tmp_store_z[ci]
		cj = 0
		for n_id in tmp_store_id:
			n_x = tmp_store_x[cj]
			n_y = tmp_store_y[cj]
			n_z = tmp_store_z[cj]
			if ci < cj:
				if(math.sqrt(math.pow(o_x - n_x, 2) + math.pow(o_y - n_y, 2) + math.pow(o_z - n_z, 2)) <= 3.0):		# check distance between C and N
					disulfide_A.append(o_id.strip())
					disulfide_B.append(n_id.strip())
			cj += 1
		ci += 1
	fin.close()

###################################################################################################
###################################################################################################
def	identifyNonBondedRestraints(filepdb, filein, saltbridge_A, saltbridge_B, hbonds_A, hbonds_B, num_hbonds, hphob_A, hphob_B, num_hphob, salt_bridge_centroids_cat, salt_bridge_centroids_an, lig_aromatic_cent, aromatic_A, aromatic_B, centroid_coordinates, all_centroids):
	residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HIE", "HID", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
	fin = open(filepdb, "r")
	# find salt bridges
	tmp_store_x = []
	tmp_store_y = []
	tmp_store_z = []
	tmp_store_x2 = []
	tmp_store_y2 = []
	tmp_store_z2 = []
	tmp_store_id = []
	tmp_store_id2 = []
	tmp_store_name = []
	tmp_store_name2 = []

	for line in fin:
		resnam = line[17:20].strip()
		atmnam = line[13:16].strip()
		flag = 0
		if (resnam == 'ASP' and atmnam == 'OD2') or (resnam == 'GLU' and atmnam == 'OE2'):
			o_x = float(line[30:38])  
			o_y = float(line[38:46])
			o_z = float(line[46:54])
			o_id = line[22:27]  
			tmp_store_x.append(o_x)
			tmp_store_y.append(o_y)
			tmp_store_z.append(o_z)
			tmp_store_id.append(o_id)
			tmp_store_name.append(resnam)

		if (resnam == 'LYS' and atmnam == 'NZ') or (resnam == 'ARG' and atmnam == 'NH1'):
			n_x = float(line[30:38])
			n_y = float(line[38:46])
			n_z = float(line[46:54])
			n_id = line[22:27]			
			tmp_store_x2.append(n_x)
			tmp_store_y2.append(n_y)
			tmp_store_z2.append(n_z)
			tmp_store_id2.append(n_id)
			tmp_store_name2.append(resnam)
	flag = 0
	ci = 0
	for o_id in tmp_store_id:
		o_x = tmp_store_x[ci]
		o_y = tmp_store_y[ci]
		o_z = tmp_store_z[ci]
		cj = 0
		for n_id in tmp_store_id2:
			n_x = tmp_store_x2[cj]
			n_y = tmp_store_y2[cj]
			n_z = tmp_store_z2[cj]
			if math.sqrt(math.pow(o_x - n_x, 2) + math.pow(o_y - n_y, 2) + math.pow(o_z - n_z, 2)) <= 9.0 :
				saltbridge_A.append(o_id.strip())
				saltbridge_B.append(n_id.strip())
			cj += 1
		ci += 1
	fin.close()
	ci = 0
	for o_id in tmp_store_id:
		o_x = tmp_store_x[ci]
		o_y = tmp_store_y[ci]
		o_z = tmp_store_z[ci]
		for n_id in salt_bridge_centroids_cat:
			n_x = centroid_coordinates[n_id][0]
			n_y = centroid_coordinates[n_id][1]
			n_z = centroid_coordinates[n_id][2]
			if math.sqrt(math.pow(o_x - n_x, 2) + math.pow(o_y - n_y, 2) + math.pow(o_z - n_z, 2)) <= 9.0 :
				saltbridge_A.append(o_id.strip())
				saltbridge_B.append(str(n_id+min(all_centroids)))
		ci += 1

	ci = 0
	for o_id in tmp_store_id2:
		o_x = tmp_store_x2[ci]
		o_y = tmp_store_y2[ci]
		o_z = tmp_store_z2[ci]
		for n_id in salt_bridge_centroids_an:
			n_x = centroid_coordinates[n_id][0] 
			n_y = centroid_coordinates[n_id][1]
			n_z = centroid_coordinates[n_id][2]
			if math.sqrt(math.pow(o_x - n_x, 2) + math.pow(o_y - n_y, 2) + math.pow(o_z - n_z, 2)) <= 12.0 :
				saltbridge_A.append(o_id.strip())
				saltbridge_B.append(str(n_id+min(all_centroids)))
		ci += 1
	ci = 0
	for o_id in salt_bridge_centroids_cat:
		o_x = centroid_coordinates[o_id][0]
		o_y = centroid_coordinates[o_id][1]
		o_z = centroid_coordinates[o_id][2]
		for n_id in salt_bridge_centroids_an:
			n_x = centroid_coordinates[n_id][0]
			n_y = centroid_coordinates[n_id][1]
			n_z = centroid_coordinates[n_id][2]
			if math.sqrt(math.pow(o_x - n_x, 2) + math.pow(o_y - n_y, 2) + math.pow(o_z - n_z, 2)) <= 5.5 :
				saltbridge_A.append(str(o_id+min(all_centroids)))
				saltbridge_B.append(str(n_id+min(all_centroids)))
		ci += 1	
	#print saltbridge_A, saltbridge_B
	#------------------------------------------------------------------
        # combine multiple salt bridges between same residues
	tmp_store_A = []
	tmp_store_B = []
	for i, j in zip(saltbridge_A, saltbridge_B):
		ii = int(i)
		jj = int(j)
		# put residue id's in order
		if ii > jj:
			ii2 = jj
			jj2 = ii
		else:
			ii2 = ii
			jj2 = jj
		# check if already assigned hydrogen bonds between same residues
		flag = 0
		ci = 0
		for k, l in zip(tmp_store_A, tmp_store_B):
			# found residue pair?
			if k == ii2 and l == jj2:
				print k, l
				flag = 1
			ci += 1
		# new pair of residues
		if flag == 0:
			tmp_store_A.append(ii2)
			tmp_store_B.append(jj2)

	print "satbridge", tmp_store_A, tmp_store_B
	saltbridge_A = []
	saltbridge_B = []
	saltbridge_A = tmp_store_A
	saltbridge_B = tmp_store_B
	fin = open(filein, "r")
	#------------------------------------------------------------------
	# find h-bonded and hydrophobic residues
	# store temporarily stable hydrogen bonds
	tmp_store_A = []
	tmp_store_B = []
	tmp_phob_A = []
	tmp_phob_B = []
	tmp_phob_C = []
        
	flag = 0 
	while 1: 
		line = fin.readline()
		if len(line) < 5:
			break
		# read hydrophobic contacts
		if flag == 2:
			resid_A = line[0:9].strip()
			resid_B = line[10:18].strip()
			resid_C = line[19:27].strip()
			tmp_phob_A.append(resid_A)
			tmp_phob_B.append(resid_B)
			tmp_phob_C.append(resid_C)
		if line.find("[ HYDROPHOBIC CLUSTERS ]") >= 0:
			flag = 2
			fin.readline()
		"""
		# read hydrogen bonds
		# tdist method of reading HBonds
		if flag == 1:
			resid_A = line[11:18].strip()
			resid_B = line[37:44].strip()
			hbond_flag = line[188:192].strip()
			if hbond_flag == "yes":
				tmp_store_A.append(resid_A)
				tmp_store_B.append(resid_B)
		if line.find("[ HBONDS ]") >= 0:
			flag = 1
			fin.readline()
		"""
	fin.close()
	
#	hbplus method of reading HBonds
	i = 0
#	for files in os.listdir():
#		if files.find(".hb2") >= 0:
	hbond_file = open("%s.hb2"%filepdb[:-4])
	print hbond_file
	for line in hbond_file:
		if i >= 8:
			resid_A = int(line[1:5])
			resid_B = int(line[15:19])
			tmp_store_A.append(resid_A)
			tmp_store_B.append(resid_B)
		i += 1			
	hbond_file.close()
	print tmp_store_A, tmp_store_B
	
	#------------------------------------------------------------------
	# combine multiple hydrogen bonds between same residues
	for i, j in zip(tmp_store_A, tmp_store_B):
		ii = int(i)
		jj = int(j)
		# put residue id's in order
		if ii > jj:
			ii2 = jj
			jj2 = ii
		else:
			ii2 = ii
			jj2 = jj
		# check if already assigned hydrogen bonds between same residues
		flag = 0
		ci = 0
		for k, l in zip(hbonds_A, hbonds_B):
			# found residue pair?
			if k == ii2 and l == jj2:
				print k, l
				flag = 1
				num_hbonds[ci] += 1
			ci += 1
		# check if already assigned salt bridges between same residues
		for k, l in zip(saltbridge_A, saltbridge_B):
			kk = int(k)
			ll = int(l)
			# put residue id's in order
			if kk > ll:
				kk2 = ll
				ll2 = kk
			else:
				kk2 = kk
				ll2 = ll
#			print kk2, ll2
			if kk2 == ii2 and ll2 == jj2:
				flag = 1
		# new pair of residues
		if flag == 0:
			hbonds_A.append(ii2)
			hbonds_B.append(jj2)
			num_hbonds.append(1)
	print hbonds_A, hbonds_B
	#------------------------------------------------------------------
	# identify residue numbers for hydrophobic atoms
	# read em.pdb file and store atom number <--> residue number correspondence
	fin = open(filepdb, "r")
	# number of atoms
	count_lines = 0
	for line in fin:
		label = line[0:6].strip()
		if label.find("ATOM") == 0:
			count_lines += 1
	# store resid[atmid]
	res_atm_list = [0 for i in range(count_lines+1)]
	fin.seek(0)
	for line in fin:
		#print line
		label = line[0:6].strip()
		if label.find("ATOM") == 0:
			resid = int(line[22:27].strip())
			atmid = int(line[6:11].strip())
			res_atm_list[atmid] = resid
#	print res_atm_list
	# combine multiple hydrophobic contacts between same residues
	hphob_A_tmp = []
	hphob_B_tmp = []
	hphob_C_tmp = []
	num_hphob_tmp = []

	for i, j, k in zip(tmp_phob_A, tmp_phob_B, tmp_phob_C):
		ii = int(i)
		jj = int(j)
		kk = int(k)
		# put residue id's in order
		if ii < jj:
			if jj < kk:		# i < j < k
				ii2 = ii
				jj2 = jj
				kk2 = kk
			else:
				if ii < kk:	# i < k < j
					ii2 = ii
					jj2 = kk
					kk2 = jj
				else:				# k < i < j
					ii2 = kk
					jj2 = ii
					kk2 = jj
		else:
			if ii < kk:		# j < i < k
				ii2 = jj
				jj2 = ii
				kk2 = kk
			else:
				if jj < kk:	# j < k < i
					ii2 = jj
					jj2 = kk
					kk2 = ii
				else:				# k < j < i
					ii2 = kk
					jj2 = jj
					kk2 = ii
#		print "% 5d % 5d % 5d" % (ii2, jj2, kk2)
		
		# check if already assigned hydrophobic contacts between same atoms
		flag = 0
		ci = 0
		for k, l, m in zip(hphob_A_tmp, hphob_B_tmp, hphob_C_tmp):
			# found atom triplet?
			if k == ii2 and l == jj2 and m == kk2:
				flag = 1
				num_hphob_tmp[ci] += 1
			ci += 1
		# new triplet of atoms
		if flag == 0:
			hphob_A_tmp.append(ii2)
			hphob_B_tmp.append(jj2)
			hphob_C_tmp.append(kk2)
			num_hphob_tmp.append(1)

	#------------------------------------------------------------------
	# translate into residue pairs and count number of hydrophobic contacts which include these two residues
	for i, j, k in zip(hphob_A_tmp, hphob_B_tmp, hphob_C_tmp):
		ii = res_atm_list[i]
		jj = res_atm_list[j]
		kk = res_atm_list[k]

		# ii-jj
		flag = 0
		ci = 0
		for k, l in zip(hphob_A, hphob_B):
			# found residue pair?
			if k == ii and l == jj:
				flag = 1
				num_hphob[ci] += 1
			ci += 1
		# new triplet of residues
		if flag == 0:
			hphob_A.append(ii)
			hphob_B.append(jj)
			num_hphob.append(1)

		# ii-kk
		flag = 0
		ci = 0
		for k, l in zip(hphob_A, hphob_B):
			# found residue pair?
			if k == ii and l == kk:
				flag = 1
				num_hphob[ci] += 1
			ci += 1
		# new triplet of residues
		if flag == 0:
			hphob_A.append(ii)
			hphob_B.append(kk)
			num_hphob.append(1)

		# jj-kk
		flag = 0
		ci = 0
		for k, l in zip(hphob_A, hphob_B):
			# found residue pair?
			if k == jj and l == kk:
				flag = 1
				num_hphob[ci] += 1
			ci += 1
		# new triplet of residues
		if flag == 0:
			hphob_A.append(jj)
			hphob_B.append(kk)
			num_hphob.append(1)

	fin.seek(0)
	# find aromatic contacts between ligand and protein residues
	tmp_store_x = []
	tmp_store_y = []
	tmp_store_z = []
	tmp_store_id = []
	tmp_store_name = []
	for line in fin:
		resnam = line[17:20].strip()
		atmnam = line[13:16].strip()
		flag = 0
		if resnam in ('HIS','TRP','TYR','PHE') and atmnam == "CA":
			o_x = float(line[30:38])
                        o_y = float(line[38:46])
                        o_z = float(line[46:54])
                        o_id = line[22:27]
                        tmp_store_x.append(o_x)
                        tmp_store_y.append(o_y)
                        tmp_store_z.append(o_z)
                        tmp_store_id.append(o_id)
                        tmp_store_name.append(resnam)
	flag = 0
	ci = 0
	for o_id in tmp_store_id:
		o_x = tmp_store_x[ci]
		o_y = tmp_store_y[ci]
		o_z = tmp_store_z[ci]
		for n_id in lig_aromatic_cent:
			n_x = centroid_coordinates[n_id][0]
			n_y = centroid_coordinates[n_id][1]
			n_z = centroid_coordinates[n_id][2]
                        if math.sqrt(math.pow(o_x - n_x, 2) + math.pow(o_y - n_y, 2) + math.pow(o_z - n_z, 2)) <= 9.0 :
                                aromatic_A.append(o_id.strip())
                                aromatic_B.append(str(n_id+min(all_centroids)))
		ci += 1
                       
	fin.close()
###################################################################################################
###################################################################################################
def	generateForceConstantMatrix(filepdb, force_const, g_width, g_V0, shake_on, shake_list, dr_list, coord_list, mass_list, centroids_mass, bonded_A, bonded_B, disulfide_A, disulfide_B, saltbridge_A, saltbridge_B, hbonds_A, hbonds_B, num_hbonds, hphob_A, hphob_B, num_hphob, aromatic_A, aromatic_B, secondary_structure, int_resn, sheet_resn, all_centroids, centroid_coordinates):
	global num_residues

	residues = ["ALA", "ARG",  "ASN",  "ASP",  "CYS",  "GLN",  "GLU",  "GLY", "HIS", "HIE", "HID", "ILE",  "LEU",  "LYS",  "MET",  "PHE",  "PRO", "SER", "THR",  "TRP",  "TYR",  "VAL"]
	masses =   [71.09, 156.19, 115.09, 114.11, 103.15, 128.14, 129.12, 57.05, 137.14, 137.14, 137.14, 113.16, 113.16, 128.17, 131.19, 147.18, 97.12, 87.08, 101.11, 186.12, 163.18, 99.14]
	fin = open(filepdb, "r")
	num_residues += len(all_centroids)
	tmp_store_x = scipy.zeros(num_residues)
	tmp_store_y = scipy.zeros(num_residues)
	tmp_store_z = scipy.zeros(num_residues)
	tmp_store_id = scipy.zeros(num_residues, int)
	fin.seek(0)	
	# store coordinates for CA atoms
#	num_residues = 0 
	ci = 0       
	for line in fin:
		label = line[0:6].strip()
		resnam = line[17:20]
		atmnam = line[13:16].strip()
		#print resnam, atmnam
		for res_list, mass in zip(residues, masses):
			if label == 'ATOM' and resnam == res_list and atmnam == 'CA':
				#print resnam, atmnam,line[22:27]
				n_x = float(line[30:38])  
				n_y = float(line[38:46])
				n_z = float(line[46:54])
				n_id = line[22:27]  
				tmp_store_x[ci] = n_x
				tmp_store_y[ci] = n_y
				tmp_store_z[ci] = n_z
				tmp_store_id[ci] = n_id
				mass_list[ci] = mass
				coord_list[3*ci  ] = n_x
				coord_list[3*ci+1] = n_y
				coord_list[3*ci+2] = n_z
				ci += 1
#				print resnam, atmnam,line[22:27],n_x, ci
#				print tmp_store_id[ci]
#				num_residues += 1
#	ci += 1
#	print ci	
	fin.close()
	# Add ligand centroids
	print len(all_centroids)
	for cent in range(len(all_centroids)):
		tmp_store_x[ci] = centroid_coordinates[cent][0]
		tmp_store_y[ci] = centroid_coordinates[cent][1]
		tmp_store_z[ci] = centroid_coordinates[cent][2]
		tmp_store_id[ci] = ci
		print "Check", tmp_store_id[ci]
		coord_list[3*ci  ] = centroid_coordinates[cent][0]
		coord_list[3*ci+1] = centroid_coordinates[cent][1]
		coord_list[3*ci+2] = centroid_coordinates[cent][2]
		mass_list[ci] = centroids_mass[cent]
		ci += 1
	"""
	for m in range(len(mass_list)):
		mass_matrix[3*m][3*m] = math.sqrt(mass_list[m])
		mass_matrix[3*m+1][3*m+1] = math.sqrt(mass_list[m])
		mass_matrix[3*m+2][3*m+2] = math.sqrt(mass_list[m])
	print mass_matrix
	"""	

	dx_list = scipy.zeros((num_residues, num_residues))
	dy_list = scipy.zeros((num_residues, num_residues))
	dz_list = scipy.zeros((num_residues, num_residues))
	dr2_list = scipy.zeros((num_residues, num_residues))

	# compute distance matrix (use C function)
	code_distance =	"""
									int			i, j;
									double	x1, y1, z1, x2, y2, z2, dx, dy, dz, dr2, dr;

									for(i = 0; i < num_residues; i++)
									{
										x1 = tmp_store_x(i);
										y1 = tmp_store_y(i);
										z1 = tmp_store_z(i);
										for(j = i+1; j < num_residues; j++)
										{
											x2 = tmp_store_x(j);
											y2 = tmp_store_y(j);
											z2 = tmp_store_z(j);
											dx = x1 - x2;
											dy = y1 - y2;
											dz = z1 - z2;
											dr2 = dx*dx + dy*dy + dz*dz;
											dr = sqrt(dr2);
											dx_list(i,j)  = dx;
											dy_list(i,j)  = dy;
											dz_list(i,j)  = dz;
											dr2_list(i,j) = dr2;
											dr_list(i,j)  = dr;
											dx_list(j,i)  = dx;
											dy_list(j,i)  = dy;
											dz_list(j,i)  = dz;
											dr2_list(j,i) = dr2;
											dr_list(j,i)  = dr;
											//printf("dr: %.3f\ti: %d\tj: %d\\n",dr,i,j);
										}
										//printf("next\\n");
									}
									"""
	weave.inline(code_distance,
							['num_residues', 'tmp_store_x', 'tmp_store_y', 'tmp_store_z', 'dx_list', 'dy_list', 'dz_list', 'dr2_list', 'dr_list'],
							type_converters = converters.blitz,
							compiler='gcc')
			
	
	#flog = open("../aaa.log", "w")
	
	# add bonds
	for i, j in zip(bonded_A, bonded_B):
		ii = int(i)-1
		jj = int(j)-1
	#	if shake_on == 0:
		if k_bond > force_const[ii][jj]:
			force_const[ii][jj] = k_bond
			force_const[jj][ii] = k_bond
	for cent in range(len(all_centroids)-1):
		ii = int(all_centroids[cent]-1)
		jj = int(all_centroids[cent])
		if k_bond > force_const[ii][jj]:
			force_const[ii][jj] = 100
			force_const[jj][ii] = 100
		

	# add 1-3 constant
	ci = 0
	for i in range(len(bonded_A)-1):
		ii  = int(bonded_A[ci])-1
		jj  = int(bonded_B[ci])-1
		jj2 = int(bonded_A[ci+1])-1
		kk  = int(bonded_B[ci+1])-1
		if jj == jj2:
			if k_bond_1_3 > force_const[ii][kk]:
				force_const[ii][kk] = k_bond_1_3
				force_const[kk][ii] = k_bond_1_3
		ci += 1

	# add 1-4 constant
	ci = 0
	for i in range(len(bonded_A)-2):
		ii  = int(bonded_A[ci])-1
		jj  = int(bonded_B[ci])-1
		jj2 = int(bonded_A[ci+1])-1
		kk  = int(bonded_B[ci+1])-1
		kk2 = int(bonded_A[ci+2])-1
		ll  = int(bonded_B[ci+2])-1
		if jj == jj2 and kk == kk2:
			if k_bond_1_4 > force_const[ii][ll]:
				force_const[ii][ll] = k_bond_1_4
				force_const[ll][ii] = k_bond_1_4
		ci += 1

	# add 1-4 attraction in helices
	ci = 0
	for i in range(len(bonded_A)-2):
		ii  = int(bonded_A[ci])-1
		if secondary_structure[ii] == "H":
			jj  = int(bonded_B[ci])-1
			jj2 = int(bonded_A[ci+1])-1
			kk  = int(bonded_B[ci+1])-1
			kk2 = int(bonded_A[ci+2])-1
			ll  = int(bonded_B[ci+2])-1
			if jj == jj2 and kk == kk2:
				if secondary_structure[jj] == "H" and secondary_structure[kk] == "H" and secondary_structure[ll] == "H":
					if k_helix_1_4 > force_const[ii][ll]:
						force_const[ii][ll] = k_helix_1_4
						force_const[ll][ii] = k_helix_1_4
		ci += 1

	# add attraction in sheets
	ci = 0
	for i in range(len(bonded_A)-2):
		ii  = int(bonded_A[ci])-1
		jj  = ii + int(int_resn[ii])
		if sheet_resn[ii] == "E":
			if k_sheet_1_4 > force_const[ii][jj]:
				force_const[ii][jj] = k_sheet_1_4
				force_const[jj][ii] = k_sheet_1_4
		ci += 1
	# add disulfides
	for i, j in zip(disulfide_A, disulfide_B):
		ii = int(i)-1
		jj = int(j)-1
	#	if shake_on == 0:
		if k_disulfide > force_const[ii][jj]:
			force_const[ii][jj] = k_disulfide
			force_const[jj][ii] = k_disulfide
	
	# add salt bridges
	for i, j in zip(saltbridge_A, saltbridge_B):
		ii = int(i)-1
		jj = int(j)-1
		if ii >= all_centroids[0]-1 or jj >= all_centroids[0]-1:
			force_const[ii][jj] = 10
			force_const[jj][ii] = 10
			continue
		if k_salt > force_const[ii][jj]:
			force_const[ii][jj] = k_salt
			force_const[jj][ii] = k_salt


	# add H-bonds
	for i, j in zip(hbonds_A, hbonds_B):
		ii = int(i)-1
		jj = int(j)-1
		if ii >= all_centroids[0]-1 or jj >= all_centroids[0]-1:
			force_const[ii][jj] = 10
			force_const[jj][ii] = 10
			continue
		if k_hbonds > force_const[ii][jj]:
			force_const[ii][jj] = k_hbonds
			force_const[jj][ii] = k_hbonds

	# add ligand and protein aromatic bonds
	for i, j in zip(aromatic_A, aromatic_B):
		ii = int(i)-1
		jj = int(j)-1
	        if k_aromatic > force_const[ii][jj]:
			force_const[ii][jj] = k_aromatic
			force_const[jj][ii] = k_aromatic

	# add hphob
	for i, j in zip(hphob_A, hphob_B):
		ii = int(i)-1
		jj = int(j)-1
		if k_hphob > force_const[ii][jj]:
			force_const[ii][jj] = k_hphob
			force_const[jj][ii] = k_hphob

	# add density
	for i in range(num_residues):
		for j in range(i+1, num_residues):
			if dr_list[i][j] <= r_cutoff:
				if k_density > force_const[i][j]:
					force_const[i][j] = k_density
					force_const[j][i] = k_density
	
#	for i in range(num_residues):
#		for j in range(num_residues):
#		print force_const[i]
#			print "% 7.2f"%force_const[i][j]
#			flog.write("% 7.2f " % force_const[i][j])
#		print "\n"
#		flog.write("\n")
		
#	flog.close()
###################################################################################################
###################################################################################################
def correctBondDeviations(num_residues, coord_list, norm_v, new_coord_list, bonded_A, bonded_B, atom_shift, ang_dev, violations, dec):

	dec = float(dec)
	#print len(bonded_A), bonded_A 
	#print len(bonded_B), bonded_B
	k = 0
	dr_list = scipy.zeros(len(bonded_A))
	for i,j in zip(bonded_A,bonded_B):
		ii = int(i)-1
		jj = int(j)-1
#		print ii, jj
		x1 = new_coord_list[3*ii]
		y1 = new_coord_list[3*ii+1]
		z1 = new_coord_list[3*ii+2]
		x2 = new_coord_list[3*jj]
		y2 = new_coord_list[3*jj+1]
		z2 = new_coord_list[3*jj+2]
		dx = x1 - x2
		dy = y1 - y2
		dz = z1 - z2
		dr2 = dx*dx + dy*dy + dz*dz
		dr = math.sqrt(dr2)
		dr_list[k] = dr
		k += 1

	#print dr_list
	# check if bond length of any pair of atoms > 5.0
	bond_devi = []
	for d in range(len(dr_list)):
		if dr_list[d] > 4.0:
			bond_devi.append(d)
	#print "bond_deviations\t", bond_devi
	
	#check if movement of atoms > 4.5 ALSO CHECK IF 4.5 CORRESPONDS TO BAD BOND LENGTHS IN DIFFERENT SYSTEMS
	tmp_atom_shift_viol = []
	for s in range(len(atom_shift)-1):
		 if atom_shift[s] > 5.0:						
			tmp_atom_shift_viol.append(s)

	# adjacent atom shift may lead to large shift in atoms
	atom_shift_viol = []
	if tmp_atom_shift_viol != []:
		for n in range(num_residues):
			if n in tmp_atom_shift_viol and n in bond_devi:
				atom_shift_viol.append(n)
			if n in tmp_atom_shift_viol and (n-1 in bond_devi or n+1 in bond_devi):
				atom_shift_viol.append(n)
		if tmp_atom_shift_viol[-1] == num_residues:
			atom_shift_viol.append(tmp_atom_shift_viol[-1])
		atom_shift_viol = list(set(atom_shift_viol))
		atom_shift_viol.sort()
		#print atom_shift_viol

	#print "atom_shift_deviations\t", tmp_atom_shift_viol
	#print "ang_deviations\t", ang_dev
	if 0 in bond_devi:
		violations.append(0)
	for i in range(1,num_residues):
		if i in bond_devi:
			violations.append(i)
			#if i in atom_shift_viol or i in ang_dev:
			#	violations.append(i)
	max_viol_allowd = len(violations)*0.05
	#print "violations", violations
	
	#return violations
	#THRID APPROACH
	"""
	# Making subsets of connected atoms which have violations
	atom_set = []
	atom_subset = []
	#for s in range(1,len(ang_dev)):
	for s in range(1,len(violations)):
		if ang_dev[s]-ang_dev[s-1] == 1:
			if ang_dev[s-1] not in atom_subset:
				atom_subset.append(ang_dev[s-1])
			if ang_dev[s] not in atom_subset:
				atom_subset.append(ang_dev[s])
		else:
			atom_set.append((atom_subset))
			atom_subset = []
			if ang_dev[s] not in atom_subset:
				atom_subset.append(ang_dev[s])
	atom_set.append((atom_subset))
	atom_subset = []
	print atom_set
	"""

	"""
	for t in range(len(atom_set)):
		atom_subset=[]
		for c in atom_set[t]:
			atom_subset.append(c)
		if atom_subset[0] == 0:
			continue
			
			next = atom_subset[-1]+1 # next = last+1 residue in the first set of atoms
			print next
			atom_subset.reverse()
			for c in range(len(atom_subset)): # fix the movements in reversed order
				i = atom_subset[c]
				print new_coord_list[3*i], new_coord_list[3*i+1], new_coord_list[3*i+2]
				new_coord_list[3*i] = coord_list[3*i] + norm_v[3*next]*math.sqrt(num_residues)*6.50
				new_coord_list[3*i+1] = coord_list[3*i+1] + norm_v[3*next+1]*math.sqrt(num_residues)*6.50
				new_coord_list[3*i+2] = coord_list[3*i+2] + norm_v[3*next+2]*math.sqrt(num_residues)*6.50
				#next = i
				print next, new_coord_list[3*i], new_coord_list[3*i+1], new_coord_list[3*i+2]
			
		elif atom_subset[-1] == num_residues-1:
			continue

		else:
			for c in range(len(atom_subset)-1):
				i = atom_subset[c]
	"""		
			

	# SECOND APPROACH
	#check which set of joined atoms have atom movements > 4.5 
	
	atom_subset=[]
	atom_set=[]
	# Making subsets of connected atoms which have violations
	"""
	for s in range(1,len(atom_shift_viol)):
		if atom_shift_viol[s]-atom_shift_viol[s-1] == 1:
			if atom_shift_viol[s-1] not in atom_subset:
				atom_subset.append(atom_shift_viol[s-1])
			if atom_shift_viol[s] not in atom_subset:
				atom_subset.append(atom_shift_viol[s])
		else:
			atom_set.append((atom_subset))
			atom_subset = []
			if atom_shift_viol[s] not in atom_subset:
				atom_subset.append(atom_shift_viol[s])
	"""
	flag = 0
	if len(violations) == 0:
		return violations, new_coord_list
	elif len(violations) == 1:
		atom_set.append((violations))					
	else:
		atom_subset.append((violations[0]))	
		for s in range(1,len(violations)): #e.g. [0,292], [0,1,2], [0,1,3,5], [0,1,3,4], [0,2,3]
			if violations[s]-violations[s-1] == 1:
				if violations[s-1] not in atom_subset:
					atom_subset.append(violations[s-1])
				if violations[s] not in atom_subset:
					atom_subset.append(violations[s])
				flag = 1
			else:
				flag = 0
				if atom_subset != []:
					atom_set.append((atom_subset))
				atom_subset = []
				#atom_subset.append(violations[s])
				if violations[s] not in atom_subset:
					atom_subset.append(violations[s])
					flag = 1
					continue
				
	if flag == 1:
		atom_set.append((atom_subset))
	atom_subset = []
	#print "atom_set", atom_set
	violations = []
	mark = 0
	end1 = 0
	end2 = 0
	if atom_set[-1][-1] == num_residues-2 and atom_set[0][0] == 0:
		# include bothe ends in violations
		mark = 2
		end1 = len(atom_set[0])
		end2 = len(atom_set[-1])
	elif atom_set[-1][-1] == num_residues-2 and atom_set[0][0] != 0:
		# include only one end in violations
		mark = 1
		end2 = len(atom_set[-1])
	elif atom_set[-1][-1] != num_residues-2 and atom_set[0][0] == 0:
		mark = 1
		end1 = len(atom_set[0])
	else:
		pass
	# For this to work, protein must have one atom that has correct bond lengths
	for t in range(len(atom_set)):
		atom_subset = []
		if atom_set[t][0] == 0:
			for c in atom_set[t]:
				atom_subset.append(c) 
			atom_subset.append(int(atom_set[t][-1])+1)
		elif atom_set[t][-1] == num_residues-2:
			atom_subset.append(int(atom_set[t][0])-1)
			for c in atom_set[t]:
				atom_subset.append(c)
		else:
			atom_subset.append(int(atom_set[t][0])-1)
			for c in atom_set[t]:
				atom_subset.append(c)
			atom_subset.append(int(atom_set[t][-1])+1)
		
		if atom_subset[0] == 0: #check if first set of atoms belong to beginning residues
			#atom_subset.append(int(atom_subset[-1])+1)
			#print "first type\n"
			#print len(atom_subset)
			check = 0
			for c in range(len(atom_subset)):
				if mark == 2:
					if end1 <= max_viol_allowd*0.5 and end2 <= max_viol_allowd*0.5 :
						violations.append(c)
					elif end1 > max_viol_allowd*0.5 and end2 <= max_viol_allowd*0.5 :
						if check < max_viol_allowd*0.5+(max_viol_allowd*0.5-end2):
							violations.append(c)				
					elif end1 <= max_viol_allowd*0.5 and end2 > max_viol_allowd*0.5 :
						violations.append(c)
					else:
						if check < max_viol_allowd*0.5:
							violations.append(c)
					check += 1
				if mark == 1:	
					if end1 != 0 and end1 <= max_viol_allowd:
						violations.append(c)
					elif end1 != 0 and end1 > max_viol_allowd:
						if check < max_viol_allowd:
							violations.append(c)
					else:
						pass
					check += 1
			#continue
			next = atom_subset[-1]+1 # next = last+1 residue in the first set of atoms
			#print next
			atom_subset.reverse()
			for c in range(len(atom_subset)-1): # fix the movements in reversed order
				i = atom_subset[c] 
				bond_len_prev = 9999.99
				bond_len = dr_list[i]
				iter1 = 0
				if bond_len <= 4.0 and bond_len > 3.6:
					#print "i", i, bond_len
					tmp_new_coord_list_x1 = new_coord_list[3*i]
					tmp_new_coord_list_y1 = new_coord_list[3*i+1]
					tmp_new_coord_list_z1 = new_coord_list[3*i+2]
					while bond_len <= 4.0 and bond_len > 3.6 and iter1 < 10:
						iter1 += 1
						iter2 = 0
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec
						dx = tmp_new_coord_list_x1 - new_coord_list[3*next]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*next+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*next+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						#if bond_len <= 4.0 and bond_len > 3.6: # else: we don't change anything about i 
						#	pass
							#print atom_subset[c+1]
						if bond_len > 4.0 or bond_len <= 3.6:
						#else:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec		#CORRECT
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
							dx = tmp_new_coord_list_x1 - new_coord_list[3*next]
							dy = tmp_new_coord_list_y1 - new_coord_list[3*next+1]
							dz = tmp_new_coord_list_z1 - new_coord_list[3*next+2]
							bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						if atom_subset[c] == atom_subset[-1]:
							continue
						else:
							#print atom_subset[c+1]
							j = atom_subset[c+1] # add condition for checking bond length later
							tmp_new_coord_list_x2 = new_coord_list[3*j]
							tmp_new_coord_list_y2 = new_coord_list[3*j+1]
							tmp_new_coord_list_z2 = new_coord_list[3*j+2]
							#bond_len2 = dr_list[j]
							bond_len2_prev = 9999.99
						 
							dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
							dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
							dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
							bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
							
							if bond_len2 <= 4.0 and bond_len2 > 3.6:
								#print "j", j, bond_len2, "break1"
								new_coord_list[3*j] = tmp_new_coord_list_x2
								new_coord_list[3*j+1] = tmp_new_coord_list_y2
								new_coord_list[3*j+2] = tmp_new_coord_list_z2
								dr_list[j] = bond_len2
								break	
							else:
								while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 10:
									iter2 += 1
									if bond_len2_prev < bond_len2:
										#print "j", j, bond_len2, "break2"
										break
									bond_len2_prev = bond_len2
									tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
									dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
									dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
									dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
									bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
									#if bond_len2 <= 4.0 and bond_len2 > 3.6: # NOT REQUIRED
									#	break
								if bond_len2_prev < bond_len2:
									new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
									dr_list[j] = bond_len2_prev
								else:
								#if bond_len2 <= 4.0 and bond_len2 > 3.6:
									new_coord_list[3*j] = tmp_new_coord_list_x2
									new_coord_list[3*j+1] = tmp_new_coord_list_y2
									new_coord_list[3*j+2] = tmp_new_coord_list_z2
									dr_list[j] = bond_len2
									break
								#else:
								#	new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
								#	new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
								#	new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
								#	dr_list[j] = bond_len2_prev	
									
					new_coord_list[3*i] = tmp_new_coord_list_x1
					new_coord_list[3*i+1] = tmp_new_coord_list_y1
					new_coord_list[3*i+2] = tmp_new_coord_list_z1
					dr_list[i] = bond_len
					next = i
					#continue
				else:
					#print "i", i, bond_len
					tmp_new_coord_list_x1 = new_coord_list[3*i]
					tmp_new_coord_list_y1 = new_coord_list[3*i+1]
					tmp_new_coord_list_z1 = new_coord_list[3*i+2]
					while (bond_len > 4.0 or bond_len <= 3.6) and iter1 < 10:
						iter1 += 1
						iter2 = 0
						#print bond_len
						if bond_len_prev < bond_len:
							#print "i", i, bond_len, "break3"
							break
						bond_len_prev = bond_len
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec
						dx = tmp_new_coord_list_x1 - new_coord_list[3*next]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*next+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*next+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						#if bond_len <= 4.0 and bond_len > 3.6:  
							#print atom_subset[c+1]
						#	pass
						#else: #CHECK AGAIN
						if bond_len_prev < bond_len:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
							#dx = tmp_new_coord_list_x1 - new_coord_list[3*next]
							#dy = tmp_new_coord_list_y1 - new_coord_list[3*next+1]
							#dz = tmp_new_coord_list_z1 - new_coord_list[3*next+2]
							#bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						if atom_subset[c] == atom_subset[-1]:
							continue
						else:
							#if atom_subset[c+1] >= atom_subset[-1]:
							j = atom_subset[c+1]
							bond_len2_prev = 9999.99
							tmp_new_coord_list_x2 = new_coord_list[3*j]
							tmp_new_coord_list_y2 = new_coord_list[3*j+1]
							tmp_new_coord_list_z2 = new_coord_list[3*j+2]
							dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
							dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
							dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
							bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
							if bond_len2 <= 4.0 and bond_len2 > 3.6:
								#print "j ",j, bond_len2, "break4"
								new_coord_list[3*j] = tmp_new_coord_list_x2
								new_coord_list[3*j+1] = tmp_new_coord_list_y2
								new_coord_list[3*j+2] = tmp_new_coord_list_z2
								dr_list[j] = bond_len2
								#break
							else:
								while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 10:
									iter2 += 1
									if bond_len2_prev < bond_len2:
										#print "j ",j, bond_len2, "break5"
										break
									bond_len2_prev = bond_len2
									tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
									dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
									dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
									dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
									bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
									#if bond_len2 <= 4.0 and bond_len2 > 3.6: # NOT REQUIRED
									#	break
								if bond_len2_prev < bond_len2:
									new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
									dr_list[j] = bond_len2_prev
								#if bond_len2 <= 4.0 and bond_len2 > 3.6:
								else:
									new_coord_list[3*j] = tmp_new_coord_list_x2
									new_coord_list[3*j+1] = tmp_new_coord_list_y2
									new_coord_list[3*j+2] = tmp_new_coord_list_z2
									dr_list[j] = bond_len2
								#else:
								#	new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
								#	new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
								#	new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
								#	dr_list[j] = bond_len2_prev
					#if bond_len <= 4.0 and bond_len > 3.6:
					if bond_len_prev < bond_len:
						dr_list[i] = bond_len_prev
					else:
						new_coord_list[3*i] = tmp_new_coord_list_x1
						new_coord_list[3*i+1] = tmp_new_coord_list_y1
						new_coord_list[3*i+2] = tmp_new_coord_list_z1
						dr_list[i] = bond_len
					next = i
					#else:				
					#	new_coord_list[3*i] = tmp_new_coord_list_x1 + norm_v[3*i]*math.sqrt(num_residues)*dec*-1.0
					#	new_coord_list[3*i+1] = tmp_new_coord_list_y1 + norm_v[3*i+1]*math.sqrt(num_residues)*dec*-1.0
					#	new_coord_list[3*i+2] = tmp_new_coord_list_z1 + norm_v[3*i+2]*math.sqrt(num_residues)*dec*-1.0
					#	dr_list[i] = bond_len_prev

		elif atom_subset[-1] == num_residues-2 :#check if the last set of atoms belong to the end terminal
			#atom_subset.append(int(atom_subset[0])-1)
			#print "Second type\n"
			#print len(atom_subset)
			check = 0
			for c in atom_subset:
				if mark == 2:
					if end1 <= max_viol_allowd*0.5 and end2 <= max_viol_allowd*0.5 :
						violations.append(c+1)
					elif end1 > max_viol_allowd*0.5 and end2 <= max_viol_allowd*0.5 :
						violations.append(c+1)
					elif end1 <= max_viol_allowd*0.5 and end2 > max_viol_allowd*0.5 :
						if check < max_viol_allowd*0.5 + (max_viol_allowd*0.5-end1):
							violations.append(c+1)
					else:
						if check < max_viol_allowd*0.5:
							violations.append(c+1)
                                        check += 1
                                if mark == 1:
                                        if end2 != 0 and end2 <= max_viol_allowd:
                                                violations.append(c+1)
                                        elif end2 != 0 and end2 > max_viol_allowd:
                                                if check < max_viol_allowd:
                                                        violations.append(c+1)
                                        else:
                                                pass
                                        check += 1
#			continue
			prev = atom_subset[0]-1 # prev = first-1 residue in the last set of atoms
			for c in range(len(atom_subset)-1): # fix the movements in reversed order
				i = atom_subset[c] 
				bond_len_prev = 9999.99
				bond_len = dr_list[prev]	#bond_len = dr_list[i]
				iter1 = 0
				if bond_len <= 4.0 and bond_len > 3.6:
					#print "i", i, bond_len
					tmp_new_coord_list_x1 = new_coord_list[3*i]
					tmp_new_coord_list_y1 = new_coord_list[3*i+1]
					tmp_new_coord_list_z1 = new_coord_list[3*i+2]
					while bond_len <= 4.0 and bond_len > 3.6 and iter1 < 10:
						iter1 += 1
						iter2 = 0
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec
						dx = tmp_new_coord_list_x1 - new_coord_list[3*prev]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*prev+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*prev+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						#if bond_len <= 4.0 and bond_len > 3.6:  
							#print atom_subset[c+1]
						#	pass
						if bond_len > 4.0 or bond_len <= 3.6:
						#else:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
							dx = tmp_new_coord_list_x1 - new_coord_list[3*prev]
							dy = tmp_new_coord_list_y1 - new_coord_list[3*prev+1]
							dz = tmp_new_coord_list_z1 - new_coord_list[3*prev+2]
							bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						if atom_subset[c] == atom_subset[-1]:
							continue
						else:
						#if atom_subset[c+1] != atom_subset[-1]:
							#print atom_subset[c+1]
							j = atom_subset[c+1] # add condition for checking bond length later
							tmp_new_coord_list_x2 = new_coord_list[3*j]
							tmp_new_coord_list_y2 = new_coord_list[3*j+1]
							tmp_new_coord_list_z2 = new_coord_list[3*j+2]
#							bond_len2 = dr_list[i]	#bond_len2 = dr_list[j]
							bond_len2_prev = 9999.99
						 
							dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
							dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
							dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
							bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
							
							if bond_len2 <= 4.0 and bond_len2 > 3.6:
									#print "j", j, bond_len2, "break1"
									new_coord_list[3*j] = tmp_new_coord_list_x2
									new_coord_list[3*j+1] = tmp_new_coord_list_y2
									new_coord_list[3*j+2] = tmp_new_coord_list_z2
									dr_list[i] = bond_len2
									break
							else:
								while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 10:
									iter2 += 1
									if bond_len2_prev < bond_len2:
										#print "j", j, bond_len2, "break2"
										break
									bond_len2_prev = bond_len2
									tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
									dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
									dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
									dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
									bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
									if bond_len2 <= 4.0 and bond_len2 > 3.6: # NOT REQUIRED
										break
								if bond_len2_prev < bond_len2:
									new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
									dr_list[i] = bond_len2_prev
								#if bond_len2 <= 4.0 and bond_len2 > 3.6:
								else:
									new_coord_list[3*j] = tmp_new_coord_list_x2
									new_coord_list[3*j+1] = tmp_new_coord_list_y2
									new_coord_list[3*j+2] = tmp_new_coord_list_z2
									dr_list[i] = bond_len2
								#else:
								#	new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
								#	new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
								#	new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
								#	dr_list[i] = bond_len2_prev	
					new_coord_list[3*i] = tmp_new_coord_list_x1
					new_coord_list[3*i+1] = tmp_new_coord_list_y1
					new_coord_list[3*i+2] = tmp_new_coord_list_z1
					dr_list[prev] = bond_len
					prev = i
					#continue
				else:
					#print "i", i, bond_len
					tmp_new_coord_list_x1 = new_coord_list[3*i]
					tmp_new_coord_list_y1 = new_coord_list[3*i+1]
					tmp_new_coord_list_z1 = new_coord_list[3*i+2]
					while (bond_len > 4.0 or bond_len <= 3.6) and iter1 < 10:
						iter1 += 1
						iter2 = 0
						#print bond_len
						if bond_len_prev < bond_len:
							#print "i", i, bond_len, "break3"
							break
						bond_len_prev = bond_len
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec
				
						dx = tmp_new_coord_list_x1 - new_coord_list[3*prev]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*prev+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*prev+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						#if bond_len <= 4.0 and bond_len > 3.6: # else: we don't change anything about i 
							#print atom_subset[c+1]
						#	pass
						#else:
						if bond_len_prev < bond_len:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
							#dx = tmp_new_coord_list_x1 - new_coord_list[3*prev]
							#dy = tmp_new_coord_list_y1 - new_coord_list[3*prev+1]
							#dz = tmp_new_coord_list_z1 - new_coord_list[3*prev+2]
							#bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						if atom_subset[c] == atom_subset[-1]:
							continue
						else:
							j = atom_subset[c+1]
							tmp_new_coord_list_x2 = new_coord_list[3*j]
							tmp_new_coord_list_y2 = new_coord_list[3*j+1]
							tmp_new_coord_list_z2 = new_coord_list[3*j+2]
							#bond_len2 = dr_list[j]
							bond_len2_prev = 9999.99
							
							dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
							dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
							dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
							bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
							if bond_len2 <= 4.0 and bond_len2 > 3.6:
								#print "j", j, bond_len2, "break4"
								new_coord_list[3*j] = tmp_new_coord_list_x2
								new_coord_list[3*j+1] = tmp_new_coord_list_y2
								new_coord_list[3*j+2] = tmp_new_coord_list_z2
								dr_list[i] = bond_len2
								#break
							else:
								while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 10:
									iter2 += 1
									if bond_len2_prev < bond_len2:
										#print "j", j, bond_len2, "break5"
										break
									bond_len2_prev = bond_len2
									tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
									dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
									dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
									dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
									bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
								if bond_len2_prev < bond_len2:
									new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
									dr_list[i] = bond_len2_prev
								#if bond_len2 <= 4.0 and bond_len2 > 3.6:		
								else:
									new_coord_list[3*j] = tmp_new_coord_list_x2
									new_coord_list[3*j+1] = tmp_new_coord_list_y2
									new_coord_list[3*j+2] = tmp_new_coord_list_z2
									dr_list[i] = bond_len2
								#else:
								#	new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
								#	new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
								#	new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
								#	dr_list[i] = bond_len2_prev			
					#if bond_len <= 4.0 and bond_len > 3.6:
					if bond_len_prev < bond_len:
						dr_list[prev] = bond_len_prev
					else:
						new_coord_list[3*i] = tmp_new_coord_list_x1
						new_coord_list[3*i+1] = tmp_new_coord_list_y1
						new_coord_list[3*i+2] = tmp_new_coord_list_z1
						dr_list[prev] = bond_len
					#else:
					#	new_coord_list[3*i] = tmp_new_coord_list_x1 + norm_v[3*i]*math.sqrt(num_residues)*dec*-1.0
					#	new_coord_list[3*i+1] = tmp_new_coord_list_y1 + norm_v[3*i+1]*math.sqrt(num_residues)*dec*-1.0
					#	new_coord_list[3*i+2] = tmp_new_coord_list_z1 + norm_v[3*i+2]*math.sqrt(num_residues)*dec*-1.0
					#	dr_list[prev] = bond_len_prev
					prev = i
		else:
			print "Third type\n"	# No rule about atom_subset of length 1
			print atom_subset
			#for xxx in range(len(atom_subset)):
			#	print dr_list[atom_subset[xxx]]
			#atom_subset.insert(0,int(atom_subset[0])-1)
			#atom_subset.append(int(atom_subset[-1])+1)
			prev = atom_subset[0]-1
			next = atom_subset[-1]+1
			if atom_shift[prev] < atom_shift[next]:
				direction = "forward" # move forwards
				direc = prev
				print "move forwards"

			else:
				direction = "backward" # move backwards
				direc = next
				atom_subset.reverse()
				print "move backwards"


			#print atom_subset
			for c in range(len(atom_subset)-1):
				i = atom_subset[c] 
				bond_len_prev = 9999.99
				tmp_new_coord_list_x1 = new_coord_list[3*i]
				tmp_new_coord_list_y1 = new_coord_list[3*i+1]
				tmp_new_coord_list_z1 = new_coord_list[3*i+2]
				dx = tmp_new_coord_list_x1 - new_coord_list[3*direc]
				dy = tmp_new_coord_list_y1 - new_coord_list[3*direc+1]
				dz = tmp_new_coord_list_z1 - new_coord_list[3*direc+2]
				bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
				#print "i",i,bond_len, tmp_new_coord_list_x1, tmp_new_coord_list_y1, tmp_new_coord_list_z1
				iter1= 0
				# if bond_len < 3.6 STOP!!!!!!!!
				#print "i", i, bond_len
				if bond_len <= 4.0 and bond_len > 3.6:
					while (bond_len <= 4.0 and bond_len > 3.6) and iter1 < 25:
						iter1 += 1
						iter2=0
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec 
						dx = tmp_new_coord_list_x1 - new_coord_list[3*direc]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*direc+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*direc+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
#						print direc,bond_len,tmp_new_coord_list_x1,tmp_new_coord_list_y1,tmp_new_coord_list_z1
						#if bond_len <= 4.0 and bond_len > 3.6:  
						#	pass
						#	print direc,bond_len,tmp_new_coord_list_x1,tmp_new_coord_list_y1,tmp_new_coord_list_z1
						if bond_len > 4.0 or bond_len <= 3.6:
						#else:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
							dx = tmp_new_coord_list_x1 - new_coord_list[3*direc]
							dy = tmp_new_coord_list_y1 - new_coord_list[3*direc+1]
							dz = tmp_new_coord_list_z1 - new_coord_list[3*direc+2]
							bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						#	print direc,bond_len,tmp_new_coord_list_x1,tmp_new_coord_list_y1,tmp_new_coord_list_z1
															
						j = atom_subset[c+1]
						bond_len2_prev = 9999.99
						tmp_new_coord_list_x2 = new_coord_list[3*j]
						tmp_new_coord_list_y2 = new_coord_list[3*j+1]
						tmp_new_coord_list_z2 = new_coord_list[3*j+2]
						dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
						dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
						dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
						bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
						#print i,bond_len2, tmp_new_coord_list_x2, tmp_new_coord_list_y2, tmp_new_coord_list_z2
						if bond_len2 <= 4.0 and bond_len2 > 3.6:
							#print "j", j, bond_len2, "break1"
							new_coord_list[3*j] = tmp_new_coord_list_x2
							new_coord_list[3*j+1] = tmp_new_coord_list_y2
							new_coord_list[3*j+2] = tmp_new_coord_list_z2
							if direction == "backward":
								dr_list[j] = bond_len2
								if j == atom_subset[-1]:
									dx = tmp_new_coord_list_x2 - new_coord_list[3*prev]
									dy = tmp_new_coord_list_y2 - new_coord_list[3*prev+1]
									dz = tmp_new_coord_list_z2 - new_coord_list[3*prev+2]
									bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
									dr_list[prev] = bond_len3
									#print "n1", prev, bond_len3
							if direction == "forward":
								dr_list[i] = bond_len2
								if j == atom_subset[-1]:
									dx = tmp_new_coord_list_x2 - new_coord_list[3*next]
									dy = tmp_new_coord_list_y2 - new_coord_list[3*next+1]
									dz = tmp_new_coord_list_z2 - new_coord_list[3*next+2]
									bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
									dr_list[j] = bond_len3
									#print "n2", next, bond_len3
							break 
						elif bond_len2 <= 3.6:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
							dx = tmp_new_coord_list_x1 - new_coord_list[3*direc]
							dy = tmp_new_coord_list_y1 - new_coord_list[3*direc+1]
							dz = tmp_new_coord_list_z1 - new_coord_list[3*direc+2]
							bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)

							dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
							dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
							dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
							bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
							if direction == "backward":
								dr_list[j] = bond_len2
								if j == atom_subset[-1]:
									dr_list[prev] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*prev])*(new_coord_list[3*j]-new_coord_list[3*prev]) + (new_coord_list[3*j+1]-new_coord_list[3*prev+1])*(new_coord_list[3*j+1]-new_coord_list[3*prev+1]) + (new_coord_list[3*j+2]-new_coord_list[3*prev+2])*(new_coord_list[3*j+2]-new_coord_list[3*prev+2]))
							if direction == "forward":
								dr_list[i] = bond_len2_prev
								if j == atom_subset[-1]:
									dr_list[j] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*next])*(new_coord_list[3*j]-new_coord_list[3*next]) + (new_coord_list[3*j+1]-new_coord_list[3*next+1])*(new_coord_list[3*j+1]-new_coord_list[3*next+1]) + (new_coord_list[3*j+2]-new_coord_list[3*next+2])*(new_coord_list[3*j+2]-new_coord_list[3*next+2]))
							
							break
						else:
#							while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 25:
							while bond_len2 > 4.0 and iter2 < 25:
								iter2 += 1
								if bond_len2_prev < bond_len2:
									#print "j", j, bond_len2, bond_len2_prev, "break2"
									#flag = 1
									break  
								tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
								tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
								tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
								bond_len2_prev = bond_len2
								dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
								dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
								dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
								bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
								#print i,bond_len2, tmp_new_coord_list_x2, tmp_new_coord_list_y2, tmp_new_coord_list_z2

							if bond_len2_prev < bond_len2:
								new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
								if direction == "backward":
									dr_list[j] = bond_len2_prev
									if j == atom_subset[-1]:
										dr_list[prev] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*prev])*(new_coord_list[3*j]-new_coord_list[3*prev]) + (new_coord_list[3*j+1]-new_coord_list[3*prev+1])*(new_coord_list[3*j+1]-new_coord_list[3*prev+1]) + (new_coord_list[3*j+2]-new_coord_list[3*prev+2])*(new_coord_list[3*j+2]-new_coord_list[3*prev+2]))
										#print "n3", prev, dr_list[prev]
								if direction == "forward":
									dr_list[i] = bond_len2_prev
									if j == atom_subset[-1]:
										dr_list[j] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*next])*(new_coord_list[3*j]-new_coord_list[3*next]) + (new_coord_list[3*j+1]-new_coord_list[3*next+1])*(new_coord_list[3*j+1]-new_coord_list[3*next+1]) + (new_coord_list[3*j+2]-new_coord_list[3*next+2])*(new_coord_list[3*j+2]-new_coord_list[3*next+2]))		
										#print "n4", next, dr_list[j]
							#if bond_len2 <= 4.0 and bond_len2 > 3.6: #bond_len3 <= 4.0 and bond_len3 > 3.6
							else:
								new_coord_list[3*j] = tmp_new_coord_list_x2
								new_coord_list[3*j+1] = tmp_new_coord_list_y2
								new_coord_list[3*j+2] = tmp_new_coord_list_z2
								if direction == "backward":
									dr_list[j] = bond_len2
									if j == atom_subset[-1]:
										dx = tmp_new_coord_list_x2 - new_coord_list[3*prev]
										dy = tmp_new_coord_list_y2 - new_coord_list[3*prev+1]
										dz = tmp_new_coord_list_z2 - new_coord_list[3*prev+2]
										bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
										dr_list[prev] = bond_len3
										#print "n5", prev, bond_len3
								if direction == "forward":
									dr_list[i] = bond_len2
									if j == atom_subset[-1]:
										dx = tmp_new_coord_list_x2 - new_coord_list[3*next]
										dy = tmp_new_coord_list_y2 - new_coord_list[3*next+1]
										dz = tmp_new_coord_list_z2 - new_coord_list[3*next+2]
										bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
										dr_list[j] = bond_len3
										#print "n6", next, bond_len3
							"""
							else:
								new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
								if direction == "backward":
									dr_list[j] = bond_len2_prev
								
									if j == atom_subset[-1]:
										dr_list[prev] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*prev])*(new_coord_list[3*j]-new_coord_list[3*prev]) + (new_coord_list[3*j+1]-new_coord_list[3*prev+1])*(new_coord_list[3*j+1]-new_coord_list[3*prev+1]) + (new_coord_list[3*j+2]-new_coord_list[3*prev+2])*(new_coord_list[3*j+2]-new_coord_list[3*prev+2]))

								if direction == "forward":
									dr_list[i] = bond_len2_prev
									if j == atom_subset[-1]:
										dr_list[j] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*next])*(new_coord_list[3*j]-new_coord_list[3*next]) + (new_coord_list[3*j+1]-new_coord_list[3*next+1])*(new_coord_list[3*j+1]-new_coord_list[3*next+1]) + (new_coord_list[3*j+2]-new_coord_list[3*next+2])*(new_coord_list[3*j+2]-new_coord_list[3*next+2]))
							"""
					new_coord_list[3*i] = tmp_new_coord_list_x1
					new_coord_list[3*i+1] = tmp_new_coord_list_y1
					new_coord_list[3*i+2] = tmp_new_coord_list_z1
					if direction == "backward":
						dr_list[i] = bond_len
					if direction == "forward":
						dr_list[direc] = bond_len
					direc = i
				elif bond_len <= 3.6:
					continue
				else:
					#print "i", i, bond_len, tmp_new_coord_list_x1, tmp_new_coord_list_y1, tmp_new_coord_list_z1
					#tmp_new_coord_list_x1 = new_coord_list[3*i]
					#tmp_new_coord_list_y1 = new_coord_list[3*i+1]
					#tmp_new_coord_list_z1 = new_coord_list[3*i+2]
					#while (bond_len > 4.0 or bond_len <= 3.6) and iter1 < 25:
					while bond_len > 4.0 and iter1 < 25:
						iter1 += 1
						iter2 = 0
						if bond_len_prev < bond_len:
							#print "i", i, bond_len, bond_len_prev, "break3"
							break
						bond_len_prev = bond_len
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec
						dx = tmp_new_coord_list_x1 - new_coord_list[3*direc]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*direc+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*direc+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						#print direc, bond_len, tmp_new_coord_list_x1, tmp_new_coord_list_y1, tmp_new_coord_list_z1
						#if bond_len <= 4.0 and bond_len > 3.6:  
							#print atom_subset[c+1]
							#pass 
						if bond_len_prev < bond_len:
						#else:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
							
															
						j = atom_subset[c+1]
						bond_len2_prev = 9999.99
						tmp_new_coord_list_x2 = new_coord_list[3*j]
						tmp_new_coord_list_y2 = new_coord_list[3*j+1]
						tmp_new_coord_list_z2 = new_coord_list[3*j+2]
						dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
						dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
						dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
						bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
						#print j, bond_len2, tmp_new_coord_list_x2, tmp_new_coord_list_y2, tmp_new_coord_list_z2
						if bond_len2 <= 4.0 and bond_len2 > 3.6:
							#print "j", j, bond_len2, "break4"
							new_coord_list[3*j] = tmp_new_coord_list_x2
							new_coord_list[3*j+1] = tmp_new_coord_list_y2
							new_coord_list[3*j+2] = tmp_new_coord_list_z2
							if direction == "backward":
								dr_list[j] = bond_len2
								if j == atom_subset[-1]:
									dx=tmp_new_coord_list_x2-new_coord_list[3*prev]
									dy=tmp_new_coord_list_y2-new_coord_list[3*prev+1]
									dz=tmp_new_coord_list_z2-new_coord_list[3*prev+2]
									bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
									dr_list[prev] = bond_len3
									#print "n7", prev, bond_len3
							if direction == "forward":
								dr_list[i] = bond_len2
								if j == atom_subset[-1]:
									dx=tmp_new_coord_list_x2-new_coord_list[3*next]
									dy=tmp_new_coord_list_y2-new_coord_list[3*next+1]
									dz=tmp_new_coord_list_z2-new_coord_list[3*next+2]
									bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
									dr_list[j] = bond_len3
									#print "n8", next, bond_len3
							#break
						else:
							while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 25:
								iter2 += 1
								if bond_len2_prev < bond_len2:
									#print "j", j, bond_len2, bond_len2_prev, "break5"
									break
								bond_len2_prev = bond_len2
								tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
								tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
								tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
								dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
								dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
								dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
								bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
								#print j,bond_len2, tmp_new_coord_list_x2, tmp_new_coord_list_y2, tmp_new_coord_list_z2
								
							if bond_len2_prev < bond_len2:
								new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
								if direction == "backward":
									dr_list[j] = bond_len2_prev
									if j == atom_subset[-1]:
										dr_list[prev] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*prev])*(new_coord_list[3*j]-new_coord_list[3*prev]) + (new_coord_list[3*j+1]-new_coord_list[3*prev+1])*(new_coord_list[3*j+1]-new_coord_list[3*prev+1]) + (new_coord_list[3*j+2]-new_coord_list[3*prev+2])*(new_coord_list[3*j+2]-new_coord_list[3*prev+2]))
										#print "n9", prev, dr_list[prev]
								if direction == "forward":
									dr_list[i] = bond_len2_prev
									if j == atom_subset[-1]:
										dr_list[j] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*next])*(new_coord_list[3*j]-new_coord_list[3*next]) + (new_coord_list[3*j+1]-new_coord_list[3*next+1])*(new_coord_list[3*j+1]-new_coord_list[3*next+1]) + (new_coord_list[3*j+2]-new_coord_list[3*next+2])*(new_coord_list[3*j+2]-new_coord_list[3*next+2]))
										#print "n10", next, dr_list[j]
							#if bond_len2 <= 4.0 and bond_len2 > 3.6:
							else:
								new_coord_list[3*j] = tmp_new_coord_list_x2
								new_coord_list[3*j+1] = tmp_new_coord_list_y2
								new_coord_list[3*j+2] = tmp_new_coord_list_z2
								if direction == "backward":
									dr_list[j] = bond_len2
									if j == atom_subset[-1]:
										dx=tmp_new_coord_list_x2-new_coord_list[3*prev]
										dy=tmp_new_coord_list_y2-new_coord_list[3*prev+1]
										dz=tmp_new_coord_list_z2-new_coord_list[3*prev+2]
										bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
										dr_list[prev] = bond_len3
										#print "n11", prev, bond_len3
								if direction == "forward":
									dr_list[i] = bond_len2
									if j == atom_subset[-1]:
										dx = tmp_new_coord_list_x2-new_coord_list[3*next]
										dy = tmp_new_coord_list_y2-new_coord_list[3*next+1]
										dz = tmp_new_coord_list_z2-new_coord_list[3*next+2]
										bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
										dr_list[j] = bond_len3				
										#print "n12", next, bond_len3
							"""
							else:
								new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
								#dr_list[j] = bond_len2_prev
								if direction == "backward":
									dr_list[j] = bond_len2_prev
									if j == atom_subset[-1]:
										dr_list[prev] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*prev])*(new_coord_list[3*j]-new_coord_list[3*prev]) + (new_coord_list[3*j+1]-new_coord_list[3*prev+1])*(new_coord_list[3*j+1]-new_coord_list[3*prev+1]) + (new_coord_list[3*j+2]-new_coord_list[3*prev+2])*(new_coord_list[3*j+2]-new_coord_list[3*prev+2]))
								if direction == "forward":
									dr_list[i] = bond_len2_prev
									if j == atom_subset[-1]:
										dr_list[j] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*next])*(new_coord_list[3*j]-new_coord_list[3*next]) + (new_coord_list[3*j+1]-new_coord_list[3*next+1])*(new_coord_list[3*j+1]-new_coord_list[3*next+1]) + (new_coord_list[3*j+2]-new_coord_list[3*next+2])*(new_coord_list[3*j+2]-new_coord_list[3*next+2]))
							"""
					#if bond_len <= 4.0 and bond_len > 3.6:
					if bond_len_prev < bond_len:
						if direction == "backward":
							dr_list[i] = bond_len_prev
						if direction == "forward":
							dr_list[direc] = bond_len_prev
					else:					
						new_coord_list[3*i] = tmp_new_coord_list_x1
						new_coord_list[3*i+1] = tmp_new_coord_list_y1
						new_coord_list[3*i+2] = tmp_new_coord_list_z1
						if direction == "backward":
							dr_list[i] = bond_len
						if direction == "forward":
							dr_list[direc] = bond_len
					direc = i
					#else:
					#	new_coord_list[3*i] = tmp_new_coord_list_x1 + norm_v[3*i]*math.sqrt(num_residues)*dec*-1.0
					#	new_coord_list[3*i+1] = tmp_new_coord_list_y1 + norm_v[3*i+1]*math.sqrt(num_residues)*dec*-1.0
					#	new_coord_list[3*i+2] = tmp_new_coord_list_z1 + norm_v[3*i+2]*math.sqrt(num_residues)*dec*-1.0
					#	if direction == "backward":
					#		dr_list[i] = bond_len_prev
					#	if direction == "forward":
					#		dr_list[direc] = bond_len_prev
					#	direc = i
					#print dr_list[i], dr_list[j]
		#for xxx in range(len(atom_subset)):
		#	print dr_list[atom_subset[xxx]]

#			print dr_list				
	
	# FIRST APPROACH
	"""				
	#for bd in range(len(bond_dev)):
	for d in reversed(range(len(dr_list)-1)):
		dr_prev = 0.0
		dr_opt = []
		if (dr_list[d] > 4.5) and (dr_list[d+1] <= 4.5 and dr_list[d+1] >= 3.6):
			dr = dr_list[d]
			increment = 0.0
			#dr = 0.0
			while (dr > 4.5 and dr <= dr_list[d]):
				#print d
				x = (coord_list[3*d] + norm_v[3*d]*increment) - new_coord_list_p[3*(d+1)]
				y = (coord_list[3*d+1] + norm_v[3*d+1]*increment) - new_coord_list_p[3*(d+1)+1]
				z = (coord_list[3*d+2] + norm_v[3*d+2]*increment) - new_coord_list_p[3*(d+1)+2]
				#print x,y,z
				dr = math.sqrt(x*x + y*y + z*z)
				dr_opt.append((dr,increment))
				#print dr, increment
				start_dist_x = coord_list[3*d] - new_coord_list_p[3*(d+1)]
				start_dist_y = coord_list[3*d+1] - new_coord_list_p[3*(d+1)+1]
				start_dist_z = coord_list[3*d+2] - new_coord_list_p[3*(d+1)+2]
				start_dist = math.sqrt(start_dist_x*start_dist_x + start_dist_y*start_dist_y + start_dist_z*start_dist_z)
				increment += 0.1
				#if dr < start_dist:
				#	continue
				#else:
				#	print "%.2f\n"%increment
			#print dr_opt
			dr_smallest = dr_opt[0][0]	
			inc_smallest = dr_opt[0][1]
			for alpha in range(1,len(dr_opt)):
				#print alpha
				if dr_opt[alpha][0] < dr_smallest:
					dr_smallest = dr_opt[alpha][0]
					inc_smallest = dr_opt[alpha][1]
			print dr_smallest, inc_smallest
			new_coord_list_p[3*d] = coord_list[3*d] + norm_v[3*d]*inc_smallest
			new_coord_list_p[3*d+1] = coord_list[3*d+1] + norm_v[3*d+1]*inc_smallest
			new_coord_list_p[3*d+2] = coord_list[3*d+2] + norm_v[3*d+2]*inc_smallest
			# Update the distances list
			dr_list[d] = dr_smallest
			print dr_list
			#d = 0
	"""
	print dr_list
	print "Bonds check complete\n"
#	sys.exit(2)
	return violations, new_coord_list
###################################################################################################
###################################################################################################
def BondDeviations(num_residues, coord_list, norm_v, new_coord_list, bonded_A, bonded_B, atom_shift, ang_dev, violations, dec):

	dec = float(dec)
	print len(bonded_A), bonded_A 
	print len(bonded_B), bonded_B
	k = 0
	dr_list = scipy.zeros(len(bonded_A))
	for i,j in zip(bonded_A,bonded_B):
		ii = int(i)-1
		jj = int(j)-1
#		print ii, jj
		x1 = new_coord_list[3*ii]
		y1 = new_coord_list[3*ii+1]
		z1 = new_coord_list[3*ii+2]
		x2 = new_coord_list[3*jj]
		y2 = new_coord_list[3*jj+1]
		z2 = new_coord_list[3*jj+2]
		dx = x1 - x2
		dy = y1 - y2
		dz = z1 - z2
		dr2 = dx*dx + dy*dy + dz*dz
		dr = math.sqrt(dr2)
		dr_list[k] = dr
		k += 1

	#print dr_list
	# check if bond length of any pair of atoms > 5.0
	bond_devi = []
	for d in range(len(dr_list)):
		if dr_list[d] > 4.0:
			bond_devi.append(d)
	print "bond_deviations\t", bond_devi
	
	#check if movement of atoms > 4.5 ALSO CHECK IF 4.5 CORRESPONDS TO BAD BOND LENGTHS IN DIFFERENT SYSTEMS
	tmp_atom_shift_viol = []
	for s in range(len(atom_shift)-1):
		 if atom_shift[s] > 5.0:						
			tmp_atom_shift_viol.append(s)

	# adjacent atom shift may lead to large shift in atoms
	atom_shift_viol = []
	if tmp_atom_shift_viol != []:
		for n in range(num_residues):
			if n in tmp_atom_shift_viol and n in bond_devi:
				atom_shift_viol.append(n)
			if n in tmp_atom_shift_viol and (n-1 in bond_devi or n+1 in bond_devi):
				atom_shift_viol.append(n)
		if tmp_atom_shift_viol[-1] == num_residues:
			atom_shift_viol.append(tmp_atom_shift_viol[-1])
		atom_shift_viol = list(set(atom_shift_viol))
		atom_shift_viol.sort()
		#print atom_shift_viol

	print "atom_shift_deviations\t", tmp_atom_shift_viol
	print "ang_deviations\t", ang_dev
	for i in range(num_residues):
		if i in bond_devi:
			violations.append(i)
			#if i in atom_shift_viol or i in ang_dev:
			#	violations.append(i)
	max_viol_allowd = len(violations)*0.05
	print "violations", violations
	
	#return violations
	#THRID APPROACH
	"""
	# Making subsets of connected atoms which have violations
	atom_set = []
	atom_subset = []
	#for s in range(1,len(ang_dev)):
	for s in range(1,len(violations)):
		if ang_dev[s]-ang_dev[s-1] == 1:
			if ang_dev[s-1] not in atom_subset:
				atom_subset.append(ang_dev[s-1])
			if ang_dev[s] not in atom_subset:
				atom_subset.append(ang_dev[s])
		else:
			atom_set.append((atom_subset))
			atom_subset = []
			if ang_dev[s] not in atom_subset:
				atom_subset.append(ang_dev[s])
	atom_set.append((atom_subset))
	atom_subset = []
	print atom_set
	"""

	"""
	for t in range(len(atom_set)):
		atom_subset=[]
		for c in atom_set[t]:
			atom_subset.append(c)
		if atom_subset[0] == 0:
			continue
			
			next = atom_subset[-1]+1 # next = last+1 residue in the first set of atoms
			print next
			atom_subset.reverse()
			for c in range(len(atom_subset)): # fix the movements in reversed order
				i = atom_subset[c]
				print new_coord_list[3*i], new_coord_list[3*i+1], new_coord_list[3*i+2]
				new_coord_list[3*i] = coord_list[3*i] + norm_v[3*next]*math.sqrt(num_residues)*6.50
				new_coord_list[3*i+1] = coord_list[3*i+1] + norm_v[3*next+1]*math.sqrt(num_residues)*6.50
				new_coord_list[3*i+2] = coord_list[3*i+2] + norm_v[3*next+2]*math.sqrt(num_residues)*6.50
				#next = i
				print next, new_coord_list[3*i], new_coord_list[3*i+1], new_coord_list[3*i+2]
			
		elif atom_subset[-1] == num_residues-1:
			continue

		else:
			for c in range(len(atom_subset)-1):
				i = atom_subset[c]
	"""		
			

	# SECOND APPROACH
	#check which set of joined atoms have atom movements > 4.5 
	
	atom_subset=[]
	atom_set=[]
	# Making subsets of connected atoms which have violations
	"""
	for s in range(1,len(atom_shift_viol)):
		if atom_shift_viol[s]-atom_shift_viol[s-1] == 1:
			if atom_shift_viol[s-1] not in atom_subset:
				atom_subset.append(atom_shift_viol[s-1])
			if atom_shift_viol[s] not in atom_subset:
				atom_subset.append(atom_shift_viol[s])
		else:
			atom_set.append((atom_subset))
			atom_subset = []
			if atom_shift_viol[s] not in atom_subset:
				atom_subset.append(atom_shift_viol[s])
	"""
	flag = 0
	if len(violations) == 0:
		return violations, new_coord_list
	elif len(violations) == 1:
		atom_set.append((violations))					
	else:
		atom_subset.append((violations[0]))	
		for s in range(1,len(violations)): #e.g. [0,292], [0,1,2], [0,1,3,5], [0,1,3,4], [0,2,3]
			if violations[s]-violations[s-1] == 1:
				if violations[s-1] not in atom_subset:
					atom_subset.append(violations[s-1])
				if violations[s] not in atom_subset:
					atom_subset.append(violations[s])
				flag = 1
			else:
				flag = 0
				if atom_subset != []:
					atom_set.append((atom_subset))
				atom_subset = []
				#atom_subset.append(violations[s])
				if violations[s] not in atom_subset:
					atom_subset.append(violations[s])
					flag = 1
					continue
				
	if flag == 1:
		atom_set.append((atom_subset))
	atom_subset = []
	print "atom_set", atom_set
	violations = []
	mark = 0
	end1 = 0
	end2 = 0
	if atom_set[-1][-1] == num_residues-2 and atom_set[0][0] == 0:
		# include bothe ends in violations
		mark = 2
		end1 = len(atom_set[0])
		end2 = len(atom_set[-1])
	elif atom_set[-1][-1] == num_residues-2 and atom_set[0][0] != 0:
		# include only one end in violations
		mark = 1
		end2 = len(atom_set[-1])
	elif atom_set[-1][-1] != num_residues-2 and atom_set[0][0] == 0:
		mark = 1
		end1 = len(atom_set[0])
	else:
		pass
	# For this to work, protein must have one atom that has correct bond lengths
	for t in range(len(atom_set)):
		atom_subset = []
		for c in atom_set[t]:
			atom_subset.append(c) 
		
		if atom_subset[0] == 0: #check if first set of atoms belong to beginning residues
			print "first type\n"
			check = 0
			for c in range(len(atom_subset)):
				if mark == 2:
					if end1 <= max_viol_allowd*0.5 and end2 <= max_viol_allowd*0.5 :
						violations.append(c)
					elif end1 > max_viol_allowd*0.5 and end2 <= max_viol_allowd*0.5 :
						if check < max_viol_allowd*0.5+(max_viol_allowd*0.5-end2):
							violations.append(c)				
					elif end1 <= max_viol_allowd*0.5 and end2 > max_viol_allowd*0.5 :
						violations.append(c)
					else:
						if check < max_viol_allowd*0.5:
							violations.append(c)
					check += 1
				if mark == 1:	
					if end1 != 0 and end1 <= max_viol_allowd:
						violations.append(c)
					elif end1 != 0 and end1 > max_viol_allowd:
						if check < max_viol_allowd:
							violations.append(c)
					else:
						pass
					check += 1
			#continue
			next = atom_subset[-1]+1 # next = last+1 residue in the first set of atoms
			print next
			atom_subset.reverse()
			for c in range(len(atom_subset)-1): # fix the movements in reversed order
				i = atom_subset[c] 
				bond_len_prev = 9999.99
				bond_len = dr_list[i]
				iter1 = 0
				if bond_len <= 4.0 and bond_len > 3.6:
					tmp_new_coord_list_x1 = new_coord_list[3*i]
					tmp_new_coord_list_y1 = new_coord_list[3*i+1]
					tmp_new_coord_list_z1 = new_coord_list[3*i+2]
					while bond_len <= 4.0 and bond_len > 3.6 and iter1 < 10:
						iter1 += 1
						iter2 = 0
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec
						dx = tmp_new_coord_list_x1 - new_coord_list[3*next]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*next+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*next+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						if bond_len <= 4.0 and bond_len > 3.6: # else: we don't change anything about i 
							pass
							#print atom_subset[c+1]
						else:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
							dx = tmp_new_coord_list_x1 - new_coord_list[3*next]
							dy = tmp_new_coord_list_y1 - new_coord_list[3*next+1]
							dz = tmp_new_coord_list_z1 - new_coord_list[3*next+2]
							bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						if atom_subset[c] == atom_subset[-1]:
							continue
						else:
							print atom_subset[c+1]
							j = atom_subset[c+1] # add condition for checking bond length later
							tmp_new_coord_list_x2 = new_coord_list[3*j]
							tmp_new_coord_list_y2 = new_coord_list[3*j+1]
							tmp_new_coord_list_z2 = new_coord_list[3*j+2]
							#bond_len2 = dr_list[j]
							bond_len2_prev = 9999.99
						 
							dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
							dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
							dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
							bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
							
							if bond_len2 <= 4.0 and bond_len2 > 3.6:
								print "break1"
								new_coord_list[3*j] = tmp_new_coord_list_x2
								new_coord_list[3*j+1] = tmp_new_coord_list_y2
								new_coord_list[3*j+2] = tmp_new_coord_list_z2
								dr_list[j] = bond_len2
								break	
							else:
								while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 10:
									iter2 += 1
									if bond_len2_prev < bond_len2:
										break
									bond_len2_prev = bond_len2
									tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
									dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
									dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
									dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
									bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
								if bond_len2 <= 4.0 and bond_len2 > 3.6:
									new_coord_list[3*j] = tmp_new_coord_list_x2
									new_coord_list[3*j+1] = tmp_new_coord_list_y2
									new_coord_list[3*j+2] = tmp_new_coord_list_z2
									dr_list[j] = bond_len2
								else:
									new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
									dr_list[j] = bond_len2_prev	
									
					new_coord_list[3*i] = tmp_new_coord_list_x1
					new_coord_list[3*i+1] = tmp_new_coord_list_y1
					new_coord_list[3*i+2] = tmp_new_coord_list_z1
					dr_list[i] = bond_len
					#continue
				else:
					print "i", i, bond_len
					tmp_new_coord_list_x1 = new_coord_list[3*i]
					tmp_new_coord_list_y1 = new_coord_list[3*i+1]
					tmp_new_coord_list_z1 = new_coord_list[3*i+2]
					while (bond_len > 4.0 or bond_len <= 3.6) and iter1 < 10:
						iter1 += 1
						iter2 = 0
						#print bond_len
						if bond_len_prev < bond_len:
							print i, bond_len, "break", iter1
							break
						bond_len_prev = bond_len
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec
						dx = tmp_new_coord_list_x1 - new_coord_list[3*next]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*next+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*next+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						if bond_len <= 4.0 and bond_len > 3.6:  
							#print atom_subset[c+1]
							pass
						else: #CHECK AGAIN
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
						if atom_subset[c] == atom_subset[-1]:
							continue
						else:
							#if atom_subset[c+1] >= atom_subset[-1]:
							j = atom_subset[c+1]
							bond_len2_prev = 9999.99
							tmp_new_coord_list_x2 = new_coord_list[3*j]
							tmp_new_coord_list_y2 = new_coord_list[3*j+1]
							tmp_new_coord_list_z2 = new_coord_list[3*j+2]
							dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
							dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
							dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
							bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
							if bond_len2 <= 4.0 and bond_len2 > 3.6:
								print "j ",j, bond_len2, "break"
								new_coord_list[3*j] = tmp_new_coord_list_x2
								new_coord_list[3*j+1] = tmp_new_coord_list_y2
								new_coord_list[3*j+2] = tmp_new_coord_list_z2
								dr_list[j] = bond_len2
								break
							else:
								while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 10:
									iter2 += 1
									if bond_len2_prev < bond_len2:
										print "j",j, bond_len2, "break5", iter2
										break
									bond_len2_prev = bond_len2
									tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
									dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
									dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
									dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
									bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
								if bond_len2 <= 4.0 and bond_len2 > 3.6:
									new_coord_list[3*j] = tmp_new_coord_list_x2
									new_coord_list[3*j+1] = tmp_new_coord_list_y2
									new_coord_list[3*j+2] = tmp_new_coord_list_z2
									dr_list[j] = bond_len2
								else:
									new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
									dr_list[j] = bond_len2_prev
					if bond_len <= 4.0 and bond_len > 3.6:
						new_coord_list[3*i] = tmp_new_coord_list_x1
						new_coord_list[3*i+1] = tmp_new_coord_list_y1
						new_coord_list[3*i+2] = tmp_new_coord_list_z1
						dr_list[i] = bond_len
					else:				
						new_coord_list[3*i] = tmp_new_coord_list_x1 + norm_v[3*i]*math.sqrt(num_residues)*dec*-1.0
						new_coord_list[3*i+1] = tmp_new_coord_list_y1 + norm_v[3*i+1]*math.sqrt(num_residues)*dec*-1.0
						new_coord_list[3*i+2] = tmp_new_coord_list_z1 + norm_v[3*i+2]*math.sqrt(num_residues)*dec*-1.0
						dr_list[i] = bond_len_prev

		elif atom_subset[-1] == num_residues-2 :#check if the last set of atoms belong to the end terminal
			print "Second type\n"
			check = 0
			for c in atom_subset:
				if mark == 2:
					if end1 <= max_viol_allowd*0.5 and end2 <= max_viol_allowd*0.5 :
						violations.append(c+1)
					elif end1 > max_viol_allowd*0.5 and end2 <= max_viol_allowd*0.5 :
						violations.append(c+1)
					elif end1 <= max_viol_allowd*0.5 and end2 > max_viol_allowd*0.5 :
						if check < max_viol_allowd*0.5 + (max_viol_allowd*0.5-end1):
							violations.append(c+1)
					else:
						if check < max_viol_allowd*0.5:
							violations.append(c+1)
                                        check += 1
                                if mark == 1:
                                        if end2 != 0 and end2 <= max_viol_allowd:
                                                violations.append(c+1)
                                        elif end2 != 0 and end2 > max_viol_allowd:
                                                if check < max_viol_allowd:
                                                        violations.append(c+1)
                                        else:
                                                pass
                                        check += 1
#			continue
			prev = atom_subset[0]-1 # prev = first-1 residue in the last set of atoms
			for c in range(len(atom_subset)-1): # fix the movements in reversed order
				i = atom_subset[c] 
				bond_len_prev = 9999.99
				bond_len = dr_list[i]
				iter1 = 0
				if bond_len <= 4.0 and bond_len > 3.6:
					tmp_new_coord_list_x1 = new_coord_list[3*i]
					tmp_new_coord_list_y1 = new_coord_list[3*i+1]
					tmp_new_coord_list_z1 = new_coord_list[3*i+2]
					while bond_len <= 4.0 and bond_len > 3.6 and iter1 < 10:
						iter1 += 1
						iter2 = 0
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec
						dx = tmp_new_coord_list_x1 - new_coord_list[3*prev]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*prev+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*prev+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						if bond_len <= 4.0 and bond_len > 3.6:  
							#print atom_subset[c+1]
							pass
						else:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
						if atom_subset[c] == atom_subset[-1]:
							continue
						else:
						#if atom_subset[c+1] != atom_subset[-1]:
							print atom_subset[c+1]
							j = atom_subset[c+1] # add condition for checking bond length later
							tmp_new_coord_list_x2 = new_coord_list[3*j]
							tmp_new_coord_list_y2 = new_coord_list[3*j+1]
							tmp_new_coord_list_z2 = new_coord_list[3*j+2]
#							bond_len2 = dr_list[j]
							bond_len2_prev = 9999.99
						 
							dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
							dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
							dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
							bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
							
							if bond_len2 <= 4.0 and bond_len2 > 3.6:
									print "break1"
									new_coord_list[3*j] = tmp_new_coord_list_x2
									new_coord_list[3*j+1] = tmp_new_coord_list_y2
									new_coord_list[3*j+2] = tmp_new_coord_list_z2
									dr_list[i] = bond_len2
									break
							else:
								while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 10:
									iter2 += 1
									if bond_len2_prev < bond_len2:
										break
									bond_len2_prev = bond_len2
									tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
									dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
									dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
									dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
									bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
								if bond_len2 <= 4.0 and bond_len2 > 3.6:
									new_coord_list_p[3*j] = tmp_new_coord_list_x2
									new_coord_list_p[3*j+1] = tmp_new_coord_list_y2
									new_coord_list_p[3*j+2] = tmp_new_coord_list_z2
									dr_list[i] = bond_len2
								else:
									new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
									dr_list[i] = bond_len2_prev	
					new_coord_list[3*i] = tmp_new_coord_list_x1
					new_coord_list[3*i+1] = tmp_new_coord_list_y1
					new_coord_list[3*i+2] = tmp_new_coord_list_z1
					dr_list[prev] = bond_len
					#continue
				else:
					print "i", i, bond_len
					tmp_new_coord_list_x1 = new_coord_list[3*i]
					tmp_new_coord_list_y1 = new_coord_list[3*i+1]
					tmp_new_coord_list_z1 = new_coord_list[3*i+2]
					while (bond_len > 4.0 or bond_len <= 3.6) and iter1 < 10:
						iter1 += 1
						iter2 = 0
						#print bond_len
						if bond_len_prev < bond_len:
							print i, bond_len
							break
						bond_len_prev = bond_len
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec
				
						dx = tmp_new_coord_list_x1 - new_coord_list[3*prev]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*prev+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*prev+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						if bond_len <= 4.0 and bond_len > 3.6: # else: we don't change anything about i 
							#print atom_subset[c+1]
							pass
						else:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
						if atom_subset[c] == atom_subset[-1]:
							continue
						else:
							j = atom_subset[c+1]
							tmp_new_coord_list_x2 = new_coord_list[3*j]
							tmp_new_coord_list_y2 = new_coord_list[3*j+1]
							tmp_new_coord_list_z2 = new_coord_list[3*j+2]
#							bond_len2 = dr_list[j]
							bond_len2_prev = 9999.99
							
							dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
							dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
							dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
							bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
							if bond_len2 <= 4.0 and bond_len2 > 3.6:
								print j, bond_len2
								new_coord_list[3*j] = tmp_new_coord_list_x2
								new_coord_list[3*j+1] = tmp_new_coord_list_y2
								new_coord_list[3*j+2] = tmp_new_coord_list_z2
								dr_list[i] = bond_len2
								break
							else:
								while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 10:
									iter2 += 1
									if bond_len2_prev < bond_len2:
										print j, bond_len2
										break
									bond_len2_prev = bond_len2
									tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
									tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
									dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
									dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
									dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
									bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
								if bond_len2 <= 4.0 and bond_len2 > 3.6:		
									new_coord_list[3*j] = tmp_new_coord_list_x2
									new_coord_list[3*j+1] = tmp_new_coord_list_y2
									new_coord_list[3*j+2] = tmp_new_coord_list_z2
									dr_list[i] = bond_len2
								else:
									new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
									new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
									dr_list[i] = bond_len2_prev			
					if bond_len <= 4.0 and bond_len > 3.6:
						new_coord_list[3*i] = tmp_new_coord_list_x1
						new_coord_list[3*i+1] = tmp_new_coord_list_y1
						new_coord_list[3*i+2] = tmp_new_coord_list_z1
						dr_list[prev] = bond_len
					else:
						new_coord_list[3*i] = tmp_new_coord_list_x1 + norm_v[3*i]*math.sqrt(num_residues)*dec*-1.0
						new_coord_list[3*i+1] = tmp_new_coord_list_y1 + norm_v[3*i+1]*math.sqrt(num_residues)*dec*-1.0
						new_coord_list[3*i+2] = tmp_new_coord_list_z1 + norm_v[3*i+2]*math.sqrt(num_residues)*dec*-1.0
						dr_list[prev] = bond_len_prev
		else:
			print "Third type\n"	# No rule about atom_subset of length 1
			prev = atom_subset[0]-1
			next = atom_subset[-1]+1
			if atom_shift[prev] < atom_shift[next]:
				direction = "forward" # move forwards
				direc = prev
				print "move forwards"

			else:
				direction = "backward" # move backwards
				direc = next
				atom_subset.reverse()
				print "move backwards"


			print atom_subset
			for c in range(len(atom_subset)-1):
				i = atom_subset[c] 
				bond_len_prev = 9999.99
				tmp_new_coord_list_x1 = new_coord_list[3*i]
				tmp_new_coord_list_y1 = new_coord_list[3*i+1]
				tmp_new_coord_list_z1 = new_coord_list[3*i+2]
				dx = tmp_new_coord_list_x1 - new_coord_list[3*direc]
				dy = tmp_new_coord_list_y1 - new_coord_list[3*direc+1]
				dz = tmp_new_coord_list_z1 - new_coord_list[3*direc+2]
				bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
				print "i",i,bond_len, tmp_new_coord_list_x1, tmp_new_coord_list_y1, tmp_new_coord_list_z1
				iter1= 0
				# if bond_len < 3.6 STOP!!!!!!!!
				if bond_len <= 4.0 and bond_len > 3.6:
					while (bond_len <= 4.0 and bond_len > 3.6) and iter1 < 10:
						iter1 += 1
						iter2=0
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec 
						dx = tmp_new_coord_list_x1 - new_coord_list[3*direc]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*direc+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*direc+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
#						print direc,bond_len,tmp_new_coord_list_x1,tmp_new_coord_list_y1,tmp_new_coord_list_z1
						if bond_len <= 4.0 and bond_len > 3.6:  
							pass
							print direc,bond_len,tmp_new_coord_list_x1,tmp_new_coord_list_y1,tmp_new_coord_list_z1
						else:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
							dx = tmp_new_coord_list_x1 - new_coord_list[3*direc]
							dy = tmp_new_coord_list_y1 - new_coord_list[3*direc+1]
							dz = tmp_new_coord_list_z1 - new_coord_list[3*direc+2]
							bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
							print direc,bond_len,tmp_new_coord_list_x1,tmp_new_coord_list_y1,tmp_new_coord_list_z1
															
						j = atom_subset[c+1]
						bond_len2_prev = 9999.99
						tmp_new_coord_list_x2 = new_coord_list[3*j]
						tmp_new_coord_list_y2 = new_coord_list[3*j+1]
						tmp_new_coord_list_z2 = new_coord_list[3*j+2]
						dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
						dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
						dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
						bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
						print i,bond_len2, tmp_new_coord_list_x2, tmp_new_coord_list_y2, tmp_new_coord_list_z2
						if bond_len2 <= 4.0 and bond_len2 > 3.6:
							print "break1"
							new_coord_list[3*j] = tmp_new_coord_list_x2
							new_coord_list[3*j+1] = tmp_new_coord_list_y2
							new_coord_list[3*j+2] = tmp_new_coord_list_z2
							if direction == "backward":
								dr_list[j] = bond_len2
								if j == atom_subset[-1]:
									dx=tmp_new_coord_list_x2-new_coord_list[3*prev]
									dy=tmp_new_coord_list_y2-new_coord_list[3*prev+1]
									dz=tmp_new_coord_list_z2-new_coord_list[3*prev+2]
									bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
									dr_list[prev] = bond_len3
							if direction == "forward":
								dr_list[i] = bond_len2
								if j == atom_subset[-1]:
									dx=tmp_new_coord_list_x2-new_coord_list[3*next]
									dy=tmp_new_coord_list_y2-new_coord_list[3*next+1]
									dz=tmp_new_coord_list_z2-new_coord_list[3*next+2]
									bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
									dr_list[j] = bond_len3
							break 
						else:
							
							while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 10:
								iter2 += 1
								if bond_len2_prev < bond_len2:
									print i, bond_len2_prev, bond_len2, "break2"
									#flag = 1
									break  
								tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
								tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
								tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
								bond_len2_prev = bond_len2
								dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
								dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
								dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
								bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
								print i,bond_len2, tmp_new_coord_list_x2, tmp_new_coord_list_y2, tmp_new_coord_list_z2
							
							if bond_len2 <= 4.0 and bond_len2 > 3.6: #bond_len3 <= 4.0 and bond_len3 > 3.6
								new_coord_list[3*j] = tmp_new_coord_list_x2
								new_coord_list[3*j+1] = tmp_new_coord_list_y2
								new_coord_list[3*j+2] = tmp_new_coord_list_z2
								if direction == "backward":
									dr_list[j] = bond_len2
									if j == atom_subset[-1]:
										dx=tmp_new_coord_list_x2-new_coord_list[3*prev]
										dy=tmp_new_coord_list_y2-new_coord_list[3*prev+1]
										dz=tmp_new_coord_list_z2-new_coord_list[3*prev+2]
										bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
										dr_list[prev] = bond_len3
								if direction == "forward":
									dr_list[i] = bond_len2
									if j == atom_subset[-1]:
										dx=tmp_new_coord_list_x2-new_coord_list[3*next]
										dy=tmp_new_coord_list_y2-new_coord_list[3*next+1]
										dz=tmp_new_coord_list_z2-new_coord_list[3*next+2]
										bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
										dr_list[j] = bond_len3
							else:
								new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
								if direction == "backward":
									dr_list[j] = bond_len2_prev
								
									if j == atom_subset[-1]:
										dr_list[prev] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*prev])*(new_coord_list[3*j]-new_coord_list[3*prev]) + (new_coord_list[3*j+1]-new_coord_list[3*prev+1])*(new_coord_list[3*j+1]-new_coord_list[3*prev+1]) + (new_coord_list[3*j+2]-new_coord_list[3*prev+2])*(new_coord_list[3*j+2]-new_coord_list[3*prev+2]))

								if direction == "forward":
									dr_list[i] = bond_len2_prev
									if j == atom_subset[-1]:
										dr_list[j] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*next])*(new_coord_list[3*j]-new_coord_list[3*next]) + (new_coord_list[3*j+1]-new_coord_list[3*next+1])*(new_coord_list[3*j+1]-new_coord_list[3*next+1]) + (new_coord_list[3*j+2]-new_coord_list[3*next+2])*(new_coord_list[3*j+2]-new_coord_list[3*next+2]))
					
					new_coord_list[3*i] = tmp_new_coord_list_x1
					new_coord_list[3*i+1] = tmp_new_coord_list_y1
					new_coord_list[3*i+2] = tmp_new_coord_list_z1
					if direction == "backward":
						dr_list[i] = bond_len
					if direction == "forward":
						dr_list[direc] = bond_len
					direc = i
				else:
					#print "i", i, bond_len, tmp_new_coord_list_x1, tmp_new_coord_list_y1, tmp_new_coord_list_z1
					#tmp_new_coord_list_x1 = new_coord_list[3*i]
					#tmp_new_coord_list_y1 = new_coord_list[3*i+1]
					#tmp_new_coord_list_z1 = new_coord_list[3*i+2]
					while (bond_len > 4.0 or bond_len <= 3.6) and iter1 < 10:
						iter1 += 1
						iter2 = 0
						if bond_len_prev < bond_len:
							print i, bond_len_prev, bond_len, "break3"
							break
						bond_len_prev = bond_len
						tmp_new_coord_list_x1 += norm_v[3*i]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_y1 += norm_v[3*i+1]*math.sqrt(num_residues)*dec
						tmp_new_coord_list_z1 += norm_v[3*i+2]*math.sqrt(num_residues)*dec
						dx = tmp_new_coord_list_x1 - new_coord_list[3*direc]
						dy = tmp_new_coord_list_y1 - new_coord_list[3*direc+1]
						dz = tmp_new_coord_list_z1 - new_coord_list[3*direc+2]
						bond_len = math.sqrt(dx*dx+dy*dy+dz*dz)
						print direc, bond_len, tmp_new_coord_list_x1, tmp_new_coord_list_y1, tmp_new_coord_list_z1
						if bond_len <= 4.0 and bond_len > 3.6:  
							#print atom_subset[c+1]
							pass 
						else:
							tmp_new_coord_list_x1 -= norm_v[3*i]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_y1 -= norm_v[3*i+1]*math.sqrt(num_residues)*dec
							tmp_new_coord_list_z1 -= norm_v[3*i+2]*math.sqrt(num_residues)*dec
															
						j = atom_subset[c+1]
						bond_len2_prev = 9999.99
						tmp_new_coord_list_x2 = new_coord_list[3*j]
						tmp_new_coord_list_y2 = new_coord_list[3*j+1]
						tmp_new_coord_list_z2 = new_coord_list[3*j+2]
						dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
						dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
						dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
						bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
						print j, bond_len2, tmp_new_coord_list_x2, tmp_new_coord_list_y2, tmp_new_coord_list_z2
						if bond_len2 <= 4.0 and bond_len2 > 3.6:
							print j, bond_len2, "break4"
							new_coord_list[3*j] = tmp_new_coord_list_x2
							new_coord_list[3*j+1] = tmp_new_coord_list_y2
							new_coord_list[3*j+2] = tmp_new_coord_list_z2
							if direction == "backward":
								dr_list[j] = bond_len2
								if j == atom_subset[-1]:
									dx=tmp_new_coord_list_x2-new_coord_list[3*prev]
									dy=tmp_new_coord_list_y2-new_coord_list[3*prev+1]
									dz=tmp_new_coord_list_z2-new_coord_list[3*prev+2]
									bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
									dr_list[prev] = bond_len3
							if direction == "forward":
								dr_list[i] = bond_len2
								if j == atom_subset[-1]:
									dx=tmp_new_coord_list_x2-new_coord_list[3*next]
									dy=tmp_new_coord_list_y2-new_coord_list[3*next+1]
									dz=tmp_new_coord_list_z2-new_coord_list[3*next+2]
									bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
									dr_list[j] = bond_len3
							break
						else:
							while (bond_len2 > 4.0 or bond_len2 <= 3.6) and iter2 < 10:
								iter2 += 1
								if bond_len2_prev < bond_len2:
									print j, bond_len2_prev, bond_len2, "break5"
									break
								bond_len2_prev = bond_len2
								tmp_new_coord_list_x2 += norm_v[3*j]*math.sqrt(num_residues)*dec
								tmp_new_coord_list_y2 += norm_v[3*j+1]*math.sqrt(num_residues)*dec
								tmp_new_coord_list_z2 += norm_v[3*j+2]*math.sqrt(num_residues)*dec
								dx = tmp_new_coord_list_x1 - tmp_new_coord_list_x2
								dy = tmp_new_coord_list_y1 - tmp_new_coord_list_y2
								dz = tmp_new_coord_list_z1 - tmp_new_coord_list_z2
								bond_len2 = math.sqrt(dx*dx+dy*dy+dz*dz)
								print j,bond_len2, tmp_new_coord_list_x2, tmp_new_coord_list_y2, tmp_new_coord_list_z2
								
							if bond_len2 <= 4.0 and bond_len2 > 3.6:
								new_coord_list[3*j] = tmp_new_coord_list_x2
								new_coord_list[3*j+1] = tmp_new_coord_list_y2
								new_coord_list[3*j+2] = tmp_new_coord_list_z2
								if direction == "backward":
									dr_list[j] = bond_len2
									if j == atom_subset[-1]:
										dx=tmp_new_coord_list_x2-new_coord_list[3*prev]
										dy=tmp_new_coord_list_y2-new_coord_list[3*prev+1]
										dz=tmp_new_coord_list_z2-new_coord_list[3*prev+2]
										bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
										dr_list[prev] = bond_len3
								if direction == "forward":
									dr_list[i] = bond_len2
									if j == atom_subset[-1]:
										dx = tmp_new_coord_list_x2-new_coord_list[3*next]
										dy = tmp_new_coord_list_y2-new_coord_list[3*next+1]
										dz = tmp_new_coord_list_z2-new_coord_list[3*next+2]
										bond_len3 = math.sqrt(dx*dx+dy*dy+dz*dz)
										dr_list[j] = bond_len3				
							else:
								new_coord_list[3*j] = tmp_new_coord_list_x2 + norm_v[3*j]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+1] = tmp_new_coord_list_y2 + norm_v[3*j+1]*math.sqrt(num_residues)*dec*-1.0
								new_coord_list[3*j+2] = tmp_new_coord_list_z2 + norm_v[3*j+2]*math.sqrt(num_residues)*dec*-1.0
								#dr_list[j] = bond_len2_prev
								if direction == "backward":
									dr_list[j] = bond_len2_prev
									if j == atom_subset[-1]:
										dr_list[prev] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*prev])*(new_coord_list[3*j]-new_coord_list[3*prev]) + (new_coord_list[3*j+1]-new_coord_list[3*prev+1])*(new_coord_list[3*j+1]-new_coord_list[3*prev+1]) + (new_coord_list[3*j+2]-new_coord_list[3*prev+2])*(new_coord_list[3*j+2]-new_coord_list[3*prev+2]))
								if direction == "forward":
									dr_list[i] = bond_len2_prev
									if j == atom_subset[-1]:
										dr_list[j] = math.sqrt((new_coord_list[3*j]-new_coord_list[3*next])*(new_coord_list[3*j]-new_coord_list[3*next]) + (new_coord_list[3*j+1]-new_coord_list[3*next+1])*(new_coord_list[3*j+1]-new_coord_list[3*next+1]) + (new_coord_list[3*j+2]-new_coord_list[3*next+2])*(new_coord_list[3*j+2]-new_coord_list[3*next+2]))
					if bond_len <= 4.0 and bond_len > 3.6:						
						new_coord_list[3*i] = tmp_new_coord_list_x1
						new_coord_list[3*i+1] = tmp_new_coord_list_y1
						new_coord_list[3*i+2] = tmp_new_coord_list_z1
						if direction == "backward":
							dr_list[i] = bond_len
						if direction == "forward":
							dr_list[direc] = bond_len
						direc = i
					else:
						new_coord_list[3*i] = tmp_new_coord_list_x1 + norm_v[3*i]*math.sqrt(num_residues)*dec*-1.0
						new_coord_list[3*i+1] = tmp_new_coord_list_y1 + norm_v[3*i+1]*math.sqrt(num_residues)*dec*-1.0
						new_coord_list[3*i+2] = tmp_new_coord_list_z1 + norm_v[3*i+2]*math.sqrt(num_residues)*dec*-1.0
						if direction == "backward":
							dr_list[i] = bond_len_prev
						if direction == "forward":
							dr_list[direc] = bond_len_prev
						direc = i
					#print dr_list[i], dr_list[j]
#			print dr_list				
	
	# FIRST APPROACH
	"""				
	#for bd in range(len(bond_dev)):
	for d in reversed(range(len(dr_list)-1)):
		dr_prev = 0.0
		dr_opt = []
		if (dr_list[d] > 4.5) and (dr_list[d+1] <= 4.5 and dr_list[d+1] >= 3.6):
			dr = dr_list[d]
			increment = 0.0
			#dr = 0.0
			while (dr > 4.5 and dr <= dr_list[d]):
				#print d
				x = (coord_list[3*d] + norm_v[3*d]*increment) - new_coord_list_p[3*(d+1)]
				y = (coord_list[3*d+1] + norm_v[3*d+1]*increment) - new_coord_list_p[3*(d+1)+1]
				z = (coord_list[3*d+2] + norm_v[3*d+2]*increment) - new_coord_list_p[3*(d+1)+2]
				#print x,y,z
				dr = math.sqrt(x*x + y*y + z*z)
				dr_opt.append((dr,increment))
				#print dr, increment
				start_dist_x = coord_list[3*d] - new_coord_list_p[3*(d+1)]
				start_dist_y = coord_list[3*d+1] - new_coord_list_p[3*(d+1)+1]
				start_dist_z = coord_list[3*d+2] - new_coord_list_p[3*(d+1)+2]
				start_dist = math.sqrt(start_dist_x*start_dist_x + start_dist_y*start_dist_y + start_dist_z*start_dist_z)
				increment += 0.1
				#if dr < start_dist:
				#	continue
				#else:
				#	print "%.2f\n"%increment
			#print dr_opt
			dr_smallest = dr_opt[0][0]	
			inc_smallest = dr_opt[0][1]
			for alpha in range(1,len(dr_opt)):
				#print alpha
				if dr_opt[alpha][0] < dr_smallest:
					dr_smallest = dr_opt[alpha][0]
					inc_smallest = dr_opt[alpha][1]
			print dr_smallest, inc_smallest
			new_coord_list_p[3*d] = coord_list[3*d] + norm_v[3*d]*inc_smallest
			new_coord_list_p[3*d+1] = coord_list[3*d+1] + norm_v[3*d+1]*inc_smallest
			new_coord_list_p[3*d+2] = coord_list[3*d+2] + norm_v[3*d+2]*inc_smallest
			# Update the distances list
			dr_list[d] = dr_smallest
			print dr_list
			#d = 0
	"""
	return violations, new_coord_list
###################################################################################################
###################################################################################################
def	generateMol2(filename, coord_list, violations):
	# open file
	fo = open(filename, "w")
	#print "violations", violations
	# find number of bonds (use distance between CA-CA < 4.5A; typically around 3.8A)
	nb = 0
	for i in range(num_residues-1):
		if i in violations: #[689]
			print i
			continue
		else:
			j = i+1
			dx = coord_list[3*i  ] - coord_list[3*j  ]
			dy = coord_list[3*i+1] - coord_list[3*j+1]
			dz = coord_list[3*i+2] - coord_list[3*j+2]
			dr = math.sqrt(dx*dx + dy*dy + dz*dz)
			if dr < 4.5:
				nb += 1
	
	# write header
	fo.write("@@<TRIPOS>MOLECULE\n****\n")
	fo.write("% 5d % 5d     0     0     0\n" % (num_residues-len(violations), nb))
	fo.write("SMALL\nUSER_CHARGES\n\n\n@@<TRIPOS>ATOM\n");

	# write coordinates
	for i in range(num_residues):
		if i in violations:
			continue
		else:
			fo.write("% 5d  CA           % 10.4f% 10.4f% 10.4f C       % 4d %3s    0.0000\n" % (i+1, coord_list[3*i  ], coord_list[3*i+1], coord_list[3*i+2], i+1, ENM_res_list[i]))

	# write bonds
	fo.write("@@<TRIPOS>BOND\n")
	nb = 0
	for i in range(num_residues-1):
		if i in violations:
			continue
		else:
			j = i+1
			dx = coord_list[3*i  ] - coord_list[3*j  ]
			dy = coord_list[3*i+1] - coord_list[3*j+1]
			dz = coord_list[3*i+2] - coord_list[3*j+2]
			dr = math.sqrt(dx*dx + dy*dy + dz*dz)
			if dr < 4.5:
				fo.write("% 5d %5d %5d    1\n" % (nb+1, i+1, j+1))
				nb += 1
	
	# close file
	fo.close()
###################################################################################################
###################################################################################################
def	generatePDB(filename, coord_list, violations, pdb_file):
	# open file
	fi = open(pdb_file)
	file_lines = []
	file_lines = fi.readlines()
	res_num_list = []
	for lines in file_lines:
		lines=lines.strip()
		if lines.split()[0] == "ATOM" and lines[12:16].find("CA") >= 0:
			res_num_list.append(int(lines[22:26]))
		if lines.split()[0] == "HETATM" :
			res_num_list.append(int(lines[22:26]))
	fi.close()

	fo = open(filename, "w")
	# write coordinates
	for i in range(num_residues):
		if i in violations:
			continue
		else:
#			fo.write("ATOM% 7d  CA  %s A% 4d    % 8.3f% 8.3f% 8.3f  1.00  0.00\n"%(i+1,ENM_res_list[i],i+1,coord_list[3*i],coord_list[3*i+1],coord_list[3*i+2]))
			#print "ATOM% 7d  CA  %s A% 4d    % 8.3f% 8.3f% 8.3f  1.00  0.00"%(i+1,ENM_res_list[i],res_num_list[i],coord_list[3*i],coord_list[3*i+1],coord_list[3*i+2])
			fo.write("ATOM% 7d  CA  %s A% 4d    % 8.3f% 8.3f% 8.3f  1.00  0.00\n"%(i+1,ENM_res_list[i],res_num_list[i],coord_list[3*i],coord_list[3*i+1],coord_list[3*i+2]))
	fo.close()

###################################################################################################
###################################################################################################
def	RemoveCOMVelocities(va, mass, num_residues):
	COM_vel = scipy.zeros((3))
	COM_mass = 0.0;

	for i in range(num_residues):
		for j in range(3):
			COM_vel[j] += va[3*i+j]*mass[i]
		COM_mass += mass[i]
#	print "COM: %f %f %f  : %f\n" % (COM_vel[0], COM_vel[1], COM_vel[2], COM_mass)

	if COM_mass > 0.00:
		for j in range(3):
			COM_vel[j] /= COM_mass;

	for i in range(num_residues):
		for j in range(3):
			va[3*i+j] -= COM_vel[j]


###################################################################################################
###################################################################################################
def	CalculateTemperature(va, mass, num_residues):
	global Temp
	vFactor = 1.2028e6

	NumConstraints = 0

	Temp = 0.0
	for i in range(num_residues):
		for j in range(3):
			Temp += mass[i]*pow(va[3*i+j], 2)
	Temp *= vFactor/(3.0*num_residues - NumConstraints)


###################################################################################################
###################################################################################################
def	calculateForces(force_const, g_width, g_V0, dr_list, coord_list, fa, num_residues):

	dr_vector   = scipy.zeros((num_residues, num_residues, 3))
	dr_list_act = scipy.zeros((num_residues, num_residues))

	# compute distance matrix (use C function)
	code_forces =	  """
									int			i, j, i3, j3;
									double	x1, y1, z1, x2, y2, z2, dx, dy, dz, dr2, dr, drinv, drr, expdrr, prefact;

									for(i = 0; i < num_residues; i++)
									{
										i3 = 3*i;
										fa(i3  ) = 0.0;
										fa(i3+1) = 0.0;
										fa(i3+2) = 0.0;
									}
										
									for(i = 0; i < num_residues; i++)
									{
										i3 = 3*i;
										x1 = coord_list(i3);
										y1 = coord_list(i3+1);
										z1 = coord_list(i3+2);
										for(j = i+1; j < num_residues; j++)
										{
											if(force_const(i,j) > 0.001 || force_const(i,j) < -0.001)
											{
												j3 = 3*j;
												x2 = coord_list(j3);
												y2 = coord_list(j3+1);
												z2 = coord_list(j3+2);
												// distance
												dx = x1 - x2;
												dy = y1 - y2;
												dz = z1 - z2;
												dr2 = dx*dx + dy*dy + dz*dz;
												dr = sqrt(dr2);
												drinv = 1/dr;
												dr_vector(i,j,0)  = dx*drinv;
												dr_vector(i,j,1)  = dy*drinv;
												dr_vector(i,j,2)  = dz*drinv;
												dr_list_act(i,j)  = dr;
												// force
												drr = dr_list_act(i,j) - dr_list(i,j);
												dx = drr*dr_vector(i,j,0);
												dy = drr*dr_vector(i,j,1);
												dz = drr*dr_vector(i,j,2);
												prefact = force_const(i,j);
												//printf("% 5d % 5d : % 10.5f % 10.5f % 10.5f : % 10.5f : % 10.5f\\n", i, j, dx, dy, dz, force_const(i,j), drr);
												fa(i3  ) -= prefact*dx;
												fa(i3+1) -= prefact*dy;
												fa(i3+2) -= prefact*dz;
												fa(j3  ) += prefact*dx;
												fa(j3+1) += prefact*dy;
												fa(j3+2) += prefact*dz;
											}
											if(g_V0(i,j) > 0.001 || g_V0(i,j) < -0.001)
											{
												j3 = 3*j;
												x2 = coord_list(j3);
												y2 = coord_list(j3+1);
												z2 = coord_list(j3+2);
												// distance
												dx = x1 - x2;
												dy = y1 - y2;
												dz = z1 - z2;
												dr2 = dx*dx + dy*dy + dz*dz;
												dr = sqrt(dr2);
												drinv = 1/dr;
												dr_vector(i,j,0)  = dx*drinv;
												dr_vector(i,j,1)  = dy*drinv;
												dr_vector(i,j,2)  = dz*drinv;
												dr_list_act(i,j)  = dr;
												// force
												drr = dr_list_act(i,j) - dr_list(i,j);
												expdrr = drr*exp(-0.5*drr*drr/(g_width(i,j)*g_width(i,j)));
												dx = expdrr*dr_vector(i,j,0);
												dy = expdrr*dr_vector(i,j,1);
												dz = expdrr*dr_vector(i,j,2);
												prefact = g_V0(i,j)/(g_width(i,j)*g_width(i,j));
												//printf("% 5d % 5d : % 10.5f % 10.5f % 10.5f : % 10.5f : % 10.5f % 10.5f % 10.5f\\n", i, j, dx, dy, dz, prefact, drr, expdrr, g_width(i,j));
												fa(i3  ) += prefact*dx;
												fa(i3+1) += prefact*dy;
												fa(i3+2) += prefact*dz;
												fa(j3  ) -= prefact*dx;
												fa(j3+1) -= prefact*dy;
												fa(j3+2) -= prefact*dz;
											}
										}
									}
									"""
	weave.inline(code_forces,
							['num_residues', 'force_const', 'g_width', 'g_V0', 'coord_list', 'dr_vector', 'dr_list', 'dr_list_act', 'fa'],
							type_converters = converters.blitz,
							compiler='gcc')



###################################################################################################
###################################################################################################
def	writeCoordinates(fcoord, xa, num_residues):
	n = 0
	for i in range(num_residues):
		for j in range(3):
			if n >= 10:
				fcoord.write("\n")
				n = 0
			fcoord.write("% 8.3f" % xa[3*i+j])
			n += 1
	if n >= 1:
		fcoord.write("\n")
		

###################################################################################################
###################################################################################################
def	initializeVelocities(va, mass_list, Temp0, num_residues):
	vFactor = math.sqrt(Temp0)*9.118243e-4;

	for i in range(num_residues):
		vPre = vFactor/math.sqrt(mass_list[i]);
		va[3*i]   = numpy.random.normal(0.0, vPre)
		va[3*i+1] = numpy.random.normal(0.0, vPre)
		va[3*i+2] = numpy.random.normal(0.0, vPre)
#		print "vel0: %d %f %f %f   %f" % (i, va[3*i], va[3*i+1], va[3*i+2], Temp0)

	RemoveCOMVelocities(va, mass_list, num_residues)
	CalculateTemperature(va, mass_list, num_residues)

	# scale velocities to achieve initial temperature
	lambd = math.sqrt(Temp0/Temp);

	for i in range(num_residues):
		for j in range(3):
			va[3*i+j] *= lambd
		#print "vel0: %d %f %f %f   %f" % (i, va[3*i], va[3*i+1], va[3*i+2], Temp0)

###################################################################################################
###################################################################################################
def	berendsenThermostat(va, mass_list, Temp0, dt, tau):

	CalculateTemperature(va, mass_list)

	lambd = (dt/tau)*(Temp0/Temp - 1);
	if lambd < -0.99:
		lambd = -0.99
	lambd = math.sqrt(1.0 + lambd)


	# scale velocities
	for i in range(num_residues):
		for j in range(3):
			va[3*i+j] *= lambd

###################################################################################################
###################################################################################################
def	shakeBonds(xa, xprev, mass_list, dt, dr2_list, bond_list_A, bond_list_B, num_residues):
	#print bond_list_A, bond_list_B
	inv_mass = scipy.zeros((num_residues))
	for i in range(num_residues):
		inv_mass[i] = 1/mass_list[i]
		#print "% 3d:% 8.3f% 8.3f% 8.3f "%(i,xa[3*i],xa[3*i+1],xa[3*i+2]) ---- CORRECT
		#sys.exit(2)
#	print xa[0], getframeinfo(currentframe()).lineno
	NumBondList = len(bond_list_A)
	for i in range(NumBondList):
		bi = bond_list_A[i]
		bxi = 3*bi
		bj = bond_list_B[i]
		bxj = 3*bj
		rij0 = xa[bxi]   - xa[bxj];
		rij1 = xa[bxi+1] - xa[bxj+1];
		rij2 = xa[bxi+2] - xa[bxj+2];
		r2 = rij0*rij0 + rij1*rij1 + rij2*rij2;
#		print "%d-%d:%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%f"%(bi,bj,xa[bxi],xa[bxi+1],xa[bxi+2],xa[bxj],xa[bxj+1],xa[bxj+2],rij0,rij1,rij2,r2)
#	print xa[0], getframeinfo(currentframe()).lineno
	redMass = scipy.zeros((NumBondList))
	for i in range(NumBondList):
		ii = bond_list_A[i]
		jj = bond_list_B[i]
		redMass[i] = 1/(inv_mass[ii] + inv_mass[jj])
#	print xa[0], getframeinfo(currentframe()).lineno
	code_shake =	  """
		//printf("%5.3f\\n",xa(51));
		double	ShakeAccuracy = 0.0001;
		//double        ShakeAccuracy = 1.0;
		int			ShakeNumber = 10000;
		double	Null = 0.00001;
		//printf("%5.3f\\n",xa(51));
    int			i, j0, j1, j2, k, n, allOK;
    int     bi, bj, bxi, bxj;
    double	*rprevij0, *rprevij1, *rprevij2, *rprev2, r2, rij0, rij1, rij2, fact, dr0, dr1, dr2, dz, dij;

    int     NumCoord = 3*num_residues;

    double	timestep = dt;
    //printf("%5.3f\\n",xa(51));
    rprevij0 = (double *) calloc(NumBondList, sizeof(double));
    rprevij1 = (double *) calloc(NumBondList, sizeof(double));
    rprevij2 = (double *) calloc(NumBondList, sizeof(double));

    rprev2 = (double *) calloc(NumBondList, sizeof(double));
    //printf("%5.3f\\n",xa(51));
    // 1. Distances of previous step
    for(i = 0; i < NumBondList; i++)
    {
            bi = bond_list_A(i);     bxi = 3*bi;
            bj = bond_list_B(i);     bxj = 3*bj;

            rprevij0[i] = xprev(bxi)   - xprev(bxj);
            rprevij1[i] = xprev(bxi+1) - xprev(bxj+1);
            rprevij2[i] = xprev(bxi+2) - xprev(bxj+2);

            rprev2[i] = rprevij0[i]*rprevij0[i] + rprevij1[i]*rprevij1[i] + rprevij2[i]*rprevij2[i];
            //rprev2[i] = dr2_list(bi, bj);
            //printf("%d:xprev0:% 5.3f\txprev1:% 5.3f\txprev2:% 5.3f\\n",bi,xprev(bxi),xprev(bxi+1),xprev(bxi+2));
            //printf("%d-%d:%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%f\\n",bi,bj,xprev(bxi),xprev(bxi+1),xprev(bxi+2),xprev(bxj),xprev(bxj+1),xprev(bxj+2),rprevij0[i],rprevij1[i],rprevij2[i],rprev2[i]);
            printf("cry_bond: %3.3f\\n",rprev2[i]);
    }
    //printf("%5.3f\\n",xa(51));
  



    // 2. Iterate
    // a. Protein + ligand
    allOK = 0;	n = 0;
    while(allOK == 0)
    {
	//printf("1: %.3f\\n",xa(51));
        allOK = 1;
        for(i = 0; i < NumBondList; i++)
        {
            		bi = bond_list_A(i);     bxi = 3*bi;
            		bj = bond_list_B(i);     bxj = 3*bj;

                rij0 = xa(bxi)   - xa(bxj);
                rij1 = xa(bxi+1) - xa(bxj+1);
                rij2 = xa(bxi+2) - xa(bxj+2);

                r2 = rij0*rij0 + rij1*rij1 + rij2*rij2;
		printf("eig_bond: %3.3f\\n",r2);
                fact = fabs((r2 - rprev2[i])/(2*r2));

		//printf("2: %.3f\\n",xa(51));

                if(fact > ShakeAccuracy)
                {
                   // printf("3: %.3f\\n",xa(51));
                    //printf("fact is greater that ShakeAccuracy\t");
                    allOK = 0;
                    if(n > ShakeNumber)
                    {
                        printf("Shake error ! Bond length between atoms %d %d exceeds shake limits:\\nActual length: % 8.4f  Equilibrium bond length: % 8.4f\\n",
                               bi+1, bj+1, sqrt(r2), sqrt(rprev2[i]));
												return 0;
                        //exit(0);
                    }
                    //printf("fact: %.3f\\n",fact);

                    dij = rprevij0[i]*rij0 + rprevij1[i]*rij1 + rprevij2[i]*rij2;
                    //printf("%d-%d\trij0: %.3f\trij1: %.3f\trij2: %.3f\\n",bxi,bxj,rij0,rij1,rij2);
                    //printf("dij: %.3f\trprevij0: %.3f\trprevij1: %.3f\trprevij2: %.3f\\n",dij,rprevij0[i],rprevij1[i],rprevij2[i]);
                    if(dij < Null)
                    {
                        dij = Null;
                    }

                    dz = redMass(i)*(rprev2[i] - r2)/(2.0*dij);
                    //printf("4: %.3f\\n",xa(51));
                    dr0 = dz*rprevij0[i];
                    dr1 = dz*rprevij1[i];
                    dr2 = dz*rprevij2[i];

                    xa(bxi)   += dr0*inv_mass(bi);
                    xa(bxi+1) += dr1*inv_mass(bi);
                    xa(bxi+2) += dr2*inv_mass(bi);
                    if(xa(bxi) >= 100.0 || xa(bxi+1) >= 100.00 || xa(bxi+2) >= 100.00)
                    {
			printf("%d: %.3f\t%.3f\\n",bxi,xa(bxi),xprev(bxi));
                        //exit(0);
                    }
                    //printf("5: %.3f\t%.3f\t%.3f\\n",xa(48),xa(51),xa(54));
                    xa(bxj)   -= dr0*inv_mass(bj);
                    xa(bxj+1) -= dr1*inv_mass(bj);
                    xa(bxj+2) -= dr2*inv_mass(bj);
                    //printf("6: %.3f\\n",xa(51));		
                }
        }
	//printf("%5.3f\\n",xa(51));
	//printf("%d\\n",n);
        n++;
	
    }
    //printf("shakeNumber:% 3d\\n",n);
    free(rprevij0);
    free(rprevij1);
    free(rprevij2);
    free(rprev2);

									"""
#	print xa[0], getframeinfo(currentframe()).lineno
	weave.inline(code_shake,
							['num_residues', 'dt', 'NumBondList', 'xa', 'xprev', 'bond_list_A', 'bond_list_B', 'dr2_list', 'inv_mass', 'redMass'],
							type_converters = converters.blitz,
							compiler='gcc')


###################################################################################################
###################################################################################################
def	runMD(dir_in, force_const, g_width, g_V0, dr_list, coord_list,  eig_coord_list, mass_list, num_steps, dT, bond_list_A, bond_list_B, shake_on, num_res, temp0, gam, freq_coors, freq_reassign_velocity):

	import scipy.weave
	import scipy.weave.converters

	print num_steps, dT, temp0, gam
	print "**************************", freq_coors
	# initialize temporary arrays
	va    = scipy.zeros((3*num_res))
	vhalf = scipy.zeros((3*num_res))
	vprev = scipy.zeros((3*num_res))
	xa    = scipy.zeros((3*num_res))
	xprev = scipy.zeros((3*num_res))
	fa    = scipy.zeros((3*num_res))
	dr2_list = scipy.zeros((num_res, num_res))
	for i in range(num_res):
		for j in range(num_res):
			dr2_list[i][j] = dr_list[i][j]*dr_list[i][j]
	

	#----------------------------------------------
	# initialize coordinates
	for i in range(num_res):
		for j in range(3):
			xa[3*i+j]    = eig_coord_list[3*i+j]
			xprev[3*i+j] = coord_list[3*i+j]
	#	print xa[3*i],xa[3*i+1],xa[3*i+2]
	# initialize velocities
	initializeVelocities(va, mass_list, temp0, num_res)

	# store velocities in array for previous time step: vprev
	for i in range(num_res):
		for j in range(3):
			vprev[3*i+j] = va[3*i+j]

	#----------------------------------------------
	# run over num_timesteps using BBK integrator
	# some constant variables
	minus_half_gamma_dt = 1 - 0.5*gam*dT
	plus_half_gamma_dt  = 1/(1 + 0.5*gam*dT)
#	random_factor = 0.00128957*math.sqrt(gamma*Temp0/dt)
	random_factor = 0.00128957*math.sqrt(3*gam*temp0*dT)
#	random_factor = 0.00128957*math.sqrt(gamma*Temp0*dt)
	force_factor = 0.00041869

#	print minus_half_gamma_dt, plus_half_gamma_dt

	# open coordinate file for output
	print os.getcwd()
	#fcoord = open("%s/output.crd"%dir_in, "w")
	fcoord = open("output.crd", "w")
	fcoord.write("test\n")
	nc = 0

	nv = 0
	for timestep in range(num_steps):
		
		# calculate forces
		#NR calculateForces(force_const, g_width, g_V0, dr_list, xa, fa, num_res)
#		for i in range(num_residues):
#			if abs(fa[3*i]) > 100.0 or abs(fa[3*i+1]) > 100.00 or abs(fa[3*i+2]) > 100.00:
#				print i, fa[3*i], fa[3*i+1], fa[3*i+2]
		#NR
		"""
		# first half kick
		for i in range(num_res):
			for j in range(3):
				ij = 3*i+j
#				rf = math.sqrt(dt/mass_list[i])*random_factor*numpy.random.normal(0.0, 1.0)
				rf = math.sqrt(1.0/mass_list[i])*random_factor*numpy.random.normal(0.0, 1.0)
				vhalf[ij] = minus_half_gamma_dt*vprev[ij] + 0.5*(force_factor*dT*fa[ij]/mass_list[i] + rf)
#				print "% 5d % 5d % 10.5f % 10.5f % 10.5f % 10.5f % 10.5f % 10.5f | % 10.5f" % (i, j, minus_half_gamma_dt, vprev[ij], force_factor*dt*fa[ij]/mass_list[i], fa[ij], rf, mass_list[i], vhalf[ij])
		# remove COM velocities
		#RemoveCOMVelocities(vhalf, mass_list, num_res)

		# constant T
#		berendsen/Thermostat(vhalf, mass_list, Temp0, dt, tau)
		CalculateTemperature(vhalf, mass_list, num_res)
#		print "Temp_a", timestep, Temp
		
		# drift
		for i in range(num_res):
			for j in range(3):
				ij = 3*i+j
				xprev[ij] = xa[ij]
#				print "% 5d % 5d % 10.5f % 10.5f" % (i, j, xa[ij], dt*vhalf[ij])
				xa[ij] += dT*vhalf[ij]
#				print "% 5d % 5d % 10.5f % 10.5f % 10.5f" % (i, j, xa[ij], dt, vhalf[ij])
		"""
#		for i in range(num_res):
#			for j in range(3):
#				print xa[3*i+j]
		# shake
		if shake_on == 1:
			shakeBonds(xa, xprev, mass_list, dT, dr2_list, bond_list_A, bond_list_B, num_res)
		#for i in range(num_res):
		#	print ("xa:% 5.3f% 5.3f% 5.3f"%(xa[3*i],xa[3*i+1],xa[3*i+2]))
		#NR	for i in range(num_res):
		#NR		for j in range(3):
		#NR			ij = 3*i+j
		#NR			vhalf[ij]	= (xa[ij] - xprev[ij])/dT
		#sys.exit(2)
		"""
		CalculateTemperature(vhalf, mass_list, num_res)
#		print "Temp_b", timestep, Temp

		# calculate forces
		calculateForces(force_const, g_width, g_V0, dr_list, xa, fa, num_res)
		# second half kick
		for i in range(num_res):
			for j in range(3):
				ij = 3*i+j
#				rf = math.sqrt(dt/mass_list[i])*random_factor*numpy.random.normal(0.0, 1.0)
				rf = math.sqrt(1.0/mass_list[i])*random_factor*numpy.random.normal(0.0, 1.0)
				va[ij] = plus_half_gamma_dt*(vhalf[ij] + 0.5*(force_factor*dT*fa[ij]/mass_list[i] + rf))
#				print "% 5d % 5d % 10.5f % 10.5f % 10.5f % 10.5f % 10.5f % 10.5f | % 10.5f" % (i, j, plus_half_gamma_dt, vhalf[ij], force_factor*dt*fa[ij]/mass_list[i], fa[ij], rf, mass_list[i], va[ij])

		CalculateTemperature(va, mass_list, num_res)
#		print "Temp_c", timestep, Temp

		"""
		# shake


		# write coordinates
		nc += 1
		if nc == freq_coors:
			print timestep
			writeCoordinates(fcoord, xa, num_res)
			nc = 0

		# reassign velocities
		nv += 1
		if nv == freq_reassign_velocity:
			initializeVelocities(va, mass_list, temp0, num_res)
			nv = 0

		# store velocities in array for previous time step: vprev
		for i in range(num_res):
			for j in range(3):
				vprev[3*i+j] = va[3*i+j]

	# close coordinate file
	fcoord.close()
		
###################################################################################################
###################################################################################################
def	readCommandFile(filename):
	global 	dt, num_timesteps, num_iterations, pdb_file, filemdp, filetdist, shake_on, initial_file, target_file
	global 	k_bond, k_bond_1_3, k_bond_1_4, k_helix_1_3, k_sheet_1_4, k_helix_1_4, k_hbonds, k_disulfide, k_hphob, k_density, r_cutoff, lambda_hbond, beta_hbond, lambda_hphob, beta_hphob, k_salt, gam_const, k_aromatic
	global	g_k_bond, g_k_helix_1_3, g_k_sheet_1_4, g_k_helix_1_4, g_k_hbonds, g_k_disulfide, g_k_hphob, g_k_density, g_r_cutoff, g_lambda_hbond, g_beta_hbond, g_lambda_hphob, g_beta_hphob
	global	g_V0_bond, g_V0_helix_1_3, g_V0_sheet_1_4, g_V0_helix_1_4, g_V0_hbonds, g_V0_disulfide, g_V0_hphob, g_V0_density
	global	tau, Temp0, gamma, freq_coord, freq_reassign_vel, min_move, max_move, move_step

	curdir = os.getcwd()

	fi = open(filename)
	for i in fi:
		line = i.strip().split()
		if len(line) > 0:
			# gamma_constant 
			if line[0] == "gamma_constant":
				gam_const = float(line[1])
			# k_bond
			if line[0] == "k_bond":
				k_bond = float(line[1])
			# k_bond_1_3
			if line[0] == "k_bond_1_3":
				k_bond_1_3 = float(line[1])
			# k_bond_1_4
			if line[0] == "k_bond_1_4":
				k_bond_1_4 = float(line[1])
			# k_helix_1_3
			if line[0] == "k_helix_1_3":
				k_helix_1_3 = float(line[1])
			# k_sheet_1_4
			if line[0] == "k_sheet_1_4":
				k_sheet_1_4 = float(line[1])
			# k_helix_1_4
			if line[0] == "k_helix_1_4":
				k_helix_1_4 = float(line[1])
			# k_hbonds
			if line[0] == "k_hbonds":
				k_hbonds = float(line[1])
			# k_disulfide
			if line[0] == "k_disulfide":
				k_disulfide = float(line[1])
			# k_salt
			if line[0] == "k_salt":
				k_salt = float(line[1])
			# k_hphob
			if line[0] == "k_hphob":
				k_hphob = float(line[1])
			# k_density
			if line[0] == "k_density":
				k_density = float(line[1])
			# k_aromatic
			if line[0] == "k_aromatic":
				k_aromatic = float(line[1])
			# r_cutoff
			if line[0] == "r_cutoff":
				r_cutoff = float(line[1])
			# lambda_hbond
			if line[0] == "lambda_hbond":
				lambda_hbond = float(line[1])
			# beta_hbond
			if line[0] == "beta_hbond":
				beta_hbond = float(line[1])
			# lambda_hphob
			if line[0] == "lambda_hphob":
				lambda_hphob = float(line[1])
			# beta_hphob
			if line[0] == "beta_hphob":
				beta_hphob = float(line[1])

	# gaussian potentials
	# width of gaussian (=1/alpha^2 in exponent)
			# g_k_bond
			if line[0] == "g_k_bond":
				g_k_bond = float(line[1])
			# g_k_helix_1_3
			if line[0] == "g_k_helix_1_3":
				g_k_helix_1_3 = float(line[1])
			# g_k_sheet_1_4
			if line[0] == "g_k_sheet_1_4":
				g_k_sheet_1_4 = float(line[1])
			# g_k_helix_1_4
			if line[0] == "g_k_helix_1_4":
				g_k_helix_1_4 = float(line[1])
			# g_k_hbonds
			if line[0] == "g_k_hbonds":
				g_k_hbonds = float(line[1])
			# g_k_disulfide
			if line[0] == "g_k_disulfide":
				g_k_disulfide = float(line[1])
			# g_k_hphob 
			if line[0] == "g_k_hphob":
				g_k_hphob = float(line[1])
			# g_k_density
			if line[0] == "g_k_density":
				g_k_density = float(line[1])
			# g_r_cutoff
			if line[0] == "g_r_cutoff":
				g_r_cutoff = float(line[1])
			# g_lambda_hbond
			if line[0] == "g_lambda_hbond":
				g_lambda_hbond = float(line[1])
			# g_beta_hbond
			if line[0] == "g_beta_hbond":
				g_beta_hbond = float(line[1])
			# g_lambda_hphob
			if line[0] == "g_lambda_hphob":
				g_lambda_hphob = float(line[1])
			# g_beta_hphob
			if line[0] == "g_beta_hphob":
				g_beta_hphob = float(line[1])
	# maximum potential (positive is repulsive)
			# g_V0_bond
			if line[0] == "g_V0_bond":
				g_V0_bond = float(line[1])
			# g_V0_helix_1_3
			if line[0] == "g_V0_helix_1_3":
				g_V0_helix_1_3 = float(line[1])
			# g_V0_sheet_1_4
			if line[0] == "g_V0_sheet_1_4":
				g_V0_sheet_1_4 = float(line[1])
			# g_V0_helix_1_4
			if line[0] == "g_V0_helix_1_4":
				g_V0_helix_1_4 = float(line[1])
			# g_V0_hbonds
			if line[0] == "g_V0_hbonds":
				g_V0_hbonds = float(line[1])
			# g_V0_disulfide
			if line[0] == "g_V0_disulfide":
				g_V0_disulfide = float(line[1])
			# g_V0_hphob
			if line[0] == "g_V0_hphob":
				g_V0_hphob = float(line[1])
			# g_V0_density
			if line[0] == "g_V0_density":
				g_V0_density = float(line[1])

			# timestep
			if line[0] == "dt":
				dt = float(line[1])
			# num_timesteps
			if line[0] == "num_timesteps":
				num_timesteps = int(line[1])
			# num_iterations
			if line[0] == "num_iterations":
				num_iterations = int(line[1])
			# pdb_file
			if line[0] == "pdb_file":
				pdb_file = line[1]
			if line[0] == "initial_file":
				initial_file = line[1]
			if line[0] == "target_file":
				target_file = line[1]
			# filemdp
			if line[0] == "filemdp":
				filemdp = line[1]
			# filetdist
			if line[0] == "filetdist":
				filetdist = line[1]
			# shake_on
			if line[0] == "shake_on":
				shake_on = int(line[1])
			# tau
			if line[0] == "tau":
				tau = float(line[1])
			# Temp0
			if line[0] == "Temp0":
				Temp0 = float(line[1])
			# gamma
			if line[0] == "gamma":
				gamma = float(line[1])
			# freq_coord
			if line[0] == "freq_coord":
				freq_coord = int(line[1])
			# freq_reassign_vel
			if line[0] == "freq_reassign_vel":
				freq_reassign_vel = int(line[1])
			# min move
			if line[0] == "min_move":
				min_move = float(line[1])
			# max move
			if line[0] == "max_move":
				max_move = float(line[1])
			# move step
			if line[0] == "move_step":
				move_step = float(line[1])
	fi.close()

	print k_bond
	print k_bond_1_3
	print k_bond_1_4
	print k_helix_1_4
	print k_sheet_1_4
	print k_hbonds
	print k_salt
	print k_disulfide
	print k_hphob
	print k_aromatic

###################################################################################################
###################################################################################################
def	alignTrajectory(dir_in):

#	os.popen("rm output_align.pdb").read()
#	os.chdir("%s/%s"%(curdir, dir_in))
	# print aligned trajectory file
	frun = open("ptraj_script", "w")
#	frun.write("babel -imol2 start_%d%d.mol2 -opdb start_%d%d.pdb &>/dev/null\n"%(int(dir_in[8:9]), int(dir_in[9:10]), int(dir_in[8:9]), int(dir_in[9:10])))
	
	frun.write("/usr/local/amber10/bin/ptraj start_%d%d.pdb &>/dev/null << EOF\n"%(int(dir_in[8:9]),int(dir_in[9:10])))
#	frun.write("trajin min_fin.pdb\n")
	frun.write("trajin output.crd\n")
	frun.write("rms first \n")
	frun.write("trajout output_align.crd\n")
	frun.write("trajout output_align.pdb PDB append\n")
	frun.write("go\n")
	frun.write("EOF\n")
	frun.close()
	os.popen("chmod a+x ptraj_script").read()
	os.popen("./ptraj_script").read()


###################################################################################################
###################################################################################################
def	computeSinglePointEnergy(md_outfile):
		# minimize PDB file
		MinimizeStructure(pdb_file, md_outfile, filemdp)

		# assign secondary structure
		secondary_structure = []
		assignSecondaryStructure("em.pdb", secondary_structure)

		# tCONCOORD analysis
		tCONCOORD_Analysis(filetdist)

		# ANM modes
		

		# read PDB file and determine number of residues
		identifyNumberOfResidues("em.pdb")

		# read PDB file and identify connecting residues (CA-CA) and disulfide bonds
		disulfide_A = []
		disulfide_B = []
		disulfide_A_chain = []
		disulfide_B_chain = []
		bonded_A = []
		bonded_B = []

		identifyBondedRestraints("em.pdb", bonded_A, bonded_B, disulfide_A, disulfide_B, disulfide_A_chain, disulfide_B_chain)
		# combine CA-CA bonds and disulfide bonds into one bonded list
		num_in_bonded_list = len(bonded_A) + len(disulfide_A)
		bond_list_A = scipy.zeros((num_in_bonded_list), int)
		bond_list_B = scipy.zeros((num_in_bonded_list), int)
		ci = 0
		for i, j in zip(bonded_A, bonded_B):
			ii = int(i)-1
			jj = int(j)-1
			bond_list_A[ci] = ii
			bond_list_B[ci] = jj
			ci += 1
		for i, j in zip(disulfide_A, disulfide_B):
			ii = int(i)-1
			jj = int(j)-1
			bond_list_A[ci] = ii
			bond_list_B[ci] = jj
			ci += 1
#		print disulfide_A, disulfide_B
#		print bonded_A
#		print bonded_B
#		print bond_list_A
#		print bond_list_B

		# read hbonds.dat file for hydrogen-bond and hyrdophobic contact identification
		hbonds_A = []
		hbonds_B = []
		num_hbonds = []
		hphob_A = []
		hphob_B = []
		num_hphob = []

		identifyNonBondedRestraints(file_pdb, "hbonds.dat", hbonds_A, hbonds_B, num_hbonds, hphob_A, hphob_B, num_hphob, saltbridge_A, saltbridge_B)
		# scan over parameter sets
		

		# generate force constant matrix
		force_const = scipy.zeros((num_residues, num_residues))
		g_width = scipy.zeros((num_residues, num_residues))
		g_V0 = scipy.zeros((num_residues, num_residues))
		shake_list = scipy.zeros((num_residues, num_residues))
		dr_list = scipy.zeros((num_residues, num_residues))
		coord_list = scipy.zeros((3*num_residues))
		mass_list = scipy.zeros((num_residues))
		generateForceConstantMatrix("em.pdb", force_const, g_width, g_V0, shake_on, shake_list, dr_list, coord_list, mass_list, bonded_A, bonded_B, disulfide_A, disulfide_B, hbonds_A, hbonds_B, num_hbonds, hphob_A, hphob_B, num_hphob, secondary_structure)

		

###################################################################################################
def	calculateNMA(new_coord_list,num_residues,force_const):
			hessian = scipy.zeros((3*num_residues,3*num_residues))
			#fout=open("%s.pdb.sparsehessian_bondCut"%pdb_file[:-4], "w")
			dr_list = scipy.zeros((num_residues,num_residues))
		
			code_hessian = """
								int jj, kk;
								double x1, x2, y1, y2, z1, z2, bx, by, bz, dr2;

								for(jj = 0; jj < num_residues; jj++)
								{
									x1 = new_coord_list(3*jj);
									y1 = new_coord_list(3*jj+1);
									z1 = new_coord_list(3*jj+2);
									for(kk = 0; kk < num_residues; kk++)
									{
										if(jj == kk)
											continue;
										// distance between atoms
										else
										{
											x2 = new_coord_list(3*kk);
											y2 = new_coord_list(3*kk+1);
											z2 = new_coord_list(3*kk+2);

											bx = x1 - x2;
											by = y1 - y2;
											bz = z1 - z2;

											dr2 = bx*bx+by*by+bz*bz;
											dr_list(jj,kk) = sqrt(dr2);
											dr_list(kk,jj) = sqrt(dr2);

											// diagonals of diagonal super elements (for j)
											hessian(3*jj,3*jj) += force_const(jj,kk)*bx*bx/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*jj+1) += force_const(jj,kk)*by*by/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*jj+2) += force_const(jj,kk)*bz*bz/(dr_list(jj,kk)*dr_list(jj,kk));

											// off-diagonals of diagonal superelements (for j)
											hessian(3*jj,3*jj+1) += force_const(jj,kk)*bx*by/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj,3*jj+2) += force_const(jj,kk)*bx*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*jj) += force_const(jj,kk)*by*bx/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*jj+2) += force_const(jj,kk)*by*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*jj) += force_const(jj,kk)*bx*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*jj+1) += force_const(jj,kk)*by*bz/(dr_list(jj,kk)*dr_list(jj,kk));

											// diagonals of off-diagonals superelements (for j&k)
											hessian(3*jj,3*kk) = -1.0*force_const(jj,kk)*bx*bx/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*kk+1) = -1.0*force_const(jj,kk)*by*by/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*kk+2) = -1.0*force_const(jj,kk)*bz*bz/(dr_list(jj,kk)*dr_list(jj,kk));
	
											// off-diagonals of off-diagonal superelements (for j&k)
											hessian(3*jj,3*kk+1) = -1.0*force_const(jj,kk)*bx*by/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj,3*kk+2) = -1.0*force_const(jj,kk)*bx*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*kk) = -1.0*force_const(jj,kk)*by*bx/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*kk+2) = -1.0*force_const(jj,kk)*by*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*kk) = -1.0*force_const(jj,kk)*bx*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*kk+1) = -1.0*force_const(jj,kk)*by*bz/(dr_list(jj,kk)*dr_list(jj,kk));
										}
									}
								}
								"""								
			scipy.weave.inline(code_hessian,
												['num_residues', 'force_const', 'new_coord_list', 'dr_list', 'hessian'],
												type_converters = scipy.weave.converters.blitz,
												compiler='gcc')
			w,v = linalg.eig(hessian)
			idx = w.argsort()
			w = w[idx]
			v = v[:,idx]
			eigen_outfile = open("eigen_values.txt","a")
			for idx in range(len(w)):
				eigen_outfile.write("%f\n"%w[idx].real)
			eigen_outfile.close()

###################################################################################################
###################################################################################################
def	runLangevinDynamics(pdb_file):

#	for it in range(num_iterations):
	global centroid_coordinates, all_centroids, func_grps
	for it in range(num_iterations,num_iterations+1):
		#print it
		start_pdbs = []
		if it == 0:
#			start_pdbs.append("%s/%s"%(curdir, pdb_file))
			start_pdbs.append("%s"% pdb_file)
		else:
			tmp_pdbs = os.listdir(curdir)
			for pdbs in tmp_pdbs:
				if pdbs.find("ENM_%d" % (int(it)-1)) >=0 and pdbs.find(".pdb") >= 0 and len(pdbs) == 10:
#					start_pdbs.append("%s/%s"%(curdir,pdbs))
#					start_pdbs.append("%s"%pdbs) 			
					start_pdbs.append("%s"%pdb_file)
		start_pdbs.sort()
		#print start_pdbs
#		for pdb_num in range(len(start_pdbs)):
		for pdb_num in range(file_num,file_num+1):

			prepareProtein4ENM(start_pdbs[pdb_num], "prepared4ENM_%d%d.pdb"%(num_iterations,file_num))
			file_pdb = "prepared4ENM_%d%d.pdb"%(num_iterations,file_num)
			# minimize PDB file
			if it >= 1:
				#MinimizeStructure(start_pdbs[pdb_num], md_outfile, filemdp)
				md_outfile = os.getcwd()
				MinimizeStructure(file_pdb, md_outfile, filemdp)

			# assign secondary structure
			int_resn = []
			sheet_resn = []
			secondary_structure = []
			# assignSecondaryStructure("em.pdb", secondary_structure)
			assignSecondaryStructure(file_pdb, secondary_structure, int_resn, sheet_resn)
			# super_secondary_structure.append((secondary_structure))

			# tCONCOORD analysis
			tCONCOORD_Analysis(filetdist, file_pdb)


			# file_pdb = start_pdbs[pdb_num]
			# HBPlus hydrogen bonds
			hbplus_Analysis(file_pdb)
			# file_pdb = "prepared4ENM_%d%d.pdb"%(num_iterations,file_num)
			#sys.exit(2)

			num_residues = 0
			# read PDB file and determine number of residues
			num_residues = identifyNumberOfResidues(file_pdb)
			num_prot_residues = num_residues
			print "Number of protein residues", num_residues
			
			# read PDB file and identify connecting residues (CA-CA) and disulfide bonds
			disulfide_A = []
			disulfide_B = []
			bonded_A = []
			bonded_B = []

			identifyBondedRestraints(file_pdb, bonded_A, bonded_B, disulfide_A, disulfide_B)
			# combine CA-CA bonds and disulfide bonds into one bonded list
			num_in_bonded_list = len(bonded_A) + len(disulfide_A)
			bond_list_A = scipy.zeros((num_in_bonded_list), int)
			bond_list_B = scipy.zeros((num_in_bonded_list), int)
			ci = 0
			for i, j in zip(bonded_A, bonded_B):
				ii = int(i)-1
				jj = int(j)-1
				bond_list_A[ci] = ii
				bond_list_B[ci] = jj
				ci += 1
			for i, j in zip(disulfide_A, disulfide_B):
				ii = int(i)-1
				jj = int(j)-1
				bond_list_A[ci] = ii
				bond_list_B[ci] = jj
				ci += 1
			
			# Prepare a 3d_points.txt file of ligand coordinates
			# Find coarse grained points (like cluster centroids)
			# Find groups in ligand 
			# Map functional groups to CG points
			salt_bridge_centroids_cat = []
			salt_bridge_centroids_an = []
			all_centroids = []
			centroid_coordinates = []
			for z in range(len(lig_names)):
				L = Ligand(lig_names[z], lig_nums[z])
				#centroid_coords, labels = L.userdefined(lig_coords, "centroid.txt")
				#num_lig_atoms = L.numAtoms
				#centroid_coords,labels = L.KMeansClusterCentroids_Ligand(file_pdb, lig_names[z], lig_nums[z], lig_coords)
				centroid_coords,labels = L.KMedoids_Ligand(file_pdb, lig_names[z], lig_nums[z], lig_coords)
				#centroid_coords,labels = L.AffinityPropagation_Ligand(lig_coords)
				#centroid_coords,labels = L.HierarchialClusterCentroids_Ligand(lig_coords, num_lig_atoms)
				#centroid_coords,labels = L.DBSCANClusterCentroids_Ligand(lig_coords)

				#centroids = MapCentroidsToLigGroups("hbonds.dat", labels, num_residues)
				centroids = MapCentroidsToLigGroups("%s.hb2"%file_pdb[:-4], file_pdb, labels, num_residues, lig_names[z], lig_nums[z]) 
				#salt_bridge_centroids_pos, salt_bridge_centroids_neg = IdentifySaltbridgesBwProtLig(file_pdb, centroids, centroid_coords, lig_coords, lig_names[z], lig_nums[z], num_residues)
				func_grps = L.findFuncGroups()
				salt_bridge_centroids_pos, salt_bridge_centroids_neg = IdentifySaltbridgesBwProtLig(func_grps, centroids, centroid_coords, lig_coords, lig_names[z], lig_nums[z], num_residues)
				for y in salt_bridge_centroids_pos:
					salt_bridge_centroids_cat.append(num_residues - num_prot_residues + y)
				for z in salt_bridge_centroids_neg:
					salt_bridge_centroids_an.append(num_residues - num_prot_residues + z)

				lig_aromatic_cent = IdentifyAromaticBwProtLig(func_grps, centroid_coords, lig_coords)

				num_residues += len(centroids)
				for y,z in zip(centroids, centroid_coords):
					print y, z
					all_centroids.append(y)
				#	centroid_coordinates.append(z.tolist())	
					centroid_coordinates.append(z)

				# Add centroids to ENM_res_list
				for l in range(len(centroids)):
					ENM_res_list.append('GLY')									
			print centroid_coordinates 
			
			
			# read hbonds.dat file for hydrogen-bond and hyrdophobic contact identification
			hbonds_A = []
			hbonds_B = []
			num_hbonds = []
			hphob_A = []
			hphob_B = []
			num_hphob = []
			saltbridge_A = [] 
			saltbridge_B = []
			aromatic_A = []
			aromatic_B = []

			#identifyNonBondedRestraints(file_pdb, "hbonds_new.dat", saltbridge_A, saltbridge_B, hbonds_A, hbonds_B, num_hbonds, hphob_A, hphob_B, num_hphob, salt_bridge_centroids_cat, salt_bridge_centroids_an, centroid_coords, centroids)
			identifyNonBondedRestraints(file_pdb, "hbonds.dat", saltbridge_A, saltbridge_B, hbonds_A, hbonds_B, num_hbonds, hphob_A, hphob_B, num_hphob, salt_bridge_centroids_cat, salt_bridge_centroids_an, lig_aromatic_cent, aromatic_A, aromatic_B, centroid_coordinates, all_centroids)
			#print hbonds_A, hbonds_B
			# generate force constant matrix
			force_const = scipy.zeros((num_residues, num_residues))
			g_width = scipy.zeros((num_residues, num_residues))
			g_V0 = scipy.zeros((num_residues, num_residues))
			shake_list = scipy.zeros((num_residues, num_residues))
			dr_list = scipy.zeros((num_residues, num_residues))
			coord_list = scipy.zeros((3*num_residues))
			mass_list = scipy.zeros((num_residues))
			generateForceConstantMatrix(file_pdb, force_const, g_width, g_V0, shake_on, shake_list, dr_list, coord_list, mass_list, centroids_mass, bonded_A, bonded_B, disulfide_A, disulfide_B, saltbridge_A, saltbridge_B, hbonds_A, hbonds_B, num_hbonds, hphob_A, hphob_B, num_hphob, aromatic_A, aromatic_B, secondary_structure, int_resn, sheet_resn, all_centroids, centroid_coordinates)
			violations = []
#			generatePDB("test.pdb", coord_list, violations, pdb_file)

			# perturb coordinates
			# initialize temporary arrays
			va    = scipy.zeros((3*num_residues))
			vhalf = scipy.zeros((3*num_residues))
			vprev = scipy.zeros((3*num_residues))
			xa    = scipy.zeros((3*num_residues))
			xprev = scipy.zeros((3*num_residues))
			xa_plus = scipy.zeros((3*num_residues))
			xa_minus = scipy.zeros((3*num_residues))
			fa    = scipy.zeros((3*num_residues))
			dr2_list = scipy.zeros((num_residues, num_residues))
			for i in range(num_residues):
				for j in range(num_residues):
					dr2_list[i][j] = dr_list[i][j]*dr_list[i][j]


			#----------------------------------------------
			# initialize coordinates
			for i in range(num_residues):
				for j in range(3):
					xa[3*i+j]    = coord_list[3*i+j]
					xprev[3*i+j] = coord_list[3*i+j]

			# initialize velocities
			initializeVelocities(va, mass_list, Temp0, num_residues)

			# store velocities in array for previous time step: vprev
			for i in range(num_residues):
				for j in range(3):
					vprev[3*i+j] = va[3*i+j]

			# some constant variables
			minus_half_gamma_dt = 1 - 0.5*gamma*dt
			plus_half_gamma_dt  = 1/(1 + 0.5*gamma*dt)
#			random_factor = 0.00128957*math.sqrt(gamma*Temp0/dt)
			random_factor = 0.00128957*math.sqrt(3*gamma*Temp0*dt)
#			random_factor = 0.00128957*math.sqrt(gamma*Temp0*dt)
			force_factor = 0.00041869

			for i in range(num_residues):
				for j in range(3):
					ij = 3*i+j
#					rf = math.sqrt(dt/mass_list[i])*random_factor*numpy.random.normal(0.0, 1.0)
					#rf = math.sqrt(1.0/mass_list[i])*random_factor*numpy.random.normal(0.0, 1.0)
#					vhalf[ij] = minus_half_gamma_dt*vprev[ij] + 0.5*(force_factor*dt*fa[ij]/mass_list[i] + rf)
					vhalf[ij] = minus_half_gamma_dt*vprev[ij]*50
				print "vel0:\t %d\t %f\t %f\t %f\t   %f" % (i, vhalf[3*i], vhalf[3*i+1], vhalf[3*i+2], Temp0)

			for i in range(num_residues):
				for j in range(3):
					ij = 3*i+j
#					print "% 5d % 5d % 10.5f % 10.5f" % (i, j, xa[ij], dt*vhalf[ij])
					xa_plus[ij] = xa[ij] + dt*vhalf[ij]
					xa_minus[ij] = xa[ij] - dt*vhalf[ij]
					
			# generateMol2
			generatePDB("originalCoords.pdb", xa, [], file_pdb)	
			generatePDB("perturbedCoords_plus.pdb", xa_plus, [], file_pdb)
			# generate hessian matrix	
			print "Printing hessian matrix"
			print num_residues, len(coord_list), len(force_const), len(dr_list)	
			hessian = scipy.zeros((3*num_residues,3*num_residues))
			fout=open("%s.pdb.sparsehessian_bondCut"%pdb_file[:-4], "w")

			for jj in range(num_residues):
				for kk in range(num_residues):
					fout.write("%10.3f\t"%(force_const[jj][kk]))
				fout.write("\n")
			
			code_hessian = """
								int jj, kk;
								double x1, x2, y1, y2, z1, z2, bx, by, bz;
								double bx_plus, by_plus, bz_plus, bx_minus, by_minus, bz_minus;

								for(jj = 0; jj < num_residues; jj++)
								{
									x1 = coord_list(3*jj);
									y1 = coord_list(3*jj+1);
									z1 = coord_list(3*jj+2);

									for(kk = 0; kk < num_residues; kk++)
									{
										if(jj == kk)
											continue;
										// distance between atoms
										else
										{
											x2 = coord_list(3*kk);
											y2 = coord_list(3*kk+1);
											z2 = coord_list(3*kk+2);
				
											bx = x1 - x2;
											by = y1 - y2;
											bz = z1 - z2;

											bx_plus = xa_plus(3*jj) - xa_plus(3*kk);
											by_plus = xa_plus(3*jj+1) - xa_plus(3*kk+1);
											bz_plus = xa_plus(3*jj+2) - xa_plus(3*kk+2);
											bx_minus = xa_minus(3*jj) - xa_minus(3*kk);
											by_minus = xa_minus(3*jj+1) - xa_minus(3*kk+1);
											bz_minus = xa_minus(3*jj+2) - xa_minus(3*kk+2);

											// diagonals of diagonal super elements (for j)
											//hessian(3*jj,3*jj) += 0.5*force_const(jj,kk)*(bx_plus*bx_plus - 2*bx*bx + bx_minus*bx_minus )/((xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2)*(xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2));
											//hessian(3*jj+1,3*jj+1) += 0.5*force_const(jj,kk)*(by_plus*by_plus - 2*by*by + by_minus*by_minus)/((xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) -y2)*(xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) -y2));
											//hessian(3*jj+2,3*jj+2) += 0.5*force_const(jj,kk)*(bz_plus*bz_plus - 2*bz*bz + bz_minus*bz_minus)/((xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2)*(xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2));

											// off-diagonals of diagonal superelements (for j)
											//hessian(3*jj,3*jj+1) +=  0.5*force_const(jj,kk)* (bx_plus*bx_plus + by_plus*by_plus - bx_plus*by_minus - bx_minus*by_plus + bx_minus*bx_minus + by_minus*by_minus )/((xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2)*(xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) - y2));
											//hessian(3*jj,3*jj+2) +=  0.5*force_const(jj,kk)* (bx_plus*bx_plus + bz_plus*bz_plus - bx_plus*bz_minus - bx_minus*bz_plus + bx_minus*bx_minus + bz_minus*bz_minus )/((xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2)*(xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2));
											//hessian(3*jj+1,3*jj) += 0.5*force_const(jj,kk)* (by_plus*by_plus + bx_plus*bx_plus - by_plus*bx_minus - by_minus*bx_plus + by_minus*by_minus + bx_minus*bx_minus )/((xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) - y2)*(xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2));
											//hessian(3*jj+1,3*jj+2) += 0.5*force_const(jj,kk)* (by_plus*by_plus + bz_plus*bz_plus - by_plus*bz_minus - by_minus*bz_plus + by_minus*by_minus + bz_minus*bz_minus )/((xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) - y2)*(xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2));
											//hessian(3*jj+2,3*jj) += 0.5*force_const(jj,kk)* (bz_plus*bz_plus + bx_plus*bx_plus - bx_plus*bz_minus - bx_minus*bz_plus  + bz_minus*bz_minus + bx_minus*bx_minus )/((xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2)*(xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2));
											//hessian(3*jj+2,3*jj+1) += 0.5*force_const(jj,kk)* (bz_plus*bz_plus + by_plus*by_plus - by_plus*bz_minus - by_minus*bz_plus + bz_minus*bz_minus + by_minus*by_minus)/((xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2)*(xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) - y2));
											
											// diagonals of off-diagonals superelements (for j&k)
											//hessian(3*jj,3*kk) = -1.0*force_const(jj,kk)*(bx_plus*bx_plus - 2*bx*bx + bx_minus*bx_minus )/((xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2)*(xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2));
											//hessian(3*jj+1,3*kk+1) = -1.0*force_const(jj,kk)*(by_plus*by_plus - 2*by*by + by_minus*by_minus)/((xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) -y2)*(xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) -y2));
											//hessian(3*jj+2,3*kk+2) = -1.0*force_const(jj,kk)*(bz_plus*bz_plus - 2*bz*bz + bz_minus*bz_minus)/((xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2)*(xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2));


											// off-diagonals of off-diagonal superelements (for j&k)
											//hessian(3*jj,3*kk+1) = -1.0*force_const(jj,kk)* (bx_plus*bx_plus + by_plus*by_plus - bx_plus*by_minus - bx_minus*by_plus + bx_minus*bx_minus + by_minus*by_minus ) / ((xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2)*(xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) - y2));
											//hessian(3*jj,3*kk+2) = -1.0*force_const(jj,kk)* (bx_plus*bx_plus + bz_plus*bz_plus - bx_plus*bz_minus - bx_minus*bz_plus + bx_minus*bx_minus + bz_minus*bz_minus ) / ((xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2)*(xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2));
											//hessian(3*jj+1,3*kk) = -1.0*force_const(jj,kk)* (by_plus*by_plus + bx_plus*bx_plus - by_plus*bx_minus - by_minus*bx_plus + by_minus*by_minus + bx_minus*bx_minus ) / ((xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) - y2)*(xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2));
											//hessian(3*jj+1,3*kk+2) += 0.5*force_const(jj,kk)* (by_plus*by_plus + bz_plus*bz_plus - by_plus*bz_minus - by_minus*bz_plus + by_minus*by_minus + bz_minus*bz_minus ) / ((xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) - y2)*(xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2));
											//hessian(3*jj+2,3*kk) += 0.5*force_const(jj,kk)* (bz_plus*bz_plus + bx_plus*bx_plus - bx_plus*bz_minus - bx_minus*bz_plus  + bz_minus*bz_minus + bx_minus*bx_minus ) / ((xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2)*(xa_plus(3*jj) - x1 + xa_plus(3*kk) - x2));
											//hessian(3*jj+2,3*kk+1) += 0.5*force_const(jj,kk)* (bz_plus*bz_plus + by_plus*by_plus - by_plus*bz_minus - by_minus*bz_plus + bz_minus*bz_minus + by_minus*by_minus) / ((xa_plus(3*jj+2) - z1 + xa_plus(3*kk+2) - z2)*(xa_plus(3*jj+1) - y1 + xa_plus(3*kk+1) - y2));

											// diagonals of diagonal super elements (for j)
											//hessian(3*jj,3*jj) += force_const(jj,kk)*(xa_plus(3*kk)-xa_plus(3*jj))*(xa_plus(3*kk)-xa_plus(3*jj))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+1,3*jj+1) += force_const(jj,kk)*(xa_plus(3*kk+1)-xa_plus(3*jj+1))*(xa_plus(3*kk+1)-xa_plus(3*jj+1))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+2,3*jj+2) += force_const(jj,kk)*(xa_plus(3*kk+2)-xa_plus(3*jj+2))*(xa_plus(3*kk+2)-xa_plus(3*jj+2))/(dr_list(jj,kk)*dr_list(jj,kk));

											// off-diagonals of diagonal superelements (for j)
											//hessian(3*jj,3*jj+1) += force_const(jj,kk)*(xa_plus(3*kk)-xa_plus(3*jj))*(xa_plus(3*kk+1)-xa_plus(3*jj+1))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj,3*jj+2) += force_const(jj,kk)*(xa_plus(3*kk)-xa_plus(3*jj))*(xa_plus(3*kk+2)-xa_plus(3*jj+2))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+1,3*jj) += force_const(jj,kk)*(xa_plus(3*kk+1)-xa_plus(3*jj+1))*(xa_plus(3*kk)-xa_plus(3*jj))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+1,3*jj+2) += force_const(jj,kk)*(xa_plus(3*kk+1)-xa_plus(3*jj+1))*(xa_plus(3*kk+2)-xa_plus(3*jj+2))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+2,3*jj) += force_const(jj,kk)*(xa_plus(3*kk+2)-xa_plus(3*jj+2))*(xa_plus(3*kk)-xa_plus(3*jj))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+2,3*jj+1) += force_const(jj,kk)*(xa_plus(3*kk+2)-xa_plus(3*jj+2))*(xa_plus(3*kk+1)-xa_plus(3*jj+1))/(dr_list(jj,kk)*dr_list(jj,kk));

											// diagonals of off-diagonals superelements (for j&k)
											//hessian(3*jj,3*kk) = -1.0*force_const(jj,kk)*(xa_plus(3*kk)-xa_plus(3*jj))*(xa_plus(3*kk)-xa_plus(3*jj))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+1,3*kk+1) = -1.0*force_const(jj,kk)*(xa_plus(3*kk+1)-xa_plus(3*jj+1))*(xa_plus(3*kk+1)-xa_plus(3*jj+1))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+2,3*kk+2) = -1.0*force_const(jj,kk)*(xa_plus(3*kk+2)-xa_plus(3*jj+2))*(xa_plus(3*kk+2)-xa_plus(3*jj+2))/(dr_list(jj,kk)*dr_list(jj,kk));

											// off-diagonals of off-diagonal superelements (for j&k)
											//hessian(3*jj,3*kk+1) = -1.0*force_const(jj,kk)*(xa_plus(3*kk)-xa_plus(3*jj))*(xa_plus(3*kk+1)-xa_plus(3*jj+1))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj,3*kk+2) = -1.0*force_const(jj,kk)*(xa_plus(3*kk)-xa_plus(3*jj))*(xa_plus(3*kk+2)-xa_plus(3*jj+2))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+1,3*kk) = -1.0*force_const(jj,kk)*(xa_plus(3*kk+1)-xa_plus(3*jj+1))*(xa_plus(3*kk)-xa_plus(3*jj))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+1,3*kk+2) = -1.0*force_const(jj,kk)*(xa_plus(3*kk+1)-xa_plus(3*jj+1))*(xa_plus(3*kk+2)-xa_plus(3*jj+2))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+2,3*kk) = -1.0*force_const(jj,kk)*(xa_plus(3*kk+2)-xa_plus(3*jj+2))*(xa_plus(3*kk)-xa_plus(3*jj))/(dr_list(jj,kk)*dr_list(jj,kk));
											//hessian(3*jj+2,3*kk+1) = -1.0*force_const(jj,kk)*(xa_plus(3*kk+2)-xa_plus(3*jj+2))*(xa_plus(3*kk+1)-xa_plus(3*jj+1))/(dr_list(jj,kk)*dr_list(jj,kk));

											// diagonals of diagonal super elements (for j)
											hessian(3*jj,3*jj) += force_const(jj,kk)*bx*bx/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*jj+1) += force_const(jj,kk)*by*by/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*jj+2) += force_const(jj,kk)*bz*bz/(dr_list(jj,kk)*dr_list(jj,kk));

											// off-diagonals of diagonal superelements (for j)
											hessian(3*jj,3*jj+1) += force_const(jj,kk)*bx*by/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj,3*jj+2) += force_const(jj,kk)*bx*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*jj) += force_const(jj,kk)*by*bx/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*jj+2) += force_const(jj,kk)*by*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*jj) += force_const(jj,kk)*bx*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*jj+1) += force_const(jj,kk)*by*bz/(dr_list(jj,kk)*dr_list(jj,kk));

											// diagonals of off-diagonals superelements (for j&k)
											hessian(3*jj,3*kk) = -1.0*force_const(jj,kk)*bx*bx/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*kk+1) = -1.0*force_const(jj,kk)*by*by/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*kk+2) = -1.0*force_const(jj,kk)*bz*bz/(dr_list(jj,kk)*dr_list(jj,kk));
	
											// off-diagonals of off-diagonal superelements (for j&k)
											hessian(3*jj,3*kk+1) = -1.0*force_const(jj,kk)*bx*by/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj,3*kk+2) = -1.0*force_const(jj,kk)*bx*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*kk) = -1.0*force_const(jj,kk)*by*bx/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+1,3*kk+2) = -1.0*force_const(jj,kk)*by*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*kk) = -1.0*force_const(jj,kk)*bx*bz/(dr_list(jj,kk)*dr_list(jj,kk));
											hessian(3*jj+2,3*kk+1) = -1.0*force_const(jj,kk)*by*bz/(dr_list(jj,kk)*dr_list(jj,kk));
										}
									}
								}
								"""								
			scipy.weave.inline(code_hessian,
												['num_residues', 'force_const', 'coord_list', 'dr_list', 'hessian', 'xa_plus', 'xa_minus'],
												type_converters = scipy.weave.converters.blitz,
												compiler='gcc')

			
			# sparse hessian matrix
			#	for jj in range(num_residues):
			#		j = jj+1
			#		fout.write("%8d%8d%25.15e\n"%(3*j-2, 3*j-2, hessian[3*j-2][1]))
			#		fout.write("%8d%8d%25.15e\n"%(3*j-2, 3*j-1, hessian[3*j-2][2]))
			#		fout.write("%8d%8d%25.15e\n"%(3*j-2, 3*j, hessian[3*j-2][3]))
			#		fout.write("%8d%8d%25.15e\n"%(3*j-1, 3*j-1, hessian[3*j-1][2]))
			#		fout.write("%8d%8d%25.15e\n"%(3*j-1, 3*j, hessian[3*j-1][3]))
			#		fout.write("%8d%8d%25.15e\n"%(3*j, 3*j, hessian[3*j][3]))
			
			# full hessian matrix
			#new_hessian = scipy.zeros((3*num_residues,3*num_residues))
			print num_residues
			for jj in range(3*num_residues):
				for kk in range(3*num_residues):
					#new_hessian[jj][kk] = hessian[jj][kk]
					fout.write("%10.3f\t%d"%(hessian[jj][kk],kk))
				fout.write("\n")
			
			fout.close()

		
			#	mass_file = open("mass_file.txt","w")
			#	for m in range(num_residues*3):
			#		for n in range(num_residues*3):
			#			mass_file.write("%7.3f"%mass_matrix[m][n])
			#		mass_file.write("\n")
			#	mass_file.close()
			
			# return hessian
		
			sys.setrecursionlimit(20000)
			job_server = pp.Server()
			ncpu = job_server.get_ncpus()
			print "number of cpus detected: %d\nwill use all of them\n" %ncpu
			job_server.set_ncpus(ncpu)
			jobs = []
			
			print "Calculating eigen vectors"
			#print hessian
			# Calculating modes and new conformation
			#w, v = linalg.eigs(hessian, k=10)
			w,v = linalg.eig(hessian)
			idx = w.argsort()
			w = w[idx]
			v = v[:,idx]
			#for idx in range(len(w)):
			#	print w[idx].real
			for idx in range(len(w)):
				if w[idx].real > 0.0000000009:
					 break
			print idx
			
			norm_v = scipy.zeros((3*(num_residues-len(all_centroids)),3*num_residues))
			print norm_v.shape
			for val in range(3*num_residues): # CHANGE NO OF EIGENMODES, CURRENT = 10
				magn_v = 0.000000
				for z in range(num_residues-len(all_centroids)):
					magn_v += v[3*z][val].real*v[3*z][val].real + v[3*z+1][val].real*v[3*z+1][val].real + v[3*z+2][val].real*v[3*z+2][val].real
				magn_v = math.sqrt(magn_v)
				for z in range(num_residues-len(all_centroids)):
					norm_v[3*z][val] = v[3*z][val].real/magn_v
					norm_v[3*z+1][val] = v[3*z+1][val].real/magn_v
					norm_v[3*z+2][val] = v[3*z+2][val].real/magn_v

			print "norm_v ready\n", norm_v.shape
			fo = open("evec_%d%d_%s-%s-%02.1f_hb2wL.mol2"%(num_iterations,file_num,initial_file[-8:-4],target_file[-8:-4],r_cutoff),"w")

			for val in range(idx,idx+10): # CHANGE NO OF EIGENMODES, CURRENT = 10
				# write header
				fo.write("@<TRIPOS>MOLECULE\n")
				fo.write("MODEL%d\n"%(val-idx))
				#fo.write("% 5d % 5d     0     0     0\n" % ((num_residues-len(all_centroids))*2, num_residues-len(all_centroids)))
				fo.write("% 5d % 5d     0     0     0\n" % (num_residues*2, num_residues))
				fo.write("SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");
				# write coordinates
				for z in range(num_residues-len(all_centroids)): # first all structures in one folder, then next folder
#					print coord_list[3*i], norm_v[3*i]
					tmp_x = coord_list[3*z] + norm_v[3*z][val]*math.sqrt(num_residues-len(all_centroids))
					tmp_y = coord_list[3*z+1] + norm_v[3*z+1][val]*math.sqrt(num_residues-len(all_centroids))
					tmp_z = coord_list[3*z+2] + norm_v[3*z+2][val]*math.sqrt(num_residues-len(all_centroids))
					fo.write("%s%s % 9.3f% 8.3f% 8.3f  O.3     % 4d %3s    0.0000\n" % (str(2*z+1).ljust(5),"CA".rjust(7), coord_list[3*z], coord_list[3*z+1], coord_list[3*z+2], 2*z+1, ENM_res_list[z]))
					fo.write("%s%s % 9.3f% 8.3f% 8.3f  C.3     % 4d %3s    0.0000\n" % (str(2*z+2).ljust(5),"CA".rjust(7), tmp_x, tmp_y, tmp_z, 2*z+2, ENM_res_list[z]))
				x = z+1
				print x
				for z in range(len(all_centroids)):
					print z
					tmp_x = centroid_coordinates[z][0] + v[3*(x+z)][val].real*math.sqrt(num_residues)
					tmp_y = centroid_coordinates[z][1] + v[3*(x+z)+1][val].real*math.sqrt(num_residues)
					tmp_z = centroid_coordinates[z][2] + v[3*(x+z)+2][val].real*math.sqrt(num_residues)
					fo.write("%s%s % 9.3f% 8.3f% 8.3f  O.3     % 4d %3s    0.0000\n" % (str(2*(x+z)+1).ljust(5),"LA".rjust(7), centroid_coordinates[z][0], centroid_coordinates[z][1], centroid_coordinates[z][2], 2*(x+z)+1, "GLY"))
					fo.write("%s%s % 9.3f% 8.3f% 8.3f  C.3     % 4d %3s    0.0000\n" % (str(2*(x+z)+2).ljust(5),"LA".rjust(7), tmp_x, tmp_y, tmp_z, 2*(x+z)+2, "GLY"))
				# write bonds
				fo.write("@<TRIPOS>BOND\n")
				nb = 0
				for z in range(num_residues-len(all_centroids)):
					fo.write("% 5d %5d %5d    1\n" % (nb+1, 2*z+1, 2*z+2))
					nb += 1
				for x in range(len(all_centroids)):
					fo.write("% 5d %5d %5d    1\n" % (nb+1, 2*z+1, 2*z+2))
					z += 1
					nb += 1
				fo.write("\n")

			fo.close()	
			"""	
			feigen = open("eigenvalues_%d%d.txt"%(num_iterations,file_num),"w")
			eigen = []
			for val in range(idx,idx+5): # CHANGE NO OF EIGENMODES, CURRENT = 10
				eigen.append(w[val].real)
				feigen.write("%f\n%f\n"%(w[val].real,w[val].real))
			feigen.close()
			print eigen
			#sys.exit(2)
			
			angle_btw_vec = scipy.zeros((5,num_residues-1)) #scipy.zeros((10,num_residues-1)) # CHANGE NO OF EIGENMODES, CURRENT = 10
			for val in range(idx,idx+5):# CHANGE NO OF EIGENMODES, CURRENT = 10
				for i in range(num_residues-1):
					norm_fac = math.sqrt(v[3*i][val-idx].real*v[3*i][val-idx].real + v[3*i+1][val-idx].real*v[3*i+1][val-idx].real + v[3*i+2][val-idx].real*v[3*i+2][val-idx].real) * math.sqrt(v[3*(i+1)][val-idx].real*v[3*(i+1)][val-idx].real + v[3*(i+1)+1][val-idx].real*v[3*(i+1)+1][val-idx].real + v[3*(i+1)+2][val-idx].real*v[3*(i+1)+2][val-idx].real)
					angle_btw_vec[val-idx][i] = (v[3*i][val-idx].real*v[3*(i+1)][val-idx].real + v[3*i+1][val-idx].real*v[3*(i+1)+1][val-idx].real + v[3*i+2][val-idx].real*v[3*(i+1)+2][val-idx].real)/norm_fac

			#for val in range(idx,idx+100):
			#	print angle_btw_vec[val-idx]
			#sys.exit(2)
			# check if angle between two adjacent eigenvectors of any pair of atoms > 10.0e-8			
			ang_mean = []
			for val in range(idx,idx+5): # CHANGE NO OF EIGENMODES, CURRENT = 10
				mean = 0.0
				for ang in range(len(angle_btw_vec[val-idx])):
					mean += angle_btw_vec[val-idx][ang]
				mean /= len(angle_btw_vec[val-idx])
				ang_mean.append(mean)
				
			ang_std = []
			ang_dev = [[]]*5 # CHANGE NO OF EIGENMODES, CURRENT = 10
			for val in range(idx,idx+5): # CHANGE NO OF EIGENMODES, CURRENT = 10
				std = 0.0
				for ang in range(len(angle_btw_vec[val-idx])):
					std += (angle_btw_vec[val-idx][ang] - ang_mean[val-idx])*(angle_btw_vec[val-idx][ang] - ang_mean[val-idx])
					if (angle_btw_vec[val-idx][ang] - ang_mean[val-idx])*(angle_btw_vec[val-idx][ang] - ang_mean[val-idx])/num_residues > 0.00000001:
						ang_dev[val-idx].append(ang)
				std /= len(angle_btw_vec[val-idx])
				ang_std.append(std)
			#print ang_dev
			ang_dev = [[]]*5 # CHANGE NO OF EIGENMODES, CURRENT = 10
			for val in range(idx,idx+5): # CHANGE NO OF EIGENMODES, CURRENT = 10
				ang_dev[val-idx] = []
				for ang in range(len(angle_btw_vec[val-idx])):
					if abs(angle_btw_vec[val-idx][ang]) < 0.7:
						#print val, idx, val-idx, ang
						ang_dev[val-idx].append(ang)
				if ang_dev[val-idx] == []:	
					ang_dev[val-idx] = []

			start_gen = time.time()
			#***********************************************************************************************************
			RMSD = 0.0
			rmsd_file0 = open("rmsd0_%d%d.txt"%(it,pdb_num),"w")
			rmsd_file1 = open("rmsd1_%d%d.txt"%(it,pdb_num),"w")
			rmsd_file2 = open("rmsd2_%d%d.txt"%(it,pdb_num),"w")
			print initial_file, target_file
			for disp in frange(min_move,max_move,move_step):
				for val in range(idx,idx+5): # CHANGE NO OF EIGENMODES, CURRENT = 10
				
					new_coord_list_m = scipy.zeros((3*num_residues))
					new_coord_list_p = scipy.zeros((3*num_residues))
					
					magn_v = 0.000000
					norm_v = scipy.zeros((3*num_residues))
					for i in range(num_residues):
						magn_v += v[3*i][val].real*v[3*i][val].real + v[3*i+1][val].real*v[3*i+1][val].real + v[3*i+2][val].real*v[3*i+2][val].real
					magn_v = math.sqrt(magn_v)
					print magn_v
					for i in range(num_residues):
				#	print v[3*i][val].real, v[3*i+1][val].real, v[3*i+2][val].real
						norm_v[3*i] = v[3*i][val].real/magn_v
						norm_v[3*i+1] = v[3*i+1][val].real/magn_v
						norm_v[3*i+2] = v[3*i+2][val].real/magn_v

				
					new = "%s/tmp_ENM_%02d_%s_%d%d%02dp"% (curdir,int(disp//1),str(disp%1)[2:],it,pdb_num,val-idx)
					os.popen("rm -rf %s/"%new).read()
				
					for root, dirs, files in os.walk(new, topdown=False):
						print files, dirs, root
						for name in files:
							os.remove(os.path.join(root, name))
						for name in dirs:
							os.rmdir(os.path.join(root, name))
						os.rmdir(root)
				
					if os.path.isdir(new) == 0:
						os.mkdir(new)

					new = "%s/tmp_ENM_%02d_%s_%d%d%02dm"% (curdir,int(disp//1),str(disp%1)[2:],it,pdb_num,val-idx)
					os.popen("rm -rf %s/"%new).read()
				
					for root, dirs, files in os.walk(new, topdown=False):
						for name in files:
							os.remove(os.path.join(root, name))
						for name in dirs:
							os.rmdir(os.path.join(root, name))
						os.rmdir(root)
				
					if os.path.isdir(new) == 0:
						os.mkdir(new)
				
					atom_shift = []
					for i in range(num_residues): # first all structures in one folder, then next folder
					#	print coord_list[3*i], norm_v[3*i]
						new_coord_list_p[3*i] = coord_list[3*i] + norm_v[3*i]*math.sqrt(num_residues)*disp
						new_coord_list_p[3*i+1] = coord_list[3*i+1] + norm_v[3*i+1]*math.sqrt(num_residues)*disp
						new_coord_list_p[3*i+2] = coord_list[3*i+2] + norm_v[3*i+2]*math.sqrt(num_residues)*disp
						atom_shift_x = norm_v[3*i]*math.sqrt(num_residues)*disp
						atom_shift_y = norm_v[3*i+1]*math.sqrt(num_residues)*disp
						atom_shift_z = norm_v[3*i+2]*math.sqrt(num_residues)*disp
						atom_shift.append(math.sqrt(atom_shift_x*atom_shift_x + atom_shift_y*atom_shift_y + atom_shift_z*atom_shift_z))

					tmp_ang_dev = []
					for i in range(len(ang_dev[val-idx])):
						tmp_ang_dev.append(ang_dev[val-idx][i])
					alist = list(set(tmp_ang_dev))
					tmp_ang_dev = alist
					tmp_ang_dev.sort()
					print "start_%02d_%s_%d%d%02dp.mol2"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx)
					violations = []
					#generatePDB("test_%02d_%s_%d%d%02dp.pdb"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx), new_coord_list_p, violations, pdb_file)
					violations, new_coord_list_p = correctBondDeviations(num_residues, coord_list, norm_v, new_coord_list_p, bonded_A, bonded_B, atom_shift, tmp_ang_dev, violations, "-0.2")
					print violations
					#calculateNMA(new_coord_list_p, num_residues,force_const)
					generatePDB("start_%02d_%s_%d%d%02dp.pdb"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx), new_coord_list_p, [], pdb_file)
					hompdb ="start_%02d_%s_%d%d%02dp.pdb"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx)
					viol = ""
					for vi in violations:
						viol += "%s "%vi
					print viol
					#RMSD = computeRMSD3.align2EachStructure(hompdb,initial_file,viol)
					#rmsd_file0.write("%s\t%.3f\n"%(hompdb,float(RMSD[1].split()[1])))
					#RMSD = computeRMSD3.align2EachStructure(hompdb,start_pdbs[pdb_num],viol)
					#rmsd_file1.write("%s\t%.3f\n"%(hompdb,float(RMSD[1].split()[1])))
					RMSD = computeRMSD3.align2EachStructure(hompdb,target_file,viol)
					rmsd_file2.write("%s\t%.3f\n"%(hompdb,float(RMSD[1].split()[1])))
					#print RMSD
					#print new_coord_list_p
					#sys.exit(2)
					generateMol2("start_%02d_%s_%d%d%02dp.mol2"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx), new_coord_list_p, violations)
					#generatePDB("start_%02d_%s_%d%d%02dp.pdb"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx), new_coord_list_p, violations, pdb_file)
					dir_in = "tmp_ENM_%02d_%s_%d%d%02dp"%(int(disp//1),str(disp%1)[2:],it,pdb_num,val-idx)
					os.popen("mv start_%02d_%s_%d%d%02dp.mol2 %s/%s"% (int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx,curdir,dir_in))
					os.popen("mv start_%02d_%s_%d%d%02dp.pdb %s/%s"% (int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx,curdir,dir_in))
					os.chdir("%s/%s"%(curdir,dir_in))
#					sys.exit(2)
#					comm = "babel output_align.pdb output_align_split.pdb -m &>/dev/null"
#					os.popen(comm)

					#hbonds_A = []
					#hbonds_B = []
					#num_hbonds = []
					#hphob_A = []
					#hphob_B = []
					#num_hphob = []
					#saltbridge_A = []
					#saltbridge_B = []
					#force_const = scipy.zeros((num_residues, num_residues))
					#g_width = scipy.zeros((num_residues, num_residues))
					#g_V0 = scipy.zeros((num_residues, num_residues))
					#shake_list = scipy.zeros((num_residues, num_residues))
					#dr_list = scipy.zeros((num_residues, num_residues))
					#eig_coord_list = scipy.zeros((3*num_residues))
					#mass_list = scipy.zeros((num_residues))

					#generateForceConstantMatrix("start_%d%d.pdb"%(it,pdb_num), force_const, g_width, g_V0, shake_on, shake_list, dr_list, eig_coord_list, mass_list, bonded_A, bonded_B, disulfide_A, disulfide_B, saltbridge_A, saltbridge_B, hbonds_A, hbonds_B, num_hbonds, hphob_A, hphob_B, num_hphob, secondary_structure)

					# conformational search: run MD
					#runMD(dir_in, force_const, g_width, g_V0, dr_list, coord_list, eig_coord_list, mass_list, num_timesteps, dt, bond_list_A, bond_list_B, shake_on, num_residues, Temp0, gamma, freq_coord, freq_reassign_vel)
				#jobs.append(job_server.submit(runMD, (dir_in, force_const, g_width, g_V0, dr_list, coord_list, eig_coord_list, mass_list, num_timesteps, dt, bond_list_A, bond_list_B, shake_on, num_residues, Temp0, gamma, freq_coord, freq_reassign_vel,), (initializeVelocities, calculateForces, CalculateTemperature, RemoveCOMVelocities, shakeBonds, writeCoordinates), ("scipy", "math","numpy",)))
					os.chdir(curdir)

				###############################################################################################
					atom_shift = []
					for i in range(num_residues):
					#	print coord_list[3*i], norm_v[3*i]
						new_coord_list_m[3*i] = coord_list[3*i] - norm_v[3*i]*math.sqrt(num_residues)*disp
						new_coord_list_m[3*i+1] = coord_list[3*i+1] - norm_v[3*i+1]*math.sqrt(num_residues)*disp
						new_coord_list_m[3*i+2] = coord_list[3*i+2] - norm_v[3*i+2]*math.sqrt(num_residues)*disp
						atom_shift_x = norm_v[3*i]*math.sqrt(num_residues)*disp
						atom_shift_y = norm_v[3*i+1]*math.sqrt(num_residues)*disp
						atom_shift_z = norm_v[3*i+2]*math.sqrt(num_residues)*disp
						atom_shift.append(math.sqrt(atom_shift_x*atom_shift_x + atom_shift_y*atom_shift_y + atom_shift_z*atom_shift_z))

					violations = []
					print "start_%02d_%s_%d%d%02dm.mol2"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx)
					#generatePDB("test_%02d_%s_%d%d%02dm.pdb"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx), new_coord_list_m, violations, pdb_file)
					violations, new_coord_list_m = correctBondDeviations(num_residues, coord_list, norm_v, new_coord_list_m, bonded_A, bonded_B, atom_shift, tmp_ang_dev, violations, "0.2")
					#calculateNMA(new_coord_list_p, num_residues,force_const)
					generatePDB("start_%02d_%s_%d%d%02dm.pdb"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx), new_coord_list_m, [], pdb_file)
					hompdb ="start_%02d_%s_%d%d%02dm.pdb"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx)
					viol = ""
					for vi in violations:
						viol += "%s "%vi
					#RMSD = computeRMSD3.align2EachStructure(hompdb,initial_file,viol)
					#rmsd_file0.write("%s\t%.3f\n"%(hompdb,float(RMSD[1].split()[1])))
					#RMSD = computeRMSD3.align2EachStructure(hompdb,start_pdbs[pdb_num],viol)
					#rmsd_file1.write("%s\t%.3f\n"%(hompdb,float(RMSD[1].split()[1])))
					RMSD = computeRMSD3.align2EachStructure(hompdb,target_file,viol)
					rmsd_file2.write("%s\t%.3f\n"%(hompdb,float(RMSD[1].split()[1])))
					print violations
					generateMol2("start_%02d_%s_%d%d%02dm.mol2"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx), new_coord_list_m, violations)
					#generatePDB("start_%02d_%s_%d%d%02dm.pdb"%(int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx), new_coord_list_m, violations, pdb_file)
					dir_in = "tmp_ENM_%02d_%s_%d%d%02dm"%(int(disp//1),str(disp%1)[2:],it,pdb_num,val-idx)
					os.popen("mv start_%02d_%s_%d%d%02dm.mol2 %s/%s"% (int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx ,curdir,dir_in))
					os.popen("mv start_%02d_%s_%d%d%02dm.pdb %s/%s"% (int(disp//1),str(disp%1)[2:],it,pdb_num, val-idx,curdir,dir_in))
#					sys.exit(2)
					os.chdir("%s/%s"%(curdir,dir_in))
#					comm = "babel output_align.pdb output_align_split.pdb -m &>/dev/null"
#					os.popen(comm)
					#hbonds_A = []
					#hbonds_B = []
					#num_hbonds = []
					#hphob_A = []
					#hphob_B = []
					#num_hphob = []
					#saltbridge_A = []
					#saltbridge_B = []
					#force_const = scipy.zeros((num_residues, num_residues))
					#g_width = scipy.zeros((num_residues, num_residues))
					#g_V0 = scipy.zeros((num_residues, num_residues))
					#shake_list = scipy.zeros((num_residues, num_residues))
					#dr_list = scipy.zeros((num_residues, num_residues))
					#eig_coord_list = scipy.zeros((3*num_residues))
					#mass_list = scipy.zeros((num_residues))
					#generateForceConstantMatrix("start_%d%d.pdb"%(it,pdb_num), force_const, g_width, g_V0, shake_on, shake_list, dr_list, eig_coord_list, mass_list, bonded_A, bonded_B, disulfide_A, disulfide_B, saltbridge_A, saltbridge_B, hbonds_A, hbonds_B, num_hbonds, hphob_A, hphob_B, num_hphob, secondary_structure)
					# conformational search: run MD
					#runMD(dir_in, force_const, g_width, g_V0, dr_list, coord_list, eig_coord_list, mass_list, num_timesteps, dt, bond_list_A, bond_list_B, shake_on, num_residues, Temp0, gamma, freq_coord, freq_reassign_vel)
					#jobs.append(job_server.submit(runMD, (dir_in, force_const, g_width, g_V0, dr_list, coord_list, eig_coord_list, mass_list, num_timesteps, dt, bond_list_A, bond_list_B, shake_on, num_residues, Temp0, gamma, freq_coord, freq_reassign_vel,), (initializeVelocities, calculateForces, CalculateTemperature, RemoveCOMVelocities, shakeBonds, writeCoordinates), ("scipy", "math","numpy",)))					
					os.chdir(curdir)
				#sys.exit(2)
			rmsd_file0.close()
			rmsd_file1.close()
			rmsd_file2.close()
		
		end_gen = time.time()
		print "Generating str time: %s"%str(datetime.timedelta(seconds=end_gen-start_gen))
"""
		#sys.exit(2)
		"""							
			for job in jobs:
				result = job()
			job_server.print_stats()
			job_server.destroy()
			fout.close()
		"""
		#sys.exit(2)	
		"""
		dir_list = []
		tmp_list = os.listdir(curdir)
		for dir in tmp_list:
			#if dir.find("tmp_ENM_%d_%s_%d%d"%(int(disp//1),str(disp%1)[2:],it,pdb_num)) >= 0:
			if dir.find("tmp_ENM_") >= 0 and dir.find("_%d%d0"%(it,pdb_num)) >= 0:
				dir_list.append(dir)
		dir_list.sort()
		print dir_list
				
		
		
		#pd2_ca2.splitPDBAndaddBackbone(dir_list)
		#scwrl_run.addSidechains(dir_list)	
		
		print curdir
		start_pulchra = time.time()
		pulchra.splitPDBAndaddBackbone(dir_list)
		end_pulchra = time.time()
		print "Adding BB time: %s"%str(datetime.timedelta(seconds=end_pulchra-start_pulchra))
		new="scwrl_%d%d"%(it,pdb_num)
		os.popen("rm -rf %s/"%new).read()

		for root, dirs, files in os.walk(new, topdown=False):
			print files, dirs, root
			for name in files:
				os.remove(os.path.join(root, name))
			for name in dirs:
				os.rmdir(os.path.join(root, name))
				os.rmdir(root)

		if os.path.isdir(new) == 0:
			os.mkdir(new)

		#os.popen("mkdir scwrl_%d%d"%(it,pdb_num))
		comm = "cp tmp_ENM_*_%d%d*/start_*.rebuilt.pdb scwrl_%d%d/"%(it,pdb_num,it,pdb_num)
		os.popen(comm)
		os.popen("cp ../scwrl_run.py scwrl_%d%d"%(it,pdb_num))
#		os.popen("cp ../computeQMEAN.py scwrl_%d%d"%(it,pdb_num))
		os.chdir("%s/scwrl_%d%d"%(curdir,it,pdb_num))
		start_scwrl = time.time()
		scwrl_run.addSidechain("scwrl_%d%d"%(it,pdb_num))
		end_scwrl = time.time()
		print "Adding SC time: %s"%str(datetime.timedelta(seconds=end_scwrl-start_scwrl))

		start_score = time.time()
		energy = []
		energy = computeJerniganPotential.computePotential("scwrl_%d%d"%(it,pdb_num))
#		energy = computeQMEAN.computeQMEANandRMSD("scwrl_%d%d"%(it,pdb_num))
#		scwrl_run.addSidechains(dir_in)
#		computeQMEAN(dir_list)
#		computeGOAP.computeGOAPandRMSD(dir_list)
		os.chdir(curdir)
		end_score = time.time()
		print "Scoring time: %s"%str(datetime.timedelta(seconds=end_score-start_score))
		"""
		"""
		for dir_in in dir_list:
			print dir_in
			os.chdir("%s/%s"%(curdir,dir_in))
			print os.getcwd()
			# align trajectory
			#alignTrajectory(dir_in)
			pd2_ca2.splitPDBAndaddBackbone(dir_in)
			scwrl_run.addSidechains(dir_in)
			#score
		"""
		"""
		frmsd=open("rmsd0_%d%d.txt"%(it,pdb_num))
		rmsd_lines = frmsd.readlines()
		initial_rmsd = []
		for r in rmsd_lines:
			initial_rmsd.append((r.split()[0].strip(),r.split()[1].strip()))
		frmsd.close()

		flag = 0
		tmp = []
		new_rmsd = []
		for values in range(len(initial_rmsd)):
			if values%10 == 0: # CHANGE NO OF EIGENMODES, CURRENT = 20
				flag = 1
			if flag == 1:
				if tmp == []:
					tmp.append(initial_rmsd[values])
					flag = 0
					continue
				else:
					new_rmsd.append((tmp))
					tmp = []
					flag = 0
			tmp.append(initial_rmsd[values])
		new_rmsd.append((tmp))
		tmp = []
		rmsd_by_mode = []
		for modes in range(len(new_rmsd[0])):
			tmp = []
			for values in range(len(new_rmsd)):
				tmp.append(new_rmsd[values][modes][1])
			rmsd_by_mode.append((tmp))
		print rmsd_by_mode
		#*****************************************************************************************************
		frmsd=open("rmsd2_%d%d.txt"%(it,pdb_num))
		rmsd_lines = frmsd.readlines()
		initial_rmsd = []
		for r in rmsd_lines:
			initial_rmsd.append((r.split()[0].strip(),r.split()[1].strip()))
		frmsd.close()
		
		flag = 0
		tmp = []
		new_rmsd = []
		for values in range(len(initial_rmsd)):
			if values%10 == 0: # CHANGE NO OF EIGENMODES, CURRENT = 10
				flag = 1
			if flag == 1:
				if tmp == []:
					tmp.append(initial_rmsd[values])
					flag = 0
					continue
				else:
					new_rmsd.append((tmp))
					tmp = []
					flag = 0
			tmp.append(initial_rmsd[values])

		new_rmsd.append((tmp))
		tmp = []
		rmsd_by_mode = []
		for modes in range(len(new_rmsd[0])):
			tmp = []
			for values in range(len(new_rmsd)):
				tmp.append(new_rmsd[values][modes][1])
			rmsd_by_mode.append((tmp))
		print rmsd_by_mode
		#*****************************************************************************************************
		frmsd=open("rmsd1_%d%d.txt"%(it,pdb_num))
		rmsd_lines = frmsd.readlines()
		initial_rmsd = []
		for r in rmsd_lines:
			initial_rmsd.append((r.split()[0].strip(),r.split()[1].strip()))
		frmsd.close()

		flag = 0
		tmp = []
		new_rmsd = []
		for values in range(len(initial_rmsd)):
			if values%10 == 0: # CHANGE NO OF EIGENMODES, CURRENT = 10
				flag = 1
			if flag == 1:
				if tmp == []:
					tmp.append(initial_rmsd[values])
					flag = 0
					continue
				else:
					new_rmsd.append((tmp))
					tmp = []
					flag = 0
			tmp.append(initial_rmsd[values])

		new_rmsd.append((tmp))
		tmp = []
		rmsd_by_mode = []
		for modes in range(len(new_rmsd[0])):
			tmp = []
			for values in range(len(new_rmsd)):
				tmp.append(new_rmsd[values][modes][1])
			rmsd_by_mode.append((tmp))
		print rmsd_by_mode
		"""
		#*****************************************************************************************************
		"""
		fenergy=open("scwrl_%d%d/qmean_score.txt"%(it,pdb_num))
		energy_lines = fenergy.readlines()
		energy = []
		for ene in energy_lines:
			energy.append((ene.split()[0].strip(),ene.split()[1].strip()))
		fenergy.close()
		"""
		"""
		flag = 0
		tmp = []
		new_energy = []
		for values in range(len(energy)):
			if values%10 == 0: # CHANGE NO OF EIGENMODES, CURRENT = 10
				flag = 1
			if flag == 1:
				if tmp == []:
					tmp.append(energy[values])
					flag = 0
					continue
				else:
					new_energy.append((tmp))
					tmp = []
					flag = 0
			tmp.append(energy[values])

		new_energy.append((tmp))
		#print new_energy
		tmp = []
		energy_by_mode = []
		for modes in range(len(new_energy[0])):
			tmp = []
			for values in range(len(new_energy)):
				tmp.append(new_energy[values][modes][1])
			energy_by_mode.append((tmp))
		print energy_by_mode
		
		# prepare pdb files of 5 structures
		# ENM_01.pdb, ENM_02.pdb, ENM_03.pdb, ENM_04.pdb, ENM_05.pdb, ENM_06.pdb
		for item in range(num_pdbs):
			for direct in range(2):
				first_min = 0
				for values in range(len(energy_by_mode[item*2+direct])-1):
					#if -1.0*float(energy_by_mode[item*2+direct][values])*(1+math.(pow(float(rmsd_by_mode[item*2+(1-direct)][values]),0.67))*math.log(eigen[item]) <= -1.0*float(energy_by_mode[item*2+direct][values+1])*(1+math.pow(float(rmsd_by_mode[item*2+(1-direct)][values+1]),0.67))*math.log(eigen[item]): # 3rd function
					#if float(energy_by_mode[item*2+direct][values])*(1+math.sqrt(float(rmsd_by_mode[item*2+(1-direct)][values]))) <= float(energy_by_mode[item*2+direct][values+1])*(1+math.sqrt(float(rmsd_by_mode[item*2+(1-direct)][values+1]))): # 2nd function
					if -1.0*float(energy_by_mode[item*2+direct][values])*float(rmsd_by_mode[item*2+(1-direct)][values])*math.log(eigen[item]) <= -1.0*float(energy_by_mode[item*2+direct][values+1])*float(rmsd_by_mode[item*2+(1-direct)][values+1])*math.log(eigen[item]): # 1st function
						first_min = values
						break
					else:
						first_min = values+1
				print	new_energy[first_min][item*2+direct][0]
#				file_list.write("%s/scwrl/%s\n"%(curdir,new_energy[first_min][item*2+direct][0]))
#				comm = "cp %s/scwrl_%d%d/%s ENM_%d%d_%d%d.pdb"%(curdir,it,pdb_num,new_energy[first_min][item*2+direct][0],it-1,pdb_num+1,it,item*2+direct+1)
				if it == 0:
#					comm = "cp /home/nrana/Downloads/energy/energy/%s ENM_%d%d.pdb"%(new_energy[first_min][item*2+direct][0],it,item*2+direct+1)
					comm = "cp /home/nrana/Downloads/energy/energy/%s ENM_%d%d.pdb"%(new_energy[first_min][item*2+direct][0],it,item*2+direct)
					os.popen(comm)
				else:
					comm = "cp /home/nrana/Downloads/energy/energy/%s ENM_%d%d_%d%d.pdb"%(new_energy[first_min][item*2+direct][0],it-1,pdb_num+1,it,item*2+direct)
					os.popen(comm)
				#f1 = open("%s/ENM_%d%d.pdb"%(curdir, it, item),"w")
				#f2 = open("%s/tmp_ENM_%s/output_align_split%03d.pdb"%(curdir,fin_energy[item][0],int(fin_energy[item][1])))
				#for f2_line in f2:
				#	f1.write(f2_line)
				#f1.close()
				#f2.close()
		"""
		#return new_hessian	
		return norm_v
###################################################################################################
###################################################################################################		
def	sortByColumn(A,*args):
	import operator
	return sorted(A,key=operator.itemgetter(*args))

###################################################################################################
def	frange(x, y, jump):
	while x <= y:
		yield x
		x += jump
###################################################################################################
def	main(argv):
	global  curdir, command_file, DSSP_exe
	global	num_residues, num_pdbs, file_num
	global  lig_name, lig_num, num_lig_atoms, num_prot_atoms
	global	hessian, new_hessian, norm_v
	global 	dt, num_timesteps, num_iterations, pdb_file, filemdp, filetdist, shake_on
	global 	k_bond, k_bond_1_3, k_bond_1_4, k_helix_1_4, k_hbonds, k_disulfide, k_salt, k_hphob, k_density, r_cutoff, lambda_hbond, beta_hbond, lambda_hphob, beta_hphob, gam_const, k_aromatic
	global	g_k_bond, g_k_helix_1_3, g_k_sheet_1_4, g_k_helix_1_4, g_k_hbonds, g_k_disulfide, g_k_hphob, g_k_density, g_r_cutoff, g_lambda_hbond, g_beta_hbond, g_lambda_hphob, g_beta_hphob
	global	g_V0_bond, g_V0_helix_1_3, g_V0_sheet_1_4, g_V0_helix_1_4, g_V0_hbonds, g_V0_disulfide, g_V0_hphob, g_V0_density
	global	tau, Temp0, gamma, freq_coord, freq_reassign_vel
	global coord_list
	global lig_names, lig_nums
	lig_names = []
	lig_nums = []

	num_iterations = 1
	shake_on = 0
	curdir = os.getcwd()

	## This code will work for only a single pdb structure. If you have multiple pdb structures, please put the appropriate sections in loop for each protein structure ##
	usage = "use: ENM_concoord_Dynamics.py\n \t-r (run type; =0: single point energy; =1: Langevin-Dynamics)\n \t-c (command file)\n \t-i (input PDB file; use full path)\n \t-n (number of iterations for conformational search; default: 100)\n \t-s (number of time steps in each iteration; default: 1000)\n \t-d (time step in fs; default: 10)\n \t-e (input mdp file for energy minimization; use full path)\n \t-t (input file for tdist, e.g. input.cpf; use full path)\n \t-k (shake for CA-CA bonds and disulfide bridges)"
	try:                                
		opts, args = getopt.getopt(argv, "hr:c:p:m:g:r:i:n:s:d:e:t:k:l:u:x:") 
	except getopt.GetoptError:          
		sys.exit(2)                     
	for opt, arg in opts:
		if opt == '-h':
			print usage                     
			sys.exit()                  
		elif opt == '-c':
			command_file = arg
		elif opt == '-p':
			num_pdbs = int(arg)
		elif opt == '-g':
			file_num = int(arg)

	readCommandFile(command_file)
	print "read command file"
	
	for opt, arg in opts:
		if opt == '-r':
			run_type = int(arg)
		if opt == '-i':
			pdb_file = arg
		elif opt == '-n':
			num_iterations = int(arg)
		elif opt == '-s':
			num_timesteps = int(arg)
		elif opt == '-d':
			dt = float(arg)
		elif opt == '-e':
			filemdp = arg
		elif opt == '-t':
			filetdist = arg
		elif opt == '-k':
			shake_on = int(arg)
		elif opt == '-l':
			#lig_name = arg
			lig_names = arg.split("_")
		elif opt == '-u':
			#lig_num = int(arg)
			lig_nums =  arg.split("_")
		elif opt == '-x':
			r_cutoff = float(arg)
			
	print pdb_file, lig_names, lig_nums	
	# create directory and delete old files (if existing)
	curdir = os.getcwd()
	"""
	try:
		for root, dirs, files in os.walk(md_outfile, topdown=False):
			for name in files:
				os.remove(os.path.join(root, name))
			for name in dirs:
				os.rmdir(os.path.join(root, name))
		if os.path.isdir(md_outfile) == 0:
			os.mkdir(md_outfile)
	except:
		pass
	"""
	#prepareProtein4ENM(pdb_file, "prepared4ENM_%d%d.pdb"%(num_iterations,file_num))
	#runENMAnalysis("prepared4ENM.pdb")
	if run_type == 1:
		norm_v = runLangevinDynamics(pdb_file)
		#hessian = runLangevinDynamics("prepared4ENM_%d%d.pdb"%(num_iterations,file_num))

	
	#print len(new_hessian)
	#print coord_list	
	# Calculating eigenmodes using ProDy: FAST, n_modes=30 modes calculated
	#anm = ANM('Using external Hessian')
	#anm.setHessian(new_hessian)
	#anm.calcModes(n_modes=30)
	#print len(list(anm[:30]))
	#os.chdir(curdir)

	# Using ProDy to calculate the deformation between open-closed conformation 1LFH-1LFG
	#reference = parsePDB('%s'%pdb_file)
	ref = "%s"%initial_file[-8:-4]
	mob = "%s"%target_file[-8:-4]
#	reference = parsePDB(ref)
#	mobile = parsePDB(mob)
	print ref, mob
	reference = parsePDB('1KX9.pdb')
	mobile = parsePDB('1N8V.pdb')
	calphas = reference.select('protein and name CA')
	calphas1 = mobile.select('protein and name CA')
	calphas1, t = superpose(calphas1, calphas)
	defvec = calcDeformVector(calphas, calphas1)
	defvecnormed = defvec.getNormed()
	arr = defvecnormed.getArray()
	# Calculating overlap from ProDy eigenvectors
	for i in range(30):
		overlap = 0.0
		evec = scipy.zeros(len(calphas)*3)
		for j in range(len(calphas)*3):
			evec[j] = norm_v[j][i]
		overlap = np.dot(evec,defvecnormed.getArray())

	#	evec = scipy.zeros(12*3)
	#	for j in range(12*3):
	#		evec[j] = norm_v[j+88*3][i]
	#		overlap += evec[j]*arr[j+88*3]

		print overlap
	#print (array(list(anm[:30])) * defvecnormed).astype(float64).round(2)
	"""
	new_coord_list_p = scipy.zeros((3*num_residues))
	num_residues = 691
	for i in range(num_residues):
	#       print coord_list[3*i], norm_v[3*i]
		new_coord_list_p[3*i] = coord_list[3*i] + arr[3*i]*math.sqrt(num_residues)*6.000
		new_coord_list_p[3*i+1] = coord_list[3*i+1] + arr[3*i+1]*math.sqrt(num_residues)*6.000
		new_coord_list_p[3*i+2] = coord_list[3*i+2] + arr[3*i+2]*math.sqrt(num_residues)*6.000

		generateMol2("start.mol2", new_coord_list_p)
		generatePDB(".pdb", new_coord_list_p)
	
	
	"""
	os.chdir(curdir)

main(sys.argv[1:])

