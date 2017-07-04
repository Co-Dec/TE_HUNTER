#!/usr/bin/python
# -*- coding: Utf-8 -*-
import re
import argparse

parser = argparse.ArgumentParser(description="TE_masker is a tool made to mask TE in a reference genome. \
It masks the genomic positions and not the sequences. It doesn't mask other regions like repeated regions. \
The input must be FASTA format with Flybase headers. \ ")

parser.add_argument("-t",help="The TE reference, FASTA, Flybase headers",required=True)
parser.add_argument("-r",help="The reference genome, FASTA format, Flybase headers",required=True)
parser.add_argument("-o",help="Output file, masked genome, default = genome.masked.fa",default="genome.masked.fa")

args = parser.parse_args()


#Fonctions : ####################################################

def recup_nom(seq): #recupere le chromosome du te à partir du champ localisation
	c=0
	r=[]
	for i in seq :
		if i == "=" :
			c=1
		if i == ":" :
			c=0
		if c == 1 and i != "=" :
			r.append(i)
	return "".join(r)
	
def recupere_positions(lire) :
	"""
	Lit l'ensemble des te (le fichier banque)
	Retourne une liste de listes contenant [Début,Longueur,Chromosome]
	"""
	ref = open(lire,"r")
	transpo=[]
	liste_des_masqueurs=[]
	for lines in ref :
		if lines[0] == ">" :
			transpo = lines.split(" ")
			depend=[]
			pos = transpo[2]
			
			depend.append(re.findall('\d+',pos.replace(str(recup_nom(pos))+":",""))[0])
			depend.append(re.findall('\d+',pos.replace(str(recup_nom(pos))+":",""))[1])
			
			
			lenghteur=transpo[6].replace("length=","")
			lenghteur=lenghteur.replace(";","")
			
			liste_des_masqueurs.append([int(depend[0]),int(lenghteur),recup_nom(transpo[2])])
	ref.close()
	
	return liste_des_masqueurs


def trie_insertions(te) :
	"""
	Fonction qui trie la liste par position de début de chaque TE
	"""
	liste_pos_a_masquer = recupere_positions(te)
	liste_pos_a_masquer.sort(key=lambda c: int(c[0]))
	return liste_pos_a_masquer

def cree_dicchr (te):
	"""
	cree dictionaire pas en dur
	"""
	
	listri = trie_insertions(te)
	chrom = []
	dichrom={}
	for ch in listri :
		if ch[2] not in chrom :
			chrom.append(ch[2])
	for ch in chrom :
		dichrom[ch]=[]
	
	for ch in listri :
		
		dichrom[ch[2]].append(ch[:2])
	
	return dichrom

def masque_genome_efficace(te):
	"""
	Cree un dictionnaire Chromosome : [Debut, longueur]
	a partir de la liste précédente
	"""
	L2=[]
	L3=[]
	R3=[]
	R2=[]
	X=[]
	c4=[]
	listriee = trie_insertions(te)
	for ch in listriee :
		
		if ch[2] == "2L" :
			L2.append(ch[:2])
		if ch[2] == "2R" :
			R2.append(ch[:2])
		if ch[2] == "3L" :
			L3.append(ch[:2])
		if ch[2] == "3R" :
			R3.append(ch[:2])
		if ch[2] == "4" :
			c4.append(ch[:2])
		if ch[2] == "X" :
			X.append(ch[:2])	
	
	
	return {"2L":L2,"3L":L3,"2R":R2,"3R":R3,"4":c4,"X":X}
	

def masque_genome_entier(liste_pos_len_chr,genom,sortie):
	"""
	Fonction qui prend : Le dictionnaire crée avant - Le génome à masquer
	Retourne : Le génome masqué
	Efficace : Lit le génome et l'écrit en parallèle en une fois
	Masque les positions quand situé sur un ET
	Quand plusieurs ET commencent à la même base : Ne prend en considération que le plus court : Problème.
	"""
	genome=open(genom,"r")
	masquedgenome=open(sortie,"w")
	decompte_ET=int()
	cptrET=int()
	cptrPOS=int()
	c=0
	WARNING = False
	for line in genome :			#Charge chaque ligne du génome indépendament
		
		if line[0] == ">" :			#Si ">" ecrit l'entête
			
			c=0
			decompte_ET=0
			chrom = line.split(" ")
			chromosome = chrom[0].replace(">","")
			if chromosome in liste_pos_len_chr :
				print "###\tMasking chromosome " + str(chromosome)+".\t###"
				masquedgenome.write(str(chrom[0])+"\n")
			cptrET= 0
			
			cptrPOS = 0
			
			NON = False
		else :	
								#Sinon on est sur de l'adn
			if chromosome not in liste_pos_len_chr and NON == False:
						WARNING = True
						NON = True
			if NON == False :
				for base in line :
					interrupteur = False
					Debut_ET = False
					Milieu_ET = False
					
					if base == "\n" :
						masquedgenome.write("\n")
						
						break
						
					else :
						
						cptrPOS += 1
						
						
							
						
						if cptrET < len(liste_pos_len_chr[chromosome]) :
							
							if cptrPOS > int(liste_pos_len_chr[chromosome][cptrET][0]):
								while cptrPOS > int(liste_pos_len_chr[chromosome][cptrET][0]) :
									cptrET+=1
							if cptrPOS == int(liste_pos_len_chr[chromosome][cptrET][0]):
								
								if int(decompte_ET) < int(liste_pos_len_chr[chromosome][cptrET][1]) :
									decompte_ET = int(liste_pos_len_chr[chromosome][cptrET][1])
								
								cptrET += 1
								
								masquedgenome.write("N")
								decompte_ET = decompte_ET -1
								Debut_ET = True
								c+=1
								
						if decompte_ET != 0 and Debut_ET == False:
							masquedgenome.write("N")
							decompte_ET =decompte_ET -1		
							Milieu_ET = True
										
						if Debut_ET == False and Milieu_ET == False:
							
							masquedgenome.write(base)
		#print c					
	genome.close()
	masquedgenome.close()
	return WARNING
##################################################################
print "###\t######################\t###"
print "###\tTE_Masker started\t###"
print "###\tC.DECHAUD\t\t###"

eff= cree_dicchr(args.t)

WARNING = masque_genome_entier(eff,args.r,args.o)
print "###\tMasking is over.\t###"
print "###\t######################\t###"
print "\n"
if WARNING :
	print "----Warning : Some chromosomes are not in your annotation, deleted from masked genome.----"
