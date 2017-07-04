#!/usr/bin/python
# -*- coding: Utf-8 -*-

"""
TE_HUNTER.py
This pipeline finds Transposable elements in a genome. Based on a NGS single-end reads.
It takes 7 arguments :
	1) Your working directory (ex : /home/me/workdir/)
	2) A prefix for all the files it will generate, so you can remove them easier.
	3) A reference sequenced genome.
	4) The masked reference genome.
	5) Your reads file (Fastq input)
	6) Your Transposable elements bank.
	7) A description file of your transposable elements
"""
#	Packages required :
try :
	import sys
	import os
	import subprocess
	import string

	import argparse
	import re
	import threading as th

except ImportError :
	print "This script requires sys, os, subprocess, string, argparse, re, threading and python 2.7 to be installed!"

#	Arguments definition :
parser = argparse.ArgumentParser(description="TE_hunter is a program looking for TE in a genome. \
Return the list of those insertions with their type : 'Shared' or 'Absent'. \
It can write new insertions in another file. DO NOT GIVE THE PATH OF YOUR FILES IN THE ARGUMENTS")


#PreSTEP
parser.add_argument("-ted",help="Input file with TE description")
parser.add_argument("-tebank",help="TE bank from flybase, FASTA, Mandatory",required=True)
parser.add_argument("-workdir",help="Directory you are working in",required=True)
parser.add_argument("-preftoremove",help="File prefix to easy remove",required=True)
#Step 1
parser.add_argument("-maskedref",help="Masked reference genome, Mandatory",required=True)
parser.add_argument("-ref",help="Reference genome, Mandatory",required=True)

#Step 2
parser.add_argument("-reads",help="Single end reads file, FASTQ format, Mandatory",required=True)
parser.add_argument("-indexref",help="Reference genome index",default="index_dmel")


parser.add_argument("-mappedreadsonref",help="Samfile from reads mapped on the reference genome",default="mapped_on_ref_reads.sam")
parser.add_argument("-unmappedreadsonref",help="FASTQ file containing reads not mapped on the reference",default="unmapped_reads_to_ref.fq")
parser.add_argument("-cores",help="The number of core you want to use, !cluster!, default = 16 ",default = 16,type=int)
parser.add_argument("-indexte",help="Output files with the index of the TE bank, prefix",default="index_te")
parser.add_argument("-unmappedreadsonboth",help="FASTQ file containing reads not mapped on the reference AND not mapped to the TE bank",default="unmapped_reads_to_both.fq")

parser.add_argument("-blastedreads",help="Output file from the blast of the unmapped reads",default="unmapped.blastn")

#Step 3
parser.add_argument("-opnf",help="Output file from blast parsed",default="parsed.txt")
parser.add_argument("-opf",help="Output file from blast parsed and filtered",default="parsed_and_filtered.txt")
parser.add_argument("-trimmed",help="Trimmed reads, FASTA format",default="trimmed.fa")

parser.add_argument("-lenmin",help="Minimum length of the trimmed reads, def = 20",default = 20,type=int)

parser.add_argument("-blastedcontigs",help="Output file from the contigs blast",default="contigs.blastn")


#Step 4
parser.add_argument("-sortdcontigs",help="Fichier de sortie en FASTA des contigs triés",default="sorted_contigs.fa")
parser.add_argument("-mappedcontigsonref",help="SAM output from contig mapping on reference",default="contigs.sam")

parser.add_argument("-rd",help="Output file woth contig description",default="contig_descr.txt")

parser.add_argument("-outputfrknownte",help="Output containing known TE",default="output.known.te_hunter.txt")
parser.add_argument("-distn",help="Distance you want to look to check the contigs, def = 30",default=5,type=int)
parser.add_argument("-rest",help="Unused contigs",default="unused_contigs.txt")
parser.add_argument("-onew",help="Output file with new insertions",default="output.unknown.te_hunter.txt")

args = parser.parse_args()

#	Functions :

def recup_fam(champ_fam): #Recup famille a partir champ famille
	name = champ_fam.replace("name=","")
	nom=[]
	for j in name :
		if j != "{" :
			nom.append(j)
		else : 
			break
	name = "".join(nom) #Type element
	return name	

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

def recup_positions(loc) :	#Recupere les positions
	

	
	cp=0
	ch=[]
	k=0
	start=[]
	end=[]
		
	for j in loc :
		if j != ":" and cp == 0 :
			ch.append(j)
		else :
			cp += 1
			
	
	chd = "".join(ch)
	chd=chd.replace("loc=","")
	a = "loc="+str(chd)	
	loc_start = loc.replace(a,"")
	
	for i in loc_start :
		
		if i.isdigit():
			if k==0 :
				
				start.append(i)
			if k==1 :
				
				end.append(i)
		if i ==".": 
			k=1
	 

	starter = "".join(start) #Debut
	ender="".join(end) #Fin
	
		
	
	return [starter,ender]

def genere_ted(te_anot,descr):
	ted = open(te_anot,"r")
	descri = open(descr,"w")
	for line in ted :
		if line[0] == ">":
			element_tr = line.split(" ")
			#print element_tr
			descri.write(element_tr[0].replace(">","")+"\t"+recup_fam(element_tr[3])+"\t"+recup_nom(element_tr[2])+"\t"+str(recup_positions(element_tr[2])[0])+"\t"+str(recup_positions(element_tr[2])[1])+"\n")
			#Recup posi spas du touuut ça que je veux
	ted.close()
	descri.close()
#	Reading blast output

def remetsurun(reads):
	subprocess.call("cat {0} | tr '\n' '|' | sed -e 's/|>/||>/g' | tr '|' '\n' > {1}".format(reads,str(reads)+str("inter")),shell=True)
	read=open(str(reads)+str("inter"),"r")
	sortie=str(reads)+"rearr"
	rearr=open(sortie,"w")
	for i in read :
		if i[0] != ">" and i[0] != "\n" :
			rearr.write(i.replace("\n",""))
		else :
			rearr.write(i.split(" ")[0]+"\n")
	rearr.close()
	read.close()

def lire_length(reads): 
	"""
	Lit un fichier fasta
	Fonction qui retourne un dictionnaire contenant :
	Clé = Nom de la séquence
	Valeur = longueur de la Séquence
	Attention : séquence sur une ligne
	"""
	remetsurun(reads)
	
	reads = str(reads) + "rearr"
	
	
	read=open(reads,"r") #Lit le fichier des reads en format FASTA
	liste=[]
	seq=[]
	noms=[]
	tot = []
	bouquin = {}
	
	for line in read :
		if line[0] != ">":
			seq.append(len(line.replace("\n",""))) #Ajoute à seq la longueur de la sequence
		if line[0] == ">":
			if sum(seq) != 0 :
				tot.append(sum(seq))
			seq=[]
			noms.append(line.split(" ")[0].replace(">","").replace("\n",""))#Puis à noms le nom
	tot.append(sum(seq))
	length = len(noms)
	for j in range(length): #Parcours les deux listes et les ajoute dans le dictionnaire
		bouquin[noms[j]]=tot[j]
	
	read.close()
	return bouquin
	
def parsage(blastn,parsed,reads_non_mappes_fasta) : # Analyse blast output
	
	"""
	Fonction qui lit la sortie de blast
	Elle ecrit de ou a ou un read doit être coupé
	En prenant le plus grand écart possible
	"""
	
	blast=open(blastn,"r")
	parse=open(parsed,"w")
	
	read=""
	pave = lire_length(reads_non_mappes_fasta) #Cree le dictionnaire contenant les reads et leur longueur
	
	
	for hit in blast :
		
		ligne = hit.split("\t")
		if ligne[0] != read and read != "" : #Si jamais c'est la premiere ligne ou que j'ai deja vu ce read, je passe. 
			parse.write(str(read)+"\t"+str(debut)+"\t"+str(fin)+"\t"+str(pave[read])+"\t"+str(ETname)+"\n") #Sinon j'écris NOM + DEBUT + FIN +  LONGUEUR du read PRECEDENT
		
		if ligne[0] != read :	#Si c'est la premiere fois que je vois ce read :
				
			read = ligne[0] 		#Je stocke son nom
			ETname = ligne[1]		#Stocke le nom de l'ET
			debut = int(ligne[6])	#Son début et sa fin
			fin = int(ligne[7])
			
		else :						#Si j'ai déja vu ce read
			if int(ligne[6]) < debut :	#Je modifie ses caractéristiques que si elles sont plus intéressantes sur ce hit
				debut = int(ligne[6])
			if int(ligne[7]) > fin : 
				fin = int(ligne[7])
	
	parse.write(str(read)+"\t"+str(debut)+"\t"+str(fin)+"\t"+str(pave[read])+"\t"+str(ETname)+"\n") #Ecrit la derniere ligne
		
				
	
	blast.close()
	parse.close()
	return

def enleve_mauvais_reads(lire,ecrire) : # Filters parsed output, gets 1-X or X-end
	"""
	Fonction qui purifie le fichier parsé de blast
	Enleve les reads qui sont entiers ou qui ne vont pas jusqu'a l'une des deux extrémités
	"""
	to_parse=open(lire,"r")
	parsed=open(ecrire,"w")
	
	for hit in to_parse :
		ligne = hit.split("\t")
		
		cpt=0
		
		if ligne[1] == "1" : #SI mon read match depuis la position 1, compteur gagne 1
			cpt+=1
			
		if str(ligne[2]) == str(ligne[3]) : #Si mon read match jusqu'au bout, gagne 1
			cpt+=1
			
		if cpt == 1 : #Je l'écris que si un des deux est vrai
			parsed.write(str(hit))
				
	
	to_parse.close()
	parsed.close()
		
#	Read trimming :
	
def recup_posi(parsed): 
	"""
	Cree un dictionnaire avec : Nom de la sequence : [début,fin] de ce qu'il faut enlever et le retourne
	"""
	
	parse=open(parsed,"r")
	livre={}
	noms=[]
	vals=[]
	
	for line in parse :
		reads=line.split("\t")
		noms.append(reads[0])				#Ajoute les noms dans une liste
		vals.append([reads[1],reads[2],reads[4]])	#Ajoute les debut et fin dans une autre liste
		
	lenth = len(noms)
	for j in range(lenth):					#En fait un dictionnaire
		livre[noms[j]]=vals[j]
	
	parse.close()
	return livre
	
def lire_reads(reads): 
	"""
	Cree un dictionnaire : 
	Clé = Nom du read
	Valeurs = séquence du read
	"""
	
	read=open(reads,"r")
	liste=[]
	seq=[]
	noms=[]
	bouquin = {}
	tot=[]
	for line in read :
		if line[0] != ">":
			seq.append(line.replace("\n",""))					
		if line[0] == ">":
			noms.append(line.replace(">","").replace("\n",""))
			if seq != [] :
				tot.append("".join(seq))
			seq=[]
	tot.append("".join(seq))
	length = len(noms)
	for j in range(length):
		bouquin[noms[j]]=tot[j]
	
	read.close()
	return bouquin	

def split_read(unmapped,parsed,ecrire,TAILLE_MIN_TRIMMED_READ): 
	"""
	Prend les reads en fasta, les réécris en fasta mais sans leur partie ET
	"""
	
	unmap = open(unmapped,"r")
	trimd=open(ecrire,"a")
	
	cpt=0
	Dictionnaire_des_sequences = lire_reads(unmapped) #Dictionnaire des sequences
	Dictionnaire_des_positions = recup_posi(parsed)	#Dictionnaire des positions
	
	
	
	for line in unmap :
		if line[0] == ">" :
			cpt+=1
			nom = line.replace(">","")
			nom = nom.replace("\n","")
			
			seq=Dictionnaire_des_sequences[nom] #Recup la sequence du read
			
			if nom in Dictionnaire_des_positions :
				pos = Dictionnaire_des_positions[nom] #Recup position a supprimer du read, si le read est conservé --par le parsage, filtrage
				
			if nom not in Dictionnaire_des_positions : #Si il n'y est pas, position prend -1 pour eviter les if suivants
				pos =-1
				
			if pos != -1 :	#Si le read est conservé :
				#print 'ok'
				trimmed_read = seq[:int(pos[0])-1]+seq[int(pos[1]):] #Le read trimmé prend les valeurs complémentaires de celles données
				
				if len(trimmed_read)>TAILLE_MIN_TRIMMED_READ : #Ecrire le read trimmé que si sa taille est supérieure à la valeur choisie
					trimd.write(">"+str(nom)+" " + str(pos[2].replace("\n",""))+"\n"+str(trimmed_read)+"\n")
					
			#if cpt%1000==0 :
				#print str(cpt) + " reads traites."
			
	trimd.close()
	unmap.close()
	

#	 Step 3

def recherche_clusters(liste_des_FBTi):
	clusters_existants = []
	for ET in liste_des_FBTi :
		if ET not in clusters_existants :
			clusters_existants.append(ET)
	return clusters_existants

def cree_les_clusters(fasta_total) :
	tot = open(fasta_total,"r")
	dictionnaire_tot={}
	for line in tot :
		if line[0] == ">" :
			fbt = line.split(" ")[1].replace("\n","")
			Name = line
		else :
			if dictionnaire_tot.has_key(fbt) :
				dictionnaire_tot[fbt].append(line.replace(" \n","").replace("\n",""))
			else :
				dictionnaire_tot[fbt] = []
				dictionnaire_tot[fbt].append(line.replace(" \n","").replace("\n",""))
			
	tot.close()		
	return dictionnaire_tot
	
def decoupe_le_dico(dico,parts):
	number=0
	pages = len(dico)
	chapitre = pages/2
	decoupage=[]
	for i in range(parts) :
		
		temp_dic = {}
		decoupage.append(temp_dic)
		
	
	for cle in dico :
			
		decoupage[number%parts][cle] = dico[cle]
		number+=1
	return decoupage

def modifie_contigs(contigs,trie) :
	a_trier = open(contigs,"r")
	tried = open(trie,"w")
	for line in a_trier :
		if line[0:3] == "FBt" :
			stock_cluster = line
		elif line[0] == ">" :
			tried.write(str(line).replace("\n","")+""+str(stock_cluster))
		else :
			tried.write(line)
	a_trier.close()
	tried.close()
		
def associe_contig_TE(cap,trimmed) :
	
	"""
	Retourne un dictionnaire qui contient : Contig = ET
	A partir de la sortie de cap et du fichier des trimmed reads
	"""
	capout = open(cap,"r")
	trim=open(trimmed,"r")
	contread={}
	cont_et={}
	readte={}
	assos=[]
	readnames=[]
	cptr=0
	for line in capout : 
		contig = line.split(" ")
		ici=0
		
		if contig[0] == "*******************" :
			assos.append("".join(contig[1:3]))
			cptr = 1
			ici = 1
		if cptr == 1 and ici == 0 :
			taille = len(line)
			readnames.append(line[:taille-2])
			
			
			cptr=0
			ici=0
		if line == "DETAILED DISPLAY OF CONTIGS\n" :
			break
	for i in range(len(assos)) :
		contread[assos[i]] = readnames[i]
	
	readnames=[]
	tenames=[]
	for line in trim :
		if line[0] == ">" :
			ret = line.split(" ")
			
			readte[ret[0].replace(">","")]=ret[1]
	
	for i in range(len(assos)) :
		
		cont_et[assos[i]] = readte[contread[assos[i]]].replace("\n","")
		
	
	return cont_et

def cb_delem(blast): 
	
	"""
	Quand on a récupéré les contigs, on les blast contre les ET, puis on vient ici,
	Cette fonction récupère la liste des contigs qu'on ne veut plus garder.
	"""
	
	bl=open(blast,"r")
	liste=[]
	
	for line in bl :
		contig = line.split("\t")
		if contig[0] not in liste :
			liste.append(contig[0])
	
	return liste
	
	bl.close()

def trie_les_bad_contigs(tous,je_les_veux,blast): 
	"""
	Lit le blast
	Ecrit les contigs que n'ont pas eu de hit, c'est a dire ceux qui sont génomiques
	"""
	
	tutti = open(tous,"r")
	tried = open(je_les_veux,"w")
	ct=0
	
	je_les_veux_pas=cb_delem(blast)
	
	for line in tutti :
		if line[0] == ">" :
			ct=0
			nom=line.replace(">","").replace("\n","")
			
			if nom not in je_les_veux_pas:
				
				ct= 1
		if ct == 1 :
			tried.write(str(line))
	
	tutti.close()
	tried.close()

def fasta_splitter(fast,parti) :
	import math
	#   To read/write FASTA. Wrap in a try/except
	try:
	    from Bio import SeqIO
	except ImportError:
	    print "This script requires BioPython to be installed!"
	
	
	#   Set these names so that it is easier to follow them in the code
	to_split = fast
	num_seqs = parti
	
	#   Make sure that the file provided is readable. 
	try:
	    open(to_split, 'r')
	except IOError:
	    print "Error! " + to_split + " does not exist, or is not readable!"
	    exit(1)
	
	#   And make sure that the second argument is a positive integer
	#   The try-Except block tests for integer, and the if tests for positive
	try:
	    int(num_seqs)
	except ValueError:
	    print "Error! Please provide a positive integer!"
	    exit(1)
	if int(num_seqs) <= 0:
	    print "Error! Please provide a positive integer!"
	    exit(1)
	
	#   Read in the FASTA file
	all_fasta = list(SeqIO.parse(to_split, 'fasta'))
	#   How many files to split into?
	num_files = int(math.ceil(len(all_fasta)/float(num_seqs)))
	#   Print a little message
	print "Will split " + to_split + " into " + str(num_files) + " files, with " + str(num_seqs) + " seqs per file."
	#   Start splitting!
	i = 1
	while i <= num_files:
	    #   Calculate the start
	    start = int(i-1) * int(num_seqs)
	    #   and the end
	    end = int(i) * int(num_seqs)
	    #   Generate the filename
	    #   strip off the extension, and get the rest
	    filename = to_split.split('.')[:-1][0] + '_' + str(i) + '.fasta'
	    #   Write the sequences into the file
	    SeqIO.write(all_fasta[start:end], filename, 'fasta')
	    #   increment the counter
	    i += 1



# Results :

def recup_carac(lire,ecrire):
	et = open(lire,"r") #Fichier a lire
	carac = open(ecrire,"w") #Fichier a ecrire
	element = []
	for line in et :
		if line[0] == ">" : #Regarde que les lignes titre
			k=0
			nom=[]
			start=[]
			end=[]
			ender=""
			starter=""
			name=""
			cp=0
			ch=[]
			element = line.split(" ") #Separe la ligne
			ID = element[0].replace(">","") #Id de la ligne
			loc=element[2] #Recup la localisation loc=3R:6158563..6167089;
			
			for j in loc :
				if j != ":" and cp == 0 :
					ch.append(j)
				else :
					cp += 1
					
			
			chd = "".join(ch)
			chd=chd.replace("loc=","")
			a = "loc="+str(chd)	
			loc_start = loc.replace(a,"")
			for i in loc_start :
				
				if i.isdigit():
					if k==0 :
						
						start.append(i)
					if k==1 :
						
						end.append(i)
				if i ==".": 
					k=1
			starter = "".join(start) #Debut
			ender="".join(end) #Fin
			name = element[3].replace("name=","")
			for j in name :
				if j != "{" :
					nom.append(j)
				else : 
					break
			name = "".join(nom) #Type element
			
			carac.write(str(ID)+"\t"+str(chd)+"\t"+str(starter)+"\t"+str(ender)+"\t"+str(name)+"\n")
	et.close()
	carac.close()

def read_sam(sam,decrire,ted) :
	"""
	Fonction qui lit un fichier sam et renvoie :
	nom du read + chr + positions génomique début + positions génomique fin
	en ajoutant la longueur a la position de début
	"""
	samm = open(sam,"r")
	dec=open(decrire,"w")
	dic = open(str(args.workdir)+str(ted),"r")
	fbtfam={}
	for linu in dic :
		fbtin = linu.split("\t")
		fbtfam[fbtin[0]]=fbtin[1]
	for line in samm :
		
		if line[0:6] == "Contig" :
			read = line.split("\t")
			#print read[0]
			if read[3] != "0" :
				dec.write(str(read[0].split(" ")[0])+"\t"+str(read[2])+"\t"+str(read[3])+"\t"+str(int(read[3])+int(re.findall('\d+',read[5])[0]))+"\t"+str(read[0][read[0].find("FBti"):])+"\t"+str(fbtfam[str(read[0][read[0].find("FBti"):])])+"\n")
			
	samm.close()
	dec.close()
	dic.close()

def regarde_a_cote(TE,contigs,SET_VALUE_AFTER_AND_BEFORE,chrs):
	"""
	Fonction qui regarde a coté d'UN TE si des contigs sont présents.
	Oui de chaque coté : Shared
	Oui + non : Doubt
	Non : Absent
	"""
	cot=open(contigs,"r")
	#TE c'est [nom,deb,fin]
	
	
	avant = int(TE[3]) - SET_VALUE_AFTER_AND_BEFORE
	apres = int(TE[4].replace("\n","")) + SET_VALUE_AFTER_AND_BEFORE
	
	post_valid="ABS"
	prec_valid="ABS"
	
	for line in cot :
		contig = line.split("\t")
		
		if contig[1] == chrs :
			deb = contig[2]
			fin = contig[3]
		
			if int(avant) > int(deb) and int(avant) < int(fin): #Regarde si avant tombe dans un contig
				if contig[5].replace("\n","")==TE[1] :
					prec_valid=contig[0]
			if int(apres) > int(deb) and int(apres) < int(fin): #Regarde après
				if contig[5].replace("\n","")==TE[1] :
					post_valid=contig[0]
	
	if prec_valid != "ABS" and post_valid != "ABS" : #Si les deux Shared 
		return ["Shared",prec_valid,post_valid]
	elif prec_valid != "ABS" or post_valid != "ABS": #Si un seul doubt
		return ["Shared",prec_valid,post_valid]
	else :
		return ["Absent",prec_valid,post_valid]	#Si aucun Absent
	
	
	cot.close

def revue_des_TE(TE,contigs,resultat_typage,valeur):
	"""
	Fonction qui regarde tout les TE et ecrit leur etat, shared etc..
	Retourne aussi les contigs qui ont été utilisés
	"""
	
	contigs_utilises=[]
	res = open(resultat_typage,"w")
	
	te = open(TE,"r")
	for line in te :
		transpo = line.split("\t")
		
		typ = regarde_a_cote(transpo,contigs,valeur,transpo[2])
		res.write(str(transpo[0])+"\t"+str(transpo[1])+"\t"+typ[0]+"\t"+typ[1]+"\t"+typ[2]+"\t"+transpo[2].replace("\n","")+"\t"+transpo[3]+"\t"+transpo[4].replace("\n","")+"\n")

		if typ[0] == "Shared" or typ[0] == "Doubt": #Ecrit le contigs utilisés par des ET validés
			if typ[1] != "ABS" :
				contigs_utilises.append(typ[1])
			if typ[2] != "ABS" :
				contigs_utilises.append(typ[2])
				
	te.close()
	
	res.close()
	return contigs_utilises
			
	

def trouve_nvelles_inser(contigs,revue,ecrire):
	"""
	Fonction qui prend les contigs utilisés et ecrit ceux inutilisés
	"""
	cont = open(contigs,"r")
	restants=open(ecrire,"w")
	for line in cont :
		contig = line.split("\t")
		if contig[0] not in revue :
			restants.write(line)
	cont.close()
	restants.close()

def compare_les_restants(remain,sortie,te_desc):
	"""
	Regarde si deux contigs restants sont collés
	A 10 bases pres
	Ecrit les contigs collés et leur positiosn génomiques
	"""
	ted = open(te_desc,"r")
	ref={}
	for line in ted :
		tete=line.split("\t")
		ref[tete[0]]=tete[1]
	rem=open(remain,"r")
	
	sor = open(sortie,"w")
	nvelles=[]
	used_contigs=[]
	#contig_to_te = associe_contig_TE(str(args.workdir)+str("cap3.out"),str(args.workdir)+str(args.preftoremove)+str(args.trimmed))
	for line in rem :
		contig=line.split("\t")
		dep=contig[2]
		fin = contig[3]
		ch = contig[1]
		FBT = contig[4]
		rem2 = open(remain,"r")
		for ligne in rem2 :
			contigue = ligne.split("\t")
			dep2=contigue[2]
			fin2=contigue[3]
			ch2=contigue[1]
			
			if int(dep) > int(fin2)-10 and int(dep) - 90 > int(dep2) and int(dep)-150 < int(fin2) and ch==ch2 and dep not in used_contigs and dep2 not in used_contigs and fin not in used_contigs and fin2 not in used_contigs and contigue[4] == FBT:
				sor.write(str(ch)+"\t"+ FBT.replace("\n","") +"\t"+ str(ref[FBT.replace("\n","")].replace("\n",""))+ "\t" + str(dep2) + "\t" + str(fin)+"\n")
				used_contigs.append(dep)
				used_contigs.append(dep2)
				used_contigs.append(fin)
				used_contigs.append(fin2)
			if int(fin) < int(dep2)+10 and int(fin) + 90 > int(dep2) and int(fin) + 150 < int(fin2) and ch == ch2 and dep not in used_contigs and dep2 not in used_contigs and fin not in used_contigs and fin2 not in used_contigs and FBT == contigue[4]:
				sor.write(str(ch)+"\t"+ FBT.replace("\n","") +"\t"+ str(ref[FBT.replace("\n","")].replace("\n",""))+ "\t" + str(dep) + "\t" + str(fin2)+"\n")
				used_contigs.append(dep)					#FBTi00555735		 fbt					#Famille ref[fbt]
				used_contigs.append(dep2)
				used_contigs.append(fin)
				used_contigs.append(fin2)
		rem2.close()
		
	
	rem.close()
	sor.close()	

#Pipeline :
#	Step 1 :

#	Recover unmapped to the reference reads :
if not os.path.exists("{0}{2}{1}".format(args.workdir,args.unmappedreadsonboth,args.preftoremove)) :
	if not os.path.exists(str(args.workdir)+str(args.indexref)+".4.bt2") :
		try :
			subprocess.call("/usr/remote/bin/bowtie2-build {0}{1} {0}{2} -q".format(args.workdir,args.maskedref,args.indexref),shell=True)
			print "Masked genome index is generated. Starting the mapper."
		
		except IOError :
			print "Your masked reference file doens't exists or is not readable."
			exit(1)
		
		except :
			print "Error while building bowtie2 index."
			exit(1)	
	else :
		print "Your index seems already done, skipping ..."
	
	if not os.path.exists(str(args.workdir)+str(args.preftoremove)+str(args.mappedreadsonref)) :	
		try :
			subprocess.call("/usr/remote/bin/bowtie2 -p {0} -x {1}{2} -r -q {1}{3} -S {1}{6}{4} --un {1}{6}{5} --end-to-end".format(args.cores,args.workdir,args.indexref,args.reads,args.mappedreadsonref,args.unmappedreadsonref,args.preftoremove),shell=True)
			print "First mapping is over."
		
		except IOError :
			print "Your reads file doesn't exists or is not readable."
			exit(1)
		
		except :
			print "Error while mapping the reads to the masked genome."
			exit(1)
	else :
		print "Reads seem already mapped to the reference. Skipping ..."	

#	Retrieve unmapped to the TE reads :
	if not os.path.exists(str(args.workdir)+str(args.indexte)) :
		try :
			subprocess.call("/usr/remote/bin/bowtie2-build {2}{0} {2}{1} -q".format(args.tebank,args.indexte,args.workdir),shell=True)
			print "TE index is ready."
		
		except IOError:
			print "Your TE bank does not exists or is not readable."
			exit(1)
		
		except :
			print "Error make the TE index."
			exit(1)
	else :
		print "TE index already exists, skipping ..."
				
	try :
		subprocess.call("/usr/remote/bin/bowtie2 -p {0} -x {5}{1} -r -q {5}{6}{2} -S {5}{6}{3} --un {5}{6}{4} --end-to-end".format(args.cores,args.indexte,args.unmappedreadsonref,args.mappedreadsonref,args.unmappedreadsonboth,args.workdir,args.preftoremove),shell=True)
		print "Reads are mapped."
	except  :
		print "Error while mapping reads to the TE index."
		exit(1)
	print "First step is over."

else :
	print "It seems like reads are already filtered. Skipping step 1 ..."

#	Step 2 :

#	Reads trimming :
#	Thread to start blast in parallel :
class Trimmage(th.Thread) :
	"""
	Thread qui va permettre de blaster en parallèle mes reads non mappés
	"""
	#LOCK for the output trimmed
	lock_trimmed = th.Lock()
	
	def __init__(self,fastq,para,treda) :
		th.Thread.__init__(self)
		self.fastq=fastq
		self.para=para
		self.treda=treda
	
	def run(self):
		
		try :
			self.lance_blast()
		except :
			print "Error running blast."
			exit(1)
		
		try :
			parsage(str(args.workdir)+str(args.preftoremove)+str(args.blastedreads)+str(".")+str(self.para),str(args.workdir)+str(args.preftoremove)+str(args.opnf)+str(".")+str(self.para),str(args.workdir)+str(args.preftoremove)+"unmapped_reads"+"_"+str(self.treda)+str(".fasta"))	
		except IOError : 
			print "Input problem to parse the blast output."
			exit(1)
		except Exception as err :
			print "Error parsing blast."
			raise err
		
		try :
			enleve_mauvais_reads(str(args.workdir)+str(args.preftoremove)+str(args.opnf)+str(".")+str(self.para),str(args.workdir)+str(args.preftoremove)+str(args.opf)+str(".")+str(self.para))
		except :
			print "Error filtering the parse output."
			exit(1)
		
		try :	
			self.lance_split()
		except :
			print "Error trimming the reads and writing them."
			exit(1)
		
		try :	
			subprocess.call("rm {0}{1}{2}{3}{4}".format(args.workdir,args.preftoremove,args.opnf,".",str(self.para)),shell=True)
			
			subprocess.call("rm {0}{1}{2}{3}{4}".format(args.workdir,args.preftoremove,args.blastedreads,".",str(self.para)),shell=True)
		
			subprocess.call("rm {0}{1}{2}{3}{4}".format(args.workdir,args.preftoremove,"unmapped_reads_",str(self.treda),".fasta"),shell=True)

			subprocess.call("rm {}{}{}{}{}".format(args.workdir,args.preftoremove,args.opf,".",str(self.para)),shell=True)
			subprocess.call("rm {}{}{}{}*".format(args.workdir,args.tebank,str(self.para),"."),shell=True)
		except :
			print "Error deleting files from the triming thread."
			exit(1)
	def lance_split(self):
		self.lock_trimmed.acquire()
		try :
			split_read(str(args.workdir)+str(args.preftoremove)+"unmapped_reads_"+str(self.treda)+str(".fastarearr"),str(args.workdir)+str(args.preftoremove)+str(args.opf)+str(".")+str(self.para),str(args.workdir)+str(args.preftoremove)+str(args.trimmed),args.lenmin)		
		except :
			print "Error while splitting the reads."
			exit(1)
		finally :
			self.lock_trimmed.release()
			
	def lance_blast(self):
		try :
			print ("blastall -p blastn -d {2}{0} -i {2}{3} -o {2}{4}{1} -m 8 -e 10e-10".format(str(args.tebank)+str(self.para),str(args.blastedreads)+str(".")+str(self.para),args.workdir,str(args.preftoremove)+str("unmapped_reads")+"_"+str(self.treda)+str(".fasta"),args.preftoremove))
			subprocess.call("blastall -p blastn -d {2}{0} -i {2}{3} -o {2}{4}{1} -m 8 -e 10e-10".format(str(args.tebank)+str(self.para),str(args.blastedreads)+str(".")+str(self.para),args.workdir,str(args.preftoremove)+str("unmapped_reads")+"_"+str(self.treda)+str(".fasta"),args.preftoremove),shell=True)
		except :
			print "Error blasting the reads."
			exit(1)

#	Start the tread :
if not os.path.exists(str(args.workdir)+str(args.preftoremove)+str(args.trimmed)):
	
	# Fastq to fasta for the unmapped reads :
	try :
		subprocess.call("awk '{0}' {1} > {2}".format('NR % 4 == 1 {gsub("@",">",$O); print} NR % 4 == 2 {print $0}',str(args.workdir)+str(args.preftoremove)+str(args.unmappedreadsonboth),str(args.workdir)+str(args.preftoremove)+str("unmapped_reads.fa")),shell=True)
	except :
		print "Error while converting reads from FastQ to fastA."
		exit(1)
		
	#	Count the number of reads :
	try :
		fd = open(str(args.workdir)+str(args.preftoremove)+str("unmapped_reads.fa"), 'r')
		japps = 0
		for linz in fd :
			if linz[0] == ">":
				japps+=1
		fd.close()
	except :
		print "Error while counting the number of unmapped reads."
		exit(1)
	
	# Make the database for the blast :
	try :	
		subprocess.call("formatdb -i {} -p F".format(str(args.workdir)+str(args.tebank)),shell=True)
		print "Database is ready for the blast, starting soon."
	except :
		print "Error while making TE database for the blast."
		exit(1)
	
	# Copy the database for the threads :
	try :	
		for numberz in range(args.cores):
			subprocess.call("cp {0}.nin {0}{1}.nin".format(str(args.workdir)+str(args.tebank),numberz),shell=True)
			subprocess.call("cp {0}.nhr {0}{1}.nhr".format(str(args.workdir)+str(args.tebank),numberz),shell=True)
			subprocess.call("cp {0}.nsq {0}{1}.nsq".format(str(args.workdir)+str(args.tebank),numberz),shell=True)
	except :
		print "Error while making a copy of the database for blast."
		exit(1)
	
	#	Split the fasta file of the unmapped reads :
	try :
		fasta_splitter(str(args.workdir)+str(args.preftoremove)+str("unmapped_reads.fa"),int(int(japps)/int(args.cores))+1)
	except :
		print "Error while splitting the unmapped reads."
		exit(1)
		
	#	Create a file for the trimmed reads :
	try :
		f=open(str(args.workdir)+str(args.preftoremove)+str(args.trimmed),"w")
		f.close()
	except :
		print "Error while creating the file for the trimmed reads, maybe bad input name."
		exit(1)
	
	fqablast=[]
	thepara=0
	nomsdesfasta = []
	
	#	Starting the thread :
	for awkdumec in range(args.cores) :
		nomsdesfasta.append(int(awkdumec)+1)
	for treadation in nomsdesfasta :
		fqablast.append(Trimmage(str(args.workdir)+str(args.preftoremove)+str("unmapped_reads")+"_"+str(treadation)+str(".fasta"),thepara,treadation))
		thepara+=1
	for Tredent in fqablast :
		Tredent.start()
	for tadam in fqablast :
		tadam.join()
	
#	Step 3 :

#	Assembly :
#	Thread for the assembly :
class Assemblage(th.Thread) :
	
	"""This thread is made for the parralel assembly of the clusters"""
	
	#LOCK OBLECT
	my_lock = th.Lock()
	
	
	def __init__(self,dico,para):
		th.Thread.__init__(self)
		self.dico = dico
		self.para= para
	
	def run(self):
		
		"""Il va faire ça quand je vias executer le Thread"""
		"""Ecrire un fasta du dico, executer CD-HIT"""
		try :
			c=0
			for read in self.dico :
				assemb = open("{0}temp_fasta{1}.fa".format(str(args.workdir)+str(args.preftoremove),read),"a")
				
				for seq in self.dico[read] :
					assemb.write(">{0} {1}\n".format(c,read))
					assemb.write(str(seq)+"\n")
					c +=1
			
				assemb.close()
				self.utilise_cap3("{0}temp_fasta{1}.fa".format(str(args.workdir)+str(args.preftoremove),read))
				self.concatene("{0}temp_fasta{1}.fa".format(str(args.workdir)+str(args.preftoremove),read),read)
				subprocess.call("rm {0}temp_fasta{1}*".format(str(args.workdir)+str(args.preftoremove),read),shell=True)
				subprocess.call("rm {0}cap3temp_fasta{1}*".format(str(args.workdir)+str(args.preftoremove),read),shell=True)
		except :
			print "Error while assembling the trimmed reads, maybe a thread error."
			exit(1)		
	def utilise_cap3(self,clust) :
		
		subprocess.call("/usr/remote/bin/cap3 {0} > {1}".format(clust,"{0}cap3{1}.out".format(str(args.workdir)+str(args.preftoremove),clust.rsplit("/")[-1].replace(args.preftoremove,""))),shell=True)
	
	def concatene(self,clust,read) :
		self.my_lock.acquire()
		try :
			subprocess.call("echo {1} >> {0}contigs.fa".format(str(args.workdir)+str(args.preftoremove),read),shell=True)

			subprocess.call("cat {0} >> {1}contigs.fa".format(clust+".cap.contigs",str(args.workdir)+str(args.preftoremove)),shell=True)
			
		except :
			print "Error while cating the contigs."
			exit(1)
		finally :
			self.my_lock.release()

#	Start the thead :
if not os.path.exists(str(args.workdir)+str(args.preftoremove)+"contigs_cluster.fa"):
	try :
		dicto = decoupe_le_dico(cree_les_clusters(str(args.workdir)+str(args.preftoremove)+str(args.trimmed)),args.cores)			
	except :
		print "Error while splitting the trimmed reads in clusters."
		exit(1)
		
	contigs_apres_cap = open(str(args.workdir)+str(args.preftoremove)+"contigs.fa","w")
	contigs_apres_cap.close()			
	
				
	#Création des threads
	threads = []
	for thread in range(args.cores) :
		threads.append(Assemblage(dicto[thread],thread))
	
	for tread in threads :
		tread.start()
	
	for td in threads :
		td.join()
	
	
	#subprocess.call("rm /pandata/dechaud/passage_cdhit/output*",shell=True)
	
	modifie_contigs(str(args.workdir)+str(args.preftoremove)+"contigs.fa",str(args.workdir)+str(args.preftoremove)+"contigs_cluster.fa")
else :
	print "Clusters already exists, skipping ..."	
#	Check if the contigs are real contigs :
if not os.path.exists(str(args.workdir)+str(args.preftoremove)+str(args.blastedcontigs)):
	try :
		subprocess.call("blastall -p blastn -d {0} -i {1} -o {2} -m 8 -e 10e-20".format(str(args.workdir)+str(args.tebank),str(args.workdir)+str(args.preftoremove)+"contigs_cluster.fa",str(args.workdir)+str(args.preftoremove)+str(args.blastedcontigs)),shell=True)
	except :
		print "Error while blasting the contigs."
		exit(1)
	print "Contigs are blasted."
else :
	print "Contigs already blasted, skipping ..."

#	Retrieve the right contigs :	
if not os.path.exists(str(args.workdir)+str(args.preftoremove)+str(args.sortdcontigs)):
	try :
		tout_les_contigs = str(args.workdir)+str(args.preftoremove)+"contigs_cluster.fa"
		trie_les_bad_contigs(tout_les_contigs,str(args.workdir)+str(args.preftoremove)+str(args.sortdcontigs),str(args.workdir)+str(args.preftoremove)+str(args.blastedcontigs))
		print "Contigs are filtered"
	except :
		print "Error while sorting the contigs."
		exit(1)
else : 
	print "Contigs aleady sorted, skipping ..."

#	Step 4 :

if args.ted :
	#	Contigs mapping :
	if not os.path.exists(str(args.workdir)+str(args.preftoremove)+str(args.mappedcontigsonref)) :
		try :
			subprocess.call("/usr/remote/bin/bowtie2-build {0} {1}index_ref_dmel -q".format(str(args.workdir)+str(args.ref),args.workdir),shell=True)
			subprocess.call("/usr/remote/bin/bowtie2 -p {0} -x {3}index_ref_dmel -r -f {1} -S {2} --end-to-end".format(args.cores,str(args.workdir)+str(args.preftoremove)+str(args.sortdcontigs),str(args.workdir)+str(args.preftoremove)+str(args.mappedcontigsonref),args.workdir),shell=True)
			print "Contigs are mapped to the genome."
		except :
			print "Error while mapping contigs to the genome."
			exit(1)
	else :
		print "Contigs are already mapped to the genome, skipping ..."		
	
	#	Write results :	
	if not os.path.exists(str(args.workdir)+str(args.preftoremove)+str(args.outputfrknownte)) :
		try :
			read_sam(str(args.workdir)+str(args.preftoremove)+str(args.mappedcontigsonref),str(args.workdir)+str(args.preftoremove)+str(args.rd),args.ted)
		except :
			print "Error reading samfile of the contigs, maybe problem with the mapper."
			exit(1)
		
		try :	
			rev = revue_des_TE(str(args.workdir)+str(args.ted),str(args.workdir)+str(args.preftoremove)+str(args.rd),str(args.workdir)+str(args.preftoremove)+str(args.outputfrknownte),args.distn)
		except :
			print "Error while finding shared TE."
			exit(1)
		try :
			
			trouve_nvelles_inser(str(args.workdir)+str(args.preftoremove)+str(args.rd),rev,str(args.workdir)+str(args.preftoremove)+str(args.rest))
			compare_les_restants(str(args.workdir)+str(args.preftoremove)+str(args.rest),str(args.workdir)+str(args.preftoremove)+str(args.onew),str(args.workdir)+str(args.ted))
			print "TE_Hunter finished. Thank you for using it."
		except:
			print "Error while finding new insertions."
			exit(1)
	else :
		print "Output files already exists, please delete output files to restart the analysis."		

else :
	
	if not os.path.exists(str(args.workdir)+str(args.preftoremove)+str(args.mappedcontigsonref)) :
		
		
		try :
			subprocess.call("/usr/remote/bin/bowtie2-build {0} {1}index_ref_dmel -q".format(str(args.workdir)+str(args.ref),args.workdir),shell=True)
			subprocess.call("/usr/remote/bin/bowtie2 -p {0} -x {3}index_ref_dmel -r -f {1} -S {2} --end-to-end".format(args.cores,str(args.workdir)+str(args.preftoremove)+str(args.sortdcontigs),str(args.workdir)+str(args.preftoremove)+str(args.mappedcontigsonref),args.workdir),shell=True)
			print "Contigs are mapped to the genome."
		except :
			print "Error while mapping contigs to the genome."
			exit(1)
	else :
		print "Contigs are already mapped to the genome, skipping ..."		
	
	#	Write results :	
	if not os.path.exists(str(args.workdir)+str(args.preftoremove)+str(args.outputfrknownte)) :
		try :
			genere_ted(args.workdir+args.tebank,args.workdir+args.preftoremove+"tedescr.txt")
			read_sam(str(args.workdir)+str(args.preftoremove)+str(args.mappedcontigsonref),str(args.workdir)+str(args.preftoremove)+str(args.rd),args.preftoremove+"tedescr.txt")
		except :
			print "Error reading samfile or creating TE description format."
			exit(1)
		
		try :	
			rev = revue_des_TE(args.workdir+args.preftoremove+"tedescr.txt",str(args.workdir)+str(args.preftoremove)+str(args.rd),str(args.workdir)+str(args.preftoremove)+str(args.outputfrknownte),args.distn)
		except :
			print "Error while finding shared TE."
			exit(1)
		try :
			
			trouve_nvelles_inser(str(args.workdir)+str(args.preftoremove)+str(args.rd),rev,str(args.workdir)+str(args.preftoremove)+str(args.rest))
			compare_les_restants(str(args.workdir)+str(args.preftoremove)+str(args.rest),str(args.workdir)+str(args.preftoremove)+str(args.onew),args.workdir+args.preftoremove+"tedescr.txt")
			print "TE_Hunter finished. Thank you for using it."
		except:
			print "Error while finding new insertions."
			exit(1)
	else :
		print "Output files already exists, please delete output files to restart the analysis."		


# 	End
