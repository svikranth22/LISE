#!/usr/bin/python3

# LISE2011: to cluster significant triangles looking for active site

# LISE2022: rewritten in Python in 2022

# open

import sys
import os
import re
import getopt
import urllib.request
from Bio.PDB.PDBIO import Select
from Bio.PDB import PDBIO
from Bio.PDB import PDBParser

HOME_DIR = os.getcwd()
PDB_DIR = os.path.join(HOME_DIR, "PDB")
RESULTS_DIR = os.path.join(HOME_DIR, "results")
DETAILS_DIR = os.path.join(HOME_DIR, "details")
TMP_DIR = os.path.join(HOME_DIR, "tmp")
ILIGS = os.path.join(HOME_DIR, "iligs")
RCSB_URL = "https://files.rcsb.org/view/"

num_ipro = 0
num_ilig = 0


def ata(nm, aa,ele):
    at = ' '
    if ele == 'C':
        if nm[2] == ' ':
            at = '0'  # atom type assignment2022# 0: C2P#
        elif nm[2] == 'A':
            at = '2'  # atom type assignment2022# 2: C3P#
        elif nm[2] == 'B':
            if aa == "SER" or aa == "THR":
                at = '2'  # atom type assignment2022# 2: C3P#
            else:
                at = '1'  # atom type assignment2022# 1: C3N#
        elif nm[2] == 'G':
            if aa == "VAL" or aa == "LEU" or aa == "ILE" or aa == "MET" or aa == "THR" or aa == "PRO" or aa == "GLU" or aa == "GLN" or aa == "ARG" or aa == "LYS":
                at = '1'  # atom type assignment2022# 1: C3N#
            elif aa == "ASP" or aa == "ASN":
                at = '0'  # atom type assignment2022# 0: C2P#
            elif aa == "PHE" or aa == "TYR" or aa == "TRP":
                at = '3'  # atom type assignment2022# 3: CRN#
            elif aa == "HIS":
                at = '3'  # atom type assignment2022# 4: CRP#
            else:#if aa == "VAL" or aa == "LEU" or aa == "ILE" or aa == "MET" or aa == "PRO" or aa == "GLU" or aa == "GLN" or aa == "ARG" or aa == "LYS":
                at = ''  # atom type assignment2022# ?: C3N#
                #input('weird atom type')
        elif nm[2] == 'D':
            if aa == "LEU" or aa == "ILE" or aa == "LYS":
                at = '1'  # atom type assignment2022# 1: C3N#
            elif aa == 'GLU' or aa == 'GLN':
                at = '0'  # atom type assignment2022# 0: C2P#
            elif aa == 'ARG' or aa == 'PRO':
                at = '2'  # atom type assignment2022# 2: C3P#
            elif aa == 'PHE' or aa == 'TYR':
                at = '3'  # atom type assignment2022# 3: CRN#
            elif aa == 'HIS':
                at = '3'  # atom type assignment2022# 4: CRP#
            elif aa == 'TRP':
                if nm[3] == '1':
                    at = '4'  # atom type assignment2022# 4: CRP#
                elif nm[3] == '2':
                    at = '3'  # atom type assignment2022# 3: CRN#
                else:
                    at = ''  # atom type assignment2022# ?: C#
                    #input('weird atom type')
            else:
                at = ''  # atom type assignment2022# ?: C#
                #input('weird atom type')
        elif nm[2] == 'E':
            if aa == 'MET':
                at = '1'  # atom type assignment2022# 1: C3N#
            elif aa == 'LYS':
                at = '2'  # atom type assignment2022# 2: C3P#
            elif aa == 'PHE' or aa == 'TYR':
                at = '3'  # atom type assignment2022# 3: CRN#
            elif aa == 'HIS':
                at = '3'  # atom type assignment2022# 4: CRP#
            elif aa == 'TRP':
                if nm[3] == '2':
                    at = '4'  # atom type assignment2022# 4: CRP#
                elif nm[3] == '3':
                    at = '3'  # atom type assignment2022# 3: CRN#
                else:
                    at = ''  # atom type assignment2022# ?: C#
                    #input('weird atom type')
            else:
                at = ''  # atom type assignment2022# ?: C#
                #input('weird atom type')
        elif nm[2] == 'Z':
            if aa == 'ARG':
                at = '2'  # atom type assignment2022# 2: C3P#
            elif aa == 'PHE' or aa == 'TRP':
                at = '3'  # atom type assignment2022# 3: CRN#
            elif aa == 'TYR':
                at = '4'  # atom type assignment2022# 4: CRP#
            else:
                at = ''  # atom type assignment2022# ?: C#
        elif nm[2] == 'H':
            at = '3'  # atom type assignment2022# 3: CRN#
        else:
            at = ''  # atom type assignment2022# ?: C#
    elif ele == 'N':
        if nm[2] == ' ':
            if aa == 'PRO':
                at = '6'  # atom type assignment2022# 14: no type/NLA#
            else:
                at = '8'  # atom type assignment2022# 8: NRD#
        elif nm[2] == 'H' or nm[2] == 'Z': #Arg or Lys
            at = '7'  # atom type assignment2022# 7: NLC#
        elif aa == 'ASN' or aa == 'GLN' or aa == 'TRP':
            at = '8'  # atom type assignment2022# 8(6?): NRD/NLB#
        elif aa == 'HIS': #or aa == 'TRP':
            at = '9'  # atom type assignment2022# 9: NRE#
        elif aa == 'ARG': #NE of ARG
            at = '6'  # atom type assignment2022# 6: NLA#
        else:
            at = ''  # atom type assignment2022# ?: N#
            #input('weird atom type')
    elif ele == 'O':
        if nm[2] == ' ':
            at = '10'  # atom type assignment2022# 10: O2A#
        elif nm[2] == 'G' or nm[2] == 'H': #SER, THR or TYR
            at = '11'  # atom type assignment2022# 11: O3B#
        elif aa == 'ASN' or aa == 'GLN':
            at = '10'  # atom type assignment2022# 10: O2A#
        elif aa == 'ASP' or aa == 'GLU':
            at = '12'  # atom type assignment2022# 12: OLC#
        elif nm[2] == 'X':
            at = '12'  # atom type assignment2022# 12: OLC#
        else:
            at = ''  # atom type assignment2022# ?: N#
            #input('weird atom type')
    elif ele == 'S':
        at = '13' #atom type assignment2022# 13: S3N#
    else:
        at = ''  # atom type assignment2022# ?: ?#
        #input('weird atom type')

    return at


def process_pdb(PDB_path, IPRO_path, ILIG_path):

    pdbfile = open(PDB_path, 'r')
    inf = pdbfile.readlines()
    pdbfile.close()

    r = {'C':'1.7', 'N':'1.55','O':'1.52','S':'1.8'}

    pno = []  # atom no.
    pnm = []  # atom name
    presnm = []  # residue name
    pchain = []  # chain name
    presno = []  # residue no.
    px = []  # x coordinate
    py = []  # y coordinate
    pz = []  # z coordinate
    ptag = []
    ptail = []
    pele = [] #element
    pat = [] #atom type
    #palter = [] #alternative location
    lno = []  # atom no.
    lnm = []  # atom name
    lresnm = []  # residue name
    lchain = []  # chain name
    lresno = []  # residue no.
    lx = []  # x coordinate
    ly = []  # y coordinate
    lz = []  # z coordinate
    ltag = []
    ltail = []
    lele = [] #element
    #lat = [] #atom type
    global num_ipro  # nof pro atom
    global num_ilig  # nof lig atom
    # write to the first intermediate file

    iprofile = open(IPRO_path, "w")
    iligfile = open(ILIG_path, "w")
    
    s = ""

    for i in inf:

        a = re.match(r'(ATOM|HETATM)', i)

        if a and a.group() == 'ATOM':
            
            if i[16] != ' ' and i[16] != 'A':

                continue

            ptag.append(i[0:6]) #ATOM card/tag
            pno.append(int(i[6:11])) #atom number
            pnm.append(i[12:16]) #atom name
            presnm.append(i[17:20]) #residue/amino acid name
            pchain.append(i[21]) #chain name
            presno.append(int(i[22:26])) #residue number
            px.append(float(i[30:38]))
            py.append(float(i[38:46]))
            pz.append(float(i[46:54]))
            ptail.append(i[54:78])
            pele.append(i[77])
            pat.append(ata(pnm[num_ipro],presnm[num_ipro],pele[num_ipro])) #atom type assignment
            
            try:
                s = presnm[num_ipro]+" "+str(presno[num_ipro])+" "+ pnm[num_ipro]+" "+str(px[num_ipro])+" "+str(py[num_ipro])+" "+str(pz[num_ipro])+" "+r[pele[num_ipro]]+" "+pat[num_ipro] +'\n'
            except:
                s = presnm[num_ipro] + " " + str(presno[num_ipro]) + " " + pnm[num_ipro] + " " + str(px[num_ipro]) + " " + str(py[num_ipro]) + " " + str(pz[num_ipro]) + " 0 0\n"

            iprofile.write(s)

            num_ipro += 1

        elif a and a.group() == 'HETATM':  # and i[17:20] == 'STL':

            if i[17:20] == 'HOH':

                continue

            elif i[17:20] == ' CL' or i[17:20] == 'SO4':

                continue

            m = re.match('(ZN|MG|CA|FE|NA|MN| K|NI|CU|CO|CD|HG|PT|MO|AL)', i[76:78])

            if m:

                ptag.append(i[0:6])  # ATOM card/tag
                pno.append(int(i[6:11]))  # atom number
                pnm.append(i[12:16])  # atom name

                presnm.append(i[17:20])  # residue/amino acid name
                pchain.append(i[21])  # chain name
                presno.append(int(i[22:26]))  # residue number
                px.append(float(i[30:38]))
                py.append(float(i[38:46]))
                pz.append(float(i[46:54]))

                ptail.append(i[54:78])
                pele.append(i[77])
                pat.append(ata(pnm[num_ipro], presnm[num_ipro], pele[num_ipro]))  # atom type assignment

                s = presnm[num_ipro] + " " + str(presno[num_ipro]) + " " + pnm[num_ipro] + " " + str(px[num_ipro]) + " " + str(py[num_ipro]) + " " + str(pz[num_ipro]) + " 1.4 5\n"

                iprofile.write(s)
                num_ipro += 1
                continue

            ltag.append(i[0:6])
            lno.append(int(i[6:11]))
            lnm.append(i[12:16])
            lresnm.append(i[17:20])
            lchain.append(i[21])
            lresno.append(int(i[22:26]))
            lx.append(float(i[30:38]))
            ly.append(float(i[38:46]))
            lz.append(float(i[46:54]))
            ltail.append(i[54:78])
            lele.append(i[76:78])
            
            s = lresnm[num_ilig]+" "+str(lresno[num_ilig])+" "+ lnm[num_ilig]+" "+str(lx[num_ilig])+" "+str(ly[num_ilig])+" "+str(lz[num_ilig])+" "+'\n'

            iligfile.write(s)
            num_ilig += 1

    iprofile.seek(0)
    iprofile.close()

    iligfile.close()


def make_dirs():
    if not os.path.exists(HOME_DIR):
        os.mkdir(HOME_DIR)
    if not os.path.exists(PDB_DIR):
        os.mkdir(PDB_DIR)
    if not os.path.exists(TMP_DIR):
        os.mkdir(TMP_DIR)
    if not os.path.exists(RESULTS_DIR):
        os.mkdir(RESULTS_DIR)
    if not os.path.exists(DETAILS_DIR):
        os.mkdir(DETAILS_DIR)
    if not os.path.exists(ILIGS):
        os.mkdir(ILIGS)
    

def download_pdb(pdb_id):
    file_name = pdb_id.upper()+'.pdb' #.upper: make string upper cases
    source = RCSB_URL + file_name
    path = os.path.join(PDB_DIR, file_name)
    try:
        urllib.request.urlretrieve(source, path)
        return True
    except(urllib.error.HTTPError or TimeoutError or urllib.error.URLError):
        print('Error: Could not download',file_name, file=sys.stderr)
        return False


def print_usage():
    print('Usage:')
    print('\tlise -i <PDBid>')
    print('\tlise -f <PDBfilepath>')


def main(argv):
    
    if len(argv) == 0:
        print('Error: No sepcified agruments', file = sys.stderr)
        print_usage()
        sys.exit(2)
     
    try:
      opts, args = getopt.getopt(argv,"hi:f:",["help", "id=", "filepath="])
    except getopt.GetoptError:
        print('Error: Invalid arguments', file = sys.stderr)
        print_usage()
        sys.exit(2)
    
    make_dirs()
    
    for opt,arg in opts:
        if opt in ("-h", "--help"):
            print_usage()
            sys.exit(1)
        elif opt in ("-i", "--id"):
            
            if not re.match('^\d\w{3}$', arg):
                print('Error: Invalid PDB ID', file = sys.stderr)
                sys.exit(2)

            PDB_id = arg.upper()
            PDB_path = os.path.join(PDB_DIR, PDB_id + ".pdb")

            if not os.path.exists(PDB_path):
                print("Fetching...")
                if not download_pdb(PDB_id):
                    sys.exit(2)

            RESULTS_path = os.path.join(RESULTS_DIR, PDB_id + "_top10.pdb")
            DETAILS_path = os.path.join(DETAILS_DIR, PDB_id + "_top3.pdb")
                
            print("RESULTS_PATH:", RESULTS_path)
            print("DETAILS_PATH:", DETAILS_path)

            if (os.path.exists(RESULTS_path) or os.path.exists(DETAILS_path)) is False:

                IPRO_path = os.path.join(TMP_DIR, PDB_id + "_ipro.txt")
                ILIG_path = os.path.join(TMP_DIR, PDB_id + "_ilig.txt")
                ILIG_pdb = os.path.join(ILIGS, PDB_id + "_temp.pdb")
                
                class IlligSelect(Select):
                    def accept_residue(self, residue):
                        if residue.id[0] != "W":
                            if residue.id[0] != " ":
                                return 1
                            else:
                                return 0

                parser = PDBParser()
                structure = parser.get_structure(PDB_id, PDB_path)
                io = PDBIO()
                io.set_structure(structure)
                io.save(ILIG_pdb, IlligSelect())

                print("IPRO_PATH:", IPRO_path)
                print("ILIG_PATH:", ILIG_path)

                process_pdb(PDB_path, IPRO_path, ILIG_path)

                print("IPRO_ATOMS:", num_ipro)
                print ("ILIG_ATOMS:", num_ilig)

            sys.exit(1)

        elif opt in ("-f", "--filepath"):
            print_usage()
            sys.exit(1)


if __name__ == "__main__":
    main(sys.argv[1:])