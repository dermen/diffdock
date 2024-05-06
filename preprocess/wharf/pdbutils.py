

"""
Code for the manipulation and correction of PDB files.

TJL, May 2022
<thomas.joseph.lane@gmail.com>
"""

import re
from wharf.logger import log


# ------------------------ REFERENCE LIBRARY ----------------------------------


# --- CRYSTALLIZATION REAGENTS
# a list of "uninteresting" crystallization reagents
# that we mostly want to remove from PDB files
CRYST_REAGENTS = ['DMS', 'IMD', 'HOH', 'CA', 'CL', 'PEG', 'S', 'EDO', 
                  'GOL', 'MPD', 'AU', 'ZN', 'SO4', 'NO3', 'SE', 'MG', 
                  'HG', 'MES', 'NI', 'FMT', 'ACT', 'IPA', 'PO4', 'ACE'
                  'H2S', 'PG4', 'MLA', 'MRD', 'DIO', 'DTZ', 'NA', 'ZN',
                  'SE', 'GLY', 'AU', 'ACY']


# --- PDB slices for ATOM and HETATM records
# https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
_RECORD         = slice( 0,  6)
_ATOM_NUMBER    = slice( 6, 11)
_ATOM_NAME      = slice(11, 16)
_ALT_LOC_IND    =       16
_RESIDUE_NAME   = slice(17, 20)
_CHAIN_ID       =       21
_RESIDUE_NUMBER = slice(22, 26)
_COORDINATES    = slice(31, 54)

_HYDROGENS = ['H', 'H1', 'H2', 'H3', 'HE', 'HE3', 'HG1', 'HH', 'HH2', 'HZ', 'HZ2', 'HZ3',
              '1H', '1HA', '1HE2', '1HH1', '1HH2', '1HZ', '2H', '2HA', '2HE2',
              '2HH1', '2HH2', '2HZ', '3H', '3HB', '3HD2', '3HE', '3HG1', '3HZ',
              'HA', 'HB2', 'HB3', 'HG', 'HA2', 'HA3', 'HD1', 'HD2', 'HD3', 'HE1',
              'HE2', 'HE3', 'HG1', 'HG2', 'HG3', 'HH11', 'HH12', 'HH21', 'HH22',
              'HB1', 'HZ1', 'HB', 'HG11', 'HG12', 'HG13', 'HG12', 'HG21', 'HG22',
              'HG23', 'HE11', 'HE21', 'HE12', 'HE22', 'HD11', 'HD12', 'HD22', 'HD21',
              'HD13', 'HD31', 'HD23']

# --- ATOM ORDERING CONVENTION FOR PDB FILES
# https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1992.pdf
# see also https://bmrb.io/ref_info/atom_nom.tbl [not used yet, but possibly useful]
_ATOM_ORDER_TABLE = {
    'ALA' : ['N', 'CA', 'C', 'O', 'CB'],
    'ARG' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
    'ASP' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2'],
    'CYS' : ['N', 'CA', 'C', 'O', 'CB', 'SG'],
    'GLU' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
    'PHE' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'GLY' : ['N', 'CA', 'C', 'O'],
    'HIS' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
    'ILE' : ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'],
    'LYS' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'],
    'LEU' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'],
    'MET' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'],
    'ASN' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2'],
    'PRO' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'],
    'GLN' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'],
    'SER' : ['N', 'CA', 'C', 'O', 'CB', 'OG'],
    'THR' : ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'],
    'VAL' : ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2'],
    'TRP' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 
             'CE3', 'CZ2', 'CZ3', 'CH2'],
    'TYR' : ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ',
             'OH'],

    # universal -- these are basically things (hydrogens!) I am not sure about :(
    #              they get appended to the end of every AA
    'UNI' : ['OXT',] + _HYDROGENS
}


# -----------------------------------------------------------------------------


def _pdb_to_dict(pdb_file_text):
    """
    turns a pdb file into a dict, one entry for each ATOM or HETATM record
      * keys are the first 22 characters
      * values are tuple: (the line position, last 22 characters)
    """

    d = {}
    for i,l in enumerate(pdb_file_text.split('\n')):
        if l[_RECORD] in ['ATOM  ', 'HETATM']:
            d[ l[11:26] ] = (i, l[:11], l[26:])
        
    return d


def _dict_to_pdb(pdb_dict):

    output = []

    for k in pdb_dict.keys():
        output.append( (pdb_dict[k][0], pdb_dict[k][1] + k + pdb_dict[k][2]) )

    txt = '\n'.join( [ x[1] for x in sorted(output) ] )

    return txt


def remove_hydrogens(pdb_file_text):
    """
    Remove ATOM and CONECT records accociated with (probably)
    uninteresting crystallization reagents.

    Paramaters
    ---------
    pdb_file_text : str
        pdb loaded as a string

    Returns
    -------
    pdb_file_text : str
        modified pdb string
    """

    output_text = []

    for l in pdb_file_text.split('\n'):

        if l[_RECORD] in ['ATOM  ', 'HETATM']:
            atom_name = l[_ATOM_NAME].strip()

            mg = re.search('([\d+|\s]H)', atom_name)

            if atom_name in _HYDROGENS:
                pass
            elif mg:
                log.debug('possible hydrogen excluded:' + atom_name)
            else:
                output_text.append(l)
        else:
            output_text.append(l)

    return '\n'.join(output_text)


def remove_cryst_reagents(pdb_file_text):
    """
    Remove ATOM and CONECT records accociated with (probably)
    uninteresting crystallization reagents.

    Paramaters
    ---------
    pdb_file_text : str
        pdb loaded as a string

    Returns
    -------
    pdb_file_text : str
        modified pdb string
    """

    output_text = []
    reagent_atoms = []

    def _is_cryst_reg(line):
        return (line[_RECORD] == 'HETATM') and \
            (line[_RESIDUE_NAME].strip() in CRYST_REAGENTS)
    
    for l in pdb_file_text.split('\n'):

        if _is_cryst_reg(l):
            atom_number = l[_ATOM_NUMBER].strip()
            reagent_atoms.append(atom_number)

    # in a second pass:
    # remove any CONECT records that involve the reagent atoms
    for l in pdb_file_text.split('\n'):

        if _is_cryst_reg(l):
            pass # remove line

        elif l[_RECORD] == 'CONECT':

            conect_atoms = l[6:].split()

            keep_conect_line = True
            for ca in conect_atoms:
                if ca in reagent_atoms:
                    keep_conect_line = False

            if keep_conect_line:
                output_text.append(l)

        else:
            output_text.append(l)

    return '\n'.join(output_text)


def split_ligands_protein(pdb_file_text):
    """
    Split the protein and any ligands into separate entries.

    Retains for all:
        - CRYST1
        - ORIGX{N}
        - SCALE{N}
    
    Retains for protein only:
        - TER

    Paramaters
    ---------
    pdb_file_text : str
        pdb loaded as a string

    Returns
    -------
    protein_ligands : list
        list of pdb_strings, with the protein as one element, and any
        individual ligands as the others
    """

    header       = []
    protein      = {}
    ligands      = {} # 3-letter code --> pdb_str
    ligand_atoms = {} # 3-letter code --> [ atom numbers ]

    header_records = ['CRYST1', 
                      'ORIGX1', 'ORIGX2', 'ORIGX3', 
                      'SCALE1', 'SCALE2', 'SCALE3']

    # 1. sort protein and individual ligands
    # right now, ligands are sorted by name only, so multiple instances
    # of the same ligand will end up all together in one entry

    for l in pdb_file_text.split('\n'):

        if l[_RECORD] in header_records:
            header.append(l)
        
        elif l[_RECORD] in ['ATOM  ', 'TER   ']:
            chain_id = l[_CHAIN_ID]
            if chain_id not in protein.keys():
                protein[chain_id] = []

            protein[chain_id].append(l)
        
        elif l[_RECORD] == 'HETATM':
            
            ligand_code = l[_RESIDUE_NAME]
            atom_number = l[_ATOM_NUMBER].strip()
            chain_id = l[_CHAIN_ID]

            if ligand_code in ligands.keys():
                if chain_id in ligands[ligand_code].keys():
                    ligands[ligand_code][chain_id].append(l)
                    ligand_atoms[ligand_code][chain_id].append(atom_number)
                else:
                    ligands[ligand_code][chain_id] = [l]
                    ligand_atoms[ligand_code][chain_id] = [atom_number]

            else: # first entry
                ligands[ligand_code]      = {}
                ligands[ligand_code][chain_id] = [l]
                ligand_atoms[ligand_code] = {}
                ligand_atoms[ligand_code][chain_id] = [atom_number]

            assert ligands.keys() == ligand_atoms.keys()

            for key in ligands.keys():
                assert ligands[key].keys() == ligand_atoms[key].keys()

    # 2. add CONECT record to ligands

    for l in pdb_file_text.split('\n'):
        if l[_RECORD] == 'CONECT':
            conect_atoms = l[6:].split()

            n_matches = 0

            # see if atoms in the CONECT record are found in our
            # previously saved list of atoms that make up each ligand

            for lig in ligand_atoms.keys():
                for chain_id in ligand_atoms[lig].keys():
                    if conect_atoms[0] in ligand_atoms[lig][chain_id]:
                        n_matches += 1
                
                        for a in conect_atoms:
                            if a not in ligand_atoms[lig][chain_id]:
                                log.warn('atom %s in CONECT not found in %s' % (a, lig))

                        ligands[lig][chain_id].append(l)

            if n_matches > 1:
                log.warn('single CONECT spans multiple ligands!')


    # 3. add header to all, format final output
    output = {}
    for chain_id in protein.keys():
        protein_here = header + protein[chain_id]
        output[chain_id] = [ '\n'.join(protein_here) ]

        for lig in ligands.keys():
            if chain_id in ligands[lig].keys():
                ligands[lig][chain_id] = header + ligands[lig][chain_id]
                output[chain_id].append( '\n'.join(ligands[lig][chain_id]) )

    return output

def enforce_cannonical_pdb_order(pdb_file_text):
    """
    Reorder atoms according to the PDB convention & renumber atom numbers
    so they are sequential. Does not change residue numbering.

    Paramaters
    ---------
    pdb_file_text : str
        pdb loaded as a string

    Returns
    -------
    pdb_file_text : str
        modified pdb string
    """


    def _atom_index(pdb_line):

        # returns an integer with the correct order of an atom in a residue
        # pdb_line argument is a single line of a PDB file
    
        a_name = pdb_line[_ATOM_NAME].strip()
        r_name = pdb_line[_RESIDUE_NAME].strip()
       
        if r_name not in _ATOM_ORDER_TABLE:
            raise KeyError('%s not valid amino acid' % r_name)

        atom_ordering = _ATOM_ORDER_TABLE[r_name] + _ATOM_ORDER_TABLE['UNI']
        if a_name not in atom_ordering:
            print(pdb_line)
            raise KeyError('%s not valid atom for AA %s' % (a_name, r_name))

        #print('***' , l[_RESIDUE_NUMBER], l)
        res_num = int(pdb_line[_RESIDUE_NUMBER])
        order_in_residue = atom_ordering.index(a_name)
        altconf = pdb_line[_ALT_LOC_IND]

        return (res_num, order_in_residue, altconf)


    def _renumber(residue_lines, atom_number):
        for i in range(len(residue_lines)):
            cl = residue_lines[i]
            residue_lines[i] = cl[:6] + str(atom_number).rjust(5) + cl[11:]
            atom_number += 1
        return atom_number


    header_lines  = []
    atom_lines    = []
    trailer_lines = []

    header = True
    text_split = pdb_file_text.split('\n')
    for l in text_split:

        if l[_RECORD].startswith('ATOM'):
            header = False # --> non-ATOM records go in trailer_lines from now on
            atom_lines.append(l)

        else: # not an ATOM record
            if header:
                header_lines.append(l)
            else:
                trailer_lines.append(l)

    atom_lines.sort(key=_atom_index)
    _renumber(atom_lines, 1)

    output = header_lines + atom_lines + trailer_lines

    assert len(output) == len(text_split), 'lines missing'

    return '\n'.join(output)


def intersection(*pdb_texts):
    """
    Returns many PDBs, but only with ATOM records that are present in all.
    """
    
    common_keys = None
    dicts = []

    for pt in pdb_texts:

        d = _pdb_to_dict(pt)

        if common_keys is None:
            common_keys = set(d.keys())
        else:
            common_keys = common_keys.intersection( set(d.keys()) )

        dicts.append(d)
    
    for d in dicts:
        to_rm = set(d.keys()) - common_keys
        for rm_key in to_rm:
            d.pop(rm_key)

        assert len(d.keys()) == len(common_keys)

    return [ _dict_to_pdb(d) for d in dicts ]


def alt_conf_residues(pdb_file_text):
    """
    Identify the residues in a PDB file that have alternative conformations
    modelled.

    Paramaters
    ---------
    pdb_file_text : str
        pdb loaded as a string

    Returns
    -------
    residues : list of ints
        A list of residue indicies (1-indexed, matching PDB) for which
        the passed PDB file contains one or more alternative conformers.
    """

    residues = []

    for l in pdb_file_text.split('\n'):

        # if there is an altconf, then
        # 16-th position is "A", "B", ...
        if l[_RECORD] == 'ATOM  ' and l[_ALT_LOC_IND] != ' ':

            chain         = l[_CHAIN_ID]
            residue_index = int( l[_RESIDUE_NUMBER] )
            residue_aa    = l[_RESIDUE_NAME]

            r = (chain, residue_aa, residue_index)

            if r not in residues:
                residues.append(r)

    return residues


def remove_altconfs(pdb_file_text):
    """
    Reorder atoms according to the PDB convention & renumber atom numbers
    so they are sequential. Does not change residue numbering.

    Paramaters
    ---------
    pdb_file_text : str
        pdb loaded as a string

    Returns
    -------
    pdb_file_text : str
        modified pdb string
    """

    subset_pdb = []

    for l in pdb_file_text.split('\n'):
        if l[_RECORD] == 'ATOM  ':

            # if there is an altconf, then
            # 16-th position is "A", "B", ...

            if l[_ALT_LOC_IND] == ' ':
                 subset_pdb.append(l)

            elif l[_ALT_LOC_IND] == 'A':
                # remove alt tag
                new_l = l[:_ALT_LOC_IND] + ' ' + l[_ALT_LOC_IND+1:]
                subset_pdb.append(new_l)

        else:
            subset_pdb.append(l)

    return '\n'.join(subset_pdb)


def split_models(pdb_file_text):
    """
    Split a multi-conformer PDB into N single model PDBs.

    Parameters
    ----------
    pdb_file_text : str
        The PDB, loaded as a single string

    Returns
    -------
    models : list
        list of strings, each a PDB model
    """

    r = []

    model       = ''
    model_index = 0

    for l in pdb_file_text.split('\n'):
        model += l + '\n'
        if l == 'ENDMDL':
            r.append(model)
            model = ''
            model_index += 1

    return r


def get_atom_coords(pdb_file_text, atom_search_terms):
    """
    Match a string to a PDB or PDBQT file and return the x/y/z coordinates
    of the atom that matches that string.

    Raise errors if we see no match or more than one match.

    Paramaters
    ---------
    pdb_file_text : str
        pdb loaded as a string

    atom_search_terms : list of str
        a list of "terms" that must all be in a single PDB file line --
        the coordinates returned will match that line

        example: ["CA", "CYS", "A", "145"]

    Returns
    -------
    residues : list of ints
        A list of residue indicies (1-indexed, matching PDB) for which
        the passed PDB file contains one or more alternative conformers.
    """

    n_match = 0

    for l in pdb_file_text.split('\n'):
        if l.startswith('ATOM'):

            # position 31 is where coordinates begin
            search = [ l[:31].find(t) for t in atom_search_terms ]
            if -1 not in search:
                coords = [ float(x) for x in l[_COORDINATES].split() ]
                n_match += 1

    if n_match == 0:
        raise RuntimeError('no atom matching: '
            '%s found' % ( str(atom_search_terms) ))
    elif n_match > 1:
        raise RuntimeError('%d atoms matching: '
            '%s found' % ( n_match, str(atom_search_terms) ))

    return coords


