'''
make a script for gga assembly that takes in
-fasta file of insert
-fasta file of plasmid
-enzyme

Have it check gibson assembly?

outputs
-full cloned assembly
-dna orf
-aa orf

also tells you
-overhangs for plasmid
-overhangs for insert
-length of construct
-clear errors
    - mismatched overhangs
    - out of frame
    - extra cut sites
'''
import sys
import glob
import getpass
import re
import argparse
import datetime; date = datetime.datetime.now().strftime('%Y%m%d')
import numpy as np
import pandas as pd
from Bio import SeqIO, PDB, SeqUtils, Seq, SeqFeature
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#!/software/containers/john_bercow.sif

# Ryan's domesticator with some minor modifications.

# ============================================
# LIBRARIES
# ============================================
import copy
from collections import Counter

# Domesticator uses dnachisel
import dnachisel
from dnachisel import DnaOptimizationProblem, NoSolutionError
from dnachisel import DEFAULT_SPECIFICATIONS_DICT
from dnachisel import Location
from dnachisel import Specification, SpecEvaluation



def reverse_complement(seq):
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    return ''.join(complement[x] for x in reversed(seq.upper()))
def is_aa(aa):
    return aa.upper() in 'ACDEFGHIKLMNPQRSTVWY'
# ============================================
# ARGUMENTS
# ============================================
def get_arguments(argv=None):
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description=" * Generates an twist-ready fasta file for ordering individual constructs from a folder of PDBs and/or a concatenated FASTA file.\n"
                        " * Appropriate overhangs for Golden Gate cloning into entry vector(s) of interest are added automatically.\n"
                        " * Reverse translation is performed with dnachisel with constraints inherited from Ryan Kibler.\n"
                        " * RECOMMENDED: check your GG assemblies at https://goldengate.neb.com/#!/\n"
                        " * Wondering why the script is called John Bercow? https://www.youtube.com/watch?v=VYycQTm2HrM&ab_channel=TheSun\n"
                        "\n"
                        " * AVAILABLE ENTRY VECTORS:\n"
                        " *** see Benchling>DBL Database>Cloning plasmids ***\n"
                        " * EXAMPLE COMMAND: python DBL_order.py input.fasta -g gg_vector.fasta \n"
                        " * If you are using cell free mix make sure to add --cell_free flag!"
            )
    # REQUIRED
    parser.add_argument("input", help='file containing name and amino acid sequence can be multiple file types:\n'
                            'fasta - Fasta format and filename must end in .fa or .fasta\n'
                            'seq   - sequences in format $Sequence $Name. Filename must end is .txt, .seq, or .tab\n'
                            'pdb   - path to a folder of pdbs. Will pull all pdbs from that folder\n'
                           ,type=str)
    parser.add_argument(
            'gg_vector',
            help='Fasta file of plasmid for Golden Gate cloning (determines the DNA adapters). Also determines the AA tags appended to the design in the FASTA output.',
            action='store',
            type=str,
            )
    parser.add_argument(
            '-p','--project',
            help='Name of the project this order is associated with. All spaces should be underscores, not case sensitive, will find if you only use substring.\n' \
                    'e.g. --project "Cancer_diffusion" or --project "cancer" will both work.\n'\
                    'List: [Cancer_diffusion","Enzyme_diffusion_Novonesis","Snakebite_diffusion","Fraunhofer_diffusion","pMHC_Binders","Migraine_binders",\
                    "Granzyme_binders","IFNAR_binders",\
                    "DIRM_TCR_modelling",\
                    "Close-loop",\
                    "ADAC",\
                    "Nanobody_design",\
                    "CD3epsilon_design",\
                    "ProteusAI",\
                    "Lipase_design",\
                    "FGFR2b_binder_design","CL_Application","Protein_de-immunizer",\
                    "Dsup_redesign","Target_characterization_module"]',
            action='store',
            type=str,
            )

    # Restriction enzyme options
    parser.add_argument(
            '--n_overhang',
            help='cutsite overhang on the n terminus',
            action='store',
            type=str,
            default='GCCG'
            )
    parser.add_argument(
            '--c_overhang',
            help='cutsite overhang on the c terminus',
            action='store',
            type=str,
            default='GAGG'
            )
    parser.add_argument(
            '--enzyme',
            help='cutsite overhang on the c terminus',
            action='store',
            type=str,
            default='BsaI'
            )

    #Reverse translation options
    parser.add_argument(
            '--avoid',
            help='sequences to avoid in the design',
            action='store',
            nargs='+',
            default=[]
            )
    parser.add_argument(
            '--species',
            help='codon optimisation will be performed for this species (e.g. e_coli, s_cerevisiae, h_sapiens, etc...)',
            action='store',
            type=str,
            default='e_coli'
            )
    parser.add_argument(
            '--cell_free',
            help="adds constraints to reverse translate with codon table for tobacco as well (for ALICE cell free mix)",
            action='store_true'
            )
    parser.add_argument(
            '--starting_kmers_weight',
            help="starting value for the kmers_weight setting of Domesticator (Default: 10). This parameter is linearly ramped (up to 100) over --n_domesticator_steps.",
            action='store',
            type=int,
            default=10
            )
    parser.add_argument(
            '--n_domesticator_steps',
            help="maximum number of Domesticator steps attempted (Default: 10). The kmers_weight parameter (which increases synthesiability of repetitive sequences) is linearly ramped up to 100 over this number of steps.",
            action='store',
            type=int,
            default=10
            )
    parser.add_argument(
            '--max_attempts',
            help="maximum number of reverse translation attempts at each Domesticator step (Default: 20). Since Domesticator is stochastic, re-running the optimisation problem with the same parameters can lead to different solutions.",
            action='store',
            type=int,
            default=20
            )
    parser.add_argument(
            '--max_length',
            help="maximum length of gene fragments. Twists's maximum is 5000 bp, but less can be specified if sequence complexity is a issue for synthesis and you want to force the generation of smaller fragments.",
            action='store',
            type=int,
            default=5000
            )
    parser.add_argument(
            '--do_not_fix_repeats',
            help=" Turn off checking for and fixing direct repeats in reverse translated dna.",
            action='store_true',
            )
    parser.add_argument(
            '--repeat_frag_size',
            help="length of the direct repeat to be avoided",
            action='store',
            type=int,
            default=12
            )
    # Name Options
    parser.add_argument(
            '--design_prefix',
            help='designs get IDs with this prefix (e.g. LM0001, LM0002, etc...)',
            action='store',
            type=str
            )
    parser.add_argument(
            '--design_id',
            help='increment design indices from this number.',
            action='store',
            type=int
            )
    parser.add_argument(
            '--verbose',
            help="increase the verbosity of th e output (recommended).",
            action='store_true'
            )
    parser.add_argument(
            '--no_adapters',
            help="adds cut site and sticky ends, but no additional adaptor sequence",
            action='store_true'
            )

    args = parser.parse_args(argv)
    return args

def read_input_sequences(args):
    seq_dict = {
    'design_name':[],
    'aa_sequence':[],
    'well_position':[],
    'readin_order':[],
    }
    filename = args.input.split('/')[-1]
    if filename == '': 
        filename = args.input.split('/')[-2]
        filetype = 'pdb'
    elif filename.split('.')[-1] == 'fa' or filename.split('.')[-1] == 'fasta':
        filetype = 'fasta'
        filename = filename.split('.')[0]
    elif filename.split('.')[-1] == 'txt' or filename.split('.')[-1] == 'seq' or filename.split('.')[-1] == 'tab':
        filetype = 'seq'
        filename = filename.split('.')[0]
    else: filetype = 'pdb'
    if filetype == 'fasta':
        fasta_sequences = list(SeqIO.parse(args.input, 'fasta'))
        for i, fasta in enumerate(fasta_sequences):
            seq_dict['design_name'].append(str(fasta.id))
            seq_dict['aa_sequence'].append(str(fasta.seq))
            seq_dict['readin_order'].append(i)
            seq_dict['well_position'].append(w96[i])
    elif filetype == 'seq':
        with open(args.input, 'r') as seqfile:
            for i, line in enumerate(seqfile):
                seq, name = line.split()
                seq_dict['design_name'].append(str(name))
                seq_dict['aa_sequence'].append(str(seq))
                seq_dict['readin_order'].append(i)
                seq_dict['well_position'].append(w96[i])
    elif filetype == 'pdb':
        pdbs = sorted(glob.glob(f'{args.input}/*.pdb'))
        print(f'Extracting sequences from {len(pdbs)} PDBs...')
        for i, pdb in enumerate(pdbs):
            pdb_name = pdb.split('/')[-1].replace('.pdb', '')
            pdb_parser = PDB.PDBParser(PERMISSIVE=1, QUIET=True)
            structure = pdb_parser.get_structure('design', pdb)
            sequences = []
            for model in structure:
                 for chain in model:
                    seq = ''
                    for residue in chain:
                        seq += SeqUtils.IUPACData.protein_letters_3to1[residue.resname.capitalize()]
                    if seq not in sequences:
                        sequences.append(str(seq))

            for j, seq in enumerate(sequences):
                if len(sequences) > 1:
                    seq_dict['design_name'].append(str(f'{pdb_name}_{i+1}'))
                    seq_dict['aa_sequence'].append(str(seq))
                    seq_dict['readin_order'].append(f"{i}_{j}")
                else:
                    seq_dict['design_name'].append(str(f'{pdb_name}'))
                    seq_dict['aa_sequence'].append(str(seq))
                    seq_dict['readin_order'].append(i)
                    seq_dict['well_position'].append(w96[i])
    df = pd.DataFrame(seq_dict)
    return df, filename

def get_binding_adapters(args, cuts):
    hard_coded_binding_adapters = ['GTTTAAAGGTCTCGGCCGT','GGAGGTGAGACCAAAGGA']# includes cut site spacers and sticky ends for bsaI. this is the normal addition for DBE
    outside_flank = ['GTTTAAA','AAAGGA']
    fw_cut, rv_cut, spacer_len, _ = cuts[args.enzyme]
    spacer = 'GTACTACGTAATGT'
    binding_adapters = [fw_cut + spacer[:spacer_len] + args.n_overhang + "T", "G" + args.c_overhang + spacer[-1*spacer_len:] + rv_cut]
    if args.no_adapters == False:
        binding_adapters = [outside_flank[0] + binding_adapters[0],binding_adapters[1] +  outside_flank[1]]
    if args.enzyme == 'BsaI' and args.no_adapters == False:
        assert binding_adapters == hard_coded_binding_adapters, f'Binding adapters are not correct for {args.enzyme}. Please check the hard-coded binding adapters in the script.'
    return binding_adapters

def check_aa_sequences(df, args, cuts):
    binding_adapters = get_binding_adapters(args, cuts)
    max_len_dna = args.max_length - len(binding_adapters[0]) - len(binding_adapters[1])
    max_length_aa = int(np.floor(max_len_dna/3)) # max len minus adapters, divided by 3 for aa
    for i, seq in enumerate(list(df['aa_sequence'])):
        if len(seq) > max_length_aa:
            print(f'  [!] Sequence {seq} is too long for the specified maximum length ({len(seq)} vs. {max_length_aa} aa). Base dna length is 5000 bp max. Did you make a mistake in the input file?')
            sys.exit("  ERROR: Sequence too long for twist synthesis. System exiting...")
        current_letter = ""
        count = 0
        for aa in seq:
            if not is_aa(aa.upper()):
                print(f'  [!] Sequence {seq} contains non-standard amino acid {aa}. Please check your input file.')
                sys.exit("  ERROR: Non-standard amino acid found. System exiting...")
            #check for repeats of the same amino acid
            if aa == current_letter:
                count += 1
                if count >= 4:
                    print(f'  [!] Sequence number {i+1} contains more than 4 repeats of the same amino acid {aa}. This will increase sequence comnplexity and possibly be bad for your protein.')
            else:
                current_letter = aa
                count = 1
            
    return

def check_for_duplicates(aa_sequences):
    all_aa_seq = len(aa_sequences)
    unique_aa_seq = len(np.unique(list(aa_sequences.values())))
    if all_aa_seq != unique_aa_seq:
        print(f'[!] Found duplicated sequences ({unique_aa_seq} unique sequences vs. {all_aa_seq} total sequences):')

        visited = set()
        dup = [x for x in aa_sequences.values() if x in visited or (visited.add(x) or False)]
        duplicates = {d:[] for d in dup}
        for k, v in aa_sequences.items():
            if v in dup:
                duplicates[v].append(k)

        for seq, ids in duplicates.items():
            for id in ids:
                print(f'>{id}')
            print(seq + '\n-----')

        print('ERROR: Duplicates. Verify your sequences and retry. System exiting...')
        sys.exit()

def get_project(args):
    # hard coded from the projects spreadsheet in the OneNote.
    projects = ["Cancer_diffusion","Enzyme_diffusion_Novonesis","Snakebite_diffusion","Fraunhofer_diffusion","pMHC_Binders","Migraine_binders",
            "Granzyme_binders","IFNAR_binders",
            "DIRM_TCR_modelling",
            "Close-loop",
            "ADAC",
            "Nanobody_design",
            "CD3epsilon_design",
            "ProteusAI",
            "Lipase_design",
            "FGFR2b_binder_design",
            "CL_Application",
            "Protein_de-immunizer",
            "Dsup_redesign",
            "Target_characterization_module",
            "Other"
            ]
    project = []
    if args.project:
        input_project = args.project
    else:
        print(f"No project was provided. Please include a project from the list with -p argument\n{projects}")
    for p in projects:
        if input_project.lower() in p.lower():
            project.append(p)
    if len(project) == 1:
        return project[0]
    if len(project) == 0:
        print(f"ERROR: Project not found\n Options are {projects}")
        sys.exit("No Matching Project")
    if len(project) > 1:
        print(f"ERROR: Project matches multiple options, pleas specify.\n Matching options are {project}")
        sys.exit("Too many matching projects")
    
# ============================================
# Domesticator
# ============================================

# Special k-mer minimasation objective from Ryan's domesticator.
class MinimizeNumKmers(Specification):
    """Minimizes a kmer score."""

    best_possible_score = 0

    def __init__(self, k=8, location=None, boost=1.0):
        self.location = location
        self.k = k
        self.boost = boost

    def initialize_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        """Return a customized kmer score for the problem's sequence"""
        sequence = self.location.extract_sequence(problem.sequence)
        all_kmers = [sequence[i : i + self.k] for i in range(len(sequence) - self.k)]
        number_of_non_unique_kmers = sum(
            [count for kmer, count in Counter(all_kmers).items() if count > 1]
        )
        score = -(float(self.k) * number_of_non_unique_kmers) / len(sequence)
        return SpecEvaluation(
            self,
            problem,
            score=score,
            locations=[self.location],
            message="Score: %.02f (%d non-unique %d-mers)"
            % (score, number_of_non_unique_kmers, self.k),
        )

    def label_parameters(self):
        return [("k", str(self.k))]

    def short_label(self):
        return f"Avoid {self.k}mers {self.boost}"

    def __str__(self):
        """String representation."""
        return "MinimizeNum%dmers" % self.k

DEFAULT_SPECIFICATIONS_DICT["MinimizeNumKmers"] = MinimizeNumKmers
def get_tobacco_table():
    codon_table = {'F': {'TTT': 0.58, 'TTC': 0.42},
            'L': {'TTA': 0.15,'TTG': 0.24,'CTT': 0.26,'CTC': 0.13,'CTA': 0.1,'CTG': 0.11},
            'I': {'ATT': 0.5, 'ATC': 0.25, 'ATA': 0.25},
            'M': {'ATG': 1.0},
            'V': {'GTT': 0.41, 'GTC': 0.17, 'GTA': 0.17, 'GTG': 0.25},
            'S': {'TCT': 0.26,'TCC': 0.13,'TCA': 0.23,'TCG': 0.07,'AGT': 0.17,'AGC': 0.13},
            'P': {'CCT': 0.37, 'CCC': 0.13, 'CCA': 0.4, 'CCG': 0.1},
            'T': {'ACT': 0.39, 'ACC': 0.19, 'ACA': 0.34, 'ACG': 0.09},
            'A': {'GCT': 0.43, 'GCC': 0.17, 'GCA': 0.32, 'GCG': 0.08},
            'Y': {'TAT': 0.57, 'TAC': 0.43},
            '*': {'TAA': 0.42, 'TAG': 0.19, 'TGA': 0.39},
            'H': {'CAT': 0.61, 'CAC': 0.39},
            'Q': {'CAA': 0.58, 'CAG': 0.42},
            'N': {'AAT': 0.61, 'AAC': 0.39},
            'K': {'AAA': 0.49, 'AAG': 0.51},
            'D': {'GAT': 0.69, 'GAC': 0.31},
            'E': {'GAA': 0.55, 'GAG': 0.45},
            'C': {'TGT': 0.58, 'TGC': 0.42},
            'W': {'TGG': 1.0},
            'R': {'CGT': 0.15,'CGC': 0.08,'CGA': 0.11,'CGG': 0.08,'AGA': 0.33,'AGG': 0.25},
            'G': {'GGT': 0.33, 'GGC': 0.17, 'GGA': 0.34, 'GGG': 0.16}}
    return codon_table

def reverse_translate(
        amino_acid_sequence,
        kmers_weight=1.0,
        cai_weight=1.0,
        hairpins_weight=1.0,
        max_tries=10,
        species='e_coli',
        avoid=['GGTCTC', 'GAGACC'],
        count=None,
        warnings={},
        cell_free=False
    ):
    '''
    Ryan's domesticator.
    '''

    # Generate a naively reverse-translated DNA sequence first.
    naive_dna_sequence = dnachisel.reverse_translate(amino_acid_sequence)

    # Sequence optimisation will happen across the whole sequence.
    location = Location.from_biopython_location(SeqFeature.FeatureLocation(0, len(amino_acid_sequence) * 3))

    # Add optimisation objectives.
    objectives = []
    objectives.append(MinimizeNumKmers(k=8, boost=kmers_weight, location=location))
    objectives.append(dnachisel.builtin_specifications.AvoidHairpins(boost=hairpins_weight, location=location))
    

    # Add optimisation constraints.
    constraints = []
    constraints.append(dnachisel.builtin_specifications.EnforceTranslation(location=location))

    # A series of sequence patterns to remove
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("AAAAA", location=location)) # terminator
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("TTTTT", location=location)) # terminator
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("CCCCCC", location=location)) # repetitions
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("GGGGGG", location=location)) # repetitions
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("ATCTGTT", location=location)) # T7/T3 RNA polymerase pausing
    constraints.append(dnachisel.builtin_specifications.AvoidPattern("GGRGGT", location=location)) # G-quadruplex?
    constraints.append(dnachisel.builtin_specifications.UniquifyAllKmers(k=20, include_reverse_complement=True, location=location)) # ensure no 20-mers are repeated in the sequence (including reverse complement)
    if cell_free:
        cf_table = get_tobacco_table()
        min_freq = 0.1
        objectives.append(dnachisel.builtin_specifications.MaximizeCAI(species=species, boost=cai_weight, location=location))
        objectives.append(dnachisel.builtin_specifications.MaximizeCAI(codon_usage_table=cf_table, boost=cai_weight, location=location))
        constraints.append(dnachisel.builtin_specifications.AvoidRareCodons(min_frequency=min_freq,codon_usage_table=cf_table,location=location,))
        constraints.append(dnachisel.builtin_specifications.AvoidRareCodons(min_frequency=min_freq,species=species,location=location,))
    else:
        objectives.append(dnachisel.builtin_specifications.MaximizeCAI(species=species, boost=cai_weight, location=location))
        
    for seq in avoid: # GG enzyme recognition site.
        constraints.append(dnachisel.builtin_specifications.AvoidPattern(seq, location=location))

    # Organism-specific constraints
    if species == 'e_coli':
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("GGAGG", location=location)) # E. coli Shine-Dalgarno
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("TAAGGAG", location=location)) # strong RBS
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("GCTGGTGG", location=location)) # Chi site in E. coli

        # Alternative start sites: G/A rich 6 nt upstream of ATG or GTG or TTG (cryptic start sites)
        # N = A or T or C or G
        # D = A or G or T
        # R = A or G
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC523762/?page=1
        # This constraint can be too harsh and lead to problems with some sequences.
        constraints_easier = copy.deepcopy(constraints)
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("RRRRRNNNNNDTG", location=location)) # 5 nt spacing
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("RRRRRNNNNNNDTG", location=location)) # 6 nt spacing
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("RRRRRNNNNNNNDTG", location=location)) # 7 nt spacing

    elif species == 'h_sapiens':
        constraints.append(dnachisel.builtin_specifications.AvoidPattern("GCCRCCATGG", location=location)) # Kozak sequence

    # Global GC content (according to TWIST).
    constraints.append(dnachisel.builtin_specifications.EnforceGCContent(mini=0.25, maxi=0.65, location=location))

    # Local GC content (according to TWIST).
    constraints.append(dnachisel.builtin_specifications.EnforceGCContent(mini=0.35, maxi=0.65, window=50, location=location))

    # Start optimisation.
    solutions = []
    solution_found = False
    for i in range(max_tries):

        if solution_found:
            break

        try:
            if species == 'e_coli' and i >= max_tries/2:
                if count == None:
                    print('  [!] Preventing alternative start sites removed from the list of optimisation constraints.')
                else:
                    if 'start_sites' not in warnings:
                        warnings['start_sites'] = []
                    warnings['start_sites'].append(count)
                    #print(f'  [!] Preventing alternative start sites removed from the list of optimisation constraints for {count}th sequence.')
                initial_problem = DnaOptimizationProblem(naive_dna_sequence, constraints=constraints_easier, objectives=objectives, logger=None)

            else:
                initial_problem = DnaOptimizationProblem(naive_dna_sequence, constraints=constraints, objectives=objectives, logger=None)

            problem = copy.deepcopy(initial_problem)
            problem.resolve_constraints_by_random_mutations()
            problem.optimize()
            problem.resolve_constraints(final_check=True)
            solutions.append(problem)
            solution_found = True

        except NoSolutionError:
            initial_problem.max_random_iters += 1000
            solution_found = False

            continue

    if len(solutions) == 0:
        raise NoSolutionError(f"No solution found for {amino_acid_sequence}", initial_problem)

    # Return the best solution according to dnachisel if multiple attemptes were made.
    scores = [solution.objectives_evaluations().scores_sum() for solution in solutions]
    best_idx = np.argmin(scores)

    best_solution = solutions[best_idx]

    return best_solution.sequence, warnings

def rev_translate_many(df, args, cuts):
        '''
        Reverse translate amino acid sequences to DNA sequences and adjust for gene fragment length requirements.
        dna_sequences: list of reverse translated DNA sequences
        dna_fragments: list of DNA sequences with adapters added and adjusted for gene fragment length requirements
        '''
        dna_sequences = []
        avoid_seqs = [cuts[args.enzyme][0], cuts[args.enzyme][1]] + args.avoid
        warnings = {}
        for i, r in df.iterrows():
            amino_acid_sequence = r['aa_sequence']
            dna_sequence,warnings = reverse_translate(
                amino_acid_sequence,
                kmers_weight=1.0,
                cai_weight=1.0,
                hairpins_weight=1.0,
                max_tries=10,
                species=args.species,
                avoid=avoid_seqs,
                count=i,
                warnings=warnings,
                cell_free=args.cell_free
                )
            dna_sequences.append(dna_sequence)
        if 'start_sites' in warnings:
            print(f'  [!] You may have alternative start sites for {len(warnings["start_sites"])} sequences.')
        df['dna_sequence'] = dna_sequences
        return df
# ============================================
# FUNCTIONS
# ============================================
def find_dna_repeats(dna_sequence,frag_size):
    '''
    Check for repeats of the same DNA sequence within a sequence and modify if necessary.
    '''
    frag_dict = {}
    dup_count = 0
    # detect repeats within sequence
    for i in range(0, len(dna_sequence) - frag_size):
        fragment = dna_sequence[i:i+frag_size]
        rv_comp_frag = reverse_complement(fragment)
        if fragment in frag_dict:
            frag_dict[fragment].append(i)
            dup_count += 1
        elif rv_comp_frag in frag_dict:
            frag_dict[rv_comp_frag].append(i)
            dup_count += 1
        else:
            frag_dict[fragment] = [i]
    return frag_dict, dup_count

def get_new_fragment(dna_frag, species, avoid, max_tries = 3):
    count = 0
    solution_found = False
    while solution_found == False and count <= max_tries:
        constraints = []
        objectives = []
        location = Location.from_biopython_location(SeqFeature.FeatureLocation(0, len(dna_frag)))
        for seq in avoid: # GG enzyme recognition site.
            constraints.append(dnachisel.builtin_specifications.AvoidPattern(seq, location=location))
        constraints.append(dnachisel.builtin_specifications.EnforceTranslation(location=location, genetic_table="Standard"))
        if count == 0:
            objectives.append(dnachisel.builtin_specifications.MaximizeCAI( species=species, boost=1.0, location=location))
            constraints.append(dnachisel.builtin_specifications.EnforceGCContent(mini=0.4, maxi=0.65, window=12, location=location))
        else:
            constraints.append(dnachisel.builtin_specifications.EnforceGCContent(mini=0.2, maxi=0.8, window=12, location=location))
        try:
            problem = DnaOptimizationProblem(dna_frag, constraints = constraints, objectives=objectives,logger=None)
            problem.resolve_constraints()
            overlap = problem.sequence
            solution_found = True
        except:
            solution_found = False
            overlap = dna_frag
        count += 1
    return overlap, solution_found

def fix_dna_repeats(df, args, max_tries=3):
    frag_size = args.repeat_frag_size
    species = args.species
    avoid = args.avoid
    dna_sequences_fixed = []
    def get_fragment_from_index(dup_index, dna_sequence):
        start = dup_index - (dup_index % 3) # start of codon before fragment
        end = dup_index + frag_size + (3 - (dup_index + frag_size) % 3) # end of codon after fragment
        temp_frag = dna_sequence[start:end]
        new_frag, solution = get_new_fragment(temp_frag, species, avoid, max_tries=3)    
        dna_sequence = dna_sequence[:start] + new_frag + dna_sequence[end:]
        return dna_sequence, solution
    
    for i, r in df.iterrows():
        dna_sequence = r['dna_sequence']
        aa_sequence = r['aa_sequence']
        for j in range(max_tries):
            frag_dict, dup_count = find_dna_repeats(dna_sequence,frag_size)
            if dup_count > 0:
                print(f"{dup_count} duplicates found in dna sequence for protein {i}, attempt to fix number {j}")
            for frag in frag_dict:
                avoid.append(frag)
                avoid.append(reverse_complement(frag))
            for frag, loc in frag_dict.items():
                if len(loc) > 1:
                    #need to be careful here to not change aa sequence, so we need it to be a synonymous change (multiples of 3 from the start of the dna sequence)
                    duplicates = loc[1:] # skip the first occurrence
                    solution_found = True
                    for dup_index in duplicates:
                        dna_sequence, solution = get_fragment_from_index(dup_index, dna_sequence)
                        if solution == False:
                            solution_found = False
                    if solution_found == False: # if any of the others couldn't be fixed, try to fix the first occurrence
                        dna_sequence, solution = get_fragment_from_index(loc[0], dna_sequence)
            assert Seq.translate(dna_sequence) == aa_sequence, f"Internal error: DNA sequence does not match amino acid sequence after fixing repeats. Please contact @sruge with the following information:\nAA sequence: {aa_sequence}\nDNA sequence: {dna_sequence}"
        dna_sequences_fixed.append(dna_sequence)
    df['dna_sequence'] = dna_sequences_fixed
    return df

def adjust_for_fragments(df, args, cuts):
    '''
    Check size of DNA sequences pad if necessary to meet minimum length requirements for ordering gene fragments from Twist Bioscience.
    '''
    max_length = args.max_length
    avoid_seqs = args.avoid
    gg_int_adapters = get_binding_adapters(args, cuts)
    dna_fragments = []
    dna_frag_len = []
    for i,r in df.iterrows():
        dna_seq = gg_int_adapters[0] + r.dna_sequence + gg_int_adapters[1]
        # Check if the sequence fits the the size limits
        if len(dna_seq) < 300:
            pad_length = ((300 - (len(dna_seq))) // 2 )
            extra = 300 - len(dna_seq) - (2 * pad_length)
            if extra < 0 :
                extra = 0
            # Pad sequence
            pad_nocut = False
            while pad_nocut == False:
                pad5prime = ''.join(np.random.choice(['A','T','C','G'], size=pad_length + extra))
                pad3prime = ''.join(np.random.choice(['A','T','C','G'], size=pad_length))
                dna_fragment =  pad5prime + dna_seq + pad3prime
                pad_nocut = True

                for a_seq in avoid_seqs:
                    # if any cut sequence is found in the padding, try again
                    if (a_seq in pad5prime) or (a_seq in pad3prime):
                        pad_nocut = False
                    elif a_seq in pad5prime + pad3prime:
                        pad_nocut = False
                    elif dna_fragment.count(a_seq) > 1:
                        pad_nocut = False
        else:
            dna_fragment = dna_seq
        if len(dna_fragment) > max_length:
            sys.exit(f"Internal error: DNA sequence is too long to be ordered as a gene fragment ({len(dna_seq)} vs. {max_length} bp).\n"\
                        "This should have been caught earlier in the script. Contact @sruge that you are seeing this error")
        dna_fragments.append(dna_fragment)
        dna_frag_len.append(len(dna_fragment))
    df['dna_fragments'] = dna_fragments
    df['length_fragments'] = dna_frag_len
    return df

def golden_gate_assembly(df, args, cuts):
    '''
    Perform Golden Gate assembly of the plasmid and insert sequences.
    Check for mismatched overhangs, out of frame issues, and extra cut sites.
    Output the full cloned assembly, DNA ORF, and AA ORF.
    '''
    # Enzyme-specific cut characteristics.
    fw, rv, n_spacer, n_sticky = cuts[args.enzyme]
    vector_fasta = SeqIO.read(args.gg_vector, 'fasta')
    entry_vector= str(vector_fasta.seq).lower()
    df['plasmid'] = vector_fasta.id

    # GG cloning.
    vector_5prime = entry_vector.find(fw.lower()) + len(fw) + n_spacer
    vector_3prime = entry_vector.find(rv.lower()) - n_spacer
    
    gg_dict = {'well_position':[], 'ORF':[], 'exp_aa_seq':[], 'cloned_plasmid_seq':[]}
    for _, r in df.iterrows():
        aa_seq = r['aa_sequence']
        insert = r['dna_fragments'].upper()

        # Find eBlock section with cut sites facing in the correct directions.
        # This should not have cut sites in the padding regions, (Necessary for cases where cut sites are accidentlly also present in the padding regions.)
        fw_locations = np.array([x.span()[0] for x in re.finditer(fw.upper(), insert.upper())])
        rv_locations = np.array([x.span()[0] for x in re.finditer(rv.upper(), insert.upper())])

        fw_idx = []
        rv_idx = []
        delta_bp = []
        for i, f_loc in enumerate(fw_locations):
            for j, rv_loc in enumerate(rv_locations):
                delta_bp.append(rv_loc - f_loc)
                fw_idx.append(i)
                rv_idx.append(j)

        delta_bp = np.array(delta_bp)

        correct_idx = np.argwhere(delta_bp==delta_bp[delta_bp>=3*len(aa_seq)].min())[0][0]

        insert_5prime = fw_locations[fw_idx[correct_idx]] \
                        + len(fw) \
                        + n_spacer \
                        + n_sticky

        insert_3prime = rv_locations[rv_idx[correct_idx]] \
                        - n_spacer \
                        - n_sticky

        assembled_plasmid = entry_vector[:vector_3prime] \
                                + insert[insert_5prime:insert_3prime] \
                                + entry_vector[vector_5prime:]
        gg_dict['cloned_plasmid_seq'].append(assembled_plasmid.lower())
        gg_dict['well_position'].append(r['well_position'])

        # Identify the ORF that contains the insert.
        # Search for the shortest START-STOP span that contains the insert sequence.
        plasmid_seq = assembled_plasmid.lower()
        starts = np.array([s.start() for s in re.finditer('atg', plasmid_seq)])
        ends = np.array(sorted([e.end() for e in re.finditer('tag', plasmid_seq)] 
                                + [e.end() for e in re.finditer('taa', plasmid_seq)] 
                                + [e.end() for e in re.finditer('tga', plasmid_seq)]))
        current_stop = 0
        ORF = None
        circular = False
        for s in starts:
            inframe_stops = ends[np.logical_and(ends>s, (ends-s)%3==0)]
            if len(inframe_stops) == 0:
                #if no in-frame stops are found, check if the insert is at the end of the plasmid and if there is an in-frame stop before the start codon
                #still need to check that this math is correct
                start = len(plasmid_seq) - s 
                inframe_stops = ends[np.logical_and(ends<s, (ends+start)%3==0)]
                circular = True
            if len(inframe_stops) > 0:
                current_stop = inframe_stops[0]
                if circular == True:
                    coding_seq = plasmid_seq[s:] + plasmid_seq[:current_stop]
                else:
                    coding_seq = plasmid_seq[s:current_stop]
                circular = False
                if insert[insert_5prime:insert_3prime].lower() in coding_seq:
                    ORF = coding_seq
                    exp_product = str(Seq.Seq(ORF).translate())
                elif ORF == None and coding_seq in insert[insert_5prime:insert_3prime].lower():
                    ORF = coding_seq
                    exp_product = str(Seq.Seq(ORF).translate())
        if ORF == None:
            print("ERROR: CANNOT FIND OPEN READING FRAME WHEN TRYING GOLDEN GATE ASSEMBLY")
            print("No open reading frame found, do the cutsites for the vector map (fasta file) and the overhangs listed in the vector map match?")
            sys.exit(1)
        gg_dict['ORF'].append(ORF)
        gg_dict['exp_aa_seq'].append(exp_product)
    gg_df = pd.DataFrame(gg_dict)
    merged_df = pd.merge(df, gg_df, on='well_position')
    return merged_df

def check_for_a280_and_coomassie(df):
    ext_coefficient = []
    mass = []
    for _,r in df.iterrows():
        seq = r['exp_aa_seq']
        seq = seq.replace("*","")
        protparam = ProteinAnalysis(seq)
        mol_weight = protparam.molecular_weight()
        mass.append(int(np.round(mol_weight)))
        ext_coefficient.append(protparam.molar_extinction_coefficient()[1])
        base_aa = r['aa_sequence']
        if 'W' not in seq and 'Y' not in seq:
            if 'R' not in seq and 'K' not in seq and 'H' not in seq and 'P' not in seq:
                print(f'  [!] Expressed sequence {r["design_name"]} does not contain tryptophan, tyrosine, arginine, lysine, histidine, or proline.\n'\
                        ' This sequence will not be detectable by A280 or Coomassie staining. Be aware of this in your assays!\n'\
                        f'{seq}')
            else:
                print(f'  [!] Expressed sequence {r["design_name"]} does not contain tryptophan or tyrosine.\n' \
                        ' This sequence will not be detectable by A280. Be aware of this in your assays!\n'\
                        f'{seq}')
        elif 'W' not in base_aa and 'Y' not in base_aa:
            if 'W' not in base_aa and 'Y' not in base_aa and 'R' not in base_aa and 'K' not in base_aa and 'H' not in base_aa and 'P' not in base_aa:
                print(f'  [!] Base design {r["design_name"]} does not contain tryptophan, tyrosine, arginine, lysine, histidine, or proline.\n'\
                    ' If you clone into a different plasmid you may not be able to see it by A280, make sure to check!.\n'\
                    f'{base_aa}')
            else:
                print(f'  [!] Base design {r["design_name"]} does not contain tryptophan or tyrosine.\n' \
                        ' If you clone into a different plasmid you may not be able to see it by A280 or Coomassie, make sure to check!.\n'\
                        f'{base_aa}')
    df['MW'] = mass
    df["ext_coefficient"] = ext_coefficient
    return

def output(df, filename):
    '''
    outputs a CSV with all information for user. columns are:
        'design_name','aa_sequence','readin_order',
        'dna_sequence','dna_fragments','length_fragments',
        'plasmid','cloned_plasmid_seq',
        'ORF','exp_aa_seq', 
    '''
    print(f'\n\nOutputting CSV and FASTA files for {len(df)} designs...')
    print("Fasta file can be inputted into Twist for ordering gene fragments. CSV file contains all information for user.")
    print("Keep the CSV, when Stacey makes more scripts that's generally one of the input files. It may also be used in the DBL database")
    #output CSV with all information for user
    mini_df = df.copy()
    mini_df = mini_df.drop(columns=['cloned_plasmid_seq'])
    with open(f'{filename}.csv', 'w') as f:
        mini_df.to_csv(f, index=False)
    #output fasta with all fragments for twist
    with open(f'{filename}.fasta', 'w') as f:
        for _, r in df.iterrows():
            f.write(f'>{r["design_name"]}\n{r["dna_fragments"]}\n')
    return

# ============================================
# HARD-CODED DEFINITIONS
# ============================================
# FW, RV, N_spacer, N_sticky -- for generating the plasmid maps
cuts = {
    'BsaI':['GGTCTC', 'GAGACC', 1, 4],
    'SapI':['GCTCTTC', 'GAAGAGC', 1, 3]
}
# For 'no_layout' option (just fill plate from A1->H12)
w96 = [str(p) + '_' + r + str(c) for p in range(1,10) for r in 'ABCDEFGH' for c in range(1,13)]
w96_single = [r + str(c) for r in 'ABCDEFGH' for c in range(1,13)]

# For 384w formatting (no layout)
w384 = [str(p) + '_' + r + str(c) for p in range(1,10) for r in 'ABCDEFGHIJKLMNOP' for c in range(1,25)]
w384_df = pd.DataFrame([r + str(c) for r in 'ABCDEFGHIJKLMNOP' for c in range(1,25)], columns=['Well Position'])
# this gives you well positions for every other well in a 384w plate (A1, A3, A5, etc.)
# for use with a multichannel pipette
w384_manual_space = []
for p in range(1, 10):
    for r in 'ABCDEFGHIJKLMNOP':
        for c in range(1,25,2):
                w384_manual_space.append(str(p) + '_' + r + str(c))
    for r in 'ABCDEFGHIJKLMNOP':
        for c in range(2,26,2):
                w384_manual_space.append(str(p) + '_' + r + str(c))

#=====================
# Main
#=====================
if __name__ == '__main__':
    user=getpass.getuser()
    args = get_arguments()
    # read in sequences
    input_df, filename = read_input_sequences(args)
    filename = f'{date}_{user}_{filename}_{args.species}_{args.enzyme}'
    #check maximum length of sequences
    check_aa_sequences(input_df, args, cuts)
    # Reverse translate sequences and add adapters
    print("Starting reverse translation and fragment generation...")
    rev_translated_df = rev_translate_many(input_df, args, cuts)
    if not args.do_not_fix_repeats: #it's a double negative, but I think I want it on by default
        print("Fixing direct DNA repeats. If this is taking too long you can remove this with --do_not_fix_repeats.\n"
              "Alternatively you can kill the script and try to see why your proteins are making so many 12+ bp repeats\n" 
              "Long times are likely due to either repeat amino acid sequences or having stretches of high or low gc content amino acids")
        rev_translate_df = fix_dna_repeats(rev_translated_df, args, max_tries=5)
    print("Adding adapters")
    adapters_df = adjust_for_fragments(rev_translated_df, args, cuts)
    gg_df = golden_gate_assembly(adapters_df, args, cuts)
    check_for_a280_and_coomassie(gg_df)
    output(gg_df, filename)