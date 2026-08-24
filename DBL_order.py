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
import os
import re
import argparse
import datetime; date = datetime.datetime.now().strftime('%Y_%m_%d')
import numpy as np
import pandas as pd
from Bio import SeqIO, PDB, SeqUtils, Seq, SeqFeature
import json

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
# ============================================
# ARGUMENTS
# ============================================
def get_arguments(argv=None):
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter,
            description=" * Generates an IDT-ready .xlsx file for ordering eBlocks from a folder of PDBs and/or a concatenated FASTA file.\n"
                        " * Appropriate overhangs for Golden Gate cloning into entry vector(s) of interest are added automatically.\n"
                        " * Reverse translation is performed with Ryan's Domesticator.\n"
                        " * RECOMMENDED: check your GG assemblies at https://goldengate.neb.com/#!/\n"
                        " * Wondering why the script is called John Bercow? https://www.youtube.com/watch?v=VYycQTm2HrM&ab_channel=TheSun\n"
                        "\n"
                        " * AVAILABLE ENTRY VECTORS:\n"
                        " *** see /net/software/lab/johnbercow/entry_vectors/ for the FULL list ***\n"
                    f"{vec_str}\n"
            )
    # REQUIRED
    parser.add_argument("input", help='file containing name and amino acid sequence can be multiple file types:\n'
                            'fasta - Fasta format and filename must end in .fa or .fasta\n'
                            'seq   - sequences in format $Sequence $Name. Filename must end is .txt, .seq, or .tab\n'
                            'pdb   - path to a folder of pdbs. Will pull all pdbs from that folder\n'
                            '* NOTE: SCRIPT TREATS ALL SEQUENCES IN ONE FILE AS A SINGLE LIBRARY, IF YOU HAVE MULTIPLE TARGETS/LIBRARIES RUN THEM AS SEPARATE FILES',type=str)
    parser.add_argument(
            '-g,--gg_vector',
            help='Fasta file of plasmid for Golden Gate cloning (determines the DNA adapters). Also determines the AA tags appended to the design in the FASTA output.',
            action='store',
            required=True,
            type=str,
            )

    # OPTIONAL
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
            '--design_prefix',
            help='designs get IDs with this prefix (e.g. LM0001, LM0002, etc...)',
            action='store',
            type=str,
            )
    parser.add_argument(
            '--design_id',
            help='increment design indices from this number.',
            action='store',
            type=int,
            )
    parser.add_argument(
            '--design_prefix',
            help='designs get IDs with this prefix (e.g. LM0001, LM0002, etc...)',
            action='store',
            type=str,
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
            '--no_layout',
            default=True,
            help="do not apply automated layout formatting.",
            action='store_true',
            )
    parser.add_argument(
            '--no_plasmids',
            help="do not generate the cloned plasmid maps.",
            action='store_true',
            )
    parser.add_argument(
            '--verbose',
            help="increase the verbosity of th e output (recommended).",
            action='store_true',
            )
    parser.add_argument(
            '--no_adapters',
            help="adds cut site and sticky ends, but no additional adaptor sequence",
            action='store_true',
            )

    args = parser.parse_args(argv)
    return args

def read_input_sequences(args):
    seq_dict = {
    'design_name':[],
    'aa_sequence':[],
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
            seq_dict['design_name'] = str(fasta.id)
            seq_dict['aa_sequence'] = str(fasta.seq)
            seq_dict['readin_order'] = i
    elif filetype == 'seq':
        with open(args.input, 'r') as seqfile:
            for i, line in enumerate(seqfile):
                seq, name = line.split()
                seq_dict['design_name'] = str(name)
                seq_dict['aa_sequence'] = str(seq)
                seq_dict['readin_order'] = i
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
                    seq_dict['design_name'] = str(f'{pdb_name}_{i+1}')
                    seq_dict['aa_sequence'] = str(seq)
                    seq_dict['readin_order'] = f"{i}_{j}"
                else:
                    seq_dict['design_name'] = str(f'{pdb_name}')
                    seq_dict['aa_sequence'] = str(seq)
                    seq_dict['readin_order'] = i
    df = pd.DataFrame(seq_dict)
    return df, filename

def get_binding_adapters(args, cuts):
    hard_coded_binding_adapters = ['GTTTAAAGGTCTCGGCCGT','GGAGGTGAGACCAAAGGA']# includes cut site spacers and sticky ends for bsaI. this is the normal addition for DBE
    outside_flank = ['GTTTAAA','AAAGGA']
    fw_cut, rv_cut, spacer_len, _ = cuts[args.enzyme]
    spacer = 'gtactacgtaatgt'
    binding_adapters = [fw_cut + spacer[:spacer_len] + args.n_overhang + "T", "G" + args.c_overhang + spacer[-1*spacer_len:] + rv_cut]
    if args.no_adapters == False:
        binding_adapters = [outside_flank[0] + binding_adapters[0],binding_adapters[1] +  outside_flank[1]]
    if args.enzyme == 'BsaI' and args.no_adapters == False:
        assert binding_adapters == hard_coded_binding_adapters, f'Binding adapters are not correct for {args.enzyme}. Please check the hard-coded binding adapters in the script.'
    return binding_adapters

def check_max_len(df, args, binding_adapters):
    max_len_dna = args.max_length - len(binding_adapters[0]) - len(binding_adapters[1])
    max_length_aa = int(np.floor(max_len_dna/3)) # max len minus adapters, divided by 3 for aa
    for seq in list(df['aa_sequence']):
        if len(seq) > max_length_aa:
            print(f'  [!] Sequence {seq} is too long for the specified maximum length ({len(seq)} vs. {max_length_aa} aa). Base dna length is 5000 bp max. Did you make a mistake in the input file?')
            sys.exit("  [!] Sequence too long for twist synthesis. System exiting...")
    return

def check_for_a280_and_coomassie(df):
    for _,r in df.iterrows():
        seq = r['aa_sequence']
        if 'W' not in seq and 'Y' not in seq:
            print(f'  [!] Design {r["design_name"]} does not contain tryptophan or tyrosine.\n' \
                ' This sequence will not be detectable by A280. Please check your input file.\n'\
                f'{seq}')
        if 'W' not in seq and 'Y' not in seq and 'R' not in seq and 'K' not in seq and 'H' not in seq and 'P' not in seq:
            print(f'  [!] Design {r["design_name"]} does not contain tryptophan, tyrosine, arginine, lysine, histidine, or proline.\n'\
                ' This sequence will not be detectable by Coomassie staining. Please check your input file.\n'\
                f'{seq}')
    return
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

def reverse_translate(
        amino_acid_sequence,
        kmers_weight=1.0,
        cai_weight=1.0,
        hairpins_weight=1.0,
        max_tries=10,
        species='e_coli',
        avoid=['GGTCTC', 'GAGACC']
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
    objectives.append(dnachisel.builtin_specifications.MaximizeCAI(species=species, boost=cai_weight, location=location))

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
                print('  [!] Preventing alternative start sites removed from the list of optimisation constraints.')
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

    return best_solution.sequence

def find_repeats(
        sequence,
        exact_repeat_min_len=20,      # Twist hard rule
        tm_repeat_min_len=11,         # shortest repeat we bother Tm-checking
        tm_threshold=60.0,            # Twist hard rule
        near_repeat_min_len=20,
        max_mismatches=2,
    ):
    """
    repeat finder created by claude, not tested yet

    Scan `sequence` (already-optimized DNA, e.g. reverse_translate() output)
    for same-strand repeats and near-repeats that could cause synthesis or
    PCR-based assembly problems.
 
    Returns a list of hit dicts sorted worst-first (exact_repeat > high_tm_repeat
    > near_repeat, longest first within each), each with:
        start1, start2 : 0-based positions of the two copies
        length          : length of the matched region (bp)
        mismatches      : 0 for exact/Tm hits, >0 for near-repeats
        tm              : predicted duplex Tm of the matched region (None if
                           not computed, e.g. for near-repeats)
        reason          : 'exact_repeat' | 'high_tm_repeat' | 'near_repeat'
        sequence        : the matched region itself (from copy 1)
    """
    sequence = sequence.upper()
    n = len(sequence)
    hits = []
 
    # --- exact repeats (>= exact_repeat_min_len) and short high-Tm repeats ---
    seed_len = min(exact_repeat_min_len, tm_repeat_min_len)
    seed_positions = _kmer_positions(sequence, seed_len)
 
    for positions in seed_positions.values():
        if len(positions) < 2:
            continue
        for i, j in combinations(positions, 2):
            if not _is_left_maximal(sequence, i, j):
                continue  # this pair is a sub-match of an earlier, longer seed
            length = _extend_exact(sequence, i, j)
            if length < seed_len:
                continue
            a, b = min(i, j), max(i, j)
            if length >= exact_repeat_min_len:
                hits.append({
                    'start1': a, 'start2': b, 'length': length, 'mismatches': 0,
                    'tm': None, 'reason': 'exact_repeat',
                    'sequence': sequence[a:a + length],
                })
            elif length >= tm_repeat_min_len:
                segment = sequence[a:a + length]
                tm = mt.Tm_NN(segment)
                if tm >= tm_threshold:
                    hits.append({
                        'start1': a, 'start2': b, 'length': length, 'mismatches': 0,
                        'tm': tm, 'reason': 'high_tm_repeat', 'sequence': segment,
                    })
 
    exact_hits = [h for h in hits if h['reason'] in ('exact_repeat', 'high_tm_repeat')]
 
    # --- near-repeats (mismatches tolerated) ---
    near_hits = []
    for i, j in _near_repeat_anchors(sequence, near_repeat_min_len, max_mismatches):
        a, b = min(i, j), max(i, j)
        length, mismatches = _extend_with_mismatches(
            sequence, a, b, max_mismatches
        )
        if length < near_repeat_min_len or mismatches == 0:
            continue
        if _contained_in_any(exact_hits + near_hits, a, b, length):
            continue
        near_hits.append({
            'start1': a, 'start2': b, 'length': length, 'mismatches': mismatches,
            'tm': None, 'reason': 'near_repeat', 'sequence': sequence[a:a + length],
        })
 
    hits.extend(near_hits)
 
    reason_rank = {'exact_repeat': 0, 'high_tm_repeat': 1, 'near_repeat': 2}
    hits.sort(key=lambda h: (reason_rank[h['reason']], -h['length']))
    return hits


# ============================================
# FUNCTIONS
# ============================================
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

def adjust_for_fragments(df, max_length, gg_int_adapters, avoid_seqs):
    '''
    Check size of DNA sequences pad if necessary to meet minimum length requirements for ordering gene fragments from Twist Bioscience.
    '''
    dna_fragments = []
    for _, r in df.iterrows():
        dna_seq = gg_int_adapters[0] + r['dna_sequence'] + gg_int_adapters[1]
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

                pad_nocut = True
                for a_seq in avoid_seqs:
                    # if any cut sequence is found in the padding, try again
                    if (a_seq in pad5prime) or (a_seq in pad3prime):
                        pad_nocut = False
                    elif a_seq in pad5prime + pad3prime:
                        pad_nocut = False
                        

            dna_seq =  pad5prime + dna_seq + pad3prime
        if len(dna_seq) > max_length:
            sys.exit(f"Internal error: DNA sequence is too long to be ordered as a gene fragment ({len(dna_seq)} vs. {max_length} bp).\n"\
                     "This should have been caught earlier in the script. Contact @sruge that you are seeing this error")
        for avoid in avoid_seqs:
            if dna_seq.count(avoid) > 1:        
                sys.exit(f"Internal error: {avoid} shows up in DNA sequence more than once.\n"\
                        "This should not be able to occur. Contact @sruge that you are seeing this error, and tell them that their internal logic is wrong")
        dna_fragments.append(dna_seq)

    return dna_fragments

def rev_translate_and_make_fragments(df, args, cuts):
        '''
        Reverse translate amino acid sequences to DNA sequences and adjust for gene fragment length requirements.
        dna_sequences: list of reverse translated DNA sequences
        dna_fragments: list of DNA sequences with adapters added and adjusted for gene fragment length requirements
        '''
        dna_sequences = []
        dna_fragments = []
        dna_frag_len = []
        gg_adapters = get_binding_adapters(args, cuts)
        avoid_seqs = [cuts[args.enzyme][0], cuts[args.enzyme][1]] + args.avoid
        for _, r in df.iterrows():
            amino_acid_sequence = r['aa_sequence']
            dna_sequence = reverse_translate(
                amino_acid_sequence,
                kmers_weight=1.0,
                cai_weight=1.0,
                hairpins_weight=1.0,
                max_tries=10,
                species=args.species,
                avoid=avoid_seqs
                )
            dna_sequences.append(dna_sequence)
            dna_fragment = adjust_for_fragments(dna_sequence, args.max_length, gg_adapters, avoid_seqs)
            dna_fragments.append(dna_fragment)
            dna_frag_len = len(dna_fragment)
        df['dna_sequence'] = dna_sequences
        df['dna_fragments'] = dna_fragments
        df['length_fragments'] = dna_frag_len
        return df

def output(df, filename):
    '''
    outputs a CSV with all information for user. columns are:
        'design_name','aa_sequence','readin_order',
        'dna_sequence','dna_fragments','length_fragments',
        'plasmid','cloned_plasmid_seq',
        'ORF','exp_aa_seq', 
    '''
    #output CSV with all information for user
    with open(f'{filename}.csv', 'w') as f:
        df.to_csv(f, index=False)
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
    user=os.getlogin()
    args = get_arguments()
   
    # read in sequences
    input_df, filename = read_input_sequences(args)
    filename = f'{date}_{user}_{filename}_{args.species}_{enzyme}'
    #check maximum length of sequences
    check_max_len(input_df, args, binding_adapters)
    check_for_a280_and_coomassie(input_df)
    # Reverse translate sequences and add adapters
    rev_translated_df = rev_translate_and_make_fragments(input_df, args, cuts)

def golden_gate_assembly(rev_translated_df, args, cuts):
    '''
    Perform Golden Gate assembly of the plasmid and insert sequences.
    Check for mismatched overhangs, out of frame issues, and extra cut sites.
    Output the full cloned assembly, DNA ORF, and AA ORF.
    '''
    # Enzyme-specific cut characteristics.
    fw, rv, n_spacer, n_sticky = cuts[args.enzyme]
    vector_fastas = SeqIO.read(args.gg_vector, 'fasta')
    if len(vector_fastas) > 1:
        print(f'Warning: {args.gg_vector} contains multiple sequences. Taking first one, but this may not be the desired sequence.')
        vector_fastas = vector_fastas[0]
    entry_vector= vector_fastas.seq.lower()
    vector_name = vector_fastas.id
    # GG cloning.
    vector_5prime = entry_vector.find(fw) + len(fw) + n_spacer
    vector_3prime = entry_vector.find(rv) - n_spacer

    for _, r in rev_translated_df.iterrows():
        gg_v = args.gg_vector.split('/')[-1].split('.')[0]
        aa_seq = r['aa_sequence']
        insert = r['dna_fragments'].lower()

        # Find eBlock section with cut sites facing in the correct directions.
        # This should not have cut sites in the padding regions, (Necessary for cases where cut sites are accidentlly also present in the padding regions.)
        fw_locations = np.array([x.span()[0] for x in re.finditer(fw, insert)])
        rv_locations = np.array([x.span()[0] for x in re.finditer(rv, insert)])

        fw_idx = []
        rv_idx = []
        delta_bp = []
        for i, f in enumerate(fw_locations):
            for j, r in enumerate(rv_locations):
                delta_bp.append(r - f)
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
        
        eblocks[f'{gg_v}_cloned_plasmid_seq'].append(assembled_plasmid.lower())
        
        
        # Identify the ORF that contains the insert.
        # Search for the shortest START-STOP span that contains the insert sequence.
        plasmid_seq = assembled_plasmid.lower()
        starts = np.array([s.start() for s in re.finditer('atg', plasmid_seq)])
        ends = np.array(sorted([e.end() for e in re.finditer('tag', plasmid_seq)] 
                                + [e.end() for e in re.finditer('taa', plasmid_seq)] 
                                + [e.end() for e in re.finditer('tga', plasmid_seq)]))
        current_stop = 0
        ORF = None
        for s in starts:
            inframe_stops = ends[np.logical_and(ends>s, (ends-s)%3==0)]
            if len(inframe_stops) == 0:
                #if no in-frame stops are found, check if the insert is at the end of the plasmid and if there is an in-frame stop before the start codon
                #still need to check that this math is correct
                start = len(plasmid_seq) - s 
                inframe_stops = ends[np.logical_and(ends<s, (start+ends)%3==0)]
            if len(inframe_stops) > 0:
                if s > current_stop:
                    current_stop = inframe_stops[0]
                    coding_seq = plasmid_seq[s:current_stop]

                    if insert[insert_5prime:insert_3prime] in coding_seq:
                        ORF = coding_seq
                        exp_product = str(Seq.Seq(ORF).translate())
                    elif ORF == None and coding_seq in insert[insert_5prime:insert_3prime]:
                        ORF = coding_seq
                        exp_product = str(Seq.Seq(ORF).translate())
        if ORF == None:
            print("ERROR: CANNOT FIND OPEN READING FRAME WHEN TRYING GOLDEN GATE ASSEMBLY")
            print("No open reading frame found, do the cutsites for the vector map (fasta file) and the overhangs listed in the vector map match?")
            print("Make sure the cloning region does not continue back to the top of the file. We don't parse the files as circular maps")
            sys.exit(1)

        eblocks[f'ORF_from_{gg_v}'].append(ORF)
        eblocks[f'exp_aa_seq_from_{gg_v}'].append(exp_product)

    return


def cut_fragment(seq, enzyme_site, overhang_length, downstream_offset):
    """
    Cut a DNA sequence at enzyme sites and return (left_overhang, insert, right_overhang).
    """
    rev_site = reverse_complement(enzyme_site)

    # Left site
    left_idx = seq.find(enzyme_site)
    if left_idx == -1:
        raise ValueError(f"Left enzyme site not found in sequence: {seq}")
    left_overhang_start = left_idx + len(enzyme_site) + downstream_offset
    left_overhang = seq[left_overhang_start:left_overhang_start + overhang_length]

    # Right site (reverse complement)
    right_idx = seq.find(rev_site)
    if right_idx == -1:
        raise ValueError(f"Right enzyme site not found in sequence: {seq}")
    right_overhang_end = right_idx - downstream_offset
    right_overhang = seq[right_overhang_end - overhang_length:right_overhang_end]
    
    # Internal insert (between overhangs)
    insert_start = left_overhang_start + overhang_length
    insert_end = right_overhang_end - overhang_length
    insert = seq[insert_start:insert_end]

    return left_overhang, insert, right_overhang

def golden_gate_assemble(backbone, fragments, enzyme_site, overhang_length, downstream_offset):
    """
    Assemble multiple DNA fragments into a backbone using Golden Gate assembly.
    Fragments can be in any order; they will be sorted automatically by overhangs.
    """
    #everything to upper case
    backbone = backbone.upper()
    fragments = [frag.upper() for frag in fragments]
    enzyme_site = enzyme_site.upper()
    
    # Cut backbone and get its overhangs
    bb_left_oh, _, bb_right_oh = cut_fragment(backbone, enzyme_site, overhang_length, downstream_offset)

    # Cut all fragments
    frag_parts = []
    for i, frag in enumerate(fragments):
        left_oh, insert, right_oh = cut_fragment(frag, enzyme_site, overhang_length, downstream_offset)
        frag_parts.append({
            "index": i,
            "left_oh": left_oh,
            "insert": insert,
            "right_oh": right_oh
        })

    # Build lookup table for fragment left overhangs
    left_lookup = {f["left_oh"]: f for f in frag_parts}

    # Start with the fragment that matches backbone's right overhang
    if bb_right_oh not in left_lookup:
        raise ValueError(f"No fragment matches backbone right overhang: {bb_right_oh}")

    ordered = []
    current_oh = bb_right_oh
    while current_oh in left_lookup:
        frag = left_lookup.pop(current_oh)
        ordered.append(frag)
        current_oh = frag["right_oh"]

    # Final check: last fragment right overhang must match backbone left
    if ordered[-1]["right_oh"] != bb_left_oh:
        raise ValueError(
            f"Assembly overhangs do not close properly. "
            f"Expected {bb_left_oh}, got {ordered[-1]['right_oh']}"
        )

    # Final assembly: backbone prefix + inserts + backbone suffix
    # Prefix = sequence before left site; suffix = after right site
    left_site_index = backbone.find(reverse_complement(enzyme_site))
    right_site_index = backbone.find(enzyme_site)
    offset = downstream_offset
    left_backbone_seq = backbone[:left_site_index - offset]
    right_backbone_seq = backbone[right_site_index + len(reverse_complement(enzyme_site)) + offset + len(bb_right_oh):]
    
    assembled = (
        left_backbone_seq +
        "".join(f["insert"] + f["right_oh"] for f in ordered) +
        right_backbone_seq
    )

    return assembled

