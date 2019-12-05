'''
Phylogenetic Tree Generation using Branch & Bound
Can process directories of fasta MSA files.
'''

import os
import copy
import random

from Bio import SeqIO
from Bio import Phylo
from Bio import AlignIO
from Bio.Phylo.BaseTree import Clade
from Bio.Phylo import TreeConstruction
from Bio.SubsMat import MatrixInfo

from typing import *

def gen_candidate_tree_helper(new_taxon_id: str, cl: Clade):
    '''
    Returns potential "replacements" for cl

    '''
    return_vals = []

    # case of branching at this clade
    return_vals.append(Clade(clades=[cl, Clade(name=new_taxon_id, branch_length=1)], branch_length=1))

    for term in return_vals[0].get_terminals():
        if term.branch_length != 1:
            print("Wrong at initial branch off, " + str(term.branch_length))
    for i in range(len(cl.clades)):
        for n_clade in gen_candidate_tree_helper(new_taxon_id, cl.clades[i]):
            n_clade_list = cl.clades[0:i] + [n_clade] + cl.clades[i+1:]
            return_vals.append(Clade(clades=n_clade_list, branch_length=1))

    return return_vals

def generate_candidate_trees(t, new_b):
    '''
    Generates candidate tree structures produced by adding the new taxon new_b
    to tree t.
    '''
    new_trees = []
    cl = t.root
    cl.branch_length = 1
    for new_clade in gen_candidate_tree_helper(new_b, cl):
        new_t = copy.copy(t)
        new_t.root = new_clade
        new_trees.append(new_t)
    return new_trees

def get_seq_list(aln, stems=False):
    '''
    Writes seed tree for alignment, cuts off names at five characters so they are
    identical between alignments (for comparability/evaluation of trees across
    alignments), and returns a sorted lsit of sequence names for expansion.
    '''
    names = []
    for rec in aln:
        if stems:
            rec.name = rec.id = rec.name[:5]
        names.append(rec.name)
    names.sort()
    with open('data' + os.path.sep + 'known.tre', 'w') as base_tree:
        base_tree.write("("+names[0]+":1, " + names[1] + ":1):1;")

    names.pop(0)
    names.pop(0)
    return names

def backup(trees, backup_dir='backup'):
    '''
    Backs up trees in trees to directory backup_dir
    '''
    if not os.path.isdir(backup_dir):
        print("Unable to find path " + str(backup_dir) + ". Attempting to create...")
        status = os.system('mkdir ' + backup_dir)
        if status != 0 or not os.path.isdir(backup_dir):
            print("Unable to create directory. Backup failed!")
            return
        else:
            print("Directory created successfully.")

    for i,tr in enumerate(trees):
        try:
            Phylo.write(tr, backup_dir + os.path.sep + str(i) + '.tre', 'newick')
        except Exception as e:
            print(f"Unable to back up tree with index {i}!")
            print(e)

def find_good_tree(trees, data_dir, lim=25):
    '''
    Evaluates a list of trees against multiple alignments, returning the best
    scoring tree across all alignments. Requires tree terminal names and
    alignment names in data_dir are the same.
    Args:
        trees: List[Tree], a list of BioPython Phylo trees to evaluate
        data_dir: str, the path to the data directory
        lim: int, number of alignments to evaluate against (more is slower)
    '''
    best_trees = []
    scorer = TreeConstruction.ParsimonyScorer()

    files = os.listdir(data_dir)
    alns = []
    for i,file in enumerate(files[:lim]):
        aln = AlignIO.read(open(data_dir + os.path.sep + file), 'fasta')
        for rec in aln:
            rec.name = rec.id = rec.name[:5]
        alns.append(aln)
    scores = {}
    max_aln = {}

    # Score trees and track highest score for each alignment to normalize later
    print(f"Processing {len(trees)} trees...")
    for i,tree in enumerate(trees):
        print(f"\t{i+1}/{len(trees)}")
        for j,aln in enumerate(alns[:lim]):
            try:
                sco = scorer.get_score(copy.deepcopy(tree), aln)
            except Exception as e:
                print(e)
                print("Scoring failed. Did you ensure that terminal names and alignment names match?")
                return None

            scores[(i,j)] = sco
            if sco > max_aln.get(j, 0): max_aln[j] = sco

    # Computes normalized scores for each tree
    fin_scores = {}
    for i in range(len(trees)):
        for j in range(len(alns[:lim])):
            m_aln = max_aln.get(j, 0)
            if m_aln <= 0: continue
            normalized_score = scores.get((i,j), 0)/m_aln
            if not normalized_score:
                print(f"Error for tree {i} and alignment {j}.")
            fin_scores[i] = fin_scores.get(i,0) + normalized_score

    # Finds final best tree
    best_tree = -1
    for key in fin_scores.keys():
        if best_tree < 0 or fin_scores[key] < fin_scores[best_tree]:
            best_tree = key
    if best_tree < 0:
        print("No best tree found.")
        return None

    return trees[best_tree]

def get_best_trees(aln, scorer, eval_limit, limited_expand):
    '''
    Finds the best trees starting from a seed tree and expanding using seq_list.
    Uses branch and bound with early stopping and (optionally) limited expansion.
    Args:
        aln: BioPython MultipleSequenceAlignment, used to score the trees
        scorer: BioPython ParsimonyScorer, used to score the trees
        eval_limit: int, maximum number of trees to evaluate
        limited_expand: bool, whether or not to use limited expansion
    Returns:
        List[Tree]: A list of best scoring BioPython trees.
    '''
    # Generate list of sequences to iteratively add to the tree
    seq_list = get_seq_list(aln, stems=True)
    n_to_add = len(seq_list)
    # Load in base tree
    tree = Phylo.read("data" + os.path.sep + "known.tre", "newick")

    # Quirks with the BioPython loading function means we need to reset branch lengths
    for term in tree.root.get_terminals():
        term.branch_length=1
    for nterm in tree.root.get_nonterminals():
        nterm.branch_length=1

    incomplete_trees = [(tree,scorer.get_score(copy.deepcopy(tree), aln), 0)]
    best_complete_score = float('inf')
    best_list = []
    ctr = 0
    best_for_size = {}

    # Limited expansion settings
    epsilon = [.5, .5, .25, .1, .1, .05, .01, .01, .01, .01, .01]
    # Ensure we will never go out of bounds
    epsilon.extend([.001]*(len(seq_list) + 2 - len(epsilon)))
    dist = 100

    for n in range(n_to_add + 1): best_for_size[n] = float('inf')
    while incomplete_trees and ctr < eval_limit:
        if ctr%1000 == 0: print("Stage: " + str(ctr))
        ctr += 1
        curr_tree,curr_score,curr_size = incomplete_trees.pop(-1)

        # Prune
        if curr_score > best_complete_score: continue
        if curr_score < best_for_size[curr_size]:
            best_for_size[curr_size] = curr_score
        elif limited_expand and curr_score > max(best_for_size[curr_size]*(1+epsilon), best_for_size[curr_size] + dist):
            continue

        if curr_size == n_to_add: # Complete Tree
            if curr_score < best_complete_score:
                best_complete_score = curr_score
                best_list = [(curr_tree, curr_score)]

            elif curr_score == best_complete_score:
                best_list.append((curr_tree, curr_score))
        else:
            expanded_trees = []
            for next_gen in generate_candidate_trees(curr_tree, seq_list[curr_size]):
                expanded_trees.append((next_gen, scorer.get_score(copy.deepcopy(next_gen), aln), curr_size + 1))

            # We want the best (lowest) scoring trees at the end, so they are
            # expanded first (increases chances of hitting a good alignment quickly)
            expanded_trees.sort(key=lambda a: a[1], reverse=True)
            incomplete_trees.extend(expanded_trees)

    return best_list

def evaluate_directory(data_dir, eval_limit=50000, lim=5, limited_expand=False):
    '''
    Processes a directory of FASTA files, generating and evaluating trees
    using branch and bound with early stopping and limited expansion (optional)
    Args:
        data_dir: str, the path to the data directory which stores fasta files
        eval_limit: int, max number of trees to evaluate for each file. Good
                    setting depends on seq length and time you're willing to wait.
        lim: int, number of files in the data directory to process
        limited_expand: bool, whether or not to use limited expansion
    Returns:
        List[Tree], a list of BioPython Phylo trees with tied best scores
    '''
    scorer = TreeConstruction.ParsimonyScorer()
    all_best = []

    files = os.listdir(data_dir)

    for i,file in enumerate(files[:lim]):
        # Load and sort file
        print(f"Processing {file} ({i+1}/{len(files[:lim])})")
        aln = AlignIO.read(open(data_dir + os.path.sep + file), 'fasta')
        aln.sort(key=lambda a: a.id)

        result_trees = get_best_trees(aln, scorer, eval_limit, limited_expand)
        print(f"Found {len(result_trees)} trees.")
        all_best.extend(result_trees)

    return [tr[0] for tr in all_best]

data = 'marker_genes_aligned'
data_dir = 'data' + os.path.sep + data
tree_dir = 'backup'

if not os.path.isdir(data_dir):
    print("Unable to find data directory" + data_dir + ". Exiting.")
    exit(0)
if not tree_dir or not os.path.isdir(tree_dir):
    trees = evaluate_directory(data_dir, eval_limit=500)
    backup(trees, backup_dir='backup')
else:
    trees = []
    for file in os.listdir(tree_dir):
        trees.append(Phylo.read(tree_dir+os.path.sep+file, 'newick'))

good = find_good_tree(trees, data_dir, lim=15)

if good:
    Phylo.write(good, 'backup' + os.path.sep + data + '.tre', 'newick')
    Phylo.draw_ascii(good)
else:
    print("Could not find any complete tree.")
