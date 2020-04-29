#! /usr/bin/env python

# File: pyvolve_sim_seqs.py
#
# Simulate alignments under a mutation-selection model using Stephanie Spielman's pyvolve package.
# For verifying non-homogeneous mutation-selection inference method.

import sys
import csv
import re
import os
import os.path
import multiprocessing as mp
import itertools as it
import operator as op
from argparse import ArgumentParser
import pyvolve
import rpy2
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from Bio import SeqIO, AlignIO, Phylo
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO

from KimuraMutSelMatrix import *


# !TODO
# Add different breakpoint patterns & fitness distributions
# Random breaks as written / different root freqs


###---------------------------------DATA DIVISION-------------------------------------###

# Parameters for simulating birth-death trees using the apTreeshape package for R
SP_BIRTH = 0.5      # lambda parameter, standard 0.66
SP_DEATH = 0.0       # mu parameter, standard 0.0066
SP_EPSILON = 0.000001 # 'Minimum age of unobserved splits'; just make it small and forget about it
SP_AGERICHNESS = 1.   # alpha parameter; 1 for a standard BDP/coalescent tree, -0.7 imb, +1.9 bal

# Standard genetic codes
MOLECULES = pyvolve.Genetics()

# Constants for drawing simulated fitnesses
AA_DIM = 19             # Number of free fitness parameters (always one left unfree )
FIT_RAD = 5e-5          # Size of interval centred on 1 from which fitnesses are drawn uniformly
FIT_DIST = 0.0002       # Currently unused - for drawing fitness vectors at fixed Euclidean distances

###------------------------------ENVIRONMENT DIVISION---------------------------------###

def make_parser():
    parser = ArgumentParser(description=
                               "Simulate mutation-selection alignments for verification.")
    parser.add_argument("--ntrees", type=int, help=
                        "Number of tree topologies to generate for simulation.", default=1)
    parser.add_argument("--ntaxa", type=int, help=
                        "Number of taxa in each tree.", default=10)
    parser.add_argument("--balance", type=float, help=
                    """Beta parameter for tree simulation. Zero for a standard birth-death tree.
                    Negative values are less balanced and positive values are more balanced.
                    """, default=0.0)
    parser.add_argument("--nsites", type=int, help=
                        "Number of CODON sites to simulate.", default=333)
    parser.add_argument("--nBackgroundCats", type=int, help=
                        "Number of background (negative selection) profiles.", default=1)
    parser.add_argument("--nModelChangeCats", type=int, help=
                        """
                        Number of different site categories allowing a change in background models at
                        each breakpoint.""", default = 0)
    parser.add_argument("--nAdaptiveCats", type=int, help=
                        """Number of different site categories allowing a new model to arise at 
                        each breakpoint.""", default = 0)
    parser.add_argument("--kappa", type=float, help=
                        "Transition-transversion ratio for the HKY mutation model", default=1.0)
    parser.add_argument("--nbreakpts", type=int, help=
                        """Number of fitness breakpoints to be randomly generated.
                        """, default=0)
    parser.add_argument("--root_freqs", help=
                        """Draw root sequences from a Dirichlet distribution over codons; 
                        Default is to use the root model to derive frequencies at the root."
                        """, action="store_true")
    parser.add_argument("--prefix", help=
                        "String to be used in output file names.", default='mutseltest')
    return parser



###-------------------------------PROCEDURE DIVISION----------------------------------###

def r_sim_trees(ntaxa, ntrees, balance=0.):
    """Simulate a list of birth-death trees using the apTreeshape library in R."""
    ape = importr("ape")
    rts = importr("apTreeshape")
    rtrees = []
    for i in range(0, ntrees):
        tre = rts.simulate_yule(SP_EPSILON, SP_AGERICHNESS, balance, ntaxa, SP_BIRTH, SP_DEATH)
        rtrees.append(tre)
    treestrs = [ape.write_tree(r_t)[0] for r_t in rtrees]
    return treestrs
    
def roll_breakpoints(treestring, nbreaks):
    """Assign breakpoints to branches on the tree, weighting different depths equally."""
    handle = StringIO(treestring)
    tree = Phylo.read(handle, "newick")
    ds = tree.depths(unit_branch_lengths=True)
    max_depth = max(ds.values())
    depth_counts = [0]*max_depth
    for d in [v for v in ds.values() if v > 0]:
        depth_counts[d-1] += 1
    depth_probs = [1.0/(float(dc) * float(max_depth)) for dc in depth_counts]
    for (clade, depth) in ds.items():
        clade.comment = str(depth_probs[depth-1])
    probarray = [float(clade.comment) for clade in list(tree.find_clades(order="postorder"))[:-1]]
    breaks = np.random.choice(len(probarray), size=nbreaks, replace=False, p=probarray) 
    return sorted(breaks)
    

def multinomial_cats(ncats, nsites, pvec=None):
    """Assign nsites sites to ncats site categories from a dirichlet-multinomial distribution.
    Categories must contain at least one site.
    """
    if pvec is None:
        pvec = np.random.dirichlet([1.]*ncats, 1)
    cats = np.random.multinomial(nsites - len(pvec), pvec.flat, 1) + np.array([1]*len(pvec))        
    return list(cats.flat)

def gen_model_assignment(nBackgroundCats, nModelChangeCats, nAdaptiveCats, nbreakpts):

    model_count = nBackgroundCats
    bg_cats = list(xrange(0, nBackgroundCats))
    mod_assigns = [[i] * (nbreakpts + 1) for i in bg_cats]
    
    for i in xrange(0, nModelChangeCats):
        mccat = [np.random.choice(bg_cats)]
        for bp in xrange(0, nbreakpts):
            mccat.append(np.random.choice(filter(lambda x: x != mccat[-1], bg_cats)))
        mod_assigns.append(mccat)
    for i in xrange(0, nAdaptiveCats):
        if nBackgroundCats == 0:
            mod_assigns.append(list(xrange(model_count, model_count + nbreakpts + 1)))
            model_count += nbreakpts + 1
        else:
            mod_assigns.append([np.random.choice(bg_cats)]
                               + list(xrange(model_count, model_count + nbreakpts)))
            model_count += nbreakpts
        
    return mod_assigns


def gen_aa_fitness_profiles(mod_assigns, nBackgroundCats):
    """Generate vectors of fitness values for each category and breakpoint of a mixed model.
    """
    bg_mods = [uniform_fitnesses() for b in xrange(0, nBackgroundCats)]
    
    fitness_vecs = []
    
    for cat in mod_assigns:
        bg_mc_vecs = [bg_mods[m] for m in cat if m < len(bg_mods)]
        av_vecs = [uniform_fitnesses() for m in cat if m >= len(bg_mods)]
        fitness_vecs.append(bg_mc_vecs + av_vecs)

    return fitness_vecs


def uniform_fitnesses():
    return np.concatenate(([1.0], np.random.uniform(1.0 - FIT_RAD, 1.0 + FIT_RAD, AA_DIM)))

# Draws fitnesses uniformly from the surface of an 18-sphere with a specified radius
def distance_fitnesses(startpoint, radius):
    ns = np.random.normal(0, 1, AA_DIM)
    return np.concatenate([1.0], startpoint + (radius * (ns / np.linalg.norm(ns))))
    

def make_full_model_set(part_fit_vecs, kappa):
    return [make_partition_model_set(vecs, kappa) for vecs in part_fit_vecs]


def make_partition_model_set(vecs, kappa):
    paramlists = [pyvolve.MutSel_Sanity("mutsel", {"fitness":vec, "kappa":kappa, "popsize":10000})()
                  for vec in vecs]
    matrices = [KimuraMutSelMatrix("mutsel", params)() for params in paramlists]
    for i in range(0, len(paramlists)):
        paramlists[i].update({"matrix":matrices[i]})
    return [pyvolve.Model("custom", plist, name=("bp%d" % i)) for (i, plist) in enumerate(paramlists)]


def make_partition_set(cat_sizes, root_freq_set, model_set, model_assignment):
    if root_freq_set is None:
        return [pyvolve.Partition(models = ms, size=nk, root_model_name="bp0")
                for (ms, nk) in it.izip(model_set, cat_sizes)]
    else:
        root_seqs = [''.join(np.random.choice(MOLECULES.codons, size=nk, p=freqs).flat)
                     for (nk, freqs)
                     in zip(cat_sizes, root_freq_set)]
        return [pyvolve.Partition(models = ms, root_sequence=root, root_model_name="bp0")
                for (ms, root) in it.izip(model_set, root_seqs)]
    

def engrave_tree(treestr, brk, nbranch, lfile):
    """Inscribe pyvolve model flags into a tree string according to the given pattern."""
    model_flags = ["_bp%d_" %(len(brk) - brk.index(i)) if i in brk else '' for i in xrange(0,nbranch)]
    branch_strs = re.split(r':',treestr)
    flagged = (re.sub(r'(\d+\.\d+[eE]?-?\d*)', r':\1' + f, b)
	       for b, f in zip(branch_strs[1:], model_flags))
    l_tree_string = branch_strs[0] + ''.join(flagged) + ';'
    lfile.write(l_tree_string + '\n')
    ltree = pyvolve.read_tree(tree=l_tree_string)
    return ltree


def tax_labels(tree_str):
    """Return a list of strings representing the tip labels for a newick tree string."""
    handle = StringIO(tree_str)
    pytree = Phylo.read(handle, "newick")
    return [node.get_name for node in pytree.get_terminals()]


def calc_bp_mixprobs(mod_assign, cat_sizes, fitness_matrix):
    "Calculate 'true' frequency parameters for bp-specific mixture models."
    nsites = sum(cat_sizes)
    cat_freqs = [float(s)/float(nsites) for s in cat_sizes]
    mods_by_bp = zip(*mod_assign)
    fits_by_bp = zip(*fitness_matrix)
    mods_probs_fits = [zip(m, cat_freqs, fits_by_bp[i]) for i, m in enumerate(mods_by_bp)]
    bp_mixprobs = []
    
    for m_p_f in mods_probs_fits:
        prob_dict = {}
        for m, p, f  in m_p_f:
            if prob_dict.get(m) is not None:
                prob_dict[m][0] += p
            else:
                prob_dict[m] = [p, f]
        bp_mixprobs.append(sorted(list(prob_dict.values()), key = lambda x: x[0], reverse=True))
        
    return bp_mixprobs
        

def make_sim_info(cat_sizes, cat_freqs, brk, root_freq_set, model_set):
    """Return a dictionary containing the simulation parameters."""
    
    info_dict = {'bp%d_node' % (len(brk) - i) : str(b) for (i, b) in enumerate(brk)}
    info_dict.update({'part%d_size' % p : nk for p, nk in enumerate(cat_sizes) })
    for (bp, mixprobs) in enumerate(cat_freqs):
        for (modIdx, p) in enumerate(mixprobs):
            prob, fits = tuple(p)
            info_dict.update({'bp%s_relproba%d_%d' % (bp, modIdx+1, bp+1) : str(prob)})
            info_dict.update({'bp%s_%d_m%d_fit%d_%d' % (bp, modIdx+1, modIdx, i, bp+1) : str(fit)
                              for (i, fit) in enumerate(fits)})
    if root_freq_set is None :
        for (cat, bp_mset) in enumerate(model_set):
            for (bp, m) in enumerate(bp_mset):
                freqs = m.params['state_freqs']
                info_dict.update({'bp%d_m%d_f%d' % (bp, cat, codon) : f
                                  for (codon, f) in enumerate(freqs)})
    else:
        for (i, rf) in enumerate(root_freq_set):
            info_dict.update({'m%d_f%d' % (i, j) : f for (j, f) in enumerate(rf)})
            
    return info_dict


def simulate(q, lock):
    while True:
        sim_no, e = q.get()
        if e is None:
            break
        e(seqfile=None, ratefile=None, infofile=None)
        seq_dict = e.get_sequences()
        seq_tuple = sorted(seq_dict.items())
        align = AlignIO.MultipleSeqAlignment(SeqRecord(Seq(seqstr, generic_dna),
                                                       id=tax,
                                                       description='',)
                                             for (tax, seqstr) in seq_tuple)
        out_align = os.path.join(out_dir,
                                 '%(pre)s_%(sim)d.fasta' % {"pre": args.prefix, "sim": sim_no})
        AlignIO.write(align, out_align, "fasta")
        lock.acquire()
        try:
            print "Process %d reporting: wrote alignment %d" % (os.getpid(), sim_no)
        finally:
            lock.release()
        q.task_done()


###------------------------Functions from itertools recipe book---------------------###

# Attrib. Python standard library docs
def flatten(iter_of_iters):
    """Flatten one level of nesting."""
    
    return it.chain.from_iterable(iter_of_iters)

# Attrib. Python standard library docs
# grouper('ABCDEFG', 3, 'x' --> ABC DEF Gxx"
def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks."""
    
    args = [iter(iterable)] * n
    return it.izip_longest(fillvalue=fillvalue, *args)

# Attrib. Python standard library docs
def repeatfunc(func, times=None, *args):
    """Repeat calls to func with specified arguments.

    Example: repeatfunc(random.random)
    """
    if times is None:
        return it.starmap(func, it.repeat(args))
    return it.starmap(func, it.repeat(args, times))


    
###----------------------------------MAIN PROCEDURE--------------------------------###

if __name__ == "__main__":

    parser = make_parser()
    args = parser.parse_args()

    nbranch = 2*args.ntaxa-2
    out_dir = os.getcwd()

    ncats= args.nBackgroundCats + args.nModelChangeCats + args.nAdaptiveCats

    # Random generation of model parameters
    topologies = r_sim_trees(args.ntaxa, args.ntrees, balance=args.balance)
    breaks = [roll_breakpoints(topologies[i], args.nbreakpts) for i in xrange(0,  args.ntrees)]
    cat_size_vecs = [multinomial_cats(ncats, args.nsites) for t in topologies]
    model_assignments = [gen_model_assignment(args.nBackgroundCats,
                                              args.nModelChangeCats,
                                              args.nAdaptiveCats,
                                              args.nbreakpts)
                         for i in xrange(0, args.ntrees)]
    fitness_matrices = [gen_aa_fitness_profiles(m_a, args.nBackgroundCats)
                        for m_a
                        in model_assignments]
    
    if not args.root_freqs:
        root_freq_sets = [None for i in xrange(0, args.ntrees)]
    else:
        bg_root_freqs = [np.random.dirichlet([1.]*len(MOLECULES.codons), args.nBackgroundCats)
                         for i in xrange(0, args.ntrees)]
        root_freq_sets = [[brf[cat[0]]
                           for cat in m_a]
                          for (brf, ma) in zip(bg_root_freqs, model_assignments)]
        
    # Construction of pyvolve data structures
    ltreefile = open(os.path.join(out_dir, args.prefix + '_labelled.trees'), 'w')
    labelled_trees = [engrave_tree(topo, brk, nbranch, ltreefile)
                      for (topo, brk) in it.izip(topologies, breaks)]
    ltreefile.close()
    model_sets = [make_full_model_set(matrix, args.kappa) for matrix in fitness_matrices]
    partition_sets = [make_partition_set(*params)
                      for params in it.izip(cat_size_vecs,
                                            root_freq_sets,
                                            model_sets,
                                            model_assignments)]   
    
    evolvers = [pyvolve.Evolver(tree=ltree, partitions=part_set)
                for (ltree, part_set) in it.izip(labelled_trees, partition_sets)]

    
    # Pack generated simulation parameters for file output
    rate_matrices = [str(list(m.extract_rate_matrix().flat)) for m in flatten(flatten(model_sets))]
    mixprob_vecs = [calc_bp_mixprobs(m_a, nk, f)
                    for m_a, nk, f
                    in zip(model_assignments, cat_size_vecs, fitness_matrices)]
    sim_info_dicts = [make_sim_info(*info)
            for info
            in it.izip(cat_size_vecs, mixprob_vecs, breaks, root_freq_sets, model_sets)]

    
    # Write topologies and simulation info to disk
    newick_file = open(os.path.join(out_dir, args.prefix + '.trees'), 'w')
    newick_file.write('\n'.join(topologies))
    newick_file.close()

    with open(os.path.join(out_dir, args.prefix + '_siminfo.csv'), 'w') as csvfile:        
        writer = csv.DictWriter(csvfile, fieldnames = sorted(sim_info_dicts[0].keys()))
        writer.writeheader()
        writer.writerows(sim_info_dicts)

    newick_file = open(os.path.join(out_dir, args.prefix + '_qmats.txt'), 'w')
    newick_file.writelines('\n'.join(rate_matrices))
    newick_file.close()

    with open(os.path.join(out_dir, args.prefix + '_brk.csv'), 'w') as brk_file:
        writer = csv.writer(brk_file, delimiter='\t')
        for brk in breaks:
            writer.writerow([str(x) for x in brk])


    # Run simulations and write sequences to disk

    lock = mp.Lock()
    q = mp.JoinableQueue(maxsize=args.ntrees)
    
    procs = [mp.Process(target=simulate, args=(q, lock)) for i in xrange(0, 6)]

    for p in procs:
       print "Starting process " + p.name + "..."
       p.start()

    for i in enumerate(evolvers):
        q.put(i)
        
    q.join()

    for i in range(0, 6):
        q.put((None,None))
    for p in procs:
        p.join()
        print "Terminating process " + p.name + "."

    sys.exit(0)


###--------------------------------END OF FILE------------------------------------###
