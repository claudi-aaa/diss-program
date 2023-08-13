
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO


def calculate_pairwise_dm(kmer_size, seed):
    """
    Takes input of kmer size and seed and calculates pairwise distance matrix
    """
    aa_filepath = f"assign_data/rand{kmer_size}_seed{seed}.phy"
    aln = AlignIO.read(open(aa_filepath), 'phylip')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    print(dm)

