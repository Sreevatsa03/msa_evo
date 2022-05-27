# EvoAlign

A library with the functionality to align 2 or more amino acid sequences using muliple modified Smith-Waterman alignment agents. The best possible alignment is generated based on a variety of fitness criteria (scoring systems). Using this library, you can not only generate alignments but also visualize tradeoffs between fitness criteria.

# Get Started

```python
from EvoAlign import EvoAlign

# read in fasta files
all = 'data/all.fasta'
 
# create alignment environment
evo_a = EvoAlign()

# load the fasta file
evo_a.read_fasta(all)

# run alignment and get visualizations
evo_a.align(gens=1000, dom=100, status=100, show=True)

# output list of fitness criteria
print("Fitness Criteria:")
evo_a.fitness_criteria()

# saves alignment to fasta
evo_a.save_alignment()
```

# Authors

[Sreevatsa Nukala](https://github.com/Sreevatsa03), [John Drohan](https://github.com/jdrohan356), [Rachel Utama](https://github.com/rootma21), [Sanjana Bhagavtula](https://github.com/bhagavatulasa)

