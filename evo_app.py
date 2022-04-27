from EvoAlign import EvoAlign

# read in fasta files
dash = 'data/dash.fasta'
dog = 'data/canis_lupus_familiaris.fasta'
human = 'data/homo_sapiens.fasta'
mouse = 'data/mus_musculus.fasta'
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