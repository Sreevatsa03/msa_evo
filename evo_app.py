from EvoAlign import EvoAlign, Evo

dash = 'data/dash.fasta'
dog = 'data/canis_lupus_familiaris.fasta'
human = 'data/homo_sapiens.fasta'
mouse = 'data/mus_musculus.fasta'
all = 'data/all.fasta'
fastas = [dog, human, mouse]
 
# create alignment environment
evo_a = EvoAlign()

# load the fasta file
evo_a.read_fasta(all)

# run alignment and get visualizations
evo_a.align(gens=2000, dom=100, status=100, show=True)

# saves alignment to fasta
evo_a.save_alignment()