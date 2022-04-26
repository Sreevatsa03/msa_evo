from EvoAlign import EvoAlign, Evo


dash = 'data/dash.fasta'
dog = 'data/canis_lupus_familiaris.fasta'
human = 'data/homo_sapiens.fasta'
mouse = 'data/mus_musculus.fasta'
all = 'data/all.fasta'
fastas = [dog, human, mouse]
 

evo_a = EvoAlign()
evo_a.read_fasta(all)
evo_a.align(gens=2000, dom=100, status=100)