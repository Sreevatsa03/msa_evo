from EvoAlign import Chart

chart = Chart()

# create dict of fasta files to use
dog = 'data_sources/P53_test_data/canis_lupus_familiaris.fasta'
human = 'data_sources/P53_test_data/homo_sapiens.fasta'
mouse = 'data_sources/P53_test_data/mus_musculus.fasta'
fasta_dict = chart.fasta_dict([dog, human, mouse])

# plot alignment chart
chart.alignment_chart('data_sources/P53_test_data/all.fasta')
