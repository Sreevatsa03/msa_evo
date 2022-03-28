from visuals import align_jd

# create dict of fasta files to use
dog = 'data_sources/P53_test_data/canis_lupus_familiaris.fasta'
human = 'data_sources/P53_test_data/homo_sapiens.fasta'
mouse = 'data_sources/P53_test_data/mus_musculus.fasta'
fasta_dict = align_jd.fasta_dict([dog, human, mouse])

# plot alignment chart
align_jd.alignment_chart('data_sources/P53_test_data/all.fasta')
