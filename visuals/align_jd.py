import urllib.request as urlreq
from dash import Dash, html
import dash_bio as dashbio
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline


def read_files(files):

    fasta_dict = {}

    for filename in files:
        species = filename.split('/')[-1][:-6]
        with open(filename, 'r', encoding='utf-8') as infile:

            for line in infile:
                line = line.replace('\n', '')

                if '[' not in line:
                    fasta_dict[species] = fasta_dict.get(species, '') + line

    return fasta_dict


app = Dash(__name__)


def main():
    dog = '/Users/bongo/Documents/msa_evo/data_sources/' \
          'P53_test_data/canis_lupus_familiaris.fasta'
    human = '/Users/bongo/Documents/msa_evo/data_sources/' \
            'P53_test_data/homo_sapiens.fasta'
    mouse = '/Users/bongo/Documents/msa_evo/data_sources/' \
            'P53_test_data/mus_musculus.fasta'
    fasta_dict = read_files([dog, human, mouse])

    #data = urlreq.urlopen('https://git.io/alignment_viewer_p53.fasta').read().decode('utf-8')
    fasta = open(r'/Users/bongo/Documents/msa_evo/'
                 r'data_sources/P53_test_data/all.fasta')

    app.layout = html.Div([
        dashbio.AlignmentChart(
            id='alignment-viewer',
            data=fasta
        ),
    ])


    fasta.close()


if __name__ == '__main__':
    main()
    #app.run_server()


