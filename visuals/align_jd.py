import urllib.request as urlreq
from dash import Dash, html
import dash_bio as dashbio
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline


def fasta_dict(files):

    fasta_dict = {}

    for filename in files:
        species = filename.split('/')[-1][:-6]
        with open(filename, 'r', encoding='utf-8') as infile:
            for line in infile:
                line = line.replace('\n', '')

                if '[' not in line:
                    fasta_dict[species] = fasta_dict.get(species, '') + line

    return fasta_dict


def alignment_chart(filename=None):
    app = Dash()

    # fasta = urlreq.urlopen('https://git.io/alignment_viewer_p53.fasta').read().decode('utf-8')

    with open(filename, 'r', encoding='utf-8') as fasta:
        app.layout = html.Div([
            dashbio.AlignmentChart(
                id='alignment-viewer',
                data=fasta
            ),
        ])

    app.run_server()