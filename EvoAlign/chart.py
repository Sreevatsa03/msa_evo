import urllib.request as urlreq
from dash import Dash, html
import dash_bio as dashbio
from Bio import SeqIO


class Chart():

    def __init__(self):
        pass

    def alignment_chart(self, filename=None):
        app = Dash()

        fasta = urlreq.urlopen('https://git.io/alignment_viewer_p53.fasta').read().decode('utf-8')

        with open(filename, 'r', encoding='utf-8') as fasta:
            app.layout = html.Div([
                dashbio.AlignmentChart(
                    id='alignment-viewer',
                    data=fasta
                ),
            ])

        app.run_server()


