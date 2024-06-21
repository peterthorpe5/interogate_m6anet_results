import os
from collections import defaultdict
import matplotlib.pyplot as plt
from interogate.parse_gtf import parse_gff_gft
from interogate.return_dict import generate_transcript_coordinates





# Example usage:
file_path = 'data/test.gtf'  # Replace with the path to your GFF or GTF file
features = parse_gff_gft(file_path)


transcript_dict = generate_transcript_coordinates(features)

for transcript, exons in transcript_dict.items():
    for exon, coordinates in exons.items():
        out_data = (f'{transcript} exon {exon}: {coordinates}')
        print(out_data)

TEST = '''
###############################
# TESTS
# Query example
query_transcript = "AT1G01010.1"
query_position = 3
exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)

if exon_number is not None:
    print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
    print(f'Transcript {query_transcript} has {total_exons} exons.')
else:
    print(f'Position {query_position} is not found in transcript {query_transcript}.')

# Query example2
query_transcript = "AT1G01020.4"
query_position = 426
exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
print("AT1G01020.4 should have 7 exons and position 426 that should be exon 4")
if exon_number is not None:
    print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
    print(f'Transcript {query_transcript} has {total_exons} exons.')
else:
    print(f'Position {query_position} is not found in transcript {query_transcript}.')


# Query example3
query_transcript = "AT1G01020.4"
query_position = 860
exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
print("AT1G01020.4 should have 7 exons and position 860 should be exon 7")
if exon_number is not None:
    print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
    print(f'Transcript {query_transcript} has {total_exons} exons.')
else:
    print(f'Position {query_position} is not found in transcript {query_transcript}.')


# Query example4
query_transcript = "AT1G01020.4"
query_position = 861
exon_number, total_exons = query_transcript_exon(transcript_dict, query_transcript, query_position)
print("AT1G01020.4 should have 7 exons and position 861 should BREAK!!!")
if exon_number is not None:
    print(f'Transcript {query_transcript} position {query_position} is in exon {exon_number}.')
    print(f'Transcript {query_transcript} has {total_exons} exons.')
else:
    print(f'Position {query_position} is not found in transcript {query_transcript}.')

    '''
