import os
import requests
import time
from Bio.Blast import NCBIXML, NCBIWWW
from Bio import SeqIO

potato_id = '4113'
this_name = __file__.split('\\')[-1]
directory = __file__.replace(this_name, '')
#file_name = directory + 'fasta/test.fasta'
database = directory + '/db/PGSC_DM_v3.4_pep_nonredundant.fasta'
accessions_file = directory + '/sequence_list/replicats_combine.csv'

def timer(func):
        def wrapper(*args, **kwargs):

                _begin  = time.perf_counter()
          
                result = func(*args, **kwargs)

                _end  = time.perf_counter()

                perf_time = round(_end - _begin, 2)
                unit = ' sec'
                if perf_time > 60.0:
                        perf_time = round(perf_time/60, 2)
                        unit = ' min'
                
                print(f'func {func.__name__} performed in {str(perf_time)}{unit}')
                return result
        return wrapper

def split_fasta(input_file, output_prefix, max_sequences=100):

    sequences = []  
    file_count = 0
    output_files =[]  

    with open(input_file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record) 
            if len(sequences) >= max_sequences:
                output_file = f"{output_prefix}_part{file_count + 1}.fasta"
                with open(output_file, "w") as output:
                    SeqIO.write(sequences, output, "fasta")
                output_files.append(output_file)
                print(f"Written {len(sequences)} sequences to {output_file}")
                sequences = [] 
                file_count += 1 

    if sequences:
        output_file = f"{output_prefix}_part{file_count + 1}.fasta"
        with open(output_file, "w") as output:
            SeqIO.write(sequences, output, "fasta")
        output_files.append(output_file)
        print(f"Written {len(sequences)} sequences to {output_file}")

    return output_files

def find_accessions(accessions_file):
        accessions =[]
        iBAQs_list =[]
        with open(accessions_file, 'r') as file:
                accessions_lines  = file.read().splitlines()
        for accessions_line in accessions_lines:
                line = accessions_line.split(',')
                accession = line.pop(0)
                accessions.append(accession)
                iBAQs_list.append(line)
        return accessions, iBAQs_list

def find_sequences(accessions, database): 

        with open(database, 'r') as db:
                delimiter = "*\n"
                protein_sequences =  [x+delimiter for x in db.read().split(delimiter) if x]
                if '>' not in protein_sequences[-1]:
                        protein_sequences.pop(-1)
        
        output_file_name = database.replace('.fasta', '_copy.fasta')

        with open(output_file_name, 'w') as file:
 
                for accession in accessions:
                        for protein_sequence in protein_sequences: 
                                if accession in protein_sequence: 
                                        file.write(protein_sequence)
       
        return output_file_name

def fetch_organism_from_uniprot(uniprot_id):

        url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}.json"
        try:
                response = requests.get(url)
        except requests.exceptions.RequestException as e:
                print(f"Request failed: {e}")
                return f"Request failed: {e}"

        if response.status_code == 200:
            try:
                # Attempt to parse the JSON response
                data = response.json()
                organism_name = data.get('organism', {}).get('scientificName', None)
                return organism_name
            except requests.exceptions.JSONDecodeError:
                print(f"Error decoding JSON response for UniProt ID: {uniprot_id}")
                return None
        elif response.status_code == 404:
            print(f"UniProt ID {uniprot_id} not found.")
            return None
        else:
            print(f"Failed to retrieve data for UniProt ID {uniprot_id}. Status code: {response.status_code}")
            return None

@timer
def get_accessions(file_name, organism_id):
        with NCBIWWW.qblast("blastp", "swissprot", open(file_name).read(), hitlist_size=1) as result_stream:

                blast_records = NCBIXML.parse(result_stream)
                accession_list=[]
                for blast_record in blast_records:

                        input_title = blast_record.query
                        # Get the first hit with the best score
                        if blast_record.alignments:
                                best_hit = blast_record.alignments[0]  # First hit (best score)
                        
                        # Access information from the hit
                                hit_title = best_hit.hit_def.replace(',',';')  
                                #hit_id = best_hit.hit_id.split("|")[1]
                                #organism_name = fetch_organism_from_uniprot(hit_id)
                                hit_accession = best_hit.accession  
                                hit_score = best_hit.hsps[0].score
                                hit_e_value = best_hit.hsps[0].expect 
                                accession_list.append((input_title, hit_accession, hit_title)) 
                        else: 
                                accession_list.append((input_title, 'NotFound', 'NotFound'))

        return tuple(accession_list)

# @timer
# def get_accessions(file_name, organism_id):
#     accession_list = []

#     # Parse the input file (assume FASTA format) and process each sequence
#     for n, record in enumerate(SeqIO.parse(file_name, "fasta")):
#         input_title = record.description  # Get the title/description of the sequence
#         input_sequence = str(record.seq)  # Get the sequence itself
        
#         # Perform BLAST query for the current sequence
#         with NCBIWWW.qblast("blastp", "swissprot", input_sequence, hitlist_size=1) as result_stream:
#             blast_records = NCBIXML.parse(result_stream)
            
#             for blast_record in blast_records:
#                 if blast_record.alignments:
#                     # Get the first hit with the best score
#                     best_hit = blast_record.alignments[0]
#                     hit_title = best_hit.hit_def.replace(',', ';')  
#                     hit_accession = best_hit.accession  
#                     accession_list.append((input_title, hit_accession, hit_title))
#                 else:
#                     accession_list.append((input_title, 'NotFound', 'NotFound'))
#         print(str(n) + '/' + str(len(SeqIO.parse(file_name, "fasta"))), end='\r')
                

#     return tuple(accession_list)

def get_uniprot_subcellular_location(accession):
    # UniProt API URL
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    
    # Send request to UniProt API
    response = requests.get(url)
    
    if response.status_code == 200:
        # Parse the JSON response
        data = response.json()
        
        # Check if subcellular location information exists
        if 'comments' in data:
            for comment in data['comments']:
                if comment['commentType'] == 'SUBCELLULAR LOCATION':
                    # Extract the subcellular location description
                    subcellular_location = comment['subcellularLocations'][0]['location']['value'].replace(',',' ')
                    return subcellular_location
        return "Subcellular location information not found."
    
    else:
        return f"Failed to retrieve data for accession {accession}, HTTP status code: {response.status_code}"

def get_uniprot_function(accession):
    # UniProt API URL
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    
    # Send request to UniProt API
    response = requests.get(url)
    
    if response.status_code == 200:
        # Parse the JSON response
        data = response.json()
        
        # Check if subcellular location information exists
        if 'comments' in data:
            for comment in data['comments']:
                if comment['commentType'] == 'FUNCTION':
                    # Extract the subcellular location description
                    function = comment['texts'][0]['value'].replace(',', ' ')
                    return function
        return "Function information not found."
    
    else:
        return f"Failed to retrieve data for accession {accession}, HTTP status code: {response.status_code}"

def get_uniprot_biological_process(accession):
    # UniProt API URL
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    
    # Send request to UniProt API
    response = requests.get(url)
    
    if response.status_code == 200:
        # Parse the JSON response
        data = response.json()
        
        # Check if subcellular location information exists
        if 'keywords' in data:
            for ref in data['keywords']:
                if ref['category'] == 'Biological process':
                    # Extract the subcellular location description
                    process = ref['name'].replace(',', ' ')
                    return process

        return "Biological process not found."
    
    else:
        return f"Failed to retrieve data for accession {accession}, HTTP status code: {response.status_code}"
'''
Script Block
'''

if __name__ == "__main__":

        accessions, iBAQs_list = find_accessions(accessions_file)

        print(len(accessions))
        file_name = find_sequences(accessions, database)

        if len(accessions) > 100:
      
                file_names = split_fasta(file_name, directory + 'db/splitted_sequences')

                for file_name in file_names:
                        print(f' begins {file_name} data treatment...' )

                        Uniprot_accessions = get_accessions(file_name, potato_id)

                        for Uniprot_accession in Uniprot_accessions:
                        
                                with open(accessions_file.replace('.csv', '_output.csv'), 'a') as file:
                                        file.write(f'{Uniprot_accession[0]},{Uniprot_accession[2]},{Uniprot_accession[1]},{get_uniprot_subcellular_location(Uniprot_accession[1])}, {get_uniprot_function(Uniprot_accession[1])}, {get_uniprot_biological_process(Uniprot_accession[1])}')
                                        for n, accession in enumerate(accessions):
                                              if accession in Uniprot_accession[0]:
                                                    for iBAQ in iBAQs_list[n]:
                                                          file.write(f',{iBAQ}')
                                        file.write('\n')
                        print(f'{file_name} processed')
                        os.remove(file_name)

                
        else: 
                Uniprot_accessions = get_accessions(file_name, potato_id)

                for Uniprot_accession in Uniprot_accessions:

                        with open(accessions_file.replace('.csv', '_output.csv'), 'a') as file:
                                file.write(f'{Uniprot_accession},{get_uniprot_subcellular_location(Uniprot_accession)}, {get_uniprot_function(Uniprot_accession[1])}')
                                for n, accession in enumerate(accessions):
                                        if accession in Uniprot_accession[0]:
                                                for iBAQ in iBAQs_list[n]:
                                                        file.write(f',{iBAQ}')
                                file.write('\n')
                        print(f'{file_name} processed')
        print('Done') 