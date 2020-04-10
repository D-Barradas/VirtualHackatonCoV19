#!/usr/bin/env python
# coding: utf-8

# In[19]:


# %load get_nextstrain_data.py
### https://github.com/tomasMasson/covid-19-Annotations-on-Structures.git
#!/usr/bin/env python3

import json
import requests
from utils.sm_annotations import Annotation
import matplotlib.colors as mcolors

def parse_json(input_file):
    ''' Wrapper to colect the data from parse_node. Currently, results are organized as a list for each isolate. First field contains isolate name, second field mutation data and third field all the metadata'''
    
    data = []

    def parse_node(input_file, array):    
        '''
        Search recursively across a json dictionary and return the keys corresponding to name, branch_attrs (mutations info), node_attrs (metadata, e.g. authors).
        ''' 
        root = input_file
        if 'history' not in root['name'] and 'NODE' not in root['name']:
            isolate = root['name']
            mutations_node = root['branch_attrs']['mutations']
            mutations = [(mutation, mutations_node[mutation]) for mutation in mutations_node if mutation != 'nuc']
            author = root['node_attrs']['author']
            gisaig = root['node_attrs']['gisaid_epi_isl']
            if mutations:
                for item in mutations:
                    array.append([isolate, [item], author, gisaig])
        if 'children' in root:
            for record in root['children']:
                parse_node(record, array)
        return array

    results = parse_node(input_file, data)
    return results

def get_nextstrain_data():
    ''' Fetch data from nextstrain raw data site.'''

    url = 'https://data.nextstrain.org/ncov.json'
    r = requests.get(url)
    nextstrain_json = r.json()
    results = parse_json(nextstrain_json['tree'])
    return results

def give_format_output(data):
    ''' Transform the raw data from nextstrain into a .csv file. ORFs names aer replaced by UniProtKB identifiers.'''

    identifiers = {'ORF1b': 'P0DTD1', 'ORF1a': 'P0DTD1', 'S': 'P0DTC2', 'ORF3a': 'P0DTC3', 'E': 'P0DTC4', 'M': 'P0DTC5', 'ORF6': 'P0DTC6', 'ORF7a': 'P0DTC7', 'ORF7b': 'P0DTD8', 'ORF8': 'P0DTC8', 'N': 'P0DTC9', 'ORF14': 'P0DTD3', 'ORF9b': 'P0DTD2', 'ORF10': 'A0A663DJA2'}

    with open('nextstrain_data.csv', 'w') as fh:
        headers = f'Protein,ORF,Mutation,Isolate,Author,GISAID\n'
        fh.write(headers)
        for record in data:
            if record[1][0][0] in identifiers:
                protein = identifiers[record[1][0][0]]
            else:
                protein = record[1][0][0]
            orf = record[1][0][0]
            if orf == 'ORF1b':
                original = record[1][0][1][0][0]
                mutant = record[1][0][1][0][-1]
                position = int(record[1][0][1][0][1:-1]) + 4401
                mutation = original + str(position) + mutant
            else:
                mutation = record[1][0][1][0]
            author = record[2]['value']
            isolate = record[0]
            gisaid = record[3]['value']
            output = f'{protein},{orf},{mutation},{isolate},{author},{gisaid}\n'
            fh.write(output)
            
def annotation_swiss(data):
    ''' Annotation swissprot '''
    annotation = Annotation()

    identifiers = {'ORF1b': 'P0DTD1', 'ORF1a': 'P0DTD1', 'S': 'P0DTC2', 'ORF3a': 'P0DTC3', 'E': 'P0DTC4',
                   'M': 'P0DTC5', 'ORF6': 'P0DTC6', 'ORF7a': 'P0DTC7', 'ORF7b': 'P0DTD8', 'ORF8': 'P0DTC8',
                   'N': 'P0DTC9', 'ORF14': 'P0DTD3', 'ORF9b': 'P0DTD2', 'ORF10': 'A0A663DJA2'}
    
    colors = {
             'ORF1b' :'#ff0000', 
             'ORF1a' :'#fffb00', 
             'S'     :'#00ff15', 
             'ORF3a' :'#6600ff', 
             'E'     :"#0000ff",
             'M'     :'#ff0066', 
             'ORF6'  : '#f0f0f0', 
             'ORF7a' :'#0804f9', 
             'ORF7b' : '#00ffdd', 
              'ORF8' : '#00ccff',
              'N'    :'#ff00ff', 
              'ORF14':'#00ff00',
              'ORF9b':'#ff00dd', 
              'ORF10':'#0f0f0f'
             } 
    for record in data:
        if record[1][0][0] in identifiers:
            protein = identifiers[record[1][0][0]]
            color = [colors[record[1][0][0]]]
        else:
            protein = record[1][0][0]
            color = [identifiers[record[1][0][0]]]
        orf = record[1][0][0]
        if orf == 'ORF1b':
            original = record[1][0][1][0][0]
            mutant = record[1][0][1][0][-1]
            position = int(record[1][0][1][0][1:-1]) + 4401
            mutation = original + str(position) + mutant
        else:
            mutation = record[1][0][1][0]
            position = record[1][0][1][0][1:-1]
        author = record[2]['value']
        isolate = record[0]
        gisaid = record[3]['value']
#         print (color[0])
        ann = 'gisaid:%s,mutation:%s,isolate:%s,author:%s'%( gisaid, mutation,isolate, author    )
        annotation.add(protein,(int(position),int(position)),color[0],ann )
    return annotation

    
def main():
    results = get_nextstrain_data()
    a= annotation_swiss(results)
    print (a)
#     give_format_output(results)

if __name__ == '__main__':
    main()

