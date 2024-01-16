import pandas as pd

'''A module contining functions for generating and filtering DNA encoding sequences'''

BASES = 'ATCG'
def generate_3mers():
    sequences = []
    for a in BASES:
        n1 = a
        for b in BASES:
            n2 = n1 + b
            for c in BASES:
                sequences.append(n2+c)
    return(sequences)

def generate_4mers():
    sequences = []
    for a in BASES:
        n1 = a
        for b in BASES:
            n2 = n1 + b
            for c in BASES:
                n3 = n2 + c
                for d in BASES:
                    sequences.append(n3+d)
    return(sequences)

def generate_5mers():
    sequences = []
    for a in BASES:
        n1 = a
        for b in BASES:
            n2 = n1 + b
            for c in BASES:
                n3 = n2 + c
                for d in BASES:
                    n4 = n3 + d
                    for e in BASES:
                        sequences.append(n4+e)
    return(sequences)

def generate_6mers():
    sequences = []
    for a in BASES:
        n1 = a
        for b in BASES:
            n2 = n1 + b
            for c in BASES:
                n3 = n2 + c
                for d in BASES:
                    n4 = n3 + d
                    for e in BASES:
                        n5 = n4 + e
                        for f in BASES:
                            sequences.append(n5+f)
    return(sequences)

def generate_7mers():
    sequences = []
    BASES = 'ATCG'
    for a in BASES:
        n1 = a
        for b in BASES:
            n2 = n1 + b
            for c in BASES:
                n3 = n2 + c
                for d in BASES:
                    n4 = n3 + d
                    for e in BASES:
                        n5 = n4 + e
                        for f in BASES:
                            n6 = n5 + f
                            for g in BASES:
                                sequences.append(n6+g)
    return(sequences)

def generate_8mers():
    sequences = []
    BASES = 'ATCG'
    for a in BASES:
        n1 = a
        for b in BASES:
            n2 = n1 + b
            for c in BASES:
                n3 = n2 + c
                for d in BASES:
                    n4 = n3 + d
                    for e in BASES:
                        n5 = n4 + e
                        for f in BASES:
                            n6 = n5 + f
                            for g in BASES:
                                n7 = n6 + g
                                for h in BASES:
                                    sequences.append(n7+h)
    return(sequences)

def get_sequences(n):
    if n==3:
        sequences=generate_3mers()
    if n==4:
        sequences=generate_4mers()
    if n==5:
        sequences=generate_5mers()
    if n==6:
        sequences=generate_6mers()
    if n==7:
        sequences=generate_7mers()
    if n==8:
        sequences=generate_8mers()
    return(sequences)


def get_filtered_dict(sequences,strict=False):
    
    #create a dictionary of enumerated base sequences
    sequence_dictA=dict(enumerate(sequences))

    #replace any sequence values containing 3-peat bases with "*"
    for k,v in sequence_dictA.items():
        if 'AAA' in v or 'TTT' in v or 'GGG' in v or 'CCC' in v:
            sequence_dictA[k] = '*'*len(sequences[0])

    #remove any sequences that are not unique at 2 or more positions
    sequence_dictB=sequence_dictA
    if strict:
        max_matches=len(sequences[0])-3
    else:
        max_matches=len(sequences[0])-2
    for k,v in sequence_dictA.items():
        if '*' in v :
            continue
        for k1,v1 in sequence_dictB.items():
            EqualBases=0 #Initiates counter to 0
            if '*' in v1:
                continue
            if v1 == v:
                continue
            for Idx, Character in enumerate(v1):
                if Character == v[Idx]:
                    EqualBases=EqualBases+1
            if EqualBases>max_matches:
                sequence_dictB[k1] = '*'*len(sequences[0])
    return(sequence_dictA)

               
def print_sequences(sequence_dict):
    for k,v in sequence_dict.items():
        if '*' not in v:
            print(v)

def output_final_sequences(sequence_dict):
    sequence_df=pd.DataFrame()
    FinalSequences=[]
    for k,v in sequence_dict.items():
        if '*' not in v:
            FinalSequences.append(v)
    sequence_df[str(len(FinalSequences[0]))+"-mer Sequences"]=FinalSequences
    sequence_df.to_csv(str(len(FinalSequences[0]))+'-mer Sequences Output.csv')
    print('Generated a list of '+ str(len(FinalSequences))+' sequences')
