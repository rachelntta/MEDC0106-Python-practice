#TASK 1

#Question 1: calculate the total number of participants 
age = [45,35,37,49,31,45,43,30,48,36,38,44,39,36,40,42,38,30,31,39]
def count_participants(age): #define the function/ action that we want to perform
  total_participant = len(age) #use the len() command to count the total number of items in a string
  return total_participant

total_count = count_participants(age)
print(total_count)

#OR

total_number_of_participant= len(age)
print(total_number_of_participant)

#Question 2: Calculate the mean age of participants
age = [45,35,37,49,31,45,43,30,48,36,38,44,39,36,40,42,38,30,31,39]
def calculate_average(age):
  mean_age = sum(age) / len(age)
  return mean_age

  average_age = calculate_average(age)
  print("Average age:", average_age)

#OR
import statistics
print(statistics.mean(age))

age_avg = calculate_average(age)
print(age_avg)

#Question 3: Calculate the proportion of male ('M') to female ('F') participants
gender = ['F','F','F','M','F','F','F','M','F','M','F','F','F','M','F','F','F','M','F','F']

male_count = gender.count('M')
female_count = gender.count('F')

print(male_count)
print(female_count)

Proportion = male_count / female_count
print(Proportion)

#Question 4: Calculate the percentage of 'yes' ('Y') answers

answer = ['Y','N','N','Y','Y','Y','N','N','N','N','Y','N','N','N','Y','N','Y','N','N','Y']

percent_of_yes = len('Y') / len(answer)
print(percent_of_yes)

#TASK 2: count the number of each nucleotide
def count_atcg(sequence):
    count = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    a_count = 0
    t_count = 0
    c_count = 0
    g_count = 0

    for x in sequence:
      if x == 'A':
        count['A'] += 1
      elif x == 'T':
        count['T'] += 1
      elif x == 'C':
        count['C'] += 1
      elif x == 'G': 
        count['G'] += 1

probe_sequence = """\
ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAG\
CCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTG\
CGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGG\
CAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCC\
TGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAA\
CTACTGCAAC"""

print('Counts for probe sequence:')
print('A:', atcg_content['A'])
print('T:', atcg_content['T'])
print('C:', atcg_content['C'])
print('G:', atcg_content['G'])

#OR <model answer>
def count_atcg(sequence):
    count = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for nucleotide in sequence:
        count[nucleotide] += 1
    return count

probe_sequence = """\
ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAG\
CCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTG\
CGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGG\
CAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCC\
TGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAA\
CTACTGCAAC"""

atcg_content = count_atcg(probe_sequence)

print('Counts for probe sequence:')
print('A:', atcg_content['A'])
print('T:', atcg_content['T'])
print('C:', atcg_content['C'])
print('G:', atcg_content['G'])

#TASK 3: Sequence translation

def translate(sequence):
    codon_table = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
  }
    amino_acid_sequence = '' #Amino acid indicates the letter inside the '' sign

  for i in range(0, len(probe_sequence), 3): #for the 'probe_sequence' string 
      codon = probe_sequence[i:i+3] 
      amino_acid_sequence += codon_table.get(codon, 'X')
    return amino_acid_sequence
  
probe_sequence = """\
ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAG\
CCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTG\
CGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGG\
CAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCC\
TGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAA\
CTACTGCAAC"""

aa_sequence = translate(probe_sequence)
print('Translated:', aa_sequence)

#TASK 4: Sequence mutation

#Mutating the sequence
def mutate_sequence(sequence, position, mutation):
    mutated_sequence = list(probe_sequence) #create a new string/ list
    mutated_sequence[position - 1] = mutation 
    return ''.join(mutated_sequence) #the mutated sequence is converted back to a string 

probe_sequence = """\
ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAG\
CCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTG\
CGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGG\
CAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCC\
TGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAA\
CTACTGCAAC"""

mutated_sequence = mutate_sequence(probe_sequence, 12, 'A')
print(mutated_sequence)

# Now use your translate function to see why this mutation may cause disease
def translate(mutated_sequence):
    codon_table = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
  }
    amino_acid_sequence = ''

    for i in range(0, len(mutated_sequence), 3):
      codon = mutated_sequence[i:i+3] 
      amino_acid_sequence += codon_table.get(codon, 'X')
    return amino_acid_sequence
  
aa_sequence = translate(mutated_sequence)
print('Translated:', aa_sequence)
