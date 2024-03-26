import Levenshtein

# Function to compute Levenshtein distance
def compute_levenshtein_distance(seq1, seq2):
    return Levenshtein.distance(seq1, seq2)

# Function to generate cyclic permutations of a sequence
def cyclic_permutations(sequence):
    rotations = []
    for i in range(len(sequence)):
        rotated_seq = sequence[i:] + sequence[:i]
        rotations.append(rotated_seq)
    return rotations

# Read sequences from the first tabular TXT file
file_path_1 = "Enriched_2kbmotifs_enriched_10kb_regions.txt"  # Replace with the path to your first input file
sequences_1 = {}
with open(file_path_1, 'r') as file:
    for line in file:
        line = line.strip().split('\t')
        if len(line) >= 2:
            sequence_name, sequence = line[0], line[1]
            sequences_1[sequence_name] = sequence
        else:
            print(f"Skipping line: {line}")

# Read sequences from the second tabular TXT file
file_path_2 = "Enriched_2kbmotifs_reverseComp_enriched_10kb_regions.txt"  # Replace with the path to your second input file
sequences_2 = {}
with open(file_path_2, 'r') as file:
    for line in file:
        line = line.strip().split('\t')
        if len(line) >= 2:
            sequence_name, sequence = line[0], line[1]
            sequences_2[sequence_name] = sequence
        else:
            print(f"Skipping line: {line}")

# Create a list to store the results
results = []

# Compute and store the shortest Levenshtein distance for each pair-wise comparison
for name_1, seq1 in sequences_1.items():
    for name_2, seq2 in sequences_2.items():
        # Calculate Levenshtein distance for all cyclic permutations of seq2
        seq2_rotations = cyclic_permutations(seq2)

        min_distance = float('inf')
        for rotated_seq2 in seq2_rotations:
            distance = compute_levenshtein_distance(seq1, rotated_seq2)
            min_distance = min(min_distance, distance)

        results.append((name_1, name_2, min_distance))

# Write the results to a tabulated TXT file
output_file_path = "levenstein_distances_FR.txt"  # Replace with your desired output file path
with open(output_file_path, 'w') as output_file:
    for result in results:
        output_file.write(f"{result[0]}\t{result[1]}\t{result[2]}\n")

print(f"Results written to {output_file_path}")

