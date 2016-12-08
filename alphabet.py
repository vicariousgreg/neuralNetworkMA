def get_nucleotides():
    return "ACTG"

def get_amino_acids():
    return "GPAVLIMCFYWHKRQNEDST"

def handle_special_letters(alphabet, counts, letter):
    if alphabet == get_nucleotides():
        if letter == "N":
            val = 1.0 / len(alphabet)
            for letter in alphabet:
                counts[letter] += val
    elif alphabet == get_amino_acids():
        if letter == "X":
            val = 1.0 / len(alphabet)
            for letter in alphabet:
                counts[letter] += val
        elif letter == "B":
            counts["N"] += 0.50
            counts["D"] += 0.50
        elif letter == "Z":
            counts["Q"] += 0.50
            counts["E"] += 0.50
