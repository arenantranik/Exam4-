import sys

class DNASequencer:
    def __init__(self, filename):
        self.filename = filename
        self.kmers_db = {}

    def generate_kmers(self, sequence, k):
        return {sequence[i:i+k]: sequence[i+1:i+1+k] for i in range(len(sequence) - k) if i+k < len(sequence)}

    def process_sequences(self, k):
        try:
            with open(self.filename, 'r') as file:
                for line in file:
                    if not line.startswith('>'):
                        sequence = line.strip()
                        kmers = self.generate_kmers(sequence, k)
                        self.update_kmers_db(kmers)
        except FileNotFoundError:
            print("File not found. Please check the file path.")
            sys.exit(1)

    def update_kmers_db(self, kmers):
        for kmer, next_kmer in kmers.items():
            if kmer in self.kmers_db:
                self.kmers_db[kmer].add(next_kmer)
            else:
                self.kmers_db[kmer] = {next_kmer}

    def find_smallest_unique_k(self):
        k = 1
        while True:
            self.kmers_db = {}
            self.process_sequences(k)
            if all(len(v) == 1 for v in self.kmers_db.values()):
                return k
            if not self.kmers_db:  # No k-mers were found, k is too large
                return -1
            k += 1

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
        return
    
    filename = sys.argv[1]
    sequencer = DNASequencer(filename)
    smallest_k = sequencer.find_smallest_unique_k()
    print(f"The smallest k with a unique subsequent k-mer for each k-mer is: {smallest_k}")

if __name__ == "__main__":
    main()
