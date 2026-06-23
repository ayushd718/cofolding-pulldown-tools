from cofolding_pulldown_tools.fasta import clean_fasta, slice_fasta

def main():
    clean_fasta('tests/data/x.fasta')
    clean_fasta('tests/data/y.fasta')
    slice_fasta('tests/data/y_cleaned.fasta', max_length=500, window=50)
    slice_fasta('tests/data/x_cleaned.fasta', max_length=500, window=100)

if __name__ == "__main__":
    main()