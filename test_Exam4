import pytest
from Exam4 import DNASequencer  

def test_generate_kmers():
    sequencer = DNASequencer("dummy_filename")  # filename won't be used in this test
    sequence = "ATCGATCG"
    k = 3
    expected_kmers = {'ATC': 'TCG', 'TCG': 'CGA', 'CGA': 'GAT', 'GAT': 'ATC'}
    assert sequencer.generate_kmers(sequence, k) == expected_kmers

def test_update_kmers_db():
    sequencer = DNASequencer("dummy_filename")
    kmers = {'ATC': 'TCG', 'TCG': 'CGA'}
    sequencer.update_kmers_db(kmers)
    assert sequencer.kmers_db == {'ATC': {'TCG'}, 'TCG': {'CGA'}}

def test_process_sequences(mocker):
    sequencer = DNASequencer("path/to/testfile.fasta")
    mocker.patch('builtins.open', mocker.mock_open(read_data=">Header\nATCGATCG\n"))
    k = 3
    sequencer.process_sequences(k)
    expected_db = {'ATC': {'TCG'}, 'TCG': {'CGA'}, 'CGA': {'GAT'}, 'GAT': {'ATC'}}
    assert sequencer.kmers_db == expected_db

def test_find_smallest_unique_k(mocker):
    sequencer = DNASequencer("path/to/testfile.fasta")
    mocker.patch('builtins.open', mocker.mock_open(read_data=">Header\nATCGATCG\nATCGATCG\n"))
    mocker.patch.object(DNASequencer, 'process_sequences', side_effect=lambda k: None)
    mocker.patch.object(DNASequencer, 'kmers_db', new_callable=mocker.PropertyMock, 
                        side_effect=[{'AT': {'TC'}, 'TC': {'CG'}, 'CG': {'GA'}, 'GA': {'AT'}}, {}])
    assert sequencer.find_smallest_unique_k() == 2

@pytest.fixture
def dna_sequencer():
    return DNASequencer("path/to/testfile.fasta")

def test_empty_file(dna_sequencer, mocker):
    mocker.patch('builtins.open', mocker.mock_open(read_data=""))
    with pytest.raises(SystemExit) as e:
        dna_sequencer.process_sequences(3)
    assert e.type == SystemExit
    assert e.value.code == 1
