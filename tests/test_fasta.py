import csv
import os
import pytest
from cofolding_pulldown_tools.fasta import (
    clean_fasta, slice_fasta, generate_bait_prey, reformat_sliced_windows,
)

# Test data defined as strings
RAW_FASTA = """>sp|A0A087X1C5|CP2D7_HUMAN Cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=1 SV=1
MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQLRRRFGDVFSLQLAWTPVVVLNGLAAVREAMVTRGEDTADRPPAPIYQVLGFGPRSQGVILSRYGPAWREQRRFSVSTLRNLGLGKKSLEQWVTEEAACLCAAFADQAGRPFRPNGLLDKAVSNVIASLTCGRRFEYDDPRFLRLLDLAQEGLKEESGFLREVLNAVPVLPHIPALAGKVLRFQKAFLTQLDELLTEHRMTWDPAQPPRDLTEAFLAKKEKAKGSPESSFNDENLRIVVGNLFLAGMVTTSTTLAWGLLLMILHLDVQRGRRVSPGCPIVGTHVCPVRVQQEIDDVIGQVRRPEMGDQAHMPCTTAVIHEVQHFGDIVPLGVTHMTSRDIEVQGFRIPKGTTLITNLSSVLKDEAVWKKPFRFHPEHFLDAQGHFVKPEAFLPFSAGRRACLGEPLARMELFLFFTSLLQHFSFSVAAGQPRPSHSRVVSFLVTPSPYELCAVPR
>sp|A0A0B4J2F0|PIOS1_HUMAN Protein PIGBOS1 OS=Homo sapiens OX=9606 GN=PIGBOS1 PE=1 SV=1
MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS"""

CLEANED_FASTA = """>A0A087X1C5_CP2D7_HUMAN
MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQLRRRFGDVFSLQLAWTPVVVLNGLAAVREAMVTRGEDTADRPPAPIYQVLGFGPRSQGVILSRYGPAWREQRRFSVSTLRNLGLGKKSLEQWVTEEAACLCAAFADQAGRPFRPNGLLDKAVSNVIASLTCGRRFEYDDPRFLRLLDLAQEGLKEESGFLREVLNAVPVLPHIPALAGKVLRFQKAFLTQLDELLTEHRMTWDPAQPPRDLTEAFLAKKEKAKGSPESSFNDENLRIVVGNLFLAGMVTTSTTLAWGLLLMILHLDVQRGRRVSPGCPIVGTHVCPVRVQQEIDDVIGQVRRPEMGDQAHMPCTTAVIHEVQHFGDIVPLGVTHMTSRDIEVQGFRIPKGTTLITNLSSVLKDEAVWKKPFRFHPEHFLDAQGHFVKPEAFLPFSAGRRACLGEPLARMELFLFFTSLLQHFSFSVAAGQPRPSHSRVVSFLVTPSPYELCAVPR
>A0A0B4J2F0_PIOS1_HUMAN
MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS"""

def test_clean_fasta(tmp_path):
    input_file = tmp_path / "x.fasta"
    input_file.write_text(RAW_FASTA)
    
    # Run the function
    output_path = clean_fasta(input_file)
    
    # Read output
    with open(output_path, "r") as f:
        output_content = f.read()
        
    assert output_content.strip() == CLEANED_FASTA.strip()
    
    # Cleanup
    os.remove(output_path)

def test_slice_fasta(tmp_path):
    input_file = tmp_path / "x_cleaned.fasta"
    input_file.write_text(CLEANED_FASTA)
    
    # Run the function (matching the parameters that likely generated the test data)
    output_path = slice_fasta(input_file, max_length=400, window=None)
    
    # Read output
    with open(output_path, "r") as f:
        output_content = f.read()
    
    # Very basic assertion: verify header exists
    assert ">A0A087X1C5_CP2D7_HUMAN,1-258" in output_content
    
    # Cleanup
    os.remove(output_path)

def test_reformat_sliced_windows_uses_real_gene_symbol(tmp_path):
    """`clean_fasta` names windows after the UniProt entry mnemonic, which can
    diverge from the real gene symbol (CAD23 is the mnemonic for CDH23).
    reformat_sliced_windows must recover the real symbol from --gene_map."""
    gene_map_file = tmp_path / "original.fasta"
    gene_map_file.write_text(
        ">sp|Q9H251|CAD23_HUMAN Cadherin-23 OS=Homo sapiens OX=9606 GN=CDH23 PE=1 SV=1\n"
        "MNCPVLSLGSGFLFQVIEMLIFAYFASISLTESRGLFPRLENVGAFKKVSIVPTQAVCGL\n"
    )

    sliced_file = tmp_path / "long_proteins_cleaned_sliced.fasta"
    sliced_file.write_text(
        ">Q9H251_CAD23_HUMAN,1-30\n"
        "MNCPVLSLGSGFLFQVIEMLIFAYFASI\n"
        ">Q9H251_CAD23_HUMAN,20-61\n"
        "FASISLTESRGLFPRLENVGAFKKVSIVPTQAVCGL\n"
    )

    fasta_path, manifest_path = reformat_sliced_windows(sliced_file, gene_map_file)

    with open(fasta_path) as f:
        content = f.read()
    assert ">Q9H251.w1-30 GN=CDH23" in content
    assert ">Q9H251.w20-61 GN=CDH23" in content
    assert "CAD23" not in content  # the mnemonic must not leak into the output

    with open(manifest_path, newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 2
    assert rows[0] == {
        "window_id": "Q9H251.w1-30", "accession": "Q9H251", "gene": "CDH23",
        "start": "1", "end": "30", "window_len": "28",
    }

    os.remove(fasta_path)
    os.remove(manifest_path)


def test_reformat_sliced_windows_isoform_accession(tmp_path):
    """Isoform accessions (e.g. Q99102-10) contain a hyphen; the window tag
    must be period-joined so downstream accession/gene parsing on the first
    underscore still isolates the accession correctly."""
    gene_map_file = tmp_path / "original.fasta"
    gene_map_file.write_text(
        ">sp|Q99102-10|MUC4_HUMAN Isoform 10 of Mucin-4 OS=Homo sapiens OX=9606 GN=MUC4\n"
        "MTPGTQSPFFLLLLLTVLTVVTG\n"
    )
    sliced_file = tmp_path / "x_cleaned_sliced.fasta"
    sliced_file.write_text(
        ">Q99102-10_MUC4_HUMAN,1-23\n"
        "MTPGTQSPFFLLLLLTVLTVVTG\n"
    )

    fasta_path, manifest_path = reformat_sliced_windows(sliced_file, gene_map_file)
    with open(fasta_path) as f:
        content = f.read()
    assert ">Q99102-10.w1-23 GN=MUC4" in content

    os.remove(fasta_path)
    os.remove(manifest_path)


def test_generate_bait_prey_validated(tmp_path):
    """Verify bait-prey pairing, validation, and format."""
    input_file = tmp_path / "x_cleaned.fasta"
    input_file.write_text(CLEANED_FASTA)
    
    # Baits must exist in CLEANED_FASTA (which has >A0A087X1C5_CP2D7_HUMAN and >A0A0B4J2F0_PIOS1_HUMAN)
    # Bait list is comma-separated on input; output pairing lines are 'prey;bait'.
    bait = "A0A087X1C5_CP2D7_HUMAN,A0A0B4J2F0_PIOS1_HUMAN"
    
    # Run the function
    output_path = generate_bait_prey(input_file, bait, double_count=False)
    
    # Verify the output file exists
    assert os.path.exists(output_path)
    
    # Read lines and verify formatting (no newline issues, strict 'prey;bait')
    with open(output_path, "r") as f:
        lines = f.readlines()
    
    assert len(lines) > 0
    # Expected format: prey;bait
    for line in lines:
        parts = line.strip().split(";")
        assert len(parts) == 2
        assert parts[0] in ["A0A087X1C5_CP2D7_HUMAN", "A0A0B4J2F0_PIOS1_HUMAN"]
        assert parts[1] in ["A0A087X1C5_CP2D7_HUMAN", "A0A0B4J2F0_PIOS1_HUMAN"]

    # Verify validation works (should raise ValueError)
    with pytest.raises(ValueError, match="not found in the FASTA file"):
        generate_bait_prey(input_file, "INVALID_BAIT", double_count=False)

    # Cleanup
    os.remove(output_path)
