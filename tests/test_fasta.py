import os
import pytest
from cofolding_pulldown_tools.fasta import clean_fasta, slice_fasta, generate_bait_prey

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

def test_generate_bait_prey_validated(tmp_path):
    """Verify bait-prey pairing, validation, and format."""
    input_file = tmp_path / "x_cleaned.fasta"
    input_file.write_text(CLEANED_FASTA)
    
    # Baits must exist in CLEANED_FASTA (which has >A0A087X1C5_CP2D7_HUMAN and >A0A0B4J2F0_PIOS1_HUMAN)
    bait = "A0A087X1C5_CP2D7_HUMAN;A0A0B4J2F0_PIOS1_HUMAN"
    
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
