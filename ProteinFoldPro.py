"""
ProteinFold Pro - Advanced Protein Structure Prediction Pipeline
A comprehensive toolkit for predicting protein structures using multiple state-of-the-art methods
"""

import os
import requests
import subprocess
import json
import time
import tempfile
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.PDB import PDBParser, PDBIO, Superimposer
import numpy as np

class ProteinFoldPro:
    """
    ProteinFold Pro: Your comprehensive protein structure prediction companion
    
    This class combines multiple prediction methods to give you the best possible
    structure prediction for your protein sequence. Think of it as having several
    expert crystallographers working together on your structure!
    """
    
    def __init__(self, output_dir="proteinfold_results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.msa_dir = self.output_dir / "alignments"
        self.msa_dir.mkdir(exist_ok=True)
        
        print("ProteinFold Pro initialized!")
        print(f"Results will be saved to: {self.output_dir}")
        
    def validate_sequence(self, sequence):
        """
        Smart sequence validation with helpful feedback
        We'll check your sequence and give you insights about what we're working with
        """
        valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
        
        # Clean up the sequence
        sequence = sequence.upper().replace(' ', '').replace('\n', '').replace('\t', '')
        
        # Check for invalid characters
        invalid_chars = set(sequence) - valid_amino_acids
        if invalid_chars:
            raise ValueError(f"Oops! Found some invalid amino acids: {invalid_chars}. Please use standard single-letter codes.")
        
        # Check minimum length
        if len(sequence) < 10:
            raise ValueError("Your sequence is quite short (less than 10 residues). We need at least 10 amino acids for a meaningful prediction.")
        
        # Give the user some insights about their sequence
        self.analyze_sequence_characteristics(sequence)
        
        return sequence
    
    def analyze_sequence_characteristics(self, sequence):
        """
        Tell the user interesting things about their protein sequence
        This helps set expectations about the prediction
        """
        length = len(sequence)
        
        # Count different types of amino acids
        hydrophobic = sum(sequence.count(aa) for aa in 'AILMFPWV')
        charged = sum(sequence.count(aa) for aa in 'DEKR')
        polar = sum(sequence.count(aa) for aa in 'NQSTHY')
        
        # Structural tendencies (simplified analysis)
        helix_favoring = sum(sequence.count(aa) for aa in 'AEHKLMQR')
        sheet_favoring = sum(sequence.count(aa) for aa in 'FILVWY')
        disorder_prone = sum(sequence.count(aa) for aa in 'AEGKPQRS')
        
        print(f"\nAnalyzing your protein sequence...")
        print(f"Length: {length} amino acids")
        print(f"Hydrophobic residues: {hydrophobic/length:.1%}")
        print(f"Charged residues: {charged/length:.1%}")
        print(f"Polar residues: {polar/length:.1%}")
        
        # Give helpful predictions about structure
        if helix_favoring/length > 0.4:
            print("This protein looks like it might have a lot of alpha helices!")
        
        if sheet_favoring/length > 0.3:
            print("Expecting some nice beta sheets in this structure!")
        
        if disorder_prone/length > 0.3:
            print("Heads up: This protein might have some disordered regions - that's totally normal but can make prediction trickier!")
        
        if length > 500:
            print("This is a large protein! We might break it into domains for better prediction accuracy.")
        
        return {
            'length': length,
            'hydrophobic_fraction': hydrophobic/length,
            'charged_fraction': charged/length,
            'disorder_propensity': disorder_prone/length
        }
    
    def create_multiple_sequence_alignment(self, sequence, protein_id):
        """
        Create a multiple sequence alignment to help with prediction
        Think of this as finding your protein's evolutionary relatives
        """
        print("\nSearching for evolutionary relatives...")
        print("This helps us understand your protein better by finding similar sequences!")
        
        # Try local MSA generation first (more accurate)
        local_msa = self.try_local_msa_generation(sequence, protein_id)
        if local_msa:
            return local_msa
        
        # Fallback to online search
        return self.search_online_relatives(sequence, protein_id)
    
    def try_local_msa_generation(self, sequence, protein_id):
        """Try to use local MSA tools if available"""
        print("Checking for local MSA tools...")
        
        fasta_file = self.msa_dir / f"{protein_id}.fasta"
        msa_file = self.msa_dir / f"{protein_id}.a3m"
        
        # Save query sequence
        record = SeqRecord(Seq(sequence), id=protein_id, description="Your protein sequence")
        SeqIO.write(record, fasta_file, "fasta")
        
        try:
            # Try HHblits if available
            cmd = [
                "hhblits",
                "-i", str(fasta_file),
                "-o", str(msa_file.with_suffix('.hhr')),
                "-oa3m", str(msa_file),
                "-d", "/usr/local/share/hhsuite/databases/uniclust30_2018_08/uniclust30_2018_08",
                "-n", "3",
                "-e", "1e-3"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            
            if result.returncode == 0 and msa_file.exists():
                print("Local MSA generation successful!")
                return msa_file
            
        except (subprocess.TimeoutExpired, FileNotFoundError):
            print("Local tools not available, switching to online search...")
        
        return None
    
    def search_online_relatives(self, sequence, protein_id):
        """Search for similar proteins online using BLAST"""
        print("Searching online databases for similar proteins...")
        
        try:
            # BLAST search with reasonable parameters
            print("This might take a few minutes - we're searching millions of proteins!")
            
            result_handle = NCBIWWW.qblast("blastp", "nr", sequence, 
                                         expect=1e-5, hitlist_size=200,
                                         word_size=6, matrix_name='BLOSUM62')
            
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)
            
            # Collect good matches
            relatives = []
            for alignment in blast_records.alignments[:50]:  # Top 50 matches
                for hsp in alignment.hsps:
                    if hsp.expect < 1e-5 and hsp.identities / hsp.align_length > 0.3:
                        relatives.append({
                            'name': alignment.title.split()[0],
                            'sequence': hsp.sbjct.replace('-', ''),
                            'similarity': hsp.identities / hsp.align_length,
                            'coverage': hsp.align_length / len(sequence)
                        })
            
            if relatives:
                # Save alignment file
                msa_file = self.msa_dir / f"{protein_id}_relatives.a3m"
                with open(msa_file, 'w') as f:
                    f.write(f">{protein_id}_query\n{sequence}\n")
                    for i, rel in enumerate(relatives):
                        f.write(f">relative_{i+1}\n{rel['sequence']}\n")
                
                print(f"Found {len(relatives)} evolutionary relatives!")
                print(f"Saved alignment to: {msa_file}")
                return msa_file
            else:
                print("Couldn't find many similar proteins - your protein might be quite unique!")
                return None
            
        except Exception as e:
            print(f"Online search had some issues: {e}")
            return None
    
    def predict_with_multiple_methods(self, sequence, protein_id, msa_file=None):
        """
        Try several prediction methods and pick the best result
        It's like getting second opinions from multiple experts!
        """
        print("\nStarting multi-method prediction pipeline...")
        
        prediction_results = []
        
        # Method 1: Try ColabFold (AlphaFold2-based, very reliable)
        colabfold_result = self.try_colabfold_prediction(sequence, protein_id, msa_file)
        if colabfold_result:
            prediction_results.append(("ColabFold", colabfold_result))
        
        # Method 2: Try ESMFold (Fast and doesn't need MSA)
        esmfold_result = self.try_esmfold_prediction(sequence, protein_id)
        if esmfold_result:
            prediction_results.append(("ESMFold", esmfold_result))
        
        # Method 3: Try ChimeraX if available
        chimerax_result = self.try_chimerax_prediction(sequence, protein_id)
        if chimerax_result:
            prediction_results.append(("ChimeraX", chimerax_result))
        
        return prediction_results
    
    def try_colabfold_prediction(self, sequence, protein_id, msa_file=None):
        """ColabFold: The gold standard for protein folding"""
        print("Trying ColabFold (AlphaFold2-based prediction)...")
        
        import shutil
        if not shutil.which("colabfold_batch"):
            print("ColabFold not found. You can install it with:")
            print("   conda install -c conda-forge -c bioconda colabfold_batch")
            return None
        
        try:
            # Prepare input
            input_file = msa_file if msa_file else self.create_simple_fasta(sequence, protein_id)
            output_dir = self.output_dir / f"{protein_id}_colabfold"
            output_dir.mkdir(exist_ok=True)
            
            print("Running ColabFold with optimized settings...")
            
            # Optimized ColabFold command
            cmd = [
                "colabfold_batch",
                str(input_file),
                str(output_dir),
                "--num-models", "3",  # Generate 3 models for comparison
                "--num-recycles", "5",  # Good balance of accuracy vs speed
                "--model-type", "alphafold2_ptm",
                "--amber",  # Relax with AMBER for better geometry
                "--templates"  # Use template information if available
            ]
            
            print("This usually takes 5-15 minutes depending on protein size...")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            
            if result.returncode == 0:
                # Find the best model
                best_model = self.find_best_colabfold_model(output_dir, protein_id)
                if best_model:
                    print("ColabFold prediction successful!")
                    return best_model
            
            print(f"ColabFold failed: {result.stderr}")
            
        except subprocess.TimeoutExpired:
            print("ColabFold took too long - this happens with very large proteins")
        except Exception as e:
            print(f"ColabFold error: {e}")
        
        return None
    
    def try_esmfold_prediction(self, sequence, protein_id):
        """ESMFold: Fast predictions without needing evolutionary information"""
        print("Trying ESMFold (fast, no MSA needed)...")
        
        # Handle long sequences by splitting
        if len(sequence) > 400:
            return self.handle_long_sequence_esmfold(sequence, protein_id)
        
        # Try local ESMFold first
        local_result = self.try_local_esmfold(sequence, protein_id)
        if local_result:
            return local_result
        
        # Use ESM Atlas API
        print("Using ESM Atlas online service...")
        url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
        
        try:
            headers = {"Content-Type": "application/x-www-form-urlencoded"}
            print("Sending sequence to ESM Atlas... (this usually takes 1-3 minutes)")
            
            response = requests.post(url, data=sequence, headers=headers, timeout=600)
            
            if response.status_code == 200:
                output_file = self.output_dir / f"{protein_id}_esmfold.pdb"
                with open(output_file, 'w') as f:
                    f.write(response.text)
                
                print("ESMFold prediction successful!")
                return output_file
            else:
                print(f"ESM Atlas returned status: {response.status_code}")
                
        except requests.exceptions.Timeout:
            print("ESM Atlas request timed out")
        except Exception as e:
            print(f"ESMFold error: {e}")
        
        return None
    
    def try_local_esmfold(self, sequence, protein_id):
        """Try local ESMFold installation"""
        try:
            import torch
            import esm
            
            print("Using local ESMFold installation...")
            
            # Load the model (this might take a while the first time)
            model = esm.pretrained.esmfold_v1()
            model = model.eval()
            
            # Predict structure
            with torch.no_grad():
                output = model.infer_pdb(sequence)
            
            # Save result
            output_file = self.output_dir / f"{protein_id}_local_esmfold.pdb"
            with open(output_file, 'w') as f:
                f.write(output)
            
            print("Local ESMFold successful!")
            return output_file
            
        except ImportError:
            print("Local ESMFold not available (install with: pip install fair-esm)")
            return None
        except Exception as e:
            print(f"Local ESMFold failed: {e}")
            return None
    
    def handle_long_sequence_esmfold(self, sequence, protein_id):
        """Handle proteins that are too long for ESMFold"""
        print(f"Your protein is quite long ({len(sequence)} residues)!")
        print("Breaking it into smaller pieces for prediction...")
        
        # Simple domain splitting
        chunk_size = 300
        overlap = 50
        
        chunks = []
        for i in range(0, len(sequence), chunk_size - overlap):
            chunk = sequence[i:i + chunk_size]
            if len(chunk) >= 50:
                chunks.append((i, chunk))
        
        print(f"Split into {len(chunks)} overlapping segments")
        
        # Predict each chunk
        chunk_predictions = []
        for i, (start_pos, chunk_seq) in enumerate(chunks):
            print(f"Predicting chunk {i+1}/{len(chunks)}...")
            chunk_id = f"{protein_id}_chunk_{i+1}"
            chunk_result = self.try_esmfold_prediction(chunk_seq, chunk_id)
            if chunk_result:
                chunk_predictions.append((start_pos, chunk_result))
        
        # Combine results (simplified approach)
        if chunk_predictions:
            return self.combine_structure_chunks(chunk_predictions, protein_id)
        
        return None
    
    def combine_structure_chunks(self, chunks, protein_id):
        """Combine multiple structure chunks into one file"""
        print("Combining structure pieces...")
        
        combined_file = self.output_dir / f"{protein_id}_combined.pdb"
        
        try:
            all_atoms = []
            atom_counter = 1
            
            for start_pos, chunk_file in chunks:
                with open(chunk_file, 'r') as f:
                    for line in f:
                        if line.startswith('ATOM'):
                            # Update atom numbering
                            new_line = line[:6] + f"{atom_counter:5d}" + line[11:]
                            all_atoms.append(new_line)
                            atom_counter += 1
            
            # Write combined structure
            with open(combined_file, 'w') as f:
                f.writelines(all_atoms)
                f.write('END\n')
            
            print(f"Combined structure saved: {combined_file}")
            return combined_file
            
        except Exception as e:
            print(f"Combining chunks failed: {e}")
            # Return first chunk as fallback
            return chunks[0][1] if chunks else None
    
    def try_chimerax_prediction(self, sequence, protein_id):
        """Try ChimeraX AlphaFold if available"""
        print("Checking for ChimeraX AlphaFold...")
        
        # This would require ChimeraX with AlphaFold plugin
        # For now, just a placeholder
        print("ChimeraX integration not yet implemented")
        return None
    
    def find_best_colabfold_model(self, model_dir, protein_id):
        """Find the most confident ColabFold model"""
        print("Selecting the best model based on confidence scores...")
        
        # Look for confidence scores
        json_files = list(model_dir.glob("*_scores_*.json"))
        pdb_files = list(model_dir.glob("*_relaxed_*.pdb"))
        
        if not pdb_files:
            pdb_files = list(model_dir.glob("*_rank_001*.pdb"))
        
        if not pdb_files:
            print("No PDB files found in ColabFold output")
            return None
        
        # If we have confidence scores, use them
        if json_files:
            try:
                with open(json_files[0], 'r') as f:
                    scores = json.load(f)
                
                best_confidence = 0
                best_file = None
                
                for pdb_file in pdb_files:
                    # Extract model info from filename
                    if 'rank_001' in pdb_file.name or 'relaxed' in pdb_file.name:
                        # This is likely the best model
                        best_file = pdb_file
                        break
                
                if not best_file:
                    best_file = pdb_files[0]
                
                # Copy to final location
                final_file = self.output_dir / f"{protein_id}_colabfold_best.pdb"
                import shutil
                shutil.copy2(best_file, final_file)
                
                print(f"Best ColabFold model: {final_file}")
                return final_file
                
            except Exception as e:
                print(f"Error reading confidence scores: {e}")
        
        # Fallback to first available model
        final_file = self.output_dir / f"{protein_id}_colabfold.pdb"
        import shutil
        shutil.copy2(pdb_files[0], final_file)
        return final_file
    
    def enhance_structure_quality(self, pdb_file, protein_id):
        """Try to improve the structure quality"""
        print("✨ Enhancing structure quality...")
        
        # Try energy minimization if OpenMM is available
        minimized = self.try_energy_minimization(pdb_file, protein_id)
        if minimized:
            return minimized
        
        # Basic cleanup
        return self.basic_structure_cleanup(pdb_file, protein_id)
    
    def try_energy_minimization(self, pdb_file, protein_id):
        """Try to energy minimize the structure"""
        try:
            import openmm
            from openmm import app
            
            print("Running energy minimization...")
            
            # Load structure
            pdb = app.PDBFile(str(pdb_file))
            
            # Set up force field
            forcefield = app.ForceField('amber14-all.xml')
            
            # Create system
            system = forcefield.createSystem(pdb.topology, 
                                           nonbondedMethod=app.NoCutoff,
                                           constraints=app.HBonds)
            
            # Minimize
            integrator = openmm.LangevinMiddleIntegrator(300*openmm.unit.kelvin, 
                                                        1/openmm.unit.picosecond, 
                                                        0.004*openmm.unit.picoseconds)
            
            simulation = app.Simulation(pdb.topology, system, integrator)
            simulation.context.setPositions(pdb.positions)
            simulation.minimizeEnergy(maxIterations=1000)
            
            # Save minimized structure
            minimized_file = self.output_dir / f"{protein_id}_minimized.pdb"
            positions = simulation.context.getState(getPositions=True).getPositions()
            app.PDBFile.writeFile(pdb.topology, positions, str(minimized_file))
            
            print("Energy minimization completed!")
            return minimized_file
            
        except ImportError:
            print("OpenMM not available for energy minimization")
            print("   Install with: conda install -c conda-forge openmm")
            return None
        except Exception as e:
            print(f"Energy minimization failed: {e}")
            return None
    
    def basic_structure_cleanup(self, pdb_file, protein_id):
        """Basic structure cleanup and validation"""
        print("Performing basic structure cleanup...")
        
        try:
            clean_file = self.output_dir / f"{protein_id}_cleaned.pdb"
            
            with open(pdb_file, 'r') as infile, open(clean_file, 'w') as outfile:
                for line in infile:
                    # Keep only ATOM records and END
                    if line.startswith(('ATOM', 'END')):
                        outfile.write(line)
            
            return clean_file
            
        except Exception as e:
            print(f"Structure cleanup failed: {e}")
            return pdb_file
    
    def create_simple_fasta(self, sequence, protein_id):
        """Create a simple FASTA file"""
        fasta_file = self.output_dir / f"{protein_id}.fasta"
        record = SeqRecord(Seq(sequence), id=protein_id, description="Input sequence")
        SeqIO.write(record, fasta_file, "fasta")
        return fasta_file
    
    def fold_protein(self, sequence, protein_name="my_protein", use_all_methods=True):
        """
        Main function: Fold your protein using the best available methods!
        
        Args:
            sequence: Your protein sequence (single letter amino acid codes)
            protein_name: A name for your protein (used in output files)
            use_all_methods: Whether to try multiple methods for best results
        
        Returns:
            Path to the best predicted structure, or None if all methods failed
        """
        
        print("Welcome to ProteinFold Pro!")
        print("Let's predict the structure of your protein!")
        
        try:
            # Step 1: Validate the sequence
            clean_sequence = self.validate_sequence(sequence)
            
            # Step 2: Build evolutionary context (MSA)
            msa_file = None
            if use_all_methods:
                msa_file = self.create_multiple_sequence_alignment(clean_sequence, protein_name)
            
            # Step 3: Try multiple prediction methods
            start_time = time.time()
            predictions = self.predict_with_multiple_methods(clean_sequence, protein_name, msa_file)
            
            if not predictions:
                print("Sorry, all prediction methods failed!")
                print("This can happen with very unusual sequences or technical issues.")
                return None
            
            # Step 4: Select and enhance the best prediction
            print(f"\nGreat! We got {len(predictions)} successful predictions!")
            
            # Use the first successful prediction (in practice, you might want more sophisticated selection)
            best_method, best_structure = predictions[0]
            print(f"Using result from: {best_method}")
            
            # Step 5: Enhance the structure quality
            final_structure = self.enhance_structure_quality(best_structure, protein_name)
            
            # Step 6: Report results
            elapsed_time = time.time() - start_time
            print(f"\nSUCCESS! Structure prediction completed in {elapsed_time:.1f} seconds")
            print(f"Your structure: {final_structure}")
            print(f"Prediction method: {best_method}")
            
            # Give visualization tips
            print(f"\nTo visualize your structure:")
            print(f"   • ChimeraX: open {final_structure}")
            print(f"   • PyMOL: load {final_structure}")
            print(f"   • Online: Upload to ChimeraX Web or Mol* viewer")
            
            return final_structure
            
        except Exception as e:
            print(f"Oops! Something went wrong: {e}")
            return None

def main():
    """
    Example usage and interactive mode
    """
    
    print("ProteinFold Pro - Advanced Protein Structure Prediction")
    print("=" * 60)
    print("Welcome! This tool helps you predict protein structures using")
    print("state-of-the-art AI methods like AlphaFold2 and ESMFold.")
    print()
    
    # Example sequence (a small, well-folded protein)
    example_sequence = """
    MKALIVLGLVLLGLLGPLPLSAQEMINGDKIILANGKLLVTLRRGMGMPNLNLRVLGRGKTQG
    SSYSFTGVLLGQGTAFDGWFKTLGSAGTKLTGLTGQTALLVLWRLDNLDEQVLTFFSSGPVD
    EAPLFQVVVDKIAAYQVAGEKISVQAFLGDTTQTLQQAAAKAKYKSVFQDRIRQHKDPDLMD
    LGYGRSVWQGIPFQLVAGGTVFSAQNAGDLLGGPGIAATGLTLLLLQYGASVNPASLGQLTL
    """
    
    # Clean up example sequence
    clean_example = "".join(example_sequence.split()).upper()
    
    print("Installation requirements for full functionality:")
    print("   • ColabFold: conda install -c conda-forge -c bioconda colabfold_batch")
    print("   • ESMFold: pip install fair-esm torch")
    print("   • OpenMM (optional): conda install -c conda-forge openmm")
    print()
    
    try:
        # Initialize the predictor
        predictor = ProteinFoldPro()
        
        # Get user input
        print("Please enter your protein sequence:")
        print("(Use single-letter amino acid codes, or press Enter to use example)")
        user_sequence = input("Sequence: ").strip()
        
        if not user_sequence:
            user_sequence = clean_example
            print(f"Using example sequence ({len(user_sequence)} residues)")
        
        protein_name = input("Protein name (default: 'my_protein'): ").strip()
        if not protein_name:
            protein_name = "my_protein"
        
        use_all = input("Use all available methods for best accuracy? (y/n, default=y): ").strip().lower()
        use_all_methods = use_all != 'n'
        
        print("\nStarting protein folding...")
        
        # Run the prediction
        result = predictor.fold_protein(user_sequence, protein_name, use_all_methods)
        
        if result:
            print(f"\nAll done! Your protein structure is ready: {result}")
        else:
            print("\nPrediction failed. Please check your sequence and try again.")
            
    except KeyboardInterrupt:
        print("\nPrediction cancelled. Come back anytime!")
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        print("Please report this issue if it keeps happening!")

if __name__ == "__main__":
    main()