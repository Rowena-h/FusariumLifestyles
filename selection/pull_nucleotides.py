# Biopython's SeqIO module handles sequence input/output
from Bio import SeqIO
import sys

def get_cds_feature_with_qualifier_value(record, name, value):
    """Function to look for CDS feature by annotation value in sequence record."""
    # Loop over the records and features
    for feature in record.features:
        if feature.type == "CDS" and value in feature.qualifiers.get(name, []):
    	    return feature

    # Could not find it
    return None

with open(sys.argv[1], "r") as file:
    proteins_id = file.read().splitlines()

with open(sys.argv[1] + "_nucl.fa", "w") as nt_output:
    for unit in proteins_id:

        sample = unit.partition(".faa_")[-3] + ".faa_"
        protein = unit.partition(".faa_")[-1]
        genome_records = SeqIO.parse("gbff_files/concat.gbff", "genbank")	
        print("Looking at " + protein)
        found_feature = False
        for record in genome_records:
            cds_feature = get_cds_feature_with_qualifier_value(record, "protein_id", protein)
            if cds_feature is None:
                continue # feature not found in this record - try the next record
	
            found_feature = True
            
            gene_sequence = cds_feature.extract(record.seq)
             
            # Output FASTA records - note \n means insert a new line.
            # This is a little lazy as it won't line wrap the sequence:
            nt_output.write(">%s%s\n%s\n" % (sample, protein, gene_sequence))

        if not found_feature:
            print("Error: could not find feature")

print("Done")
