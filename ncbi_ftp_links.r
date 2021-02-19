#Download and read in ncbi genome data (< 3 MB file)
download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt", destfile=paste0(Sys.Date(), "_eukaryotes.txt"))
ncbi <- read.csv(paste0(Sys.Date(), "_eukaryotes.txt"), header=TRUE, sep="\t")
#Filter for assemblies with annotated proteins
ncbi <- ncbi[ncbi$Proteins != "-",]
#Filter for Fusarium taxa
ncbi.filtered <- ncbi[grep("Fusarium", ncbi$X.Organism.Name),]
ncbi.filtered <- ncbi.filtered[order(ncbi.filtered$X.Organism.Name),]
#Select taxa
ncbi.filtered <- ncbi.filtered[c(1:19, 34:36, 42, 44:54, 62, 65, 70, 76, 83, 84, 87:90, 92:96, 100:104, 108:109),]

#Download and read in file with ftp links to assemblies (< 300 MB file)
download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt", destfile=paste0(Sys.Date(), "_assembly_summary_genbank.txt"))
assembly.sum <- read.csv(paste0(Sys.Date(), "_assembly_summary_genbank.txt"), skip=1, header=TRUE, sep="\t", quote="")
#Match to ncbi data
assembly.sum <- assembly.sum[match(ncbi.filtered$Assembly.Accession, assembly.sum$X..assembly_accession),]

#Get FTP links and write to file
genomic.ftp <- paste0(assembly.sum$ftp_path, "/", assembly.sum$X..assembly_accession, "_", assembly.sum$asm_name, "_genomic.gbff.gz")
genomic.ftp <- gsub(" ", "_", genomic.ftp)
write(genomic.ftp, file="fus_ncbi_genomic")

protein.ftp <- paste0(assembly.sum$ftp_path, "/", assembly.sum$X..assembly_accession, "_", assembly.sum$asm_name, "_protein.faa.gz")
protein.ftp <- gsub(" ", "_", protein.ftp)
write(protein.ftp, file="fus_ncbi_proteins")