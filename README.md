# MutationFinder

web application to identify and visualize single nucleotide polymorphisms (SNPs) in a sample whole genome sequencing dataset.

Before Starting of the Project Let's dicuss some of the tools that I have used to complete this project. I have list steps for mac

## How to Install and Use Tools

### Install Homebrew

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

## Install BWA.

### BWA is used for aligning short reads to the reference genome:

```bash
brew install bwa
```

## Samtools is required for manipulating and processing the SAM/BAM files:

```bash
brew install samtools
```

## Download the Reference Genome

### Navigate to your working directory

```bash cd /path/to/your/working/directory````

### Download the reference genome

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/571/405/GCA_002571405.2_ASM257140v2/GCA_002571405.2_ASM257140v2_genomic.fna.gz
```

### Download the Short Read Data

### You can use the SRA toolkit to download the short read data, but an easier method is to use the fastq-dump command from the SRA toolkit:

### Install the SRA toolkit

```bash
brew install sratoolkit
```

### Download the short read data

```bash
fastq-dump --split-files --gzip SRR5892450
```

### Align the short reads to the reference genome using BWA:

### Index the reference genome

```bash
bwa index GCA_002571405.2_ASM257140v2_genomic.fna
```

### Align the short reads to the reference genome

```bash
bwa mem GCA_002571405.2_ASM257140v2_genomic.fna SRR5892450_1.fastq.gz SRR5892450_2.fastq.gz > alignment.sam
```

### Convert the SAM file to BAM format and sort it using Samtools:

### Convert SAM to BAM

```bash
samtools view -S -b alignment.sam > alignment.bam
```

### Sort the BAM file

```bash
samtools sort alignment.bam -o alignment_sorted.bam
```

### Index the BAM file

```bash
samtools index alignment_sorted.bam

Mutation-finder/
│
├── backend/
│ ├── app.py # Flask application entry point containing routes
│ ├── dataset/ # Directory to store genome data and outputs
│ │ ├── GCA_002571405.2_ASM257140v2_genomic.fna # Reference genome file
│ │ ├── alignment_reads_sorted.bam # Aligned BAM file (output from alignment step)
│ │ ├── output_variants.vcf # VCF file containing identified SNPs (output from variant caller)
│ │
│ ├── requirements.txt # Python dependencies for the Flask backend
│ └── README.md # Instructions for setting up and running the backend
```

### To run the Backend Server Follow these Steps.

## Running the Backend Server

Follow these steps to run the backend server:

- First, move into the `MutationFinder` directory.
- Install all the dependencies from `requirements.txt`.
- Start the application using the command:
  ```bash
    python3 app.py
  ```

```bash
frontend
│
├── snp-visualization/
│ ├── react-app/ # React application directory (if using React)
│ │ ├── src/
│ │ │ ├── App.js # Main React component
| | | |-- components
│ │ │ | ├── GeonmeViewer.js # Component to visualize SNPs
│ │ │ | ├── ReferenceSequenceViewer.js # Component to display the reference sequence
│ │ │ ├── SNPSummaryTable.js # Component to show SNPs in a table format
│ │ | |-- VarianceBoxPlot.js #Displays the QUAL value and position
│ │ | |-- ZoomControls.js # Zoom in and Zoom out options.
│ │ ├── public/
│ │ ├── package.json # React dependencies
│ │ └── README.md # Instructions for setting up and running the

```
