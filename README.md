# Automated dbGaP Data Retrieval
FireCloud methods that Retrieves DbGaP SRA files and converts them to BAM files.

Check in what format the SRA files were initially submitted to dbGaP. If you want to check whether a run is aligned you can run sra-stat on it:

	$ sra-stat -qx SRR7228168

If there is a <AlignInfo> section then it was submitted as aligned, if not then it was not submitted with alignment info.
There are 3 different methods on FireCloud to account for 3 different scenarios and data storage considerations:

#### ekofman/dbgap_download_fastqs_and_keep_only_bams
* This relies on the SRA file having been originally submitted as paired read fastqs (without alignment info) and not as an aligned bam.
* This stores and returns only bams, generated using bwa 

#### ekofman/dbgapDownloads_unaligned
* This relies on the SRA file having been originally submitted as paired read fastqs (without alignment info) and not as an aligned bam.
* This stores and returns both the fastqs and the bams, which means way more in storage costs. This could be useful if you're looking at RNA data for which you need the raw fastqs to do your own alignment or pseudoalignment. 

#### ekofman/dbgapDownloads
* This relies on the SRA file having been originally submitted as a bam (with alignment info) and not as fastqs
* This stores and returns only bams

# Docker
The docker docker image for all of these methods is built automatically from the Dockerfile at the top level of this repository: https://hub.docker.com/r/vanallenlab/dbgap/. Never remove or move this Dockerfile from the top level of this repository or the docker build will break. This docker image contains bwa, the SRA Toolkit, Aspera, ascp, and samtools.

# Example
Launch docker with mounted ngc file:


	docker run -it -v /path/to/mount/1000_genomes_test_data.ngc:example.ngc dbgap/dbgap:1.0


Input credentials using vdb config:


	vdb-config --import example.ngc
	output:
		Repository directory is: '/root/ncbi/dbGaP-0'.


cd to the appropriate directory (prefetch only works from within the directory):


    cd /root/ncbi/dbGaP-0/

    
Run prefetch to fetch SRA files by their SRR ID


    prefetch -t ascp -a "/usr/bin/ascp|/home/data/.aspera/connect/etc/asperaweb_id_dsa.openssh" SRR390728

Or, to download all files indicated in a .krt file:


	prefetch -t ascp -a "/usr/bin/ascp|/home/data/.aspera/connect/etc/asperaweb_id_dsa.openssh" example.krt



