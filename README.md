# automated_DbGaP_DataRetrieval
FireCloud methods that Retrieves DbGaP SRA files and converts them to BAM files (https://portal.firecloud.org/#methods/ekofman/dbgapDownloads/9)

# Example
Launch docker with mounted ngc file:


	docker run -it -v /path/to/mount/1000_genomes_test_data.ngc:example.ngc vanallenlab/aspera:1.0


Input credentials using vdb config:


	vdb-config --import mount/example.ngc
	output:
		Repository directory is: '/root/ncbi/dbGaP-0'.


cd to the appropraite directory:


    cd /root/ncbi/dbGaP-0/

    
Run prefetch to fetch SRA files by their SRR ID


    prefetch -t ascp -a "/usr/bin/ascp|/home/data/.aspera/connect/etc/asperaweb_id_dsa.openssh" SRR390728

Or, to download all files indicated in a .krt file:


	prefetch -t ascp -a "/usr/bin/ascp|/home/data/.aspera/connect/etc/asperaweb_id_dsa.openssh" example.krt


