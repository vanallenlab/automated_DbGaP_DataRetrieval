workflow downloadDbgapFiles {
    File ngcFilePath
    String projectFolderName
    String srrId
    Int diskSpace
    String memory

    call downloadDbgapFilesTask {
    	input:
        	ngcFilePath=ngcFilePath,
        	projectFolderName=projectFolderName,
        	srrId=srrId,
        	diskSpace=diskSpace,
        	memory=memory
    }

    call BwaAlignment {
    	input:
    		firstEndFastq=downloadDbgapFilesTask.fastq1,
			secondEndFastq=downloadDbgapFilesTask.fastq2
    }

    output {
        downloadDbgapFilesTask.fastq1
        downloadDbgapFilesTask.fastq2
        BwaAlignment.aligned_bwa_bam
        BwaAlignment.aligned_bwa_bam_index
    }
}

task downloadDbgapFilesTask {
    File ngcFilePath
    String projectFolderName
    String srrId
    Int diskSpace
    String memory

    command <<<
    	echo "Current directory and contents:"
        pwd
        ls -lh

        cd /

        echo "vdb-config --import ${ngcFilePath}"
        vdb-config --import ${ngcFilePath}

        echo "cd /root/ncbi/${projectFolderName}"
        cd /root/ncbi/${projectFolderName}

        echo "prefetch -t ascp -a /usr/bin/ascp|/home/data/.aspera/connect/etc/asperaweb_id_dsa.openssh ${srrId}"
        prefetch -t ascp -a "/usr/bin/ascp|/home/data/.aspera/connect/etc/asperaweb_id_dsa.openssh" ${srrId}

        echo "fastq-dump --split-files ${srrId}"
        fastq-dump -F --split-files ${srrId}

		echo "ls -lh"
        ls -lh

        echo "mv ${srrId}_1.fastq /cromwell_root/"
		mv ${srrId}_1.fastq /cromwell_root/
        echo "mv ${srrId}_2.fastq /cromwell_root/"
        mv ${srrId}_2.fastq /cromwell_root/

        echo "ls -lh /cromwell_root/"
        ls -lh /cromwell_root/
    >>>

    output {
        File fastq1="${srrId}_1.fastq"
        File fastq2="${srrId}_2.fastq"
    }

    runtime {
    	docker: "vanallenlab/dbgap:1.0"
        memory: "${memory}"
        disks: "local-disk ${diskSpace} HDD"
        bootDiskSizeGb: 60
        preemptible: 5
    }
}

task BwaAlignment {
	File refFasta
	File firstEndFastq
	File secondEndFastq
	String sampleName
	Int memoryGb
	Int diskSpaceGb
	Int cpu

	command <<<
    	# log resource usage for debugging purposes
       	function runtimeInfo() {
        	echo [$(date)]
        	echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        	echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        	echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
        	do runtimeInfo;
           	sleep 15;
       	done &

		/usr/gitc/bwa index ${refFasta}

		/usr/gitc/bwa aln ${refFasta} -q 5 -l 32 -k 2 -t ${cpu} -o 1 -f ${sampleName}.Homo_sapiens_assembly19.1.sai ${firstEndFastq}

		/usr/gitc/bwa aln ${refFasta} -q 5 -l 32 -k 2 -t ${cpu} -o 1 -f ${sampleName}.Homo_sapiens_assembly19.2.sai ${secondEndFastq}

		/usr/gitc/bwa sampe -P -f ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.sam ${refFasta} ${sampleName}.Homo_sapiens_assembly19.1.sai ${sampleName}.Homo_sapiens_assembly19.2.sai ${firstEndFastq} ${secondEndFastq}

    	echo "samtools view -S -b ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.sam > ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.bam"
    	samtools view -S -b ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.sam > ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.bam

        echo "samtools sort -o ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.sorted.bam ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.bam"
    	samtools sort -o ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.sorted.bam ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.bam

        echo "samtools index ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.bam"
        samtools index ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.sorted.bam

        echo "ls -lh"
        ls -lh
    >>>

	output {
		File aligned_bwa_bam = "${sampleName}.Homo_sapiens_assembly19.aligned_bwa.sorted.bam"
        File aligned_bwa_bam_index = "${sampleName}.Homo_sapiens_assembly19.aligned_bwa.sorted.bam.bai"
	}

	runtime {
		docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
		memory: "${memoryGb} GB"
		cpu: "${cpu}"
		disks: "local-disk ${diskSpaceGb} HDD"
	}
}
