workflow downloadDbgapFiles {
    File ngcFilePath
    String projectFolderName
    String srrId
    String sampleName
    Int diskSpace
    Int memory
    Int cpu
    File refFasta

    call downloadDbgapFilesTask {
    	input:
        	ngcFilePath=ngcFilePath,
        	projectFolderName=projectFolderName,
        	srrId=srrId,
        	sampleName=sampleName,
        	diskSpace=diskSpace,
        	memory=memory,
        	cpu=cpu,
        	refFasta=refFasta
    }

    output {
        BwaAlignment.aligned_bwa_bam
        BwaAlignment.aligned_bwa_bam_index
    }
}

task downloadDbgapFilesTask {
    File ngcFilePath
    String projectFolderName
    String srrId
    String sampleName
    Int diskSpace
    Int cpu
    File refFasta
    Int memory

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

		# Now we have ${srrId}_1.fastq
        # And we also have ${srrId}_2.fastq

        # Use bwa to align paired end fastq files
        /usr/gitc/bwa index ${refFasta}

		/usr/gitc/bwa aln ${refFasta} -q 5 -l 32 -k 2 -t ${cpu} -o 1 -f ${sampleName}.Homo_sapiens_assembly19.1.sai ${srrId}_1.fastq

		/usr/gitc/bwa aln ${refFasta} -q 5 -l 32 -k 2 -t ${cpu} -o 1 -f ${sampleName}.Homo_sapiens_assembly19.2.sai ${srrId}_2.fastq

		/usr/gitc/bwa sampe -P -f ${sampleName}.Homo_sapiens_assembly19.aligned_bwa.sam ${refFasta} ${sampleName}.Homo_sapiens_assembly19.1.sai ${sampleName}.Homo_sapiens_assembly19.2.sai ${srrId}_1.fastq ${srrId}_2.fastq

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
    	docker: "vanallenlab/dbgap:2.0"
        memory: "${memory} GB"
        cpu: "${cpu}"
        disks: "local-disk ${diskSpace} HDD"
        bootDiskSizeGb: 60
        preemptible: 5
    }
}
