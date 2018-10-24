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

    output {
        downloadDbgapFilesTask.bam
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

        echo "sam-dump ${srrId} | samtools view -bS - > ${srrId}.bam"
        sam-dump ${srrId} | samtools view -bS - > ${srrId}.bam

		echo "ls -lh"
        ls -lh

		echo "mv ${srrId}.bam /cromwell_root/"
		mv ${srrId}.bam /cromwell_root/

        echo "ls -lh /cromwell_root/"
        ls -lh /cromwell_root/
    >>>

    output {
        File bam="${srrId}.bam"
    }

    runtime {
    	docker: "vanallenlab/dbgap:1.0"
        memory: "${memory}"
        disks: "local-disk ${diskSpace} HDD"
        bootDiskSizeGb: 60
        preemptible: 5
    }
}