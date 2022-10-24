import os
import subprocess


nThreads = 6
topHatAnchor = 6
StarIndex = "/data/bioinformatics/referenceGenome/Homo_sapiens/UCSC/hg38/Sequence/starIndex.v2.7.3"
gtf = "/data/bioinformatics/referenceGenome/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"


def getstatusoutput(cmd):
    """ behave like commands.getstatusoutput which moved to
    subprocess.getstatusoutput in Python 3.
    Implmented with subprocess.check_output which is available
    in both Python 2 and 3.
    """
    status = 0
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        status = e.returncode
        output = e.output

    output = output.decode()
    if output[-1] == '\n':
        output = output[:-1]

    return status, output


def doSTARMapping(fastq):
    bams = []
    base_path = "/data/bioinformatics/projects/biohub/marco2022/"
    print("STAR Mapping...")
    for i, each in enumerate(fastq):
        prefix = each.split("/")[-1].split(".")[0]
        map_folder = base_path + "star_results_IFN/" + prefix
        if not os.path.exists(map_folder):
            # if os.path.isdir(map_folder):
            #     os.rmdir(map_folder)
            # else:
            #     os.unlink(map_folder)
            os.makedirs(map_folder)

        cmd = 'STAR --twopassMode Basic '
        cmd += ' --chimSegmentMin 2 --outFilterMismatchNmax 3'
        cmd += ' --runThreadN ' + str(max([4, nThreads]))
        cmd += ' --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate '
        cmd += '--alignSJDBoverhangMin ' + str(topHatAnchor)
        cmd += ' --alignIntronMax 299999 --genomeDir ' + StarIndex
        cmd += ' --sjdbGTFfile ' + gtf
        cmd += ' --outFileNamePrefix ' + map_folder + '/ --readFilesIn '
        # if not allow_clipping:
        #     cmd += ' --alignEndsType EndToEnd'
        cmd += each
        if each.endswith('.gz'):
            cmd += ' --readFilesCommand zcat'

        print("mapping sample %d / %d, STAR Command: \n\t%s" % (i + 1, len(fastq), cmd))
        status, output = getstatusoutput(cmd)
        print("mapping %s is done with status %s" % (each, status))

        if int(status) != 0:
            print("\n*********** error in mapping sample_%d, %s *************" % (i, each))
            print("exit status: ", status)
            print("error detail: \n%s\n" % output)
            # raise Exception()
        else:
            print(output)
        bams.append(os.path.join(map_folder, 'Aligned.sortedByCoord.out.bam'))

    b1 = open(base_path + "input/b1_IFN.txt", "w")
    b2 = open(base_path + "input/b2_IFN.txt", "w")
    for i, b in enumerate(bams):
        if i < 6:
            b1.write(b)
            if i != 5:
                b1.write(",")
        else:
            b2.write(b)
            if i != 11:
                b1.write(",")
        print(b)
    b1.close()
    b2.close()


if __name__ == '__main__':
    raw_data_path = "/data/bioinformatics/projects/biohub/marco2022/raw/"
    s1_path = raw_data_path + "s1_IFN.txt"
    s2_path = raw_data_path + "s2_IFN.txt"
    s1 = open(s1_path).readline().strip().split(",")
    s2 = open(s2_path).readline().strip().split(",")
    Fastq = s1 + s2

    doSTARMapping(Fastq)
