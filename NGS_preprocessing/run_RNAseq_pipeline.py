#!/usr/bin/env python3

#############
# LIBRARIES #
#############

import os
import subprocess
import sys
import time

#############
# FUNCTIONS #
#############

def run_command(bashCommand):
    
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    return output.strip().decode('utf-8')


def write_PBS_header(f, email, cpu, mem, directory, jobname, queue, hpc_account):
    
    f.write('#!/usr/bin/env bash\n')
    f.write("#PBS -m a\n")
    f.write("#PBS -M {}\n".format(email))
    f.write("#PBS -d {}\n".format(directory))
    f.write("#PBS -l nodes=1:ppn={},mem={}g\n".format(cpu, mem))
    f.write("#PBS -N {}\n".format(jobname))
    f.write("#PBS -o {}/{}.o.txt\n".format(outputdirs['jobdir'], jobname))
    f.write("#PBS -e {}/{}.e.txt\n".format(outputdirs['jobdir'], jobname))
    f.write("#PBS -V\n")
    f.write("#PBS -q {}\n".format(queue))
    f.write("#PBS -A {}\n".format(hpc_account))
    f.write("echo 'Running on : ' `hostname`\n")
    f.write("echo 'Start Time : ' `date`\n")
    f.write("echo 'Command:'\n")
    f.write("echo '========'\n")
    

def write_PBS_command(f, command):
    
    f.write("echo '{}'\n".format(command))
    f.write('{}\n'.format(command))


def writenextjobsubmit(f, scriptsub):
    
    f.write("\n\n## submit next step\n")
    f.write("JID=`qsub {}`\n".format(scriptsub))
    f.write("while  [[ ! $JID =~ ^[[:digit:]] ]] ; do\n")
    f.write("  sleep 5\n")
    f.write("  tmpFile=$(mktemp)\n")
    f.write("  TRASH=`qstat -x 2>$tmpFile | grep {}`\n".format(scriptsub))
    f.write("  # can be non-zero on grep-not-found and qstat-failure. only resumbit on valid qstat check.\n")
    f.write("  if  [[ $? -ne 0 ]] ; then\n")
    f.write("     if [ ! -s $tmpFile ] ; then\n")
    f.write("         JID=`qsub {}`\n".format(scriptsub))
    f.write("     fi\n")
    f.write("  fi\n")
    f.write("  rm $tmpFile\n")
    f.write("done\n")


if __name__ == '__main__':
    
    ##########################
    # COMMAND LINE ARGUMENTS #
    ##########################
    
    # read arguments from command line
    args = sys.argv[1:]
    
    # check if argument passed is an existing directory
    # only keep first argument and remove last forward slash if present
    run_directory = args[0]
    if run_directory.endswith('/'):
        run_directory = run_directory[:-1]
    if os.path.isdir(run_directory) != True:
        raise OSError('{} is not a directory.'.format(run_directory))
    
    print('Run directory: {}'.format(run_directory))
    
    
    #####################
    # LOCATION BINARIES #
    #####################
    
    binaries = {'java': "/opt/software/sun-jre-8/bin/java",
                'R': "/opt/software/R/3.2.1/bin/R",
                'python': "/opt/software/python3.5.2/bin/python3.5",
                'trimmomatic': "/opt/software/Trimmomatic-0.36/trimmomatic-0.36.jar",
                'hisat': "/opt/software/hisat2-2.0.4/hisat2"}
    
    files = {'GRCh38_index_dir': "/opt/NGS/References/GRCh38/HISAT2_index/",
             'hisat2_base_indexname': "GRCh38_noalt_h2",
             'annotationfile': "/opt/NGS/References/GRCh38/GCF_000001405.31_GRCh38.p5_genomic.gff",
             'countfile': "readcounts.txt",
             'statsfile': "Statistics.csv"}
    
    
    ###############
    # OUTPUT DIRS #
    ###############
    
    outputdirs = {'jobdir': "{}/Job_output".format(run_directory),
                  'tmpfilesdir': "{}/tmp_files".format(run_directory),
                  'tmpbindir': "{}/tmp_binaries".format(run_directory),
                  'resultsdir': "{}/Results".format(run_directory)}
    
    
    ###################
    # OTHER VARIABLES #
    ###################
    
    others = {'polyAscript':'/home/shared_data_immuno/bin/polyA_removal.pl',
              'countscript':'/home/shared_data_immuno/bin/CountTable.py',
              'countjob':'{}/Counts.sh'.format(outputdirs['tmpbindir']),
              'email': 'nicolas.deneuter@uantwerpen.be',
              'emailfile':'{}/finishmail.txt'.format(outputdirs['tmpfilesdir'])}
    
    
    ######################
    # MAKE DIR STRUCTURE #
    ######################
    
    # make outputdirs if they don't exist already
    for name, dirpath in outputdirs.items():
        print('Creating {} directory at {}'.format(name, dirpath))
        try:
            os.mkdir(dirpath)
        except OSError:
            print('{} already exist'.format(dirpath))
            
    # find sample files
    data_files = []
    
    for main_dir, subdirs, dir_files in os.walk(run_directory):
        for file in dir_files:
            if file.endswith('.fastq.gz') and main_dir not in outputdirs.values():
                data_files.append('{}/{}'.format(main_dir, file))
    
    print('{} files found for processing'.format(len(data_files)))
    
    ##################
    # START ANALYSIS #
    ##################
    
    """
    Per Sample/input file
        1. Trim with Trimmomatic
        2. Trim polyA
        3. Map with HISAT2
    Wait for all to finish
    On all files
        4. Do counts
    """
        
    sample_names = []
    samples_done = []
    
    for file in data_files:
        
        sample_name = file.split('/')[-1].replace('.fastq.gz', '')
        sample_names.append(sample_name)
        print('Processing {}'.format(sample_name))
        
        # intermediate file names
        polyA_fastq = sample_name+"_polyA.fastq.gz"
        trimmomatic_fastq = sample_name+"_Trimmomatic.fastq.gz"
        samfile = sample_name+".sam"
        sampledone = outputdirs['tmpfilesdir']+'/{}.done'.format(sample_name)
        
        # job script names
        trimmomaticjob = outputdirs['tmpbindir']+'/{}.Trimmomatic.sh'.format(sample_name)
        polyAtrimjob = outputdirs['tmpbindir']+'/{}.polyAtrim.sh'.format(sample_name)
        hisatjob = outputdirs['tmpbindir']+'/{}.HISAT2.sh'.format(sample_name)
    
        samples_done.append(sampledone)
        
        # 1. Make Trimmomatic trim
        with open(trimmomaticjob, 'w') as f:
            write_PBS_header(f, others['email'], 6, 4, run_directory, 'Immuno.{}.Trimmomatic'.format(sample_name), 'batch', 'Immuno')
            write_PBS_command(f, "{} -jar {} SE -threads 6 -phred33 {} /tmp/{} HEADCROP:20 SLIDINGWINDOW:4:15 MINLEN:30".format(
                binaries['java'], binaries['trimmomatic'], file, trimmomatic_fastq))
            write_PBS_command(f, "cp /tmp/{} {}".format(trimmomatic_fastq, outputdirs['tmpfilesdir']))
            write_PBS_command(f, "rm /tmp/{}".format(trimmomatic_fastq))
            writenextjobsubmit(f, polyAtrimjob)
        run_command('chmod +x {}'.format(trimmomaticjob))
    
        # 2. Make PolyA trim job file
        with open(polyAtrimjob, 'w') as f:
            write_PBS_header(f, others['email'], 1, 2, run_directory, 'Immuno.{}.polyAtrim'.format(sample_name), 'batch', 'Immuno')
            write_PBS_command(f, "perl {} -f {}/{} -o /tmp/{}".format(
                others['polyAscript'], outputdirs['tmpfilesdir'], trimmomatic_fastq, polyA_fastq))
            write_PBS_command(f, "cp /tmp/{} {}".format(polyA_fastq, outputdirs['tmpfilesdir']))
            write_PBS_command(f, "rm /tmp/{}".format(polyA_fastq))
            writenextjobsubmit(f, hisatjob)
        run_command('chmod +x {}'.format(polyAtrimjob))
        
        # 3. Make HISAT2 mapping job file
        with open(hisatjob, 'w') as f:
            write_PBS_header(f, others['email'], 6, 8, run_directory, 'Immuno.{}.HISAT2'.format(sample_name), 'batch', 'Immuno')
            write_PBS_command(f, "export HISAT2_INDEXES={}".format(files['GRCh38_index_dir']))
            write_PBS_command(f, "{} -p 6 -x {} -U {}/{} -S /tmp/{} --time".format(
                binaries['hisat'], files['hisat2_base_indexname'], outputdirs['tmpfilesdir'], polyA_fastq, samfile))
            write_PBS_command(f, "cp /tmp/{} {}".format(samfile, outputdirs['tmpfilesdir']))
            write_PBS_command(f, "rm /tmp/{}".format(samfile))
            # create a file to indicate that the pipeline for this sample has finished
            write_PBS_command(f, "echo 1 > {}".format(sampledone))
        run_command('chmod +x {}'.format(hisatjob))
    
        # submit first step
        print("Executing: qsub {}".format(trimmomaticjob))
        print(run_command('qsub {}'.format(trimmomaticjob)))
    
    #############################
    ## WAIT FOR JOBS TO FINISH ##
    #############################
    
    print(samples_done)
    
    ## torque/pbs does not have the -sync options, hence this workaround:
    # check if last file created in pipeline is present for each sample
    finished = 0
    while finished == 0:
        finished = 1
        for checkifdone in samples_done:
            if os.path.isfile(checkifdone) != True:
                finished = 0
                time.sleep(60)
                break
    
    print('Proceed to count\n')
    
    # 4. Make Counting job
    with open(others['countjob'], 'w') as f:
        write_PBS_header(f, others['email'], 1, 6, run_directory, 'Immuno.Count', 'batch', 'Immuno')
        write_PBS_command(f, "{} {} -i {} -a {} -o {}".format(
            binaries['python'], others['countscript'], outputdirs['tmpfilesdir'], files['annotationfile'], files['countfile']))
        write_PBS_command(f, "cp /tmp/{} {}".format(files['countfile'], outputdirs['resultsdir']))
        write_PBS_command(f, "rm /tmp/{}".format(files['countfile']))
        write_PBS_command(f, "echo 1 > {}/analysis.done".format(outputdirs['resultsdir']))
    run_command('chmod +x {}'.format(others['countjob']))
    
    ## All files made, submit count job and wait for analysis to finish
    submitcommand = 'qsub {}'.format(others['countjob'])
    print("Executing: {}\n".format(submitcommand))
    run_command(submitcommand)
    while os.path.isfile('{}/analysis.done'.format(outputdirs['resultsdir'])) != True:
        time.sleep(60)
    print('Analysis done\n')
    
    with open('{}'.format(others['emailfile']), 'w') as f:
        f.write('To: {}\r\n'.format(others['email']))
        f.write('subject: Immuno Analysis Done\r\n')
        f.write('from: Immuno.Pipeline\r\n')
        f.write('The analysis is done, the data is in {}:\r\n'.format(outputdirs['resultsdir']))
    
    subprocess.Popen("sendmail -t < {}".format(others['emailfile']).split(), stdout=subprocess.PIPE)
    
    #DONE