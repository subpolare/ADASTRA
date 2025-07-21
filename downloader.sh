download=/home/subpolare/gtrd/to_download.txt
sshpass -p 'gMO74&]j' rsync -avz --ignore-existing --files-from=${download} --progress -e 'ssh -p 60011' autosome@85.118.228.170:/ /home/subpolare/gtrd/reads
