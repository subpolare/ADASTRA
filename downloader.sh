download=/home/v.nachatoy/gtrd/batch${1}_download.txt
sshpass -p 'gMO74&]j' rsync -avz --ignore-existing --files-from=${download} --progress -e 'ssh -p 60011' autosome@85.118.228.170:/ /mnt/data/v.nachatoy/gtrd/reads
