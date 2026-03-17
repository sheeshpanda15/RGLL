#!/bin/bash -l 
#SBATCH --output=/users/k21181837/RGSS/D/CASE1.txt
#SBATCH --job-name=CASE1D
#SBATCH --nodes=1              # 任务数
#SBATCH --cpus-per-task=8   
#SBATCH --time=48:00:00
#SBATCH --mem=100G 
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=k21181837@kcl.ac.uk
echo "Hello, World! From $HOSTNAME"

# 初始化环境（不同集群路径可能不同，两个都写最稳）
source /etc/profile 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
module load r/4.3.0-gcc-13.2.0-withx-rmath-standalone-python-3.11.6
R_DIR=/users/k21181837/RGSS/D
cd "$R_DIR"

cd /users/k21181837/RGSS/D
Rscript /users/k21181837/RGSS/D/CASE1.R