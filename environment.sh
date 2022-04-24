# /bin/sh
# To use: `source environment.sh`
module load cuda
export PATH=/p/cs458/$(whoami)/local/bin:/usr/local/cuda/bin/:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
export CPATH=/usr/local/cuda/targets/x86_64-linux/include:$CPATH
