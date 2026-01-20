# bash ./testsr.sh
source activate base
conda activate ~/.conda/envs/py_reos_dev

DIR="./tests"  

echo "going into: $DIR"
which python


files=("$DIR"/*.py)

if [ ${#files[@]} -eq 0 ]; then
    echo "no .py file found."
    exit 0
fi

for f in ${files[@]}; do
    echo "=============================="
    echo "Running: $f"
    echo "=============================="
    python3 "$f"
    echo
done
