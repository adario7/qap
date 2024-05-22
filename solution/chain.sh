set -e

sh sync.sh
python extract.py
python process.py
python stats.py
