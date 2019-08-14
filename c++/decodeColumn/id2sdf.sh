i=$(grep -n $1 id.txt | cut -d: -f1)
is=$((i-2))
ie=$((i-1))
cs=$(tail -c +$((1+8*is)) mconfs.u64 | head -c 8 | decodeColumn-u64)
ce=$(tail -c +$((1+8*ie)) mconfs.u64 | head -c 8 | decodeColumn-u64)
os=$(tail -c +$((-7+8*cs)) conformers.ftr | head -c 8 | decodeColumn-u64)
oe=$(tail -c +$((-7+8*ce)) conformers.ftr | head -c 8 | decodeColumn-u64)
tail -c +$((1+os)) conformers.sdf | head -c $((oe-os)) > $1.sdf
