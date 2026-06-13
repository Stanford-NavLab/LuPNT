SO="/Users/keidaiiiyama/anaconda3/envs/lupnt/lib/python3.12/site-packages/tecsimpy/_tecsimpy.cpython-312-darwin.so"

DYLIB=$(dyld_info -dependents "$SO" 2>/dev/null | grep libtec_sim.dylib || true)
echo "$DYLIB"

# brute-force locate likely copies
python -c "import sys,glob; import site;
paths=set(site.getsitepackages()+[site.getusersitepackages()]);
hits=[];
import os
for p in paths:
    hits += glob.glob(os.path.join(p,'**','libtec_sim.dylib'), recursive=True)
print('\n'.join(hits) if hits else 'No libtec_sim.dylib under site-packages')"
